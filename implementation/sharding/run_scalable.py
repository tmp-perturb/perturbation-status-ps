#!/usr/bin/env python3
"""Run pooled-global-union PS sharding as one reproducible workflow."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import subprocess
import shutil
import sys
import time
from pathlib import Path

import yaml


HERE = Path(__file__).resolve().parent
RSCRIPT = shutil.which("Rscript")
if RSCRIPT is None:
    raise FileNotFoundError("Required executable is not on PATH: Rscript")


def run(command: list[str], log: Path) -> dict[str, object]:
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[name] = "1"
    with log.open("w", encoding="utf-8") as stream:
        process = subprocess.run(["/usr/bin/time", "-v", *command], stdout=stream, stderr=subprocess.STDOUT, env=environment)
    text = log.read_text(encoding="utf-8")
    match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", text)
    if process.returncode != 0:
        raise RuntimeError(f"Command failed; see {log}")
    return {"wall_seconds": time.time() - started, "peak_rss_kib": int(match.group(1)) if match else None, "log": str(log)}


def parallel(commands: list[tuple[list[str], Path]], workers: int) -> list[dict[str, object]]:
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(run, command, log) for command, log in commands]
        return [future.result() for future in futures]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    work = Path(config["work_dir"])
    output = Path(config["output_dir"])
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"Output directory is not empty: {output}")
    work.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    logs = work / "logs"
    shard_count = int(config["shards"])
    workers = int(config["max_workers"])
    started = time.time()
    stages: dict[str, object] = {}

    partition_outputs = [work / "shard_manifest.json", work / "target_manifest.tsv"]
    if args.resume and all(path.exists() for path in partition_outputs):
        stages["partition"] = {"status": "REUSED"}
    else:
        stages["partition"] = run([
            sys.executable, str(HERE / "prepare_validation.py"),
            "--barcode", config["barcode_table"], "--output-dir", str(work), "--shards", str(shard_count),
            "--min-target-cells", str(config.get("min_target_cells", 3))
        ], logs / "partition.log")
    input_outputs = [work / "inputs" / f"shard_{shard:03d}.rds" for shard in range(1, shard_count + 1)]
    if args.resume and all(path.exists() for path in input_outputs):
        stages["input_shards"] = {"status": "REUSED"}
    else:
        stages["input_shards"] = run([
            RSCRIPT, str(HERE / "build_input_shards.R"), config["input_rds"],
            config["barcode_table"], str(work), str(work / "inputs"), str(shard_count)
        ], logs / "input_shards.log")

    selection_commands = []
    for shard in range(1, shard_count + 1):
        selection_output = work / f"shard_{shard:03d}_response_genes.tsv"
        if args.resume and selection_output.exists() and selection_output.stat().st_size > 0:
            continue
        selection_commands.append(([
            RSCRIPT, str(HERE / "select_response_genes.R"),
            str(work / "inputs" / f"shard_{shard:03d}.rds"), config["barcode_table"],
            str(work / f"shard_{shard:03d}_targets.txt"), str(work / f"shard_{shard:03d}_response_genes.tsv"),
            config["non_target_control"], str(config["parameters"]["logfc_threshold"])
        ], logs / f"selection_{shard:03d}.log"))
    stages["response_gene_selection"] = parallel(selection_commands, workers) if selection_commands else {"status": "REUSED"}
    stages["response_gene_union"] = run([
        sys.executable, str(HERE / "merge_response_genes.py"), "--selection-dir", str(work),
        "--output-dir", str(work), "--expected-targets", str(json.loads((work / "shard_manifest.json").read_text())["targets"])
    ], logs / "response_gene_union.log")
    stages["global_gene_filter"] = run([
        RSCRIPT, str(HERE / "freeze_modeled_genes.R"), config["input_rds"],
        str(work / "global_response_genes.txt"), str(work / "global_modeled_genes.txt")
    ], logs / "global_gene_filter.log")

    score_commands = []
    for shard in range(1, shard_count + 1):
        shard_output = work / "scores" / f"shard_{shard:03d}"
        score_commands.append(([
            RSCRIPT, str(HERE / "run_shard.R"),
            str(work / "inputs" / f"shard_{shard:03d}.rds"), config["barcode_table"],
            str(work / f"shard_{shard:03d}_targets.txt"), str(work / "global_modeled_genes.txt"),
            str(shard_output), str(int(config["seed"]) + shard)
        ], logs / f"score_{shard:03d}.log"))
    stages["scoring"] = parallel(score_commands, workers)
    merge_command = [
        sys.executable, str(HERE / "merge_validate.py"), "--shard-dir", str(work / "scores"),
        "--barcode", config["barcode_table"], "--output-dir", str(output),
        "--excluded-targets", str(work / "shard_manifest.json")
    ]
    if config.get("reference_result"):
        merge_command.extend(["--canonical", config["reference_result"]])
    stages["merge"] = run(merge_command, logs / "merge.log")

    manifest = {
        "module": config["module"], "method": config["method"], "dataset": config["dataset"],
        "execution_mode": config["execution_mode"], "shards": shard_count, "max_workers": workers,
        "numerical_threads_per_worker": 1, "wall_seconds": time.time() - started, "stages": stages,
        "status": "PASS"
    }
    temporary = output / "run_manifest.json.partial"
    temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    temporary.rename(output / "run_manifest.json")
    print(json.dumps({"dataset": config["dataset"], "wall_seconds": manifest["wall_seconds"], "status": "PASS"}))


if __name__ == "__main__":
    main()
