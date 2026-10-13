#!/usr/bin/env python3
"""Run PS as an independent perturbation-status Omnibenchmark module."""

from __future__ import annotations

import argparse
import gzip
import json
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

from common import dump_yaml, load_yaml, make_archive, rscript, run, write_json

HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE / "implementation"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--data.h5ad", required=True)
    parser.add_argument("--data.spec", required=True)
    parser.add_argument("--shared_signature.bundle", required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    # The status-estimation stage has a shared input contract. PS deliberately
    # accepts but does not read the Mixscape/Mixscale signature artifact.
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    spec = json.loads(Path(getattr(args, "data.spec")).read_text(encoding="utf-8"))
    h5ad = Path(getattr(args, "data.h5ad")).resolve()
    started = time.monotonic()

    with tempfile.TemporaryDirectory(prefix=f"pse-ps-{args.name}-") as temp_text:
        work = Path(temp_text)
        input_config = load_yaml(IMPLEMENTATION / "config" / f'{spec["dataset"]}.yaml')
        input_dir = work / "ps_input"
        input_config.update({
            "dataset": spec["dataset"], "source_h5ad": str(h5ad),
            "input_dir": str(input_dir), "shared_sparse_input_dir": "",
            "seurat_rds": str(input_dir / "ps_input.rds"),
        })
        input_path = work / "ps_input.yaml"
        dump_yaml(input_config, input_path)
        run([sys.executable, IMPLEMENTATION / "adapt_ps_input.py", "--config", input_path])
        run(
            [rscript(), IMPLEMENTATION / "build_ps_seurat.R", "--config", input_path],
        )
        config = load_yaml(IMPLEMENTATION / "config" / f'{spec["dataset"]}_scalable.yaml')
        ps_output = work / "ps"
        config.update({
            "dataset": spec["dataset"], "input_rds": str(input_dir / "ps_input.rds"),
            "barcode_table": str(input_dir / "barcode.tsv.gz"), "output_dir": str(ps_output),
            "work_dir": str(work / "ps_work"), "max_workers": args.workers,
        })
        config.pop("reference_result", None)
        config_path = work / "ps_run.yaml"
        dump_yaml(config, config_path)
        run(
            [sys.executable, IMPLEMENTATION / "sharding" / "run_scalable.py", "--config", config_path],
        )
        result = ps_output / "ps_metadata.tsv"
        if not result.is_file():
            raise FileNotFoundError(result)
        with result.open("rb") as source, gzip.open(output / "status_result.tsv.gz", "wb") as destination:
            shutil.copyfileobj(source, destination)
        make_archive(ps_output, output / "method_aux.tar.gz", "ps")

    write_json({
        "dataset": spec["dataset"], "method": "ps", "status": "PASS",
        "wall_seconds": time.monotonic() - started, "workers": args.workers,
        "peak_child_rss_kib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "result_file": "status_result.tsv.gz", "aux_archive": "method_aux.tar.gz",
    }, output / "method_manifest.json")


if __name__ == "__main__":
    main()
