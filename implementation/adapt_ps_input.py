#!/usr/bin/env python3
"""Build deterministic scMAGeCK PS input files from a canonical H5AD."""

from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml


def decode(value: object) -> object:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_array(dataset: h5py.Dataset) -> np.ndarray:
    if h5py.check_dtype(vlen=dataset.dtype) is not None or dataset.dtype.kind == "O":
        return np.asarray(dataset.asstr()[:], dtype=object)
    values = dataset[:]
    if values.dtype.kind == "S":
        return np.asarray([value.decode("utf-8") for value in values], dtype=object)
    return values


def read_column(item: h5py.Dataset | h5py.Group) -> np.ndarray:
    if isinstance(item, h5py.Group) and item.attrs.get("encoding-type", "") == "categorical":
        categories = np.asarray([decode(value) for value in read_array(item["categories"])], dtype=object)
        codes = np.asarray(item["codes"][:], dtype=np.int64)
        return np.asarray([None if code < 0 else categories[code] for code in codes], dtype=object)
    return np.asarray([decode(value) for value in read_array(item)], dtype=object)


def read_index(group: h5py.Group) -> np.ndarray:
    index_name = str(decode(group.attrs["_index"]))
    return np.asarray([str(value) for value in read_array(group[index_name])], dtype=object)


def as_bool(values: np.ndarray) -> np.ndarray:
    return np.asarray([value is True or str(value).lower() == "true" for value in values], dtype=bool)


def write_counts(
    source: Path,
    destination: Path,
    selected_rows: np.ndarray,
    n_features: int,
    chunk_size: int,
) -> int:
    with h5py.File(source, "r") as handle:
        matrix = handle["layers/counts"]
        if matrix.attrs.get("encoding-type", "") != "csr_matrix":
            raise ValueError("layers/counts must be CSR")
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        nnz = int(np.sum(indptr[selected_rows + 1] - indptr[selected_rows]))
        with gzip.open(destination, "wt", encoding="utf-8", compresslevel=1) as stream:
            stream.write("%%MatrixMarket matrix coordinate integer general\n")
            stream.write("% M05 PS counts from canonical H5AD layers/counts.\n")
            stream.write(f"{n_features} {selected_rows.size} {nnz}\n")
            for output_start in range(0, selected_rows.size, chunk_size):
                rows = selected_rows[output_start : output_start + chunk_size]
                for local_index, source_row in enumerate(rows):
                    start = int(indptr[source_row])
                    end = int(indptr[source_row + 1])
                    indices = np.asarray(matrix["indices"][start:end], dtype=np.int64)
                    values = np.asarray(matrix["data"][start:end])
                    output_cell = output_start + local_index + 1
                    for feature, value in zip(indices, values):
                        stream.write(f"{feature + 1} {output_cell} {int(value)}\n")
    return nnz


def link_validated_counts(
    source_h5ad: Path,
    shared_dir: Path,
    destination: Path,
    expected_cells: int,
    allowed_structures: set[str],
) -> tuple[int, int]:
    manifest = json.loads((shared_dir / "adapter_manifest.json").read_text(encoding="utf-8"))
    if Path(manifest["source_h5ad"]).resolve() != source_h5ad.resolve():
        raise ValueError("Shared sparse counts were generated from a different H5AD")
    if int(manifest["output_shape"][1]) != expected_cells or manifest["selected_cells"] != expected_cells:
        raise ValueError("Shared sparse counts have a different selected cohort")
    shared_allowed = set(manifest["selection"]["eligible_assignment_structures"])
    if shared_allowed != allowed_structures:
        raise ValueError("Shared sparse counts use a different eligibility contract")
    source_counts = shared_dir / "counts.mtx.gz"
    if destination.exists():
        destination.unlink()
    os.link(source_counts, destination)
    return int(manifest["counts_nnz"]), int(manifest["output_shape"][0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    source = Path(config["source_h5ad"])
    output = Path(config["input_dir"])
    allowed = set(config["selection"]["eligible_assignment_structures"])
    guide_field = str(config.get("guide_field", "guide_id"))

    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        required = {"target_label", "is_ntc", "assignment_structure", "batch_id", guide_field}
        missing = required.difference(obs.keys())
        if missing:
            raise ValueError(f"Missing required obs columns: {sorted(missing)}")
        n_cells, n_features = map(int, handle["X"].attrs["shape"])
        cell_ids = read_index(obs)
        feature_ids = read_index(handle["var"])
        labels = read_column(obs["target_label"])
        structures = read_column(obs["assignment_structure"])
        ntc = as_bool(read_column(obs["is_ntc"]))
        guides = read_column(obs[guide_field])
        batches = read_column(obs["batch_id"])
        expected_ntc = labels == config["non_target_control"]
        if not np.array_equal(ntc, expected_ntc):
            raise ValueError("is_ntc is inconsistent with target_label")
        eligible = np.isin(structures, list(allowed)) & pd.notna(labels) & pd.notna(guides)
        selected = np.flatnonzero(eligible)
        if selected.size == 0:
            raise ValueError("No eligible cells")

    output.mkdir(parents=True, exist_ok=True)
    selected_cells = cell_ids[selected].astype(str)
    metadata = pd.DataFrame(
        {
            "target_label": labels[selected].astype(str),
            "gene": labels[selected].astype(str),
            "is_ntc": ntc[selected],
            "guide_id": guides[selected].astype(str),
            "batch_id": batches[selected].astype(str),
            "assignment_structure": structures[selected].astype(str),
        },
        index=selected_cells,
    )
    metadata.index.name = "cell_id"
    metadata.to_csv(output / "metadata.tsv.gz", sep="\t", compression="gzip")

    barcode = pd.DataFrame(
        {
            "cell": selected_cells,
            "barcode": guides[selected].astype(str),
            "gene": labels[selected].astype(str),
        }
    )
    barcode.to_csv(output / "barcode.tsv.gz", sep="\t", index=False, compression="gzip")
    shared_dir_value = config.get("shared_sparse_input_dir")
    if shared_dir_value:
        shared_dir = Path(shared_dir_value)
        nnz, output_features = link_validated_counts(
            source,
            shared_dir,
            output / "counts.mtx.gz",
            int(selected.size),
            allowed,
        )
        destination_features = output / "features.tsv.gz"
        if destination_features.exists():
            destination_features.unlink()
        os.link(shared_dir / "features.tsv.gz", destination_features)
        counts_materialization = "hardlink to validated shared canonical sparse counts"
    else:
        nnz = write_counts(source, output / "counts.mtx.gz", selected, n_features, int(config["row_chunk_size"]))
        output_features = n_features
        pd.DataFrame({"gene_id": feature_ids.astype(str)}).to_csv(
            output / "features.tsv.gz", sep="\t", index=False, header=False, compression="gzip"
        )
        counts_materialization = "direct export from canonical H5AD"

    manifest = {
        "module": "M05",
        "method": "Perturbation-response Score",
        "dataset": config["dataset"],
        "source_h5ad": str(source.resolve()),
        "source_shape": [n_cells, n_features],
        "output_shape": [output_features, int(selected.size)],
        "counts_nnz": nnz,
        "counts_source": "layers/counts",
        "counts_materialization": counts_materialization,
        "normalization": {"method": "LogNormalize", "scale_factor": 10000},
        "barcode_columns": ["cell", "barcode", "gene"],
        "barcode_semantics": "one target-level assignment per eligible cell",
        "guide_field": guide_field,
        "eligible_assignment_structures": sorted(allowed),
        "selected_cells": int(selected.size),
        "selected_ntc_cells": int(ntc[selected].sum()),
        "selected_target_cells": int((~ntc[selected]).sum()),
        "selected_targets": int(len(set(labels[selected][~ntc[selected]].astype(str)))),
        "cell_order": "canonical H5AD row order",
        "subset_rds_required": False,
    }
    (output / "adapter_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
