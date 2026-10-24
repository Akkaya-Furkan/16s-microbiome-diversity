"""Read only the BIOM payload from a QIIME artifact, without installing QIIME."""
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from .stats import validate_counts


def load_inputs(root):
    root = Path(root)
    sources = json.loads((root / "data/sources.json").read_text())
    for name, source in sources.items():
        if hashlib.sha256((root / "data/raw" / name).read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Input hash mismatch: {name}")
    with ZipFile(root / "data/raw/table-dada2.qza") as archive:
        names = [n for n in archive.namelist() if n.endswith("/data/feature-table.biom")]
        if len(names) != 1:
            raise ValueError("Expected exactly one BIOM payload")
        with h5py.File(io.BytesIO(archive.read(names[0])), "r") as biom:
            samples = biom["sample/ids"].asstr()[:]
            features = biom["observation/ids"].asstr()[:]
            matrix = biom["observation/matrix"]
            counts = csr_matrix((matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                                shape=(len(features), len(samples))).toarray().T
    counts = validate_counts(counts)
    metadata = pd.read_csv(root / "data/raw/sample-metadata.tsv", sep="\t", comment="#",
                           dtype={"sample-id": str}).set_index("sample-id")
    if not metadata.index.is_unique or len(set(samples)) != len(samples):
        raise ValueError("Duplicate sample identifiers")
    if set(samples) != set(metadata.index):
        raise ValueError("Table and metadata IDs must match exactly; no silent deletion")
    metadata = metadata.loc[samples]
    if metadata[["subject", "body-site", "days-since-experiment-start"]].isna().any().any():
        raise ValueError("Missing design metadata")
    return counts, metadata, features


def complete_visit_mask(metadata):
    sites = set(metadata["body-site"])
    keys = ["subject", "days-since-experiment-start"]
    valid = metadata.groupby(keys)["body-site"].agg(lambda x: len(x) == len(sites) and set(x) == sites)
    return np.array([bool(valid.loc[tuple(row)]) for row in metadata[keys].to_numpy()])
