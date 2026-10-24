import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from microdiv.io import load_inputs, complete_visit_mask
from microdiv.pipeline import run

ROOT = Path(__file__).resolve().parents[1]


def test_input_alignment_and_complete_visits():
    counts, metadata, features = load_inputs(ROOT)
    assert counts.shape == (len(metadata), len(features))
    selected = metadata.loc[complete_visit_mask(metadata)]
    groups = selected.groupby(["subject", "days-since-experiment-start"])
    assert (groups.size() == 4).all()
    assert (groups["body-site"].nunique() == 4).all()


def test_offline_end_to_end_and_committed_results(tmp_path):
    summary = run(ROOT, tmp_path)
    report = json.loads((tmp_path / "acceptance.json").read_text())
    assert report["passed"] == report["expected"] == 10
    assert summary == json.loads((ROOT / "results/summary.json").read_text())
    for name in ["alpha_diversity.tsv", "statistical_tests.tsv", "bray_curtis.tsv", "permutation_null.tsv", "rarefied_counts.tsv"]:
        actual = pd.read_csv(tmp_path / name, sep="\t", index_col=0)
        expected = pd.read_csv(ROOT / "results" / name, sep="\t", index_col=0)
        pd.testing.assert_frame_equal(actual, expected, check_exact=False, atol=1e-9, rtol=1e-9)
    tests = pd.read_csv(tmp_path / "statistical_tests.tsv", sep="\t", index_col=0)
    assert (tests.permutations_including_identity == 576).all()
    assert tests.q_bh_alpha_family.notna().sum() == 2
    null = pd.read_csv(tmp_path / "permutation_null.tsv", sep="\t", index_col=0)
    for name in tests.index:
        expected = np.mean(null[name] >= tests.loc[name, "statistic"] - 1e-8)
        assert tests.loc[name, "p_exact"] == pytest.approx(expected)


def test_committed_manifest():
    manifest = json.loads((ROOT / "results/manifest.json").read_text())
    for name, details in manifest.items():
        content = (ROOT / "results" / name).read_bytes()
        assert len(content) == details["bytes"]
        assert hashlib.sha256(content).hexdigest() == details["sha256"]


def test_download_rejects_corruption_without_partial_mutation(tmp_path):
    spec = importlib.util.spec_from_file_location("download_data", ROOT / "scripts/download_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    bad = tmp_path / "sample-metadata.tsv"
    bad.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="SHA-256"):
        module.download(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == [bad.name]
    assert bad.read_bytes() == b"corrupt"
