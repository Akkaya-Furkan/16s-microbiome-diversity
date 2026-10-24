"""Offline reproducible analysis and artifact-level acceptance checks."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .io import load_inputs, complete_visit_mask
from .stats import (alpha, rarefy, bray_curtis, pcoa, permanova_stat,
                    trajectory_permutations, alpha_site_stat, dispersion_stat,
                    exact_pvalue, bh_adjust)

METRICS = ["shannon", "observed_asvs"]
COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]


def save_table(frame, path):
    frame.to_csv(path, sep="\t", float_format="%.12g")


def run(root, output):
    root, output = Path(root), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    x, metadata, features = load_inputs(root)
    depth = int(x.sum(axis=1).min())
    normalized = rarefy(x, depth, seed=26)
    a = pd.DataFrame(alpha(normalized), index=metadata.index, columns=METRICS)
    raw_a = pd.DataFrame(alpha(x), index=metadata.index, columns=METRICS)
    mask = complete_visit_mask(metadata)
    design = metadata.loc[mask]
    labels, subjects = design["body-site"].to_numpy(), design["subject"].to_numpy()
    sites = sorted(metadata["body-site"].unique())
    distance = bray_curtis(normalized)
    coordinates, eigenvalues = pcoa(distance)
    infer_distance = distance[np.ix_(mask, mask)]
    infer_coords, _ = pcoa(infer_distance)
    f, r2 = permanova_stat(infer_distance, labels)
    observed = [f, dispersion_stat(infer_coords, labels)] + [
        alpha_site_stat(a.loc[mask, metric], labels, subjects) for metric in METRICS]
    null = []
    for permuted in trajectory_permutations(labels, subjects):
        null.append([permanova_stat(infer_distance, permuted)[0],
                     dispersion_stat(infer_coords, permuted)] + [
                         alpha_site_stat(a.loc[mask, metric], permuted, subjects) for metric in METRICS])
    null = np.array(null)
    pvalues = [exact_pvalue(obs, null[:, i]) for i, obs in enumerate(observed)]
    names = ["permanova", "dispersion_positive_pcoa", "alpha_shannon", "alpha_observed_asvs"]
    tests = pd.DataFrame({"statistic": observed, "p_exact": pvalues,
                          "permutations_including_identity": len(null)}, index=names)
    tests["r_squared"] = [r2, np.nan, np.nan, np.nan]
    tests["q_bh_alpha_family"] = [np.nan, np.nan, *bh_adjust(pvalues[2:])]
    tests.index.name = "test"
    audit = metadata[["body-site", "subject", "days-since-experiment-start"]].copy()
    audit["raw_reads"] = x.sum(axis=1)
    audit["normalized_reads"] = normalized.sum(axis=1)
    audit["complete_visit_inference"] = mask
    save_table(a.join(metadata[["body-site", "subject"]]), output / "alpha_diversity.tsv")
    save_table(raw_a, output / "alpha_unrarefied.tsv")
    save_table(a.join(metadata[["body-site"]]).groupby("body-site")[METRICS].agg(["count", "mean", "median", "std"]), output / "alpha_site_summary.tsv")
    subject_alpha = a.loc[mask].join(design[["subject", "body-site"]]).groupby(["subject", "body-site"])[METRICS].mean()
    save_table(subject_alpha, output / "alpha_subject_site_means.tsv")
    save_table(pd.DataFrame(distance, index=metadata.index, columns=metadata.index), output / "bray_curtis.tsv")
    axes = [f"PC{i+1}" for i in range(coordinates.shape[1])]
    save_table(pd.DataFrame(coordinates, index=metadata.index, columns=axes).join(metadata[["body-site", "subject"]]), output / "pcoa_coordinates.tsv")
    positive_total = eigenvalues[eigenvalues > 1e-12].sum()
    save_table(pd.DataFrame({"eigenvalue": eigenvalues, "fraction_positive_inertia": np.where(eigenvalues > 1e-12, eigenvalues / positive_total, 0)}), output / "pcoa_eigenvalues.tsv")
    save_table(tests, output / "statistical_tests.tsv")
    save_table(pd.DataFrame(null, columns=names).rename_axis("permutation"), output / "permutation_null.tsv")
    save_table(audit, output / "sample_audit.tsv")
    save_table(pd.DataFrame(normalized, index=metadata.index, columns=features), output / "rarefied_counts.tsv")
    sensitivity = pd.DataFrame({"rarefied_mean": a.mean(), "unrarefied_mean": raw_a.mean(),
                                "spearman_rank_correlation": [a[c].corr(raw_a[c], method="spearman") for c in METRICS]})
    save_table(sensitivity, output / "alpha_depth_sensitivity.tsv")
    blocks = design.groupby(["subject", "days-since-experiment-start"]).ngroups
    summary = {
        "samples": len(metadata), "asvs": len(features), "subjects": metadata.subject.nunique(),
        "body_sites": len(sites), "site_sample_counts": metadata["body-site"].value_counts().sort_index().to_dict(),
        "raw_reads": int(x.sum()), "raw_library_min": depth, "raw_library_max": int(x.sum(axis=1).max()),
        "rarefaction_depth": depth, "seed": 26, "samples_dropped_normalization": 0,
        "complete_visits": blocks, "inference_samples": int(mask.sum()),
        "excluded_from_inference": metadata.index[~mask].tolist(),
        "permutations_per_test": len(null), "permanova_f": f, "permanova_r2": r2,
        "permanova_p": pvalues[0], "dispersion_p": pvalues[1],
        "alpha_p": dict(zip(METRICS, pvalues[2:])), "alpha_q": dict(zip(METRICS, bh_adjust(pvalues[2:]))),
        "negative_eigenvalues": int((eigenvalues < -1e-12).sum()),
        "negative_inertia_fraction_absolute": float(-eigenvalues[eigenvalues < -1e-12].sum() / np.abs(eigenvalues).sum()),
        "pcoa_axes": coordinates.shape[1], "pc1_positive_inertia_fraction": float(eigenvalues[0] / positive_total),
        "pc2_positive_inertia_fraction": float(eigenvalues[1] / positive_total),
        "versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "matplotlib", "h5py"]},
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures(output, a, audit, sites, coordinates, eigenvalues / positive_total)
    files = sorted(p for p in output.iterdir() if p.suffix in {".tsv", ".json", ".png"} and p.name not in {"acceptance.json", "manifest.json"})
    manifest = {p.name: {"bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    machine_count = len([p for p in files if p.suffix in {".tsv", ".json"}])
    png_count = len([p for p in files if p.suffix == ".png"])
    checks = {
        "two_hashed_sources": len(json.loads((root / "data/sources.json").read_text())) == 2,
        "design_minimum": len(metadata) >= 20 and len(sites) == 4 and summary["subjects"] == 2 and blocks >= 5,
        "normalization": bool((normalized.sum(axis=1) == depth).all()),
        "two_alpha_metrics": a.shape == (len(metadata), 2) and bool(np.isfinite(a).all().all()),
        "distance_valid": bool(np.allclose(distance, distance.T) and np.allclose(np.diag(distance), 0) and (distance >= 0).all() and (distance <= 1).all()),
        "pcoa_axes": coordinates.shape[1] >= 2,
        "four_exact_tests": null.shape == (576, 4),
        "two_bh_values": len(summary["alpha_q"]) == 2 and bool(((tests.q_bh_alpha_family.dropna() >= 0) & (tests.q_bh_alpha_family.dropna() <= 1)).all()),
        "machine_readable_outputs": machine_count >= 10,
        "three_figures": png_count == 3 and all(p.stat().st_size > 0 for p in files),
    }
    acceptance = {"checks": checks, "passed": sum(checks.values()), "expected": len(checks),
                  "machine_readable_files": machine_count, "png_figures": png_count}
    (output / "acceptance.json").write_text(json.dumps(acceptance, indent=2) + "\n")
    if not all(checks.values()):
        raise RuntimeError(f"Acceptance checks failed: {acceptance}")
    print(json.dumps({"summary": summary, "acceptance": acceptance}, indent=2))
    return summary


def figures(output, a, audit, sites, coords, fractions):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for ax, metric, title in zip(axes, METRICS, ["Shannon diversity (natural log)", "Observed ASVs"]):
        for i, (site, color) in enumerate(zip(sites, COLORS)):
            values = a.loc[audit["body-site"] == site, metric].to_numpy()
            ax.scatter(i + np.linspace(-.12, .12, len(values)), values, color=color, alpha=.8, s=28)
            ax.plot([i-.2, i+.2], [np.median(values)] * 2, color="black", lw=2)
        ax.set(xticks=range(len(sites)), xticklabels=sites, title=title, xlabel="Body site")
    fig.suptitle("Within-sample diversity · fixed-depth rarefaction\nPoints are repeated samples; bars show medians")
    fig.savefig(output / "alpha_diversity.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    for site, color in zip(sites, COLORS):
        for subject, marker in zip(sorted(audit.subject.unique()), ["o", "^"]):
            mask = ((audit["body-site"] == site) & (audit.subject == subject)).to_numpy()
            ax.scatter(coords[mask, 0], coords[mask, 1], c=color, marker=marker, s=55,
                       label=f"{site} / {subject}", alpha=.85, edgecolors="white", linewidths=.4)
    ax.set(xlabel=f"PC1 ({fractions[0]:.1%} of positive inertia)", ylabel=f"PC2 ({fractions[1]:.1%} of positive inertia)",
           title="Bray–Curtis PCoA · all samples\nColor: body site · shape: subject")
    ax.legend(fontsize=8, loc="best")
    fig.savefig(output / "pcoa_bray_curtis.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    for i, (site, color) in enumerate(zip(sites, COLORS)):
        vals = audit.loc[audit["body-site"] == site, "raw_reads"].to_numpy()
        ax.scatter(i + np.linspace(-.12, .12, len(vals)), vals, c=color, s=30)
    ax.axhline(audit.normalized_reads.iloc[0], c="black", linestyle="--", label="Rarefaction depth")
    ax.set(xticks=range(len(sites)), xticklabels=sites, ylabel="Reads (log scale)", yscale="log", title="Library sizes before normalization")
    ax.legend()
    fig.savefig(output / "library_sizes.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.root, args.output)


if __name__ == "__main__":
    main()
