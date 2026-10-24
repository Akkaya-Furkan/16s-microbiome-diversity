# One-day scope — 2026-10-24

## Question and design
Do community compositions and within-sample diversity differ among body sites in the public QIIME 2 Moving Pictures teaching dataset?

Use the precomputed DADA2 ASV frequency table and sample metadata from the archived QIIME 2 2024.10 tutorial. No raw-read processing is planned. Analyze all matched samples descriptively; use only subject/visit blocks containing all four body sites for inference. Preserve repeated visits by permuting body-site labels consistently over each subject's entire trajectory, independently between subjects. This is a small, exploratory demonstration, not population-level evidence.

## Analysis fixed before implementation
- Download two source files with SHA-256 verification and retain provenance.
- Rarefy once without replacement to the smallest observed library size, fixed seed 26; report depth and excluded samples. Retain unrarefied alpha values as a depth-sensitivity comparison.
- Shannon entropy (natural log), observed ASVs, relative-abundance Bray–Curtis distance, classical PCoA; export all coordinates, eigenvalues, and negative-eigenvalue diagnostics.
- Body-site PERMANOVA on complete visits. Exhaustively enumerate (4!)^2 = 576 subject-restricted trajectory label permutations, including the identity. Use the fraction of statistics at least as large as observed, without Monte Carlo plus-one correction. Report pseudo-F, R-squared, and permutation count.
- Alpha group tests use between-site variation of subject-centered, visit-averaged alpha diversity. Use the same exact trajectory label space; adjust the two omnibus alpha p-values with Benjamini–Hochberg. Report site means and subject-level values. No feature-level differential-abundance claim.
- Report a similarly restricted dispersion test based on distances to group centroids in positive PCoA coordinates; document approximation from non-Euclidean distances.
- UniFrac, taxonomy, raw-read denoising, and independent-cohort validation: not done in this scope.

## Numerical acceptance criteria
1. Exactly 2 downloaded inputs with recorded URLs and SHA-256 hashes; table/metadata sample IDs align with no silent deletion.
2. At least 20 samples, 4 body sites, 2 subjects, and 5 complete subject/visit blocks available for inference; no zero-count samples after normalization.
3. Exactly 2 alpha metrics per sample, one symmetric Bray–Curtis matrix with zero diagonal and values in [0, 1], and at least 2 PCoA coordinate axes.
4. Exactly 576 permutations for each of 4 tests (PERMANOVA, dispersion, 2 alpha tests); exactly 2 BH-adjusted alpha p-values in [0, 1]. Statistical significance is not an acceptance criterion.
5. At least 10 exported machine-readable analysis tables/JSON files and exactly 3 nonempty PNG figures; a machine-generated acceptance report verifies these requirements.
6. At least 12 passing pytest tests covering numerical reference cases, normalization, invalid inputs, permutation restrictions, and end-to-end artifact consistency; CI runs pytest and reproduces outputs offline from committed source data.
7. English README includes Question, Approach, Results, Limitations, source attribution, reproducible commands, and actual computed results. Environment file, download script, and GitHub Actions workflow exist.
8. Commit and push main; compare local HEAD to `git ls-remote origin main` and read the resulting GitHub Actions conclusion.
