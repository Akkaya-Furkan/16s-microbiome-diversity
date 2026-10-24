# Statistical methods and interpretation

## Input and normalization

The QIIME 2 artifact's embedded BIOM HDF5 payload is read as a sample-by-ASV matrix. Counts must be finite, nonnegative integers and every library must be nonempty. Duplicate identifiers, unequal sample ID sets, and missing design metadata are errors. The input loader verifies both pinned SHA-256 hashes on every run. The original artifact, including its upstream DADA2 provenance, is retained.

The rarefaction depth is the minimum library size computed from the input, rather than selected after looking at p-values. A multivariate hypergeometric draw samples reads without replacement, using NumPy's generator with the recorded seed. No samples are discarded for sequencing depth. Unrarefied alpha values and rank correlations provide a descriptive depth-sensitivity comparison; no claim of robustness to alternative rarefaction draws is made.

## Alpha diversity and its group tests

Shannon entropy is `-sum(p * ln(p))` over positive ASV frequencies. Observed richness counts ASVs with nonzero reads; it is not an estimator of unobserved richness.

For inference, retain subject/visit blocks with exactly one sample at each of the observed body sites. Compute each subject's mean alpha diversity per site over their complete visits, subtract that subject's mean across sites, then calculate the sum of squared site means across subjects. Larger values indicate stronger site differences. This omnibus statistic gives each subject equal weight and is invariant to a subject-wide diversity offset. It is calibrated by the exact restricted permutations below. Benjamini–Hochberg correction applies to the prespecified family of two alpha outcomes only. PERMANOVA is a separate primary composition test; dispersion is a diagnostic. No pairwise site or taxon tests are performed.

## Bray–Curtis and ordination

Counts are divided by each library total. Bray–Curtis is `sum(abs(p_i - p_j)) / sum(p_i + p_j)`. Classical PCoA eigendecomposes the double-centered squared-distance matrix `-H D² H / 2`. All eigenvalues are exported; plotted axes and percentage denominators use positive inertia only. Negative eigenvalues are reported rather than silently described as explained variance. Axis signs are fixed by making the largest absolute coordinate positive.

## PERMANOVA and restricted permutation space

For n samples in k sites, total sums of squares are `sum(D²)/(2n)`. Within-group sums of squares are the sum over sites of `sum(D_group²)/(2n_group)`. The pseudo-F statistic is `(SS_between/(k-1))/(SS_within/(n-k))`; R-squared is `SS_between/SS_total`. These nominal degrees of freedom define the statistic, not an independent-sample F reference distribution.

Every subject receives a bijection of body-site labels, applied consistently to that subject's entire trajectory. Cartesian products of subject-specific bijections enumerate the full `(k!)^subjects` reference space. Time trajectories remain intact, and complete visits make the sampling schedule equal across sites within each subject. The observed labeling is included. Exact p-values are the fraction of all enumerated statistics at least as large as observed, allowing a numerical tie tolerance. There is no Monte Carlo sampling and no plus-one correction. Global relabelings yield duplicate statistics by design and must remain in this uniform labeling space.

This calibration assumes that whole site trajectories would be exchangeable under the null within a person. Body sites were not randomly assigned. Conditioning on subject and preserving time structure avoids naive sample-level shuffling, but does not establish population representativeness or causal identification. Only the actual number of subjects provides independent biological replication.

## Dispersion diagnostic

For each candidate labeling, recompute group centroids in the inference subset's positive PCoA coordinates, calculate distances to those centroids, and calculate a one-way F statistic of those distances. Calibrate with the same restricted labelings. This is an approximate centroid-based dispersion diagnostic because negative PCoA axes are omitted. It is not a full signed-axis PERMDISP implementation. A nonsignificant result does not establish equal dispersion or prove that a PERMANOVA effect is purely a location effect.

## Reproducibility and verification

The committed source inputs make the analysis and CI independent of source-server availability. The separate downloader checks all payload hashes before writing any new file and rejects corrupt cached files. Tests include a failed-download mutation check. Scientific numerical tests include known entropy values, exact Bray–Curtis examples, Euclidean distance reconstruction by PCoA, equivalence of Euclidean one-factor PERMANOVA to univariate ANOVA, BH reference values, and whole-trajectory permutation invariants. End-to-end tests compare regenerated statistics, counts, and distances to committed results with numerical tolerances. PNG byte identity across platforms is not required; the committed manifest verifies the delivered files.
