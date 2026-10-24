"""Numerical methods; all count matrices have samples in rows."""
from itertools import permutations, product
import numpy as np
from scipy.spatial.distance import pdist, squareform


def validate_counts(counts):
    x = np.asarray(counts)
    if x.ndim != 2 or min(x.shape) == 0:
        raise ValueError("Expected a nonempty sample-by-feature matrix")
    if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
        raise ValueError("Counts must be finite nonnegative integers")
    if (x.sum(axis=1) == 0).any():
        raise ValueError("Empty samples are not supported")
    return x.astype(np.int64)


def rarefy(counts, depth, seed=26):
    x = validate_counts(counts)
    if not isinstance(depth, (int, np.integer)) or depth < 1 or depth > x.sum(axis=1).min():
        raise ValueError("Depth must be a positive integer no larger than every library")
    rng = np.random.default_rng(seed)
    return np.array([rng.multivariate_hypergeometric(row, depth) for row in x])


def alpha(counts):
    x = validate_counts(counts)
    p = x / x.sum(axis=1, keepdims=True)
    logp = np.zeros_like(p)
    np.log(p, out=logp, where=p > 0)
    return np.column_stack((-(p * logp).sum(axis=1), (x > 0).sum(axis=1)))


def bray_curtis(counts):
    x = validate_counts(counts)
    return squareform(pdist(x / x.sum(axis=1, keepdims=True), metric="braycurtis"))


def pcoa(distance):
    d = np.asarray(distance, dtype=float)
    if d.ndim != 2 or d.shape[0] != d.shape[1] or not np.isfinite(d).all():
        raise ValueError("Expected finite square distances")
    if (d < 0).any() or not np.allclose(d, d.T) or not np.allclose(np.diag(d), 0):
        raise ValueError("Invalid distance matrix")
    h = np.eye(len(d)) - np.ones_like(d) / len(d)
    eig, vec = np.linalg.eigh(-0.5 * h @ (d * d) @ h)
    order = np.argsort(eig)[::-1]
    eig, vec = eig[order], vec[:, order]
    positive = eig > 1e-12
    coordinates = vec[:, positive] * np.sqrt(eig[positive])
    # Fix arbitrary signs to stabilize exported coordinates across LAPACK versions.
    for j in range(coordinates.shape[1]):
        anchor = np.argmax(np.abs(coordinates[:, j]))
        if coordinates[anchor, j] < 0:
            coordinates[:, j] *= -1
    return coordinates, eig


def permanova_stat(distance, labels):
    """Anderson one-factor pseudo-F and R² from squared distances."""
    d = np.asarray(distance)
    labels = np.asarray(labels)
    groups = np.unique(labels)
    n, k = len(labels), len(groups)
    if d.shape != (n, n) or not 1 < k < n:
        raise ValueError("Need aligned distances and between 2 and n-1 groups")
    total = np.sum(d ** 2) / (2 * n)
    within = sum(np.sum(d[np.ix_(labels == g, labels == g)] ** 2) /
                 (2 * np.sum(labels == g)) for g in groups)
    between = total - within
    if total <= 0 or within <= 0:
        raise ValueError("Degenerate distance variation")
    return (between / (k - 1)) / (within / (n - k)), between / total


def trajectory_permutations(labels, subjects):
    """Each subject receives one label bijection applied to all of their visits."""
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    if labels.shape != subjects.shape or labels.ndim != 1:
        raise ValueError("Labels and subjects must align")
    sites, people = np.unique(labels), np.unique(subjects)
    for person in people:
        if set(labels[subjects == person]) != set(sites):
            raise ValueError("Every subject must contain every site")
    for maps in product(list(permutations(sites)), repeat=len(people)):
        result = labels.copy()
        for person, mapping in zip(people, maps):
            lookup = dict(zip(sites, mapping))
            mask = subjects == person
            result[mask] = [lookup[label] for label in labels[mask]]
        yield result


def alpha_site_stat(values, labels, subjects):
    """Between-site sum of squares after centering each subject's site means."""
    values, labels, subjects = np.asarray(values), np.asarray(labels), np.asarray(subjects)
    means = np.array([[values[(subjects == s) & (labels == g)].mean()
                       for g in np.unique(labels)] for s in np.unique(subjects)])
    centered = means - means.mean(axis=1, keepdims=True)
    return np.sum(centered.mean(axis=0) ** 2)


def dispersion_stat(coords, labels):
    """One-way F of distances to positive-PCoA group centroids (approximate)."""
    labels = np.asarray(labels)
    z = np.empty(len(labels))
    for group in np.unique(labels):
        mask = labels == group
        z[mask] = np.linalg.norm(coords[mask] - coords[mask].mean(axis=0), axis=1)
    grand = z.mean()
    between = sum(np.sum(labels == g) * (z[labels == g].mean() - grand) ** 2
                  for g in np.unique(labels))
    within = sum(np.sum((z[labels == g] - z[labels == g].mean()) ** 2)
                 for g in np.unique(labels))
    return (between / (len(np.unique(labels)) - 1)) / (within / (len(labels) - len(np.unique(labels))))


def exact_pvalue(observed, null):
    null = np.asarray(null, dtype=float)
    if not len(null) or not np.isfinite(null).all() or not np.isfinite(observed):
        raise ValueError("Finite statistics are required")
    return float(np.mean(null >= observed - 1e-10 * max(1, abs(observed))))


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not len(p) or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("Invalid p-values")
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(q)
    result[order] = np.minimum(q, 1)
    return result
