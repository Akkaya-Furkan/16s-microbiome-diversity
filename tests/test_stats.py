from itertools import permutations, product
import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway
from microdiv.stats import (alpha, rarefy, bray_curtis, pcoa, permanova_stat,
                           trajectory_permutations, alpha_site_stat, bh_adjust,
                           exact_pvalue, validate_counts)


def test_alpha_reference():
    np.testing.assert_allclose(alpha([[1, 1, 0], [0, 0, 9]]), [[np.log(2), 2], [0, 1]])


@pytest.mark.parametrize("counts", [[[0, 0]], [[-1, 2]], [[.5, 2]], [[np.nan, 1]], [], [1, 2]])
def test_reject_invalid_counts(counts):
    with pytest.raises(ValueError):
        validate_counts(counts)


def test_rarefaction_conserves_counts_and_seed():
    x = np.array([[30, 20, 50], [20, 80, 0]])
    y = rarefy(x, 40)
    np.testing.assert_array_equal(y.sum(axis=1), [40, 40])
    assert (y <= x).all() and (y >= 0).all()
    np.testing.assert_array_equal(y, rarefy(x, 40))


@pytest.mark.parametrize("depth", [0, -1, 101, 1.5])
def test_invalid_depth(depth):
    with pytest.raises(ValueError):
        rarefy([[30, 70]], depth)


def test_bray_reference_and_depth_invariance():
    x = [[1, 0], [0, 10], [5, 5], [2, 0]]
    d = bray_curtis(x)
    assert d[0, 1] == 1 and d[0, 2] == .5 and d[0, 3] == 0
    np.testing.assert_allclose(d, d.T)
    np.testing.assert_array_equal(np.diag(d), 0)


def test_pcoa_recovers_euclidean_distances():
    points = np.array([[0, 0], [1, 0], [0, 2], [2, 3]])
    d = squareform(pdist(points))
    coords, eig = pcoa(d)
    np.testing.assert_allclose(squareform(pdist(coords)), d, atol=1e-10)
    assert sum(eig > 1e-12) == 2


def test_pcoa_retains_negative_eigenvalue_diagnostics():
    _, eig = pcoa([[0, 1, 1], [1, 0, 3], [1, 3, 0]])
    assert min(eig) < 0


def test_pcoa_rejects_asymmetric_matrix():
    with pytest.raises(ValueError):
        pcoa([[0, 1], [2, 0]])


def test_permanova_matches_univariate_anova():
    values = np.array([0, 1, 2, 6, 7, 9], dtype=float)
    f, r2 = permanova_stat(squareform(pdist(values[:, None])), [0, 0, 0, 1, 1, 1])
    assert f == pytest.approx(f_oneway(values[:3], values[3:]).statistic)
    assert r2 == pytest.approx(f / (f + 4))


def test_trajectory_permutations_preserve_repeated_visits():
    labels = np.tile(["a", "b", "c", "d"], 4)
    subjects = np.repeat(["s1", "s2"], 8)
    permuted = list(trajectory_permutations(labels, subjects))
    assert len(permuted) == 576
    assert len({tuple(p) for p in permuted}) == 576
    for p in permuted:
        for start in [0, 8]:
            np.testing.assert_array_equal(p[start:start+4], p[start+4:start+8])
    assert any(np.array_equal(labels, p) for p in permuted)


def test_exact_pvalue_includes_ties():
    assert exact_pvalue(2, [1, 2, 3, 2]) == .75


def test_bh_reference_and_order():
    np.testing.assert_allclose(bh_adjust([.04, .001, .03, .9]), [.053333333333, .004, .053333333333, .9])


@pytest.mark.parametrize("p", [[-.1], [1.1], [np.nan], []])
def test_bh_rejects_invalid_values(p):
    with pytest.raises(ValueError):
        bh_adjust(p)


def test_alpha_stat_invariant_to_subject_baseline():
    labels = np.tile(["a", "b", "c", "d"], 2)
    subjects = np.repeat(["s1", "s2"], 4)
    values = np.tile([1., 2., 3., 4.], 2)
    assert alpha_site_stat(values, labels, subjects) == pytest.approx(alpha_site_stat(values + np.repeat([0, 100], 4), labels, subjects))
