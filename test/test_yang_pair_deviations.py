"""Tests of analyses/yang_pair_deviations.py on small made-up matrices (no shared directory needed)."""
from itertools import combinations

import numpy as np
import pytest

from analyses import yang_pair_deviations as Y

ODORS = ['a', 'b', 'c', 'd', 'e', 'f']
PAIRS = list(combinations(ODORS, 2))


def test_monotone_transform_gives_zero_residuals():
    rng = np.random.default_rng(0)
    T = rng.uniform(0.0, 0.8, size=(5, len(PAIRS)))
    st = T.mean(0)
    sc = 0.2 + 0.5 * st ** 2   # an increasing function of the TNT similarities
    assert np.allclose(Y.iso_resid(sc, st), 0, atol=1e-12)
    # the same through the per-pair table: control recordings a monotone transform of the TNT mean, any noise per
    # recording averaging out to it
    C = np.tile(sc, (4, 1)) + np.array([[0.01], [-0.01], [0.02], [-0.02]])
    block = np.array([p == ('a', 'b') for p in PAIRS])
    tab, _ = Y.pair_table(PAIRS, C, T, block, n=50, seed=1)
    assert np.allclose(tab.iso_resid, 0, atol=1e-12)
    assert np.allclose(tab.apl_effect, tab.sim_tnt - tab.sim_ctrl)
    assert np.allclose(tab.d_rank, tab.rank_ctrl - tab.rank_tnt) and tab.rank_ctrl.min() == 1


def test_per_odor_effect_recovered():
    X = Y.odor_design(PAIRS, ODORS)
    assert X.shape == (15, 6) and (X.sum(1) == 2).all() and (X.sum(0) == 5).all()
    effect = np.array([-0.05, 0.02, 0.0, 0.03, -0.01, 0.01])
    res = X @ effect
    a, e = Y.odor_fit(res, X)
    assert np.allclose(a, effect) and np.allclose(e, 0)
    # fitted outside a set of pairs, the effects still come back, and what is left on those pairs is theirs alone
    special = np.array([p in {('a', 'b'), ('a', 'c'), ('b', 'd')} for p in PAIRS])
    extra = np.where(special, -0.04, 0.0)
    a3, rem = Y.odor_fit(res + extra, X, ~special)
    assert np.allclose(a3, effect) and np.allclose(rem, extra)
    # the mean residual over each odor's pairs
    M = X / X.sum(0)
    assert np.allclose(res @ M, [(res[[k for k, p in enumerate(PAIRS) if o in p]]).mean() for o in ODORS])


def test_bh_against_hand_computation():
    # sorted 0.01, 0.03, 0.04, 0.2 -> 0.04, 0.06, 0.0533, 0.2 -> running minimum from the top: 0.04, 0.0533, 0.0533, 0.2
    q = Y.bh([0.01, 0.04, 0.03, 0.2])
    assert np.allclose(q, [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.2])
    assert np.allclose(Y.bh([0.5, 0.9]), [0.9, 0.9]) and Y.bh([1.0, 1.0]).max() == 1.0


def test_two_sided_p():
    assert np.allclose(Y.two_sided_p(np.ones((9, 1))), [0.2])
    X = np.array([[-1.0], [1.0], [2.0], [3.0]])   # 1 of 4 at or below 0: tail (1 + 1) / 5
    assert np.allclose(Y.two_sided_p(X), [0.8])
    assert Y.two_sided_p(np.zeros((3, 1)))[0] == 1.0


def test_side_counts():
    below = np.array([3, 8, 0, 11])
    full = np.array([-0.1, 0.2, 0.0, -0.3])
    assert list(Y.side_counts(below, 11, full)) == [3, 3, 11, 11]


def test_single_recordings_and_resamples():
    rng = np.random.default_rng(3)
    T = rng.uniform(0.1, 0.7, size=(4, len(PAIRS)))
    C = T * 0.5 + rng.normal(0, 0.01, size=T.shape)
    n_c, n_t = Y.single_recordings(C, T)
    assert n_c.shape == (len(PAIRS),) and n_c.max() <= len(C) and n_t.max() <= len(T)
    # the draws: control first, then TNT, from one generator
    (c, t), = list(Y.resamples(C, T, 1, 7))
    rng = np.random.default_rng(7)
    ic, it = rng.integers(0, 4, 4), rng.integers(0, 4, 4)
    assert np.allclose(c, C[ic].mean(0)) and np.allclose(t, T[it].mean(0))


def test_odor_tables_shapes(monkeypatch):
    monkeypatch.setattr(Y, 'RESAMPLES', dict(pairs=20, odor_effects=20, outside_block=20))
    rng = np.random.default_rng(5)
    T = rng.uniform(0.1, 0.7, size=(4, len(PAIRS)))
    C = T * 0.5 + rng.normal(0, 0.02, size=T.shape)
    block = np.array([p in {('a', 'd'), ('b', 'e')} for p in PAIRS])
    odors, block_pairs, rows = Y.odor_tables(PAIRS, C, T, block, odors=ODORS)
    assert list(odors.odor) == ODORS and len(block_pairs) == 2
    assert np.allclose(block_pairs.iso_resid, block_pairs.predicted_from_odors + block_pairs.pair_specific)
    assert (odors.mean_lo <= odors.mean_hi).all() and set(rows) == {'block_predicted_from_odors', 'block_pair_specific',
        'odor_effect_variance_share'}
    assert 0 <= rows['odor_effect_variance_share'][0] <= 1


if __name__ == '__main__':
    pytest.main([__file__])
