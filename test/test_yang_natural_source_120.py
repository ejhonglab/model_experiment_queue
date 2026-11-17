"""Tests of analyses/yang_natural_source_120.py on small made-up matrices (no shared directory needed)."""
from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from analyses import yang_natural_source_120 as Y

ODORS = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
PAIRS = list(combinations(ODORS, 2))


def made_up_matrix(n, seed=0):
    rng = np.random.default_rng(seed)
    pts = rng.normal(size=(n, 3))
    return np.linalg.norm(pts[:, None] - pts[None], axis=-1)


def test_resid_uncorrelated_with_x():
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 1, 50)
    y = 0.3 + 0.8 * x + rng.normal(0, 0.1, 50)
    r = Y.resid(y, x)
    assert abs(np.corrcoef(r, x)[0, 1]) < 1e-10 and abs(r.mean()) < 1e-12


def test_linear_kc_gives_zero_residual_and_undefined_rho():
    D = made_up_matrix(len(ODORS))
    d4 = Y.upper(D)
    orn = Y.upper(made_up_matrix(len(ODORS), seed=2))
    is66 = np.array([a in ODORS[:5] and b in ODORS[:5] for a, b in PAIRS])
    stats = Y.make_stats(d4, orn, is66, orn[is66], d4[is66])
    s = stats(0.2 + 0.5 * orn)
    assert np.isnan(s['res_d4_120']) and np.isnan(s['res_d4_66']) and np.isnan(s['res_d3_66'])
    assert s['kc_orn_120'] == pytest.approx(1.0) and s['kc_d3_66'] == pytest.approx(Y.rho(orn[is66], d4[is66]))
    assert np.isnan(Y.rho(np.zeros(5), np.arange(5))) and Y.rho(np.arange(5), np.arange(5) ** 2) == pytest.approx(1.0)


def test_permutation_of_a_matrix_against_itself():
    D = made_up_matrix(8, seed=3)
    perms = Y.permutations(8, 500, 7)
    obs, p, sd, q95 = Y.perm_test(Y.upper(D), D, perms)
    assert obs == pytest.approx(1.0) and p == pytest.approx(1 / 501) and sd > 0 and q95 < 1
    # the fast null is Spearman's ρ of the permuted matrix
    fixed = np.random.default_rng(4).normal(size=28)
    from scipy.stats import spearmanr
    iu = np.triu_indices(8, 1)
    null = [spearmanr(fixed, D[np.ix_(q, q)][iu])[0] for q in perms[:50]]
    obs, p, sd, q95 = Y.perm_test(fixed, D, perms[:50])
    assert p == pytest.approx((np.sum(np.array(null) >= obs - 1e-12) + 1) / 51)
    assert sd == pytest.approx(np.std(null)) and q95 == pytest.approx(np.percentile(null, 95))


def test_lookup_either_orientation():
    s = pd.Series([1.0, 2.0], index=pd.MultiIndex.from_tuples([('a', 'b'), ('c', 'a')]))
    assert list(Y.lookup(s, [('a', 'b'), ('b', 'a'), ('a', 'c'), ('b', 'c')])[:3]) == [1.0, 1.0, 2.0]
    assert np.isnan(Y.lookup(s, [('b', 'c')])[0])


def test_odor_keys_and_matrix():
    cols = ['1S/X/c1', '1S/Y/c2/h1', 'other']
    keys = Y.odor_keys(cols, {'x': 'InChI=1S/X/c1', 'y': 'InChI=1S/Y/c2/h1/t8-,9-/m0/s1'})
    assert keys == {'x': ('1S/X/c1', 'exact'), 'y': ('1S/Y/c2/h1', 'stereo layers dropped')}
    with pytest.raises(KeyError):
        Y.odor_keys(cols, {'z': 'InChI=1S/Z'})
    emb = pd.DataFrame([[0, 1, 2], [1, 0, 3], [2, 3, 0]], index=cols, columns=cols, dtype=float)
    assert Y.odor_matrix(emb, keys, ['y', 'x']).tolist() == [[0, 1], [1, 0]]
    assert list(Y.upper(np.arange(16).reshape(4, 4))) == [1, 2, 3, 6, 7, 11]


def test_bootstrap_draws_control_first():
    X = {'A': np.arange(12.0).reshape(3, 4), 'B': np.arange(8.0).reshape(2, 4)}
    out = Y.bootstrap(X, lambda v: dict(m=float(v.mean())), n=4, seed=5)
    rng = np.random.default_rng(5)
    a = [X['A'][rng.integers(0, 3, 3)].mean() for _ in range(4)]
    b = [X['B'][rng.integers(0, 2, 2)].mean() for _ in range(4)]
    assert np.allclose(out['A'].m, a) and np.allclose(out['B'].m, b)
