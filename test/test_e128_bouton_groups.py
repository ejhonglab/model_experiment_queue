"""Tests of E128's reader (experiments/E128_bouton_groups.py); synthetic data only, no shared-directory reads."""
import numpy as np
import pandas as pd
import pytest

from experiments import E128_bouton_groups as E


def _fly(seed=0):
    rng = np.random.default_rng(seed)
    # three tunings: one dominant group (large on odor 0), two broad ones; plus unresponsive ROIs
    t = np.array([[20.0] + [1.0] * 11, [4.0, 5, 6, 5, 4, 5, 6, 5, 4, 5, 6, 5], [6.0, 5, 4, 5, 6, 5, 4, 5, 6, 5, 4, 5]])
    M = np.vstack([t[0] + rng.normal(0, 0.05, (30, 12)), t[1] + rng.normal(0, 0.05, (10, 12)),
                   t[2] + rng.normal(0, 0.05, (10, 12)), rng.normal(0, 0.3, (5, 12))])
    return M


def test_groups_and_shares():
    M = _fly()
    lab = E.fly_groups(M)
    assert (lab[-5:] == 0).all() and len(np.unique(lab[lab > 0])) == 3
    sh = E.variance_shares(M, lab)
    assert sh.sum() == pytest.approx(1.0)
    big = lab[0]
    assert sh.drop(0).idxmax() == big and sh[big] > 0.9
    # by hand: cov(group sum, total) / var(total)
    T = M.sum(axis=0)
    g = M[lab == big].sum(axis=0)
    assert sh[big] == pytest.approx(np.cov(g, T, bias=True)[0, 1] / T.var())


def test_readings():
    M = _fly()
    lab = E.fly_groups(M)
    r = E.readings(M, lab)
    assert r['measured'] == pytest.approx(E.cv(M.mean(axis=0)))
    # leaving out the dominant group leaves two broad, complementary groups plus the remainder: nearly flat
    assert r['leave1'] < 0.1 < r['measured']
    keep = lab != lab[0]
    assert r['leave1'] == pytest.approx(E.cv(M[keep].sum(axis=0)))
    ids = [g for g in np.unique(lab)]
    means = [M[lab == g].mean(axis=0) for g in ids]
    assert r['equal'] == pytest.approx(E.cv(np.mean(means, axis=0)))
    assert r['equal_no_rem'] == pytest.approx(E.cv(np.mean(means[1:], axis=0)))
    assert E.class_min(r) == min(r[k] for k in E.CLASS_READINGS)
    # a given order leaves out those groups, the remainder never
    r2 = E.readings(M, lab, order=[0, lab[-10]])
    assert r2['leave1'] == pytest.approx(E.cv(M[lab != lab[-10]].sum(axis=0)))


def test_null_deal_keeps_sizes_and_remainder():
    lab = E.fly_groups(_fly())
    d = E.null_deal(lab, np.random.default_rng(0))
    assert (d[lab == 0] == 0).all()
    assert sorted(np.bincount(d)) == sorted(np.bincount(lab))


def test_three_stats_matches_e119():
    from experiments import E119_stages_per_roi as E119
    rng = np.random.default_rng(1)
    X = rng.gamma(1.0, 2.0, (300, 12)) - 0.5
    st = E119.roi_stats(pd.DataFrame(X, columns=list(E119.ODORS)), 2.0)
    mine = E.three_stats(X, 2.0)
    for k in E.THREE:
        assert mine[k] == pytest.approx(st[k], abs=1e-12)


def test_deal_population_and_pass():
    U = np.arange(40, dtype=float).reshape(10, 4)
    gl = np.array(['a'] * 5 + ['b'] * 5)
    P, w = E.deal_population(U, gl, [3, 4, 2], np.random.default_rng(0))
    assert P.shape == (9, 4) and sum(w.values()) == 9
    rows = {tuple(r) for r in U}
    assert all(tuple(r) in rows for r in P)
    sp = dict(totals_cv=(0.3, 0.4), resp_cv=(0.4, 0.6), top1=(0.06, 0.07))
    assert E.model_pass(dict(totals_cv=0.35, resp_cv=0.5, top1=0.07), sp)
    assert not E.model_pass(dict(totals_cv=0.35, resp_cv=0.61, top1=0.05), sp)
    assert not E.model_pass(dict(totals_cv=0.29, resp_cv=0.5, top1=0.05), sp)


def test_classes_letters():
    A = pd.DataFrame(dict(acq=list('abcde'), fly_passes=[True] * 4 + [False],
                          **{k: [0.2] * 5 for k in E.CLASS_READINGS + ('equal_big', 'leave5')}))
    R = pd.DataFrame(dict(acq=list('aabbccddee'), below=[True] * 10))
    W = pd.DataFrame(dict(variant=['all'] * 8 + ['big'] * 8, frac_pass=[0.1] * 6 + [0.0] * 2 + [0.0] * 8))
    assert E.classes(A, R, W)['outcome'] == 'a'
    W.loc[:7, 'frac_pass'] = 0.01
    assert E.classes(A, R, W)['outcome'] == 'b'
    R['below'] = False
    assert E.classes(A, R, W)['outcome'] == 'd'
    A['fly_passes'] = False
    assert E.classes(A, R, W)['outcome'] == 'c'
