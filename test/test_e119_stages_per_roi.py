"""Tests of E119 (experiments/E119_stages_per_roi.py) and its part 0 (data_prep/remy_per_roi_data.py); synthetic data
only, no shared-directory reads."""
import numpy as np
import pandas as pd
import pytest

from data_prep import remy_per_roi_data as P
from experiments import E119_stages_per_roi as E


def test_odor_and_repeats():
    assert P.odor_of_stim('ms @ -3.0') == 'ms'
    assert P.odor_of_stim('2-but @ -3.0') == '2-but'
    assert list(P.trial_repeats(['a', 'b', 'a', 'c', 'b', 'a'])) == [0, 0, 1, 0, 1, 2]


def test_window_means_and_noise_sd():
    t = np.round(np.arange(-5.0, 19.81, 0.2), 6)
    rng = np.random.default_rng(0)
    x = rng.standard_normal((4, 3, len(t))) * 2.0 + 10.0
    wm = P.window_means(x, t, {'w': (0.5, 2.5)})['w']
    m = (t >= 0.5) & (t <= 2.5)
    assert m.sum() == 10
    np.testing.assert_allclose(wm, x[..., m].mean(axis=-1))
    sd = P.noise_sd(x, t)
    assert sd.shape == (3,)
    np.testing.assert_allclose(sd, 2.0, rtol=0.15)
    # a constant offset per trial does not change the noise SD
    np.testing.assert_allclose(P.noise_sd(x + np.arange(4)[:, None, None], t), sd)


def test_stim_rdm_and_frame_rdms_agree():
    rng = np.random.default_rng(1)
    odors_t = ['a', 'b', 'c'] * 3
    x = rng.standard_normal((9, 40, 5))
    F = P.frame_rdms(x, odors_t, ['a', 'b', 'c'])
    assert F.shape == (5, 3, 3)
    for k in range(5):
        np.testing.assert_allclose(F[k], P.stim_rdm(x[:, :, k], odors_t, ['a', 'b', 'c']), atol=1e-12)
    np.testing.assert_allclose(np.diagonal(F, axis1=1, axis2=2), 0, atol=1e-12)


def test_match_flies_one_to_one():
    rng = np.random.default_rng(2)
    n = 6
    th = []
    for _ in range(3):
        a = rng.random((n, n))
        th.append((a + a.T) / 2)
    th = np.stack(th)
    ours = {'x': th[2] + 0.001 * rng.random((n, n)), 'y': th[0], 'z': rng.random((n, n))}
    m = P.match_flies({k: (v + v.T) / 2 for k, v in ours.items()}, th).set_index('acq')
    assert m.loc['x', 'dash_fly'] == 3 and m.loc['y', 'dash_fly'] == 1
    assert m.loc['x', 'r'] > 0.99


def test_cut_trials_interpolates_and_subtracts_baseline():
    ft = np.arange(0, 100, 0.3)
    tr = np.vstack([np.sin(ft), ft])
    t = np.round(np.arange(-5.0, 19.81, 0.2), 6)
    out = P.cut_trials(tr, ft, [20.0, 50.0], t)
    assert out.shape == (2, 2, len(t))
    b = P.frame_mask(t, (-5.0, 0.0))
    np.testing.assert_allclose(np.median(out[:, :, b], axis=-1), 0, atol=1e-12)
    # the linear trace's cut is t minus its baseline median
    np.testing.assert_allclose(out[0, 1], t - np.median(t[b]), atol=1e-9)
    with pytest.raises(ValueError):
        P.cut_trials(tr, ft, [1.0], t)


def test_trial_pair_rdms():
    rng = np.random.default_rng(3)
    base = rng.standard_normal((2, 50, 1))
    odors_t = ['a', 'b', 'a', 'b', 'a', 'b']
    x = np.stack([base[0 if o == 'a' else 1] for o in odors_t]) + 0.0 * rng.standard_normal((6, 50, 1))
    D = P.trial_pair_rdms(x, odors_t, ['a', 'b'])
    np.testing.assert_allclose(D[0, 0, 0], 0, atol=1e-12)          # noise-free: perfectly reliable
    r = np.corrcoef(base[0, :, 0], base[1, :, 0])[0, 1]
    np.testing.assert_allclose(D[0, 0, 1], 1 - r, atol=1e-12)


def test_opr_top_share_and_roi_stats():
    np.testing.assert_allclose(E.opr_distribution(np.array([0, 1, 2, 3, 5, 12])), [1 / 6, 1 / 6, 2 / 6, 1 / 6, 1 / 6])
    V = np.zeros((100, 2))
    V[0, 0] = 10.0
    V[1:, 0] = 0.0
    V[:, 1] = 1.0
    assert E.top_share(V) == pytest.approx((1.0 + 0.01) / 2)
    odors = list(E.ODORS)
    X = pd.DataFrame(np.zeros((10, len(odors))), columns=odors)
    X.iloc[:5, 0] = 3.0
    X.iloc[0, 1] = -3.0
    s = E.roi_stats(X, 2.0)
    assert s['rate'] == pytest.approx(5 / (10 * len(odors)))
    assert s['per_odor'].iloc[0] == pytest.approx(0.5)
    assert s['opr'][0] == pytest.approx(0.5) and s['opr'][1] == pytest.approx(0.5)
    assert s['inhibited'] == pytest.approx(1 / (10 * len(odors)))
    assert s['resp_cv'] == pytest.approx(0.0)


def test_matched_threshold():
    odors = list(E.ODORS)
    rng = np.random.default_rng(4)
    V = pd.DataFrame(rng.random((200, len(odors))), columns=odors)
    thr = E.matched_threshold(V, 0.1)
    assert E.roi_stats(V, thr)['rate'] == pytest.approx(0.1, abs=1 / V.size + 1e-9)


def test_profile_band_and_rule():
    m, b = E.profile_band([np.array([0.0, 1.0]), np.array([1.0, 0.0]), np.array([0.5, 0.5])], 'opr')
    np.testing.assert_allclose(m, [0.5, 0.5])
    assert b == pytest.approx(0.5 * (abs(0 - 0.75) + abs(1 - 0.25)))
    assert E.scalar_differs([1, 1, 1, 1, 1, 1, 0, -1]) == (True, 1)
    assert E.scalar_differs([1, 1, 1, 1, 1, -1, -1, 0]) == (False, 0)


def test_within_sides():
    odors = list(E.ODORS)
    stats = [E.roi_stats(pd.DataFrame(np.random.default_rng(i).random((50, len(odors))), columns=odors), 0.8)
             for i in range(4)]
    sp = E.fly_spread({str(i): s for i, s in enumerate(stats)})
    w = E.within(stats[0], sp)
    assert all(w[f'{s}_side'] == 0 for s in E.SCALARS)
    assert isinstance(w['opr_in'], bool)


def test_pool_units():
    U = pd.DataFrame(np.arange(20, dtype=float).reshape(10, 2), columns=['a', 'b'])
    assert E.pool_units(U, '1', 0).shape == (10, 2)
    P3 = E.pool_units(U, 'r3', 0)
    assert P3.shape == (3, 2)
    np.testing.assert_allclose(P3['b'] - P3['a'], 1.0)
    g = np.array(['x'] * 6 + ['y'] * 4)
    G = E.pool_units(U, 'g3', 0, g)
    # x: two ROIs of 3; y: one of 3 and a remainder of 1 (< 3 / 2, dropped)
    assert G.shape == (3, 2)
    assert ((G['a'] <= 10) | (G['a'] >= 12)).all()   # every ROI holds units of one glomerulus (x: a 0..10, y: 12..18)
    assert (G['a'] >= 12).sum() == 1 and (G['a'] <= 10).sum() == 2
    assert G['a'][G['a'] <= 10].sum() == pytest.approx(np.arange(0, 12, 2).sum() / 3)
    with pytest.raises(ValueError):
        E.pool_units(U, 'g3', 0)


def test_noise_removed_noiseless_equals_trial_mean():
    rng = np.random.default_rng(5)
    X = rng.random((40, len(E.ODORS)))
    tr = np.repeat(X[None], 3, axis=0)
    for ms in E.MEASURES:
        a = E.noise_removed(tr, ms)
        b = E.trial_mean_distances(tr, ms)
        np.testing.assert_allclose(a.to_numpy(), b.to_numpy(), atol=1e-10)


def test_noisy_trials_noise_level():
    odors = list(E.ODORS)
    rng = np.random.default_rng(6)
    Z = pd.DataFrame(rng.standard_normal((4000, len(odors))), columns=odors)
    rel = pd.Series(0.5, index=odors)
    T = E.noisy_trials(Z, rel, np.random.default_rng(0))
    r = np.corrcoef(T[0, :, 0], T[1, :, 0])[0, 1]
    assert r == pytest.approx(0.5, abs=0.05)
