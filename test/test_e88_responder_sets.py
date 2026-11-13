#!/usr/bin/env python3
"""Tests for E88 (the common error split into who responds and how strongly; experiments/E88_responder_sets.py), parts
A to C (no simulations; part D is not implemented): R * L * G reproduces the thresholded cosine similarity of a matrix
exactly and R * L one minus the binarized distance, with L and G computed from their definitions; a cell with every
factor from the flies is the flies' effect and with none the model's own; at zero noise the estimators return the
truth; a per-odor cut gives the asked fractions; the noise-free factors are `H.noise_free_pair_distances`'s; and, on the
shared data, the six typical recordings are the eight without E77's STRONG_TNT, the flies' reading 1 has the
registry's thresholded mean effect, and reading 1 of one wiring reproduces E74's `read_zcell` thresholded mean effect at
the same seed. The synthetic tests need no shared data; the others are skipped without the shared directory."""

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E88_responder_sets')

ODORS = list(H.APL16_ODORS)


def _sparse(rng, n_kc=300, p=0.12, scale=3.0):
    """A KC x odor matrix with sub-threshold values 0 (about a share p responding, amplitudes from 1)."""
    x = 1.0 + rng.exponential(scale, (n_kc, len(ODORS)))
    return np.where(rng.random((n_kc, len(ODORS))) < p, x, 0.0)


def _counts(rng, n_kc=400, extra=0.0):
    """Synthetic spike counts: each KC with a few odors it fires for (more with `extra`, as with APL off)."""
    lam = np.where(rng.random((n_kc, len(ODORS))) < 0.10 + extra, rng.uniform(2, 12, (n_kc, len(ODORS))), 0.05)
    return pd.DataFrame(rng.poisson(lam).astype(float), columns=ODORS)


CAL = dict(K=12.0, noise=0.06, noise_ratio=0.55, q=0.13, drift=1.0)


def test_factors_multiply_to_the_similarity():
    rng = np.random.default_rng(1)
    x = _sparse(rng)
    s, O, R = X.factors_of_matrix(x)
    F = X.Factors(s, O, R)
    i, j = X._iu()
    N = x.shape[0]
    b = x > 0
    n = b.sum(axis=0).astype(float)
    nab = (b[:, i] & b[:, j]).sum(axis=0).astype(float)
    R_def = np.sqrt(n[i] * n[j]) / N
    L_def = nab * N / (n[i] * n[j])
    rms = np.sqrt((x ** 2).sum(axis=0) / n)
    G_def = ((x[:, i] * x[:, j]).sum(axis=0) / nab) / (rms[i] * rms[j])
    cos = (x.T @ x)[i, j] / np.sqrt((x ** 2).sum(axis=0)[i] * (x ** 2).sum(axis=0)[j])
    ok = nab > 0
    assert np.allclose((R_def * L_def * G_def)[ok], cos[ok], rtol=0, atol=1e-12)
    bf = b.astype(float)
    assert np.allclose((R_def * L_def), (bf.T @ bf)[i, j] / np.sqrt(n[i] * n[j]), rtol=0, atol=1e-12)
    assert np.allclose(R, R_def, atol=1e-14) and np.allclose(F.L[ok], L_def[ok], rtol=1e-12)
    assert np.allclose(F.G[ok], G_def[ok], rtol=1e-12)
    assert np.allclose((F.R * F.L * F.G)[ok], s[ok], rtol=0, atol=1e-12)
    assert np.allclose(s, 1 - X.E87.cos_dist(x), atol=1e-15) and np.allclose(O, 1 - X.E87.cos_dist(bf), atol=1e-15)


def test_cells_at_the_ends():
    rng = np.random.default_rng(2)
    fl = {c: X.rec_factors([_sparse(rng, p=0.07 if c == 'on' else 0.15) for _ in range(4)]).mean() for c in X.CONDS}
    mo = {c: X.rec_factors([_sparse(rng, p=0.08 if c == 'on' else 0.13) for _ in range(3)]).mean() for c in X.CONDS}
    assert np.allclose(X.cell_effect(mo, fl, ('R', 'L', 'G')), X.effect(fl), atol=1e-12, equal_nan=True)
    assert np.allclose(X.cell_effect(mo, fl, ()), X.effect(mo), atol=1e-12, equal_nan=True)
    # the binarized effect is the cell of R and L with G set to 1
    ones = {c: X.Factors(fl[c].O, fl[c].O, fl[c].R) for c in X.CONDS}   # G = s / O = 1
    assert np.allclose(X.cell_effect(ones, fl, ('R', 'L')), X.effect(fl, 'O'), atol=1e-12, equal_nan=True)


def test_estimators_return_the_truth_at_zero_noise():
    rng = np.random.default_rng(3)
    o = {'on': (None, _counts(rng)), 'off': (None, _counts(np.random.default_rng(4), extra=0.08))}
    m = X.ModelReading(o, cal=CAL, noise_scale=0.0)
    m.mult = {'on': 1.0, 'off': 1.0}
    assert m.factors('r1')['on'].n == X.N_PSEUDO
    assert 0.02 < m.rate['r1']['on'] < 0.3
    for rd in X.NOISY_READINGS:
        est, tru = m.mean_factors(rd), m.truth(rd)
        for c in X.CONDS:
            for f in X.FACTORS:
                assert np.allclose(est[c].get(f), tru[c].get(f), rtol=0, atol=1e-12, equal_nan=True), (rd, c, f)
    # with noise the readings depart from the truth, and the estimator check sees it
    n = X.ModelReading(o, cal=CAL)
    ec = X.estimator_check([n], ('r1', 'r2'))
    assert set(ec.reading) == {'r1', 'r2'} and ec['err'].abs().max() > 0


def test_reliability_target_is_met():
    rng = np.random.default_rng(5)
    o = {'on': (None, _counts(rng)), 'off': (None, _counts(np.random.default_rng(6), extra=0.08))}
    m = X.ModelReading(o, cal=CAL, rel_target={'on': 0.5, 'off': 0.5})
    mult = m.multipliers()
    for c in X.CONDS:
        if not m.mult_info[c]['at_bound']:
            assert abs(m.mult_info[c]['rel_r1b'] - 0.5) < 0.005
    m.factors('r1b')
    for c in X.CONDS:   # 1b keeps reading 1's response rates
        assert abs(m.rate['r1b'][c] - m.rate['r1'][c]) < 1e-3
    assert all(v > 0 for v in mult.values())


def test_per_odor_cut_gives_the_asked_fractions():
    rng = np.random.default_rng(7)
    x = _sparse(rng, n_kc=500, p=0.2)
    have = (x > 0).mean(axis=0)
    fr = have * rng.uniform(0.2, 0.9, len(ODORS))
    cut = X.E87.per_odor_oracle(x, fr, 'divisive', 0.5)
    assert np.array_equal((cut > 0).sum(axis=0), np.round(fr * x.shape[0]).astype(int))
    assert np.all((cut > 0) <= (x > 0))


def test_noise_free_factors_are_the_harness_estimator():
    rng = np.random.default_rng(8)
    tr = rng.normal(0, 1, (3, 200, len(ODORS))) + 2.5 * (rng.random((1, 200, len(ODORS))) < 0.12)
    thr = 1.2
    s, O, R, rel = X.factors_of_trials(X.E87.thresholded(tr, thr))
    keys = [(ODORS[a], ODORS[b]) for a, b in zip(*X._iu())]

    def tri(sr):
        sr.index = [tuple(k) for k in sr.index]
        return np.array([sr.get(k, sr.get(k[::-1], np.nan)) for k in keys])
    g = tri(H.noise_free_pair_distances(tr, ODORS, values='thresholded', thresh=thr))
    b = tri(H.noise_free_pair_distances(tr, ODORS, values='binary', thresh=thr))
    assert np.allclose(s, 1 - g, atol=1e-12, equal_nan=True) and np.allclose(O, 1 - b, atol=1e-12, equal_nan=True)
    assert np.all((rel > 0) & (rel <= 1 + 1e-12))


def test_group_masks():
    m = X.group_masks()
    assert m['ester'].sum() == 15 and m['alcohol'].sum() == 6 and m['other'].sum() == 99
    assert m['nearest4'].sum() == 4 and m['six'].sum() == 6 and not (m['nearest4'] & ~m['ester']).any()


# the shared data

def _shared():
    try:
        H.data_targets_path().read_bytes()
    except (FileNotFoundError, OSError):
        pytest.skip('the shared directory (data-target registry) is not available')


def test_six_typical_recordings():
    _shared()
    from experiments import E77_input_side_common_error as E77
    allt = list(X.E87.trialmean_mats('z')['TNT'])
    six = X.six_typical()
    assert len(allt) == 8 and len(six) == 6
    assert set(six) == set(allt) - {tuple(k) for k in E77.STRONG_TNT}


def test_flies_reading_1_is_the_registry_value():
    _shared()
    rf = X.flies_factors('z', 'r1')
    F = {c: rf[c].mean() for c in X.CONDS}
    assert np.isclose(np.nanmean(X.effect(F)), H.data_target('apl_effect_thresholded_mean'), atol=1e-12)


def test_reading_1_reproduces_e74():
    _shared()
    E74 = pytest.importorskip('experiments.E74_apl_effect_by_odors_per_kc')
    E87 = X.E87
    try:
        w, f, t = E87.model_families()['b2'][0]
        o = E87._load(f, t)
    except Exception as e:   # the saved fits are not on this host
        pytest.skip(f'no saved breadth-2 fit: {e}')
    m = X.ModelReading(o, label=w)
    blk = E74.read_zcell(m.counts['on'], m.counts['off'], m.cal, noise=True)
    mine = float(np.nanmean(X.effect(m.mean_factors('r1'))))
    assert abs(mine - blk['mean_thresholded']) < 1e-12
    assert abs(m.rate['r1']['on'] - blk['rate_on']) < 1e-15 and abs(m.rate['r1']['off'] - blk['rate_off']) < 1e-15


@pytest.mark.skipif(__import__('os').environ.get('QUICK') == '1', reason="D1's base identity: a pool of one-odor simulators (minutes)")
def test_d1_base_reproduces_e87_d2_real_wiring():
    """D1's base on the real wiring (per-odor thresholds re-solved at the flies' fractions with APL output blocked, APL
    off) reproduces E87's D2 APL-off run exactly."""
    E87 = pytest.importorskip('experiments.E87_inhibition_form')
    b = X.base_thresholds(H.results_dir(X.EXPERIMENT), 'real')
    pool = E87.WiringPool('real', n_workers=16)
    try:
        cols = list(pool.cols)
        n_kc = int(E87._load_d2('D2_b3_thr_noise_real')['off'][0].shape[0])
        r, s = pool.run(np.tile(np.array([b['thr'][c] for c in cols]), (n_kc, 1)))
    finally:
        pool.close()
    off = E87._load_d2('D2_b3_thr_noise_real')['off']
    assert np.array_equal(s.to_numpy(dtype=float), off[1].to_numpy(dtype=float))
