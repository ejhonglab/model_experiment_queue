#!/usr/bin/env python3
"""Tests for E89 (the co-response read with a checked estimate and the breadth weighting inside the models' APL;
experiments/E89_apl_breadth_weighting.py): the weightings (exponent 0 gives factors of 1, F at exponent 2 gives
`breadth_weights` exactly, S and N are 0 for a KC without a response of two or more spikes and N is at most 16, a
permuted control keeps the factors as a set within each KC type); the candidates at zero noise (each returns the truth,
the geometric mean equals the two estimates where they agree); noise level 3 at its middle share gives the flies'
same-repeat gap and repeatability; the choice and the gate on made-up tables; and C0's two identities (slow, skipped
with QUICK=1). The synthetic tests need no shared data; the others are skipped without the shared directory."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E89_apl_breadth_weighting')

ODORS = list(H.APL16_ODORS)
QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


def _off_counts(rng, n_kc=200):
    lam = np.where(rng.random((n_kc, len(ODORS))) < 0.12, rng.uniform(2, 12, (n_kc, len(ODORS))), 0.1)
    idx = pd.MultiIndex.from_arrays([np.arange(1000, 1000 + n_kc), np.where(np.arange(n_kc) % 3 == 0, "ab", np.where(np.arange(n_kc) % 3 == 1, "g", "a'b'"))],
                                    names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    return pd.DataFrame(rng.poisson(lam).astype(float), index=idx, columns=ODORS)


def _claw_weights(rng, kc_ids, claws_per_kc=4):
    ids = np.repeat(np.asarray(kc_ids), claws_per_kc)
    idx = pd.MultiIndex.from_arrays([ids, np.arange(len(ids))], names=[H.KC_ID_LEVEL, 'claw_id'])
    return pd.Series(rng.uniform(0.2, 2.0, len(ids)), index=idx)


def test_weightings_are_twelve_and_the_reference():
    ws = X.weightings()
    assert len(ws) == 13 and ws[0] == X.NO_WEIGHTING and ws[-1] == X.REFERENCE
    assert len(set(X.wtag(w) for w in ws)) == 13
    assert sorted(k for m, k in ws if m == 'F' and k > 0) == [0.25, 0.5, 1.0, 2.0]


def test_exponent_0_gives_factors_of_1_and_f2_gives_breadth_weights():
    from experiments import E10_joint_fit as J
    rng = np.random.default_rng(3)
    off = _off_counts(rng)
    ids = off.index.get_level_values(H.KC_ID_LEVEL)
    br = pd.Series(rng.exponential(5.0, len(ids)).round(), index=pd.Index(ids, name=H.KC_ID_LEVEL))
    wn = _claw_weights(rng, ids)
    for m in X.MEASURES:
        mv = X.measure_values(m, off, br)
        f0 = X.weighting_factors(mv, 0.0)
        assert (f0.to_numpy() == 1.0).all()
        assert np.array_equal(X.apl_kc_weights(wn, mv, 0.0).to_numpy(), wn.to_numpy())
    w2 = X.apl_kc_weights(wn, X.measure_values('F', off, br), 2.0)
    assert np.array_equal(w2.to_numpy(), J.breadth_weights(wn, br, 2.0).to_numpy())


def test_s_and_n_of_a_kc_without_responses():
    rng = np.random.default_rng(4)
    off = _off_counts(rng)
    off.iloc[:5] = 1.0          # one spike to every odor: no response
    off.iloc[5] = 2.0           # two spikes to every odor: 16 responses
    S = X.measure_values('S', off, None)
    N = X.measure_values('N', off, None)
    assert (S.iloc[:5] == 0).all() and (N.iloc[:5] == 0).all()
    assert N.iloc[5] == len(ODORS) and S.iloc[5] == 2.0 * len(ODORS)
    assert N.max() <= 16
    c = off.to_numpy()
    assert np.allclose(S.to_numpy(), np.where(c >= 2, c, 0).sum(axis=1))


def test_permuted_control_keeps_the_factors_within_type():
    from experiments import E80_screen as S80
    rng = np.random.default_rng(5)
    off = _off_counts(rng)
    fac = X.weighting_factors(X.measure_values('S', off, None), 1.0)
    types = off.index.get_level_values(H.KC_TYPE_LEVEL)
    fi = pd.MultiIndex.from_arrays([fac.index, types], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    for seed in X.CONTROL_SEEDS:
        p = S80.permute_within_type(pd.Series(fac.to_numpy(), index=fi), seed)
        assert not np.array_equal(p.to_numpy(), fac.to_numpy())
        for t in set(types):
            m = np.asarray(types == t)
            assert np.array_equal(np.sort(p.to_numpy()[m]), np.sort(fac.to_numpy()[m]))


@pytest.mark.skipif(QUICK or not _shared(), reason="C0's identity on the fits as they are: three simulations (minutes)")
def test_c0_f2_reproduces_the_saved_fit_real_wiring():
    from experiments import E83_cross_wiring_fit as E
    model = E.WiringModel('real', check=True)
    parts = E.wiring_parts('real')
    wk = X.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], 2.0)
    assert np.array_equal(wk.to_numpy(), model._kc0.ravel())
    X.set_weights(model, X.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], 0.5))
    assert not np.array_equal(model._kc.ravel(), model._kc0.ravel())
    X.set_weights(model, wk)
    on = model.run_on(1.0)
    X.set_weights(model, None)
    assert E._same({'on': on, 'off': model.saved['off']}, model.saved) == (0, 0)


@pytest.mark.skipif(QUICK or not _shared(), reason="C0's identity on the base: a pool of one-odor simulators (minutes)")
def test_c0_base_reproduces_e87_d2_real_wiring():
    from experiments import E80_screen as S80
    from experiments import E83_cross_wiring_fit as E
    b = X.base_thresholds('real')
    parts = E.wiring_parts('real')
    fit = E.saved_fit('real')
    pool = X.C2Pool('real', 16)
    try:
        thr = [float(b['thr'][c]) for c in pool.cols]
        tg = S80.b3_targets(X.C2_UNITS, pool.cols, int(fit['saved']['on'][0].shape[0]))
        wk = X.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], 2.0)
        outs, recs = pool.run_task(dict(thr=thr, target_on=tg['target_on'], label='test'), wk)
    finally:
        pool.close()
    d2 = X.E87._load_d2(f'D2_b3_thrapl_{X.E87.D2_UNITS}_real')
    assert E._same({'on': outs['on'], 'off': outs['off']}, d2) == (0, 0)


def _zero_noise_recording(rng, n_kc=400, p=0.12):
    """A noise-free recording (three identical repeats) with no value within 0.3 of the criterion 1."""
    C = np.clip(rng.normal(0.0, 0.15, (n_kc, len(ODORS))), -0.6, 0.6)
    resp = rng.random((n_kc, len(ODORS))) < p
    C[resp] = 1.5 + rng.exponential(2.0, int(resp.sum()))
    return C, np.repeat(C[None], 3, axis=0)


def test_candidates_return_the_truth_at_zero_noise():
    E87, E88 = X.E87, X.E88
    rng = np.random.default_rng(7)
    C, trials = _zero_noise_recording(rng)
    v = X.recording_values(C, trials, 1.0, 1.0, {k: X.LC_EDGES[k] for k in ('crit', 'coarse')})
    s, O, R = E88.factors_of_matrix(E87.thresholded(C, 1.0))
    for cand, lc, tol in (('tm', None, 1e-12), ('nr', None, 1e-12), ('gm', None, 1e-12), ('hf', None, 1e-6),
                          ('lc', 'crit', 1e-6), ('lc', 'coarse', 1e-6)):
        F = X.cand_factors(v, cand, lc)
        assert np.allclose(F.O, O, atol=tol, rtol=0, equal_nan=True), cand
        assert np.allclose(F.R, R, atol=1e-12, rtol=0, equal_nan=True) and np.allclose(F.s, s, atol=1e-12, rtol=0, equal_nan=True)
    # the geometric mean equals the two estimates where they agree
    v2 = dict(v, O_nr=v['O_tm'].copy())
    assert np.allclose(X.cand_factors(v2, 'gm').O, v['O_tm'], atol=1e-14, equal_nan=True)


def test_hf_overlap_recovers_a_known_mixture():
    rng = np.random.default_rng(8)
    n_kc, n_rep = 20000, 3
    z = rng.random((n_kc, len(ODORS))) < 0.15
    z[:, 1] = z[:, 0] & (rng.random(n_kc) < 0.6) | (~z[:, 0] & (rng.random(n_kc) < 0.05))
    h, e = 0.8, 0.03
    x = np.where(z[None], rng.random((n_rep, n_kc, len(ODORS))) < h, rng.random((n_rep, n_kc, len(ODORS))) < e).astype(float)
    O, f, prm = X.hf_overlap(x)
    zf = z.astype(float)
    O_true = (zf.T @ zf / n_kc)[X._iu()] / np.sqrt(np.outer(zf.mean(axis=0), zf.mean(axis=0))[X._iu()])
    assert np.allclose(f, zf.mean(axis=0), atol=0.01)
    assert np.abs(O - O_true).max() < 0.03


def test_same_repeat_gap_is_the_harness_statistic():
    rng = np.random.default_rng(9)
    tr = np.where(rng.random((3, 300, len(ODORS))) < 0.1, rng.exponential(2.0, (3, 300, len(ODORS))), 0.0)
    tr += rng.normal(0, 0.3, (3, 300, 1)) * (tr > 0)
    assert np.isclose(X.same_repeat_gap(tr), H.same_repeat_gap(H.distances_by_kind(tr, ODORS, values='graded')), atol=1e-12)


@pytest.mark.skipif(not _shared(), reason='needs the shared directory (a saved fit and the flies)')
def test_level3_middle_share_gives_the_flies_gap_and_repeatability():
    E88 = X.E88
    gap = X.flies_gap_target('z')
    rel = E88.flies_reliability('z')
    w, key, K = X.a_runs()['b2'][0]
    mr = X.ModelReading89(X._load_run(key), label=w, seed=X.SEEDS[1], K=K, rel_target=rel,
                          gap_target={c: v['thresholded'] for c, v in gap.items()})
    for c in X.CONDS:
        _, info = mr.level_values('L3', c)
        assert abs(info['gap'] - gap[c]['thresholded']) <= 0.005, (c, info)
        assert abs(info['reliability'] - rel[c]) <= 0.01, (c, info)
        assert abs(info['rate_tm'] - mr.rate['r1'][c]) < 0.002 and abs(info['rate_st'] - mr.rate['r1'][c]) < 0.002
    # level 1 is E88's readings 1 and 2
    v1, _ = mr.level_values('L1', 'on')
    f1, f2 = mr.factors('r1')['on'].mean(), mr.factors('r2')['on'].mean()
    assert np.allclose(v1['O_tm'], f1.O, equal_nan=True) and np.allclose(v1['O_nr'], f2.O, equal_nan=True)
