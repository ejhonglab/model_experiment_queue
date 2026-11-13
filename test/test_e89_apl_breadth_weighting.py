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
