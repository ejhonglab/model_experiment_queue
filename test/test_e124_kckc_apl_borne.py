#!/usr/bin/env python3
"""Tests for E124 (the inhibition among KCs as a part of APL; experiments/E124_kckc_apl_borne.py) and E10_joint_fit.py's
--off-kckc-scale: the matrices (the type-structured draw's totals, partners and weights; the seeds of the equal ladder),
the dose arms' APL-off KC-KC scale, the planned points, part B's reading, the corner rule, the fits' commands and gate,
and the option's identity when off (the call to optimize_mb_model_params unchanged; on pyolfsysm's real optimizer with a
made-up simulator, as test_e111_residual_apl.py). No model simulation."""

import inspect
from typing import List

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
X = pytest.importorskip('experiments.E124_kckc_apl_borne')
E10 = pytest.importorskip('experiments.E10_joint_fit')
E123 = pytest.importorskip('experiments.E123_kckc_form')


def _toy(seed=0, n_g=60, n_apbp=30, n_ab=80):
    rng = np.random.default_rng(seed)
    types = np.array(['g'] * n_g + ["a'b'"] * n_apbp + ['ab'] * n_ab, dtype=object)
    n = len(types)
    ref = rng.poisson(0.5, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.15)
    # more weight within type, as in the hemibrain
    same = types[:, None] == types[None, :]
    ref[same] *= 3.0
    np.fill_diagonal(ref, 0.0)
    ref *= n / ref.sum()
    return ref, types


@pytest.mark.parametrize('seed', [0, 1])
def test_type_structured_partners(seed):
    ref, types = _toy()
    M = X.type_structured_partners(ref, types, 16, seed)
    assert np.isclose(M.sum(), ref.sum())
    assert np.allclose(np.diag(M), 0.0)
    assert ((M > 0).sum(axis=1) == 16).all()
    for a in X.TYPES:
        rows = types == a
        # equal weights within a postsynaptic type; each KC the type's mean amount
        vals = np.unique(M[rows][M[rows] > 0])
        assert len(vals) == 1
        assert np.allclose(M[rows].sum(axis=1), ref[rows].sum() / rows.sum())
        for b in X.TYPES:
            t = ref[np.ix_(rows, types == b)].sum()
            # each block's total up to the rounding of the counts (less than one entry per block, largest remainders)
            assert abs(M[np.ix_(rows, types == b)].sum() - t) < vals[0]
    other = X.type_structured_partners(ref, types, 16, 1 - seed)
    assert not np.array_equal(M, other)
    assert np.array_equal(M, X.type_structured_partners(ref, types, 16, seed))


def test_block_shares_kept_by_types_not_by_ladder():
    ref, types = _toy()
    hb = X.block_shares(ref, types)['within_type_share']
    ts = X.block_shares(X.arm_array('TS16s0', ref, types), types)['within_type_share']
    k16 = X.block_shares(X.arm_array('K16s1', ref, types), types)['within_type_share']
    assert abs(ts - hb) < 0.02
    assert abs(k16 - hb) > 0.1


def test_arms_matrices():
    ref, types = _toy()
    assert np.array_equal(X.arm_array('HBoff25', ref, types), ref)
    assert np.array_equal(X.arm_array('HB', ref, types), ref)
    assert np.array_equal(X.arm_array('K4s0', ref, types), E123.arm_array('K4', ref, types))
    assert np.array_equal(X.arm_array('K16s0', ref, types), E123.arm_array('K16', ref, types))
    assert np.array_equal(X.arm_array('SH0', ref, types), E123.arm_array('SH0', ref, types))
    for k in (4, 16):
        a, b, c = (X.arm_array(f'K{k}s{s}', ref, types) for s in (0, 1, 2))
        assert not np.array_equal(a, b) and not np.array_equal(b, c)
        for m in (a, b, c):
            assert np.isclose(m.sum(), ref.sum()) and ((m > 0).sum(axis=1) == k).all()
    with pytest.raises(ValueError):
        X.arm_array('nope', ref, types)


def test_doses_and_plan():
    assert X.A1_ARMS == ('HBoff25', 'HBoff50', 'HBoff75')
    assert [X.off_kckc_fraction(a) for a in X.A1_ARMS] == [0.25, 0.5, 0.75]
    assert X.off_kckc_fraction('HBoff0') == 0.0 and X.off_kckc_fraction('TS16s0') is None
    p = dict(wKCKC_scale=40.0)
    assert X.off_extra_of('HBoff50', p) == dict(wKCKC_scale=20.0)
    assert X.off_extra_of('HBoff0', p) == dict(wKCKC_scale=0.0)
    assert X.off_extra_of('K16s1', p) is None
    pl = X.planned()
    assert len(pl) == 72 and len(set(pl)) == 72
    assert [a for a, _ in pl[:24]] == [a for a in X.A1_ARMS for _ in X.WIRINGS]
    assert X.planned(wirings=['P3']) == [(a, 'P3') for a in X.RUN_ORDER]


def _cls(kg):
    return {a: dict(complete=True, keeps_gain=v, cls='keeps the gain' if v else 'neither') for a, v in kg.items()}


def _rd_for(arms):
    rows = [dict(kind=a, unit=u, readout='cal', all_pattern=0.05) for a in arms for u in X.POPS]
    return pd.DataFrame(rows)


def test_b_reading():
    need = ['TS16s0', 'TS16s1'] + list(X.K_DRAWS[4]) + list(X.K_DRAWS[16])
    rd = _rd_for(need)
    kg = {a: False for a in need}
    kg.update(TS16s0=True, TS16s1=True)
    r = X.b_reading(_cls(kg), rd)
    assert r['complete'] and r['types_make_up_the_rest'] and not r['gap_of_draws']
    assert r['reading'] == 'the types make up the rest'
    kg.update(K16s1=True)
    r = X.b_reading(_cls(kg), rd)
    assert not r['types_make_up_the_rest'] and r['gap_of_draws'] and r['straddle'] == {4: False, 16: True}
    assert r['reading'] == 'a gap of draws'
    kg = {a: False for a in need}
    kg.update(K4s0=True)
    r = X.b_reading(_cls(kg), rd)
    assert r['reading'] == 'a gap of draws' and r['straddle'][4]
    kg = {a: False for a in need}
    assert X.b_reading(_cls(kg), rd)['reading'] == 'neither'
    c = _cls(kg)
    del c['K16s2']
    assert not X.b_reading(c, rd)['complete']


def test_corner_rule():
    start = dict(wAPLKC_scale=0.1)
    row = pd.Series(dict(thr_factor=1.5, wAPLKC_scale=0.05, aplpn_factor_vs_e3=1.2))
    c = X.corner_of(row, start)
    assert not c['thr_at_bound'] and not c['aplkc_at_bound'] and not c['aplpn_at_cap']
    assert X.corner_of(pd.Series(dict(thr_factor=3.99, wAPLKC_scale=0.05)), start)['thr_at_bound']
    assert X.corner_of(pd.Series(dict(thr_factor=0.251, wAPLKC_scale=0.05)), start)['thr_at_bound']
    assert X.corner_of(pd.Series(dict(thr_factor=1.0, wAPLKC_scale=0.1 / 30)), start)['aplkc_at_bound']
    assert X.corner_of(pd.Series(dict(thr_factor=1.0, wAPLKC_scale=0.1, aplpn_factor_vs_e3=3.0)), start)['aplpn_at_cap']


def test_fit_commands_and_gate():
    c = X.fit_command('A2', 'real', python='py')
    assert c[:3] == ['py', '-u', 'experiments/E10_joint_fit.py']
    s = ' '.join(c)
    assert '--off-kckc-scale 0.0' in s and '--maxfev 150' in s and '--aplpn-max 3' in s and '--loss apl' in s
    assert '--experiment E124_kckc_apl_borne/fits/A2_real' in s
    assert '--init-params results/E124_kckc_apl_borne/starts.json:A2_real' in s
    assert '--kckc-matrix' not in s and '--kckc-fixed' not in s and '--wiring' not in s
    c3 = ' '.join(X.fit_command('A3', 'P3', python='py'))
    assert '--kckc-matrix ' in c3 and 'kckc_SH0.parquet' in c3 and '--wiring claw_permutation_within_type' in c3
    assert X.a3_allowed_from({}) [0] is False
    assert X.a3_allowed_from(dict(A2=dict(complete=True, corner=True, keeps_gain=True)))[0] is False
    assert X.a3_allowed_from(dict(A2=dict(complete=True, corner=False, keeps_gain=False, cls='neither')))[0] is False
    assert X.a3_allowed_from(dict(A2=dict(complete=True, corner=False, keeps_gain=True, cls='equal')))[0] is True


# ---------------------------------------------------------------------------
# E10_joint_fit.py's --off-kckc-scale

def test_off_kckc_overrides_and_tag():
    o = {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0}}
    assert E10.off_kckc_overrides(o, None) is o
    q = E10.off_kckc_overrides(o, 0)
    assert q == {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0, 'wKCKC_scale': 0.0}}
    assert o == {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0}}
    assert type(q['off']['wKCKC_scale']) is float
    for bad in (-1.0, float('nan'), float('inf')):
        with pytest.raises(ValueError):
            E10.off_kckc_overrides(o, bad)
    assert E10.off_kckc_tag(None) == '' and E10.off_kckc_tag(0.0) == 'kckcoff0_' and E10.off_kckc_tag(12.5) == 'kckcoff12.5_'
    assert inspect.signature(E10.fit_one).parameters['off_kckc'].default is None


def _fake_base(n=12):
    idx = pd.Index([f'kc{i}' for i in range(n)], name=H.KC_ID_LEVEL)
    wnorm = {k: pd.Series(np.linspace(1, 2, n), index=idx) for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')}
    breadth = pd.Series(np.arange(n, dtype=float), index=idx)
    initial = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.1, wPNAPL_scale=1e-3)
    return initial, wnorm, breadth


def test_fit_one_call_unchanged_without_the_option(monkeypatch):
    import pyolfsysm.mb_model as mb
    got = []

    def fake(orn, **k):
        got.append(k)
        return dict(history=pd.DataFrame(), params=dict(k['initial_params']), loss=0.0)
    monkeypatch.setattr(mb, 'optimize_mb_model_params', fake)
    initial, wnorm, breadth = _fake_base()
    kw = dict(kckc0=60.0, maxfev=1, method='Nelder-Mead', thr0_factor=0.7, verbose=False)

    def loss_fn(o):
        return {}
    E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw)
    E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw, off_kckc=None)
    E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw, off_kckc=0.0)
    a, b, c = got
    assert set(a) == set(b) == set(c)
    for k, v in a.items():
        assert b[k] is v or (isinstance(v, (pd.Series, pd.DataFrame)) and b[k].equals(v)) or b[k] == v, k
        if k != 'condition_param_overrides':
            assert c[k] is v or (isinstance(v, (pd.Series, pd.DataFrame)) and c[k].equals(v)) or c[k] == v, k
    assert a['condition_param_overrides'] == {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0}}
    assert c['condition_param_overrides'] == {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0,
                                                                'wKCKC_scale': 0.0}}
    assert 'wKCKC_scale' in c['params_to_optimize']
    with pytest.raises(NotImplementedError):
        E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw, off_kckc=0.0,
                    apl_residual=E10.apl_residual_spec(0.0))


def _fake_simulator(calls: List[dict], n_kcs: int = 60):
    odors = list(H.APL16_ODORS)

    def make(orn, initial_params=None, **kws):
        base = np.random.default_rng(0).gamma(1.0, 1.0, (n_kcs, len(odors)))
        cols = pd.Index(odors, name='odor')

        def simulate(params, kc_thr_factor=None, spike_times=False):
            calls.append(dict(params))
            drive = (base * 40.0 / params['fixed_thr'] - 3.0 * params['wAPLKC_scale'] - 2.0 * params['wAPLPN_scale']
                     - 0.002 * params.get('wKCKC_scale', 0.0))
            s = pd.DataFrame(np.floor(np.clip(drive, 0.0, None) * 3.0), columns=cols)
            return (s >= 2).astype(float), s
        return simulate, dict(initial_params)
    return make


def _fit(monkeypatch, **kw):
    import pyolfsysm.mb_model as mb
    calls: List[dict] = []
    monkeypatch.setattr(mb, 'make_fixed_param_kc_simulator', _fake_simulator(calls))
    initial, wnorm, breadth = _fake_base()

    def loss_fn(outputs):
        n_on = float(outputs['on'][0].to_numpy().sum())
        n_off = float(outputs['off'][0].to_numpy().sum())
        return dict(loss=((n_on - 60.0) / 10.0) ** 2 + ((n_off - 130.0) / 10.0) ** 2)
    res, init, bounds = E10.fit_one(2.0, pd.DataFrame(), initial, wnorm, breadth, loss_fn, kckc0=20.0, maxfev=12,
                                    method='Nelder-Mead', thr0_factor=0.7, verbose=False, **kw)
    return res, calls


def test_real_optimizer_off_condition(monkeypatch):
    r0, c0 = _fit(monkeypatch)
    r1, c1 = _fit(monkeypatch, off_kckc=None)
    assert c0 == c1 and r0['loss'] == r1['loss']
    r2, c2 = _fit(monkeypatch, off_kckc=0.0)
    on, off = c2[0::2], c2[1::2]
    assert len(on) == len(off) > 4
    assert all(c['wKCKC_scale'] == 0.0 and c['wAPLKC_scale'] == 0.0 for c in off)
    assert on[0]['wKCKC_scale'] == 20.0 and any(c['wKCKC_scale'] != 0.0 for c in on)
    # each blocked run is its evaluation's APL-on parameters with the overrides
    for a, b in zip(on, off):
        assert {k: v for k, v in a.items() if k not in ('wAPLKC_scale', 'wAPLPN_scale', 'wKCKC_scale')} == \
            {k: v for k, v in b.items() if k not in ('wAPLKC_scale', 'wAPLPN_scale', 'wKCKC_scale')}
    # the first evaluation's APL-on run is the same with and without the option
    assert c2[0] == c0[0]
