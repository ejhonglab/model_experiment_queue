#!/usr/bin/env python3
"""Tests for E111 (a residual APL output in the blocked flies in place of KC-KC inhibition; model_experiment_queue):
E10_joint_fit.py's --apl-residual (`apl_residual_spec`, `apl_residual_tag`, `apl_residual_params`, `fit_one`'s
apl_residual, `with_apl_residual`) and experiments/E111_residual_apl.py's run lists, starts, readings and verdict.

The option's two identities are checked on pyolfsysm's real optimizer (`optimize_mb_model_params`) with a made-up
simulator in place of olfsysm's (`make_fixed_param_kc_simulator` monkeypatched): at f = 0 every simulation, the history,
the losses and the best outputs equal those of the fit without the option, exactly; at f = 1 every blocked run equals
the APL-on run of its evaluation. Part 0 of E111 checks both on the model itself (the hemibrain wiring of the family).
No simulations of the model here. Needs the shared directory (E10 reads the data-target registry at import); skipped
without it."""

import inspect
from typing import Dict, List

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
E10 = pytest.importorskip('experiments.E10_joint_fit')

ODORS = list(H.APL16_ODORS)


# ---------------------------------------------------------------------------
# E10_joint_fit.py's --apl-residual

def _fake_base(n=12):
    idx = pd.Index([f'kc{i}' for i in range(n)], name=H.KC_ID_LEVEL)
    wnorm = {k: pd.Series(np.linspace(1, 2, n), index=idx) for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')}
    breadth = pd.Series(np.arange(n, dtype=float), index=idx)
    initial = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.1, wPNAPL_scale=1e-3)
    return initial, wnorm, breadth


def _fake_simulator(calls: List[dict], n_kcs: int = 60):
    """A stand-in for `make_fixed_param_kc_simulator`: a deterministic simulate(params) whose spike counts fall with the
    APL output scales and the KC-KC scale and rise as the threshold falls; every call's parameters go to `calls`."""
    def make(orn, initial_params=None, **kws):
        base = np.random.default_rng(0).gamma(1.0, 1.0, (n_kcs, len(ODORS)))
        cols = pd.Index(ODORS, name='odor')

        def simulate(params, kc_thr_factor=None, spike_times=False):
            calls.append(dict(params))
            drive = (base * 40.0 / params['fixed_thr'] - 3.0 * params['wAPLKC_scale'] - 2.0 * params['wAPLPN_scale']
                     - 0.002 * params.get('wKCKC_scale', 0.0))
            s = pd.DataFrame(np.floor(np.clip(drive, 0.0, None) * 3.0), columns=cols)
            return (s >= 2).astype(float), s
        return simulate, dict(initial_params)
    return make


def _loss(seen: List[dict]):
    """A loss that moves with the parameters: the APL-on and APL-off response counts against made-up targets."""
    def loss_fn(outputs) -> Dict[str, float]:
        seen.append({c: tuple(x.copy() for x in outputs[c][:2]) for c in ('on', 'off')})
        n_on = float(outputs['on'][0].to_numpy().sum())
        n_off = float(outputs['off'][0].to_numpy().sum())
        return dict(loss=((n_on - 60.0) / 10.0) ** 2 + ((n_off - 130.0) / 10.0) ** 2, n_on=n_on, n_off=n_off)
    return loss_fn


def _fit(monkeypatch, maxfev=12, **kw):
    import pyolfsysm.mb_model as mb
    calls: List[dict] = []
    seen: List[dict] = []
    monkeypatch.setattr(mb, 'make_fixed_param_kc_simulator', _fake_simulator(calls))
    initial, wnorm, breadth = _fake_base()
    res, init, bounds = E10.fit_one(2.0, pd.DataFrame(), initial, wnorm, breadth, _loss(seen), kckc0=20.0, maxfev=maxfev,
        method='Nelder-Mead', thr0_factor=0.7, verbose=False, **kw)
    return res, calls, seen


def _same_frames(a, b) -> bool:
    return all(x.equals(y) for x, y in zip(a, b))


def test_apl_residual_spec_and_tag():
    assert E10.apl_residual_spec(None) is None and E10.apl_residual_tag(None) == ''
    s = E10.apl_residual_spec(0.05)
    assert s == dict(value=0.05, fitted=False, bounds=None) and E10.apl_residual_tag(s) == 'aplres0.05_'
    assert E10.apl_residual_tag(E10.apl_residual_spec(0.0)) == 'aplres0_'
    assert E10.apl_residual_tag(E10.apl_residual_spec(1.0)) == 'aplres1_'
    s = E10.apl_residual_spec(0.2, fitted=True)
    assert s == dict(value=0.2, fitted=True, bounds=E10.APL_RESIDUAL_BOUNDS) and E10.apl_residual_tag(s) == 'aplresfit0.2_'
    assert E10.apl_residual_tag(E10.apl_residual_spec(0.2, True, (0.0, 0.3))) == 'aplresfit0.2b0-0.3_'
    for bad in (dict(f=-0.1), dict(f=1.5), dict(f=float('nan')), dict(f=0.6, fitted=True),
                dict(f=None, fitted=True), dict(f=0.1, bounds=(0.0, 0.5)), dict(f=0.1, fitted=True, bounds=(0.3, 0.2))):
        with pytest.raises(ValueError):
            E10.apl_residual_spec(bad['f'], bad.get('fitted', False), bad.get('bounds'))


def test_apl_residual_params():
    p = dict(fixed_thr=150.0, wAPLKC_scale=0.03, wKCAPL_scale=2e-3, wAPLPN_scale=0.4, wPNAPL_scale=1e-2, wKCKC_scale=26.0)
    q0 = E10.apl_residual_params(p, 0.0)
    assert q0 == dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0) and list(q0) == list(p)
    assert all(type(q0[k]) is float for k in ('wAPLKC_scale', 'wAPLPN_scale'))
    assert E10.apl_residual_params(p, 1.0) == p
    q = E10.apl_residual_params(p, 0.1)
    assert q['wAPLKC_scale'] == 0.1 * 0.03 and q['wAPLPN_scale'] == 0.1 * 0.4
    assert {k: v for k, v in q.items() if k not in ('wAPLKC_scale', 'wAPLPN_scale')} == {
        k: v for k, v in p.items() if k not in ('wAPLKC_scale', 'wAPLPN_scale')}


def test_fit_one_call_unchanged_without_the_option(monkeypatch):
    """apl_residual=None: optimize_mb_model_params gets exactly the call without the kwarg."""
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
    E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw, apl_residual=None)
    a, b = got
    assert set(a) == set(b)
    for k, v in a.items():
        assert b[k] is v or (isinstance(v, (pd.Series, pd.DataFrame)) and b[k].equals(v)) or b[k] == v, k
    assert a['condition_param_overrides'] == {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0}}
    assert 'condition_hooks' not in a and 'extra_params' not in a
    assert inspect.signature(E10.fit_one).parameters['apl_residual'].default is None


def test_identity_f0_equals_the_blocked_condition(monkeypatch):
    """f = 0: every simulation, the history, the losses and the best outputs as without the option, exactly."""
    r0, c0, s0 = _fit(monkeypatch)
    r1, c1, s1 = _fit(monkeypatch, apl_residual=E10.apl_residual_spec(0.0))
    assert len(c0) == len(c1) > 4 and c0 == c1
    # the blocked condition's calls: the static override's 0.0
    assert all(c['wAPLKC_scale'] == 0.0 and c['wAPLPN_scale'] == 0.0 for c in c1[1::2])
    h0, h1 = r0['history'].drop(columns='time_s'), r1['history'].drop(columns='time_s')
    assert (h1[E10.APL_RESIDUAL_KEY] == 0.0).all()
    pd.testing.assert_frame_equal(h0, h1.drop(columns=E10.APL_RESIDUAL_KEY), check_exact=True)
    assert r0['loss'] == r1['loss'] and r0['losses'] == r1['losses'] and r0['params'] == r1['params']
    for c in ('on', 'off'):
        assert _same_frames(r0['condition_outputs'][c], r1['condition_outputs'][c])
    assert r1['condition_params']['off'] == r0['condition_params']['off']
    assert all(_same_frames(a[c], b[c]) for a, b in zip(s0, s1) for c in ('on', 'off'))
    assert r1['apl_residual'] == dict(value=0.0, fitted=False, initial=0.0, bounds=None, off_wAPLKC_scale=0.0,
                                      off_wAPLPN_scale=0.0)


def test_identity_f1_blocked_equals_apl_on(monkeypatch):
    """f = 1: every evaluation's blocked run has the APL-on parameters and equals the APL-on run."""
    r, calls, seen = _fit(monkeypatch, apl_residual=E10.apl_residual_spec(1.0))
    assert len(calls) == 2 * len(seen) and all(a == b for a, b in zip(calls[0::2], calls[1::2]))
    assert all(_same_frames(o['on'], o['off']) for o in seen)
    assert _same_frames(r['condition_outputs']['on'], r['condition_outputs']['off'])
    assert r['condition_params']['off'] == {k: v for k, v in r['params'].items()}


def test_fitted_residue(monkeypatch):
    """A fitted f: an extra parameter within its bounds, in the history; each blocked run at that evaluation's f."""
    spec = E10.apl_residual_spec(0.05, fitted=True)
    r, calls, seen = _fit(monkeypatch, maxfev=20, apl_residual=spec)
    h = r['history']
    f = h[E10.APL_RESIDUAL_KEY].to_numpy()
    assert f[0] == 0.05 and len(np.unique(f)) > 1 and ((f >= 0.0) & (f <= 0.5)).all()
    on, off = calls[0::2], calls[1::2]
    assert len(on) == len(off) == len(h)
    for a, b, fi in zip(on, off, f):
        assert 'apl_residual' not in a and 'apl_residual' not in b
        assert b == E10.apl_residual_params(a, fi)
    best = h.loc[h['loss'].idxmin()]
    assert r['apl_residual']['value'] == float(best[E10.APL_RESIDUAL_KEY]) == r['params'][E10.APL_RESIDUAL_KEY]
    assert r['apl_residual']['fitted'] and r['apl_residual']['bounds'] == [0.0, 0.5] and r['apl_residual']['initial'] == 0.05
    model_p = {k: v for k, v in r['params'].items() if k != E10.APL_RESIDUAL_KEY}
    assert r['condition_params']['off'] == E10.apl_residual_params(model_p, r['apl_residual']['value'])
    assert r['apl_residual']['off_wAPLKC_scale'] == r['condition_params']['off']['wAPLKC_scale']


def test_fitted_residue_beside_kckc_held(monkeypatch):
    r, calls, _ = _fit(monkeypatch, maxfev=10, kckc_fixed=0.0, apl_residual=E10.apl_residual_spec(0.2, fitted=True))
    assert all(c['wKCKC_scale'] == 0.0 for c in calls)
    assert 'wKCKC_scale' not in r['params_to_optimize'] and r['params_to_optimize'][-1] == E10.APL_RESIDUAL_KEY


def test_refusals(monkeypatch):
    initial, wnorm, breadth = _fake_base()
    kw = dict(kckc0=60.0, maxfev=1, method='Nelder-Mead', thr0_factor=0.7, verbose=False,
              apl_residual=E10.apl_residual_spec(0.1))
    proxy = dict(sigma=0.1, thr_ratio=1.0, n_recordings=1, n_trials=1, seed=0, passes=1, relax=False)
    for extra in (dict(kc_thr_factors=np.ones((2, 12))), dict(block_proxy=proxy), dict(kc_thr_factor_fixed=np.ones(12)),
                  dict(spike_times=True)):
        with pytest.raises(NotImplementedError):
            E10.fit_one(2.0, None, initial, wnorm, breadth, lambda o: {}, **kw, **extra)
    with pytest.raises(ValueError):
        E10.fit_one(2.0, None, initial, wnorm, breadth, lambda o: {}, **dict(kw, apl_residual=dict(value=0.1)))


def test_main_and_the_breadth_runs():
    """main passes the option to fit_one only; the breadth run and its pre-passes never see it (APL's output a true 0)."""
    src = inspect.getsource(E10.main)
    for s in ("'--apl-residual'", "'--fit-apl-residual'", "'--apl-residual-bounds'", 'apl_residual_tag(residual)',
              "**({} if residual is None else dict(apl_residual=residual))"):
        assert s in src, s
    for f in (E10.bouton_base, E10.input_scaled_prepass, E10.transfer_prepass):
        assert 'residual' not in inspect.getsource(f)
    assert "off = dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0)" in inspect.getsource(E10.input_scaled_prepass)
    assert "off = dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0)" in inspect.getsource(E10.transfer_prepass)


# ---------------------------------------------------------------------------
# experiments/E111_residual_apl.py

E = pytest.importorskip('experiments.E111_residual_apl')
E110 = pytest.importorskip('experiments.E110_no_kckc')
W = list(E.WIRINGS)


def _changes(**pop_and_wirings):
    """A made-up `changes_table` (index E110's measures, columns 'pop', 'mean_paired' and one per wiring): every change
    0 unless given as key=(population change, change on each of the first n wirings, n)."""
    t = pd.DataFrame(0.0, index=[m.key for m in E110.MEASURES], columns=['pop', 'mean_paired'] + W)
    for key, (pop, per, n) in pop_and_wirings.items():
        t.loc[key, 'pop'] = pop
        t.loc[key, W[:n]] = per
        t.loc[key, 'mean_paired'] = float(t.loc[key, W].mean())
    return t


def _bands(outside: int = 0):
    b = pd.DataFrame(True, index=W, columns=[x.key for x in E110.BANDS])
    b.iloc[:outside, 0] = False
    return b


def _f(value=0.05, **per):
    f = {w: value for w in W}
    f.update(per)
    return f


def test_verdict_replaces():
    v = E.verdict(_changes(), _changes(), _bands(), _f(), 'needed')
    assert v['cls'] == E.CLASSES[0] and v['margins_hold'] and v['f_inside'] and v['n_f_at_bound'] == 0
    # small changes inside every margin still replace
    v = E.verdict(_changes(pattern_cal=(0.002, 0.0025, 8), rho_on=(-0.004, -0.004, 8)), _changes(), _bands(), _f(),
                  'not determined')
    assert v['cls'] == E.CLASSES[0]


def test_verdict_helps_but_does_not_replace():
    # the pattern part worse than the family's on the population and every wiring; better than E110's on 7 of 8
    gain = dict(pattern_cal=(-0.004, -0.004, 7))
    v = E.verdict(_changes(pattern_cal=(0.01, 0.01, 8)), _changes(**gain), _bands(), _f(), 'needed')
    assert v['cls'] == E.CLASSES[1] and 'pattern_cal' in v['margins']['not_needed_fails'] and 'pattern_cal' in v['gains_vs_b']
    # a ρ higher by 0.010 against E110's on 6 of 8 is a gain too
    v = E.verdict(_changes(rho_off=(-0.02, -0.02, 8)), _changes(rho_on=(0.012, 0.012, 6)), _bands(), _f(), 'needed')
    assert v['cls'] == E.CLASSES[1]
    # margins hold but f sits at a bound on 3 of 8 (not "replaces"); the gain against E110's makes it "helps"
    v = E.verdict(_changes(), _changes(**gain), _bands(), _f(0.05, real=0.5, P1=0.4999, P2=0.0), 'needed')
    assert v['n_f_at_bound'] == 3 and not v['f_inside'] and v['cls'] == E.CLASSES[1]
    # a band outside on one run breaks "replaces"
    v = E.verdict(_changes(), _changes(**gain), _bands(outside=1), _f(), 'needed')
    assert v['cls'] == E.CLASSES[1]


def test_verdict_does_nothing():
    # a margin fails and no gain against E110's on enough wirings (5 of 8)
    v = E.verdict(_changes(loss_raw=(0.02, 0.02, 8)), _changes(pattern_cal=(-0.004, -0.004, 5)), _bands(), _f(), 'needed')
    assert v['cls'] == E.CLASSES[2] and not v['gains_vs_b']
    # f below 0.01 on 6 of 8 decides it whatever the margins
    small = {w: 0.005 for w in W[:6]}
    v = E.verdict(_changes(), _changes(pattern_cal=(-0.01, -0.01, 8)), _bands(), _f(0.05, **small), 'needed')
    assert v['cls'] == E.CLASSES[2] and v['n_f_below'] == 6
    # 0.01 itself is not below
    v = E.verdict(_changes(), _changes(), _bands(), _f(0.05, **{w: 0.01 for w in W[:6]}), 'needed')
    assert v['cls'] == E.CLASSES[0]


def test_verdict_without_a_class():
    for cls in ('not needed', 'better without', None):
        v = E.verdict(_changes(), _changes(), _bands(), _f(), cls)
        assert v['cls'] is None and v['margins_hold']


def _values(fam_vals: dict, c_vals: dict):
    rows = []
    for fam, vals in (('family', fam_vals), ('C', c_vals)):
        for i, w in enumerate(W):
            rows.append(dict(family=fam, unit=w, **{k: (v[i] if isinstance(v, list) else v) for k, v in vals.items()}))
    return pd.DataFrame(rows)


IV = dict(block=dict(lo=0.035, hi=0.075), block_ps=dict(lo=-0.043, hi=-0.003), rate_off=dict(lo=0.148, hi=0.164),
          opk_off_tv=dict(lo=0.0, hi=0.22))


def test_told_apart():
    inside = dict(block=0.05, block_ps=-0.02, rate_off=0.155, opk_off_tv=0.1)
    # the residue family's block outside on 6 wirings, the family's inside everywhere: told apart, the family favoured
    c = dict(inside, block=[0.01] * 6 + [0.05] * 2)
    t = E.told_apart(_values(inside, c), IV)
    assert t['told'] and t['cls'] == E.TOLD[0] and t['measures_apart'] == ['block'] and t['counts']['block'] == 6
    assert t['favoured'] == {'block': 'family'} and t['n_needed'] == 6 and t['n_wirings'] == 8
    # on 5 only: not told apart
    c = dict(inside, block=[0.01] * 5 + [0.05] * 3)
    t = E.told_apart(_values(inside, c), IV)
    assert not t['told'] and t['cls'] == E.TOLD[1] and t['counts']['block'] == 5
    # both outside (or both inside) tells nothing; the residue family inside and the family outside favours it
    t = E.told_apart(_values(dict(inside, rate_off=0.13), dict(inside, rate_off=0.12)), IV)
    assert not t['told']
    t = E.told_apart(_values(dict(inside, opk_off_tv=0.3), inside), IV)
    assert t['told'] and t['favoured'] == {'opk_off_tv': 'C'}


def test_f_at_bound():
    assert E.f_at_bound(0.0) and E.f_at_bound(0.005) and E.f_at_bound(0.4951) and E.f_at_bound(0.5)
    assert not E.f_at_bound(0.0051) and not E.f_at_bound(0.2) and not E.f_at_bound(0.494)


def test_run_lists_and_commands():
    c, c8, d = E.fit_specs('C'), E.fit_specs('C', with_t4=False), E.fit_specs('D')
    assert len(c) == 12 and len(c8) == 8 and len(d) == 4
    assert {s.wiring for s in c8} == set(E.S4) and {s.name for s in c if s.part == 'C2'} == {f'C2_{w}' for w in E.S4}
    assert all(s.kckc_fixed == 0.0 and s.start == ('params', f'C_{s.wiring}') for s in c)
    assert {s.f0 for s in c if s.part == 'C1'} == {0.05} and {s.f0 for s in c if s.part == 'C2'} == {0.2}
    assert all(s.kckc_fixed is None and s.f0 == 0.02 and s.start == ('exact',) for s in d)
    s = c[0]
    assert E.fit_tag(s) == 'rate0.068_aplboth_kckc0_aplresfit0.05_k2_cap3_C1'
    cmd = E.fit_command(s)
    i = cmd.index('--apl-residual')
    assert cmd[i + 1] == '0.05' and cmd[i + 2] == '--fit-apl-residual' and '--kckc-fixed' in cmd
    assert cmd[cmd.index('--init-params') + 1] == 'results/E111_residual_apl/starts.json:C_real'
    dc = E.fit_command(d[0])
    assert '--kckc-fixed' not in dc and '--init-from-exact' in dc and dc[dc.index('--apl-residual') + 1] == '0.02'
    assert E.fit_tag(d[0]) == 'rate0.068_aplboth_aplresfit0.02_k2_cap3_D'
    ch = E.check_specs()
    assert ch['f0_real'][0] == 'rate0.068_aplboth_aplres0_k2_cap3_check0'
    assert ch['f1_real'][0] == 'rate0.068_aplboth_aplres1_k2_cap3_check1'
    for name, f in (('f0_real', '0.0'), ('f1_real', '1.0')):
        cmd = ch[name][1]
        assert cmd[cmd.index('--maxfev') + 1] == '1' and cmd[cmd.index('--apl-residual') + 1] == f
        assert '--init-from-exact' in cmd and '--kckc-fixed' not in cmd and '--fit-apl-residual' not in cmd


def test_e110_class_and_the_gate(monkeypatch):
    def with_verdict(cls, final):
        monkeypatch.setattr(E, 'e110_verdict', lambda: dict(cls=cls, final=final))
    t4 = {s.name for s in E.fit_specs('C') if s.wiring in E.T4}
    for cls, final, want_class, run_t4, waits in (('needed', False, 'needed', True, False),
                                                  ('not determined', True, 'not determined', True, False),
                                                  ('not needed', False, None, False, True),
                                                  ('not needed', True, 'not needed', False, False),
                                                  ('better without', True, 'better without', False, False),
                                                  (None, False, None, False, True)):
        with_verdict(cls, final)
        assert E.e110_class() == want_class
        add, why = E.c_gate()()
        assert ({s.name for s in add} == t4) if run_t4 else not add
        assert (why is not None) == waits
        assert E.c_wirings() == (E.S4 if (final and cls in E.E110_NO_CLASS) else E.WIRINGS)
    with_verdict('not needed', True)
    assert {s.name for s in E.c_gate(force=True)()[0]} == t4


def test_block_pair_specific_matches_the_block():
    idx = pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(ODORS) for b in ODORS[i + 1:]], names=['odor_a', 'odor_b'])
    rng = np.random.default_rng(3)
    off = pd.Series(rng.uniform(0.3, 0.8, len(idx)), index=idx)
    on = off + 0.15 + rng.normal(0, 0.03, len(idx))
    bp = E.block_pair_specific(on, off)
    st, _ = H.compare_on_off(E.E83._cs(on), E.E83._cs(off))
    assert abs(bp['block_resid'] + float(st['block_fruity_small_iso'])) < 1e-12 and np.isfinite(bp['block_ps'])


def test_same_outputs_and_row_check():
    s = pd.DataFrame(np.arange(12.0).reshape(4, 3))
    r = (s >= 2).astype(float)
    a = {'on': (r, s), 'off': (r, s)}
    assert E.same_outputs(a, a, 'on', 'off')['exact']
    s2 = s.copy()
    s2.iloc[0, 0] = 1.0
    assert E.same_outputs(a, {'on': (r, s2)}, 'on', 'on')['n_counts_differ'] == 1
    row = dict(tag='x', maxfev=1, loss=1.5, fit_rate_on=0.07, rate_off=0.15, apl_residual=0.0)
    ref = dict(tag='y', maxfev=150, loss=1.5, fit_rate_on=0.07, rate_off=0.15)
    rc = E.row_check(row, ref)
    assert rc['exact'] and rc['fit_exact'] and rc['n_columns'] == 3 and rc['n_fit_columns'] == 2
    rc = E.row_check(dict(row, rate_off=0.16), ref)
    assert not rc['exact'] and rc['fit_exact'] and list(rc['differs']) == ['rate_off']
