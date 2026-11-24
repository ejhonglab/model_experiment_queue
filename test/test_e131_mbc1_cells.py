#!/usr/bin/env python3
"""Tests for E131 (MB-C1's two cells beside the inhibition among KCs; experiments/E131_mbc1_cells.py), on made-up data:
the dose arithmetic (reading 1, and on hemibrain KC IDs with olfsysm's loader, skipped without it), the arms' planned
tags and lanes, the runner's builds (the option off is E83's simulator unchanged; the plain cells pass mbc1_inh), the
cells' activity by olfsysm's recursion (reading 4) and its pulse integrals, step 0's two numbers and the pooled gate
(readings 5 to 7), the corner rule (reading 8), the classes on synthetic tables (reading 9), the flags, outcomes and
step 2's gate (readings 10, 11). No model simulation."""

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E131_mbc1_cells')


# ---------------------------------------------------------------------------
# reading 1: the dose

def test_mbc1_scale_and_arm_params():
    p = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.4, wPNAPL_scale=1e-3,
             wKCKC_scale=26.0)
    want = {'B25': (26.0, 6.5, {}), 'B100': (26.0, 26.0, {}), 'B400': (26.0, 104.0, {}), 'R50': (13.0, 13.0, {}),
            'H50': (13.0, 0.0, {}), 'X50': (13.0, 13.0, {'wKCKC_scale': 0.0}), 'D100': (26.0, 26.0, {}),
            'P100': (26.0, 26.0, {}), 'L100': (26.0, 26.0 * 3.3, {})}
    for arm, (kk, mb, off) in want.items():
        q, o = X.arm_params(arm, p, late_ratio=3.3 if arm == 'L100' else None)
        assert q['wKCKC_scale'] == pytest.approx(kk) and q['wMBC1_scale'] == pytest.approx(mb), arm
        assert o == off, arm
        assert {k: q[k] for k in p if k != 'wKCKC_scale'} == {k: v for k, v in p.items() if k != 'wKCKC_scale'}
    with pytest.raises(ValueError):
        X.mbc1_scale('L100', 26.0)
    with pytest.raises(ValueError):
        X.mbc1_scale('L100', 26.0, late_ratio=float('nan'))
    assert set(X.ARMS) == {'B25', 'B100', 'B400', 'R50', 'H50', 'X50', 'L100', 'D100', 'P100'}


def test_form_kws():
    st = dict(threshold={'real': 3.7}, late_ratio={'real': 3.2})
    assert X.form_kws('plain', 'real') == {} and X.form_kws('off', 'real') == {}
    assert X.form_kws('late', 'real') == dict(mbc1_tau=1.0)
    assert X.form_kws('pooled', 'real') == dict(mbc1_pooled=True)
    assert X.form_kws('strong', 'real', st) == dict(mbc1_thr=3.7)
    with pytest.raises(ValueError):
        X.form_kws('strong', 'real')
    with pytest.raises(ValueError):
        X.form_kws('other', 'real')


def test_uniform_activity_totals():
    """Reading 1 on made-up weights normalized as olfsysm's loader does: each cell's input weights sum to 1 and the
    equivalent matrix (output x input, summed over the cells) has total N, so at m 1 the cells' total output for uniform
    activity is the KC-KC matrix's (mean row sum 1)."""
    rng = np.random.default_rng(0)
    n, c = 40, 2
    out_c = rng.poisson(3.0, (n, c)).astype(float)
    in_c = rng.poisson(2.0, (n, c)).astype(float)
    tot = float((out_c.sum(0) * in_c.sum(0)).sum())
    w_in = in_c / in_c.sum(0)
    w_out = out_c * (in_c.sum(0) / (tot / n))
    E = w_out @ w_in.T
    assert E.sum() == pytest.approx(n)
    assert np.allclose(E, out_c @ in_c.T / (tot / n))
    hb = rng.random((n, n))
    hb *= n / hb.sum()
    a = 7.0
    s, m = 26.0, 1.0
    assert (m * s * w_out @ (w_in.T @ np.full(n, a))).sum() == pytest.approx((s * hb @ np.full(n, a)).sum())


def _mbc1_or_skip():
    mb = pytest.importorskip('pyolfsysm.mb_model')
    if not (hasattr(mb, 'connectome_MBC1_cell_weights') and all(p.exists() for p in mb.mbc1_synapse_paths.values())):
        pytest.skip('pyolfsysm without the MB-C1 cells (olfsysm >= fa3543d) or the hemibrain MB-C1 tables')
    return mb


def _hemibrain_index(mb, n=300, seed=0):
    w_out, w_in, _, _ = mb.connectome_MBC1_weights()
    with_mb = w_out.index.union(w_in.index)
    kk = mb.connectome_KCKC_weights(normalize=False).index
    rng = np.random.default_rng(seed)
    ids = np.unique(np.concatenate([with_mb.to_numpy()[:n], kk.difference(with_mb).to_numpy()[:30]]))
    ids = ids[rng.permutation(ids.size)]
    return pd.MultiIndex.from_arrays([ids, rng.choice(list(X.TYPES), ids.size)], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])


def test_dose_check_hemibrain():
    mb = _mbc1_or_skip()
    idx = _hemibrain_index(mb)
    d = X.dose_check(idx)
    assert d['ok'] and d['cells']['ok'] and d['pooled']['ok']
    assert d['cells']['n_cells'] == 2 and d['pooled']['n_cells'] == 1
    assert d['cells']['equiv_total'] == pytest.approx(len(idx)) and d['hb_total'] == pytest.approx(len(idx))
    o, i = X.cell_weights(idx)
    assert np.array_equal(o.index.to_numpy(), idx.get_level_values(H.KC_ID_LEVEL).to_numpy())
    oc, ic = X.cell_weights(idx, normalize=False)
    assert (oc.to_numpy() == np.round(oc.to_numpy())).all() and (ic.to_numpy() >= 0).all()


def test_dose_arithmetic():
    syn = pd.DataFrame([
        dict(kc_type='all', n_kcs=10, kckc_in=300.0, mbc1_out=120.0, mbc1_in=50.0, mbc1_out_c1=60.0, mbc1_in_c1=20.0,
             mbc1_out_c2=60.0, mbc1_in_c2=30.0, apl_out_claws=150.0, apl_out_calyx=250.0, mbc1_out_per_kc=12.0,
             apl_out_claws_per_kc=15.0, apl_out_calyx_per_kc=25.0, mbc1_out_share=1.0, mbc1_in_share=1.0),
        dict(kc_type='g', n_kcs=4, kckc_in=100.0, mbc1_out=100.0, mbc1_in=40.0, mbc1_out_c1=50.0, mbc1_in_c1=16.0,
             mbc1_out_c2=50.0, mbc1_in_c2=24.0, apl_out_claws=60.0, apl_out_calyx=100.0, mbc1_out_per_kc=25.0,
             apl_out_claws_per_kc=15.0, apl_out_calyx_per_kc=25.0, mbc1_out_share=100 / 120, mbc1_in_share=0.8)])
    a = X.dose_arithmetic(syn)
    assert a['m_equal_synapse'] == pytest.approx(120 / 300)
    pairs = 60 * 20 + 60 * 30
    assert a['m_equal_synapse_per_cell'] == pytest.approx({'c1': pairs / (20 * 300), 'c2': pairs / (30 * 300)})
    assert a['g_share_mbc1_vs_apl_claws'] == pytest.approx(100 / 160)
    assert a['g_share_mbc1_vs_apl_calyx'] == pytest.approx(100 / 200)
    # the synapse-weighted mean of the per-cell doses' reciprocals is the average's
    w = np.array([60.0, 60.0]) / 120
    inv = np.array([1 / a['m_equal_synapse_per_cell'][c] for c in ('c1', 'c2')])
    assert 1 / (w @ inv) == pytest.approx(a['m_equal_synapse'])


# ---------------------------------------------------------------------------
# the runner's builds

def test_simulator_builds(monkeypatch):
    seen = []

    def fake(w, params=None, extra_kws=None, **kw):
        seen.append((w, dict(params), None if extra_kws is None else dict(extra_kws)))
        return (lambda q: None), {}

    monkeypatch.setattr(X.E83, 'build_simulator', fake)
    p = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.4, wPNAPL_scale=1e-3,
             wKCKC_scale=26.0, wMBC1_scale=6.5)
    X.simulator('real', 'off', p)
    X.simulator('real', 'plain', p)
    X.simulator('real', 'late', p)
    X.simulator('P1', 'strong', p, dict(threshold={'P1': 3.3}))
    X.simulator('real', 'pooled', p)
    off, plain, late, strong, pooled = seen
    assert off[2] is None and 'wMBC1_scale' not in off[1]          # the option off: E83's simulator as it is
    assert plain[2] == dict(mbc1_inh=6.5) and plain[1]['wMBC1_scale'] == 6.5
    assert late[2] == dict(mbc1_inh=6.5, mbc1_tau=1.0)
    assert strong[2] == dict(mbc1_inh=6.5, mbc1_thr=3.3) and strong[0] == 'P1'
    assert pooled[2] == dict(mbc1_inh=6.5, mbc1_pooled=True)


def test_build_runner_cache(monkeypatch):
    n = []

    def fake(w, form, params, step0=None):
        n.append((w, form))
        return (lambda q: None), dict(w=w, form=form)

    monkeypatch.setattr(X, 'simulator', fake)
    X._CACHE.pop('runner', None)
    p = {}
    X.build_runner('B25', 'real', p)
    X.build_runner('B400', 'real', p)
    X.build_runner('H50', 'real', p)
    X.build_runner('X50', 'real', p)
    assert n == [('real', 'plain')]
    X.build_runner('B25', 'P1', p)
    X.build_runner('P100', 'P1', p)
    X.build_runner('B25', 'P1', p)
    assert n == [('real', 'plain'), ('P1', 'plain'), ('P1', 'pooled'), ('P1', 'plain')]
    X._CACHE.pop('runner', None)


def test_ndiff_identity():
    idx = pd.Index(range(5))
    r = pd.DataFrame(np.eye(5)[:, :3], index=idx)
    s = pd.DataFrame(np.arange(15.0).reshape(5, 3), index=idx)
    assert X._ndiff((r, s), (r.copy(), s.copy())) == dict(counts=0, responses=0)
    s2 = s.copy()
    s2.iloc[1, 2] += 1
    assert X._ndiff((r, s2), (r, s)) == dict(counts=1, responses=0)


# ---------------------------------------------------------------------------
# reading 4: the cells' activity

def _spikes(rng, n_kcs=20, n_odors=3, n=60, steps=None):
    steps = X.n_steps() if steps is None else steps
    p0, p1 = X.pulse_steps()
    return pd.DataFrame(dict(kc=rng.integers(0, n_kcs, n), odor=rng.integers(0, n_odors, n),
                             step=rng.integers(p0 - 50, min(p1 + 50, steps - 1), n)))


def test_activity_trace_is_olfsysms_recursion():
    rng = np.random.default_rng(0)
    spk = _spikes(rng)
    w_in = rng.random((2, 20))
    w_in /= w_in.sum(axis=1, keepdims=True)
    tr = X.activity_trace(spk, w_in, 3, 0.05)
    assert tr.shape == (2, 3, X.n_steps()) and not tr[:, :, 0].any()
    # by hand for one cell and odor
    c, o = 1, 2
    x = np.zeros(X.n_steps())
    sel = spk[spk['odor'] == o]
    np.add.at(x, sel['step'].to_numpy(), w_in[c, sel['kc'].to_numpy()])
    m, d = 0.0, 1 - X.TIME_DT / 0.05
    want = np.zeros(X.n_steps())
    for t in range(1, X.n_steps()):
        m = d * m + x[t - 1] / 0.05
        want[t] = m
    assert np.allclose(tr[c, o], want, rtol=1e-12, atol=0)


def test_cell_activity_and_pulse_integral():
    """With thr 0 the pulse output of a cell is its input weights times each KC's pulse integral (`pulse_integral`),
    so the cells' linear output is the equivalent matrix applied to the KCs' filtered activity (reading 1)."""
    rng = np.random.default_rng(1)
    spk = _spikes(rng, n=200)
    w_in = rng.random((2, 20))
    a = X.cell_activity(spk, w_in, 3, 0.05)
    I = X.pulse_integral(spk, 20, 3, 0.05)
    assert np.allclose(a['pulse_output'], w_in @ I, rtol=1e-9)
    p0, p1 = X.pulse_steps()
    assert np.allclose(a['pulse_mean'], a['pulse_output'] / ((p1 - p0) * X.TIME_DT))
    # a threshold above the peak silences the output, one below it does not
    hi = X.cell_activity(spk, w_in, 3, 0.05, thr=float(a['peak'].max()))
    assert not hi['pulse_output'].any()
    lo = X.cell_activity(spk, w_in, 3, 0.05, thr=0.5 * float(a['peak'].max()))
    assert lo['pulse_output'].sum() > 0 and (lo['pulse_output'] <= a['pulse_output'] + 1e-12).all()
    # one spike well inside the pulse integrates to ~1 over the pulse
    one = pd.DataFrame(dict(kc=[0], odor=[0], step=[p0 + 10]))
    assert X.pulse_integral(one, 1, 1, 0.05)[0, 0] == pytest.approx(1.0, abs=1e-4)
    assert X.pulse_integral(one, 1, 1, 1.0)[0, 0] < 0.5


def test_check_spikes():
    rng = np.random.default_rng(2)
    spk = _spikes(rng)
    spk['time_s'] = X.TIME_PRE_START + spk['step'] * X.TIME_DT
    C = np.zeros((20, 3))
    np.add.at(C, (spk['kc'], spk['odor']), 1)
    X.check_spikes(spk, pd.DataFrame(C))
    C[0, 0] += 1
    with pytest.raises(RuntimeError):
        X.check_spikes(spk, pd.DataFrame(C))
    C[0, 0] -= 1
    bad = spk.copy()
    bad['time_s'] += 0.01
    with pytest.raises(RuntimeError):
        X.check_spikes(bad, pd.DataFrame(C))


# ---------------------------------------------------------------------------
# readings 5 to 7: step 0's numbers and the pooled gate

def test_late_ratio():
    a50 = dict(pulse_mean=np.array([[2.0, 4.0], [1.0, 3.0]]))
    a1 = dict(pulse_mean=np.array([[1.0, 1.0], [0.5, 0.5]]))
    out_w = np.array([2.0, 1.0])
    assert X.late_ratio(a50, a1, out_w) == pytest.approx((2 * 6 + 1 * 4) / (2 * 2 + 1 * 1))
    # the late cells at the ratio deliver the plain cells' summed pulse output
    r = X.late_ratio(a50, a1, out_w)
    assert r * out_w @ a1['pulse_mean'].sum(1) == pytest.approx(out_w @ a50['pulse_mean'].sum(1))
    with pytest.raises(ValueError):
        X.late_ratio(a50, dict(pulse_mean=np.zeros((2, 2))), out_w)


def test_cells_agreement_and_pooled_gate():
    x = np.arange(16.0)
    assert X.cells_agreement(np.vstack([x, 2 * x + 1])) == pytest.approx(1.0)
    assert np.isnan(X.cells_agreement(x[None]))
    g = X.pooled_gate({'real|on': 0.99, 'real|off': 0.96, 'P1|on': 0.951})
    assert not g['open'] and g['below'] == {}
    g = X.pooled_gate({'real|on': 0.99, 'P3|off': 0.949})
    assert g['open'] and list(g['below']) == ['P3|off'] and g['min_r'] == pytest.approx(0.949)
    assert X.pooled_gate({'real|on': float('nan')})['open']


def test_gamma_input():
    idx = pd.MultiIndex.from_arrays([[1, 2, 3, 4], ['g', 'ab', 'g', "a'b'"]], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    C = pd.DataFrame([[1.0, 0.0], [2.0, 1.0], [0.0, 3.0], [1.0, 1.0]], index=idx)
    E = np.full((4, 4), 0.25)
    hb = np.eye(4)[::-1]
    g = X.gamma_input(C, 2.0, E, hb)
    assert list(g['kc_id']) == [1, 3]
    assert np.allclose(g['mbc1_input'], 2.0 * 0.25 * C.to_numpy().sum() * np.ones(2))
    assert np.allclose(g['kckc_input'], [2.0 * 2.0, 2.0 * 3.0])


# ---------------------------------------------------------------------------
# the arms' tags and lanes

def test_planned_tags():
    pts = X.planned()
    assert len(pts) == 72 and len(set(pts)) == 72
    tags = [X.tag_of(a, w) for a, w in pts]
    assert tags[:3] == ['S1_B25_real', 'S1_B100_real', 'S1_B400_real']
    assert tags[24:26] == ['S1_R50_real', 'S1_H50_real'] and tags[40] == 'S1_X50_real'
    assert tags[48:50] == ['S1_L100_real', 'S1_L100_P1'] and tags[-1] == 'S1_P100_M3'
    groups = [X.ARMS[a]['group'] for a, _ in pts]
    assert groups == sorted(groups, key=X.GROUPS.index)
    t = X.planned(X.T_WIRINGS, ['B25', 'X50'])
    assert t == [('B25', 'P3'), ('B25', 'P4'), ('B25', 'M2'), ('B25', 'M3'), ('X50', 'P3'), ('X50', 'P4'),
                 ('X50', 'M2'), ('X50', 'M3')]
    assert sum(1 for a, _ in pts if a != 'P100') == 64


def test_runnable_arms():
    base = ('B25', 'B100', 'B400', 'R50', 'H50', 'X50')
    assert X.runnable_arms(None) == base
    st = dict(pooled_gate=dict(open=False), threshold={}, late_ratio={})
    assert X.runnable_arms(st) == base + ('L100',)
    assert X.runnable_arms(st, set(X.WIRINGS)) == base + ('L100', 'D100')
    assert X.runnable_arms(st, set(X.WIRINGS) - {'M3'}) == base + ('L100',)
    st['pooled_gate']['open'] = True
    assert X.runnable_arms(st, set(X.WIRINGS)) == base + ('L100', 'D100', 'P100')


def test_claims(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    assert X.claim('S1_B25_real', 'lane')
    assert not X.claim('S1_B25_real', 'other')
    assert X.claim('S1_B25_real', 'other', reclaim=True)
    X.release('S1_B25_real')
    assert X.claim('S1_B25_real', 'lane')


def test_runs_refused_without_checks(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    with pytest.raises(SystemExit):
        X.part_runs(['real'])


# ---------------------------------------------------------------------------
# readings 8 to 11: corners, classes, flags, outcomes, the gate

def test_held_ok_and_corner_rule():
    assert X.held_ok(3540, 3540) and X.held_ok(3547, 3540) and not X.held_ok(3548, 3540)
    runs = pd.DataFrame([dict(kind='a', wiring=w, corner=w in ('real', 'P1')) for w in X.WIRINGS]
                        + [dict(kind='b', wiring=w, corner=w in ('real', 'P1', 'M3')) for w in X.WIRINGS])
    assert X.corners(runs, 'a') == dict(corner_wirings=['P1', 'real'], corner=False)
    assert X.corners(runs, 'b')['corner']


def _rd_runs(values, out_bands=None, corner=None, mean=0.178, raw=None, seed=0):
    """Made-up readings: per kind a calibrated pattern part (every target, each wiring and population), its mean APL
    effect under the readout and its raw-count loss; runs with the bands and corners given."""
    rng = np.random.default_rng(seed)
    rows, runs = [], []
    raw = raw or {}
    for kind, base in values.items():
        for unit in list(X.WIRINGS) + list(X.POPS) + ['all8']:
            r = dict(kind=kind, unit=unit, readout='cal', mean_effect=mean if kind != 'family' else 0.178)
            for t in X.targets():
                r[f'{t}_pattern'] = base + rng.normal(0, 1e-4)
            for k in ('amplitude', 'shape', 'sd_ratio', 'pearson_effect', 'rho_effect'):
                r[f'all_{k}'] = 0.5
            rows.append(r)
            rows.append(dict(kind=kind, unit=unit, readout='raw', mean_effect=0.12, all_loss=raw.get(kind, 0.1)))
        for w in X.WIRINGS:
            runs.append(dict(kind=kind, wiring=w, rate_rho_on=0.5, in_bands=not (out_bands and w in out_bands.get(kind, ())),
                             corner=bool(corner and w in corner.get(kind, ())), fixed_thr=100.0, thr_below_60=False))
    return pd.DataFrame(rows), pd.DataFrame(runs)


def test_classes():
    pytest.importorskip('experiments.E123_kckc_form')
    fam, z = 0.048, 0.063
    rd, runs = _rd_runs({'family': fam, 'E110_A': z, 'eq': 0.050, 'lo': 0.044, 'hi': 0.052, 'mid': 0.0505,
                         'cor': 0.050, 'lo_mean': 0.044, 'lo_raw': 0.044, 'lo_band': 0.044},
                        out_bands={'lo_band': ('P2',), 'hi': ('real', 'P1')}, corner={'cor': ('real', 'P1', 'P3')},
                        raw={'lo_raw': 0.106})
    c = {k: X.classify(rd, runs, k, 0.178) for k in ('eq', 'lo', 'hi', 'mid', 'cor', 'lo_raw', 'lo_band')}
    assert c['eq']['cls'] == 'equal' and c['eq']['keeps_gain']
    assert c['lo']['cls'] == 'below the family' and c['lo']['below']['below']
    assert c['hi']['cls'] == 'costs' and c['hi']['costs'] and c['hi']['runs_outside_bands'] == ['P1', 'real']
    assert c['mid']['cls'] == 'equal'
    assert c['cor']['cls'] == 'corner (no class)' and c['cor']['cls_before_corner'] == 'equal'
    assert c['lo_raw']['cls'] == 'equal' and not c['lo_raw']['below']['raw_loss_ok']
    assert c['lo_band']['cls'] == 'neither' and c['lo_band']['below']['pattern_lower'] and not c['lo_band']['equal']
    rd2, runs2 = _rd_runs({'family': fam, 'E110_A': z, 'lo': 0.044}, mean=0.21)
    c2 = X.classify(rd2, runs2, 'lo', 0.178)
    assert c2['cls'] == 'equal' and not c2['below']['mean_effect_ok']
    # "costs" against all recordings only, on both populations
    rd3, runs3 = _rd_runs({'family': fam, 'E110_A': z, 'k': 0.0515})
    rd3.loc[(rd3['kind'] == 'k') & (rd3['unit'] == 'T4') & (rd3['readout'] == 'cal'), 'all_pattern'] = 0.050
    assert X.classify(rd3, runs3, 'k', 0.178)['cls'] == 'neither'


def test_flags_outcomes_and_gate():
    def c(cls, keeps=False, equal=None, corner=False):
        return dict(complete=True, cls=cls, keeps_gain=keeps, equal=cls == 'equal' if equal is None else equal,
                    corner=corner, corner_wirings=[], runs_outside_bands=[])

    runs = pd.DataFrame([dict(kind=a, wiring=w, in_bands=True) for a in X.ARMS for w in X.WIRINGS])
    # (a): every B dose costs, R50 does not stand in, X50 does not keep the gain
    cl = {'B25': c('costs'), 'B100': c('costs'), 'B400': c('costs'), 'R50': c('equal'), 'H50': c('equal'),
          'X50': c('neither')}
    comps = {'R50|H50': dict(verdict='not told apart')}
    fl = X.flags(cl, comps)
    assert fl['R50'] == dict(keeps_gain=False, stands_in=False)
    o = X.outcomes(cl, comps, fl, runs)
    assert o['a'] is True and o['b'] is False and o['c'] is False and o['d_b_below'] == [] and not o['f']
    assert not X.step2_gate(cl, fl)['open']
    # (b): B25 equal, B100 and B400 cost; (c) R50 stands in
    cl.update(B25=c('equal'), R50=c('equal', keeps=True))
    comps = {'R50|H50': dict(verdict='R50 lower')}
    fl = X.flags(cl, comps)
    assert fl['R50']['stands_in']
    o = X.outcomes(cl, comps, fl, runs)
    assert o['a'] is False and o['b'] is True and o['largest_tolerated_dose'] == 0.25 and o['c'] is True
    g = X.step2_gate(cl, fl)
    assert g['open'] and g['opened_by'] == ['R50 stands in']
    # B25 equal, B100 neither: not (b)
    cl.update(B100=c('neither'))
    assert X.outcomes(cl, comps, X.flags(cl, comps), runs)['b'] is False
    # (d) and the gate: a dose below the family; X50 keeps the gain; a corner keeps nothing
    cl.update(B100=c('below the family'), X50=c('neither', keeps=True), R50=c('equal', keeps=True, corner=True))
    fl = X.flags(cl, comps)
    assert not fl['R50']['stands_in'] and fl['X50']['keeps_gain']
    o = X.outcomes(cl, comps, fl, runs)
    assert o['d_b_below'] == ['B100'] and o['d'] is None and o['c'] is True
    assert X.step2_gate(cl, fl)['opened_by'] == ['B100', 'X50 keeps the gain']
    # (e): a form against the plain cells
    comps = {'L100|B100': dict(verdict='L100 lower'), 'D100|B100': dict(verdict='not told apart')}
    cl.update(L100=c('equal'), D100=c('equal'))
    o = X.outcomes(cl, comps, X.flags(cl, comps), runs)
    assert o['e'] == {'L100': 'L100 lower', 'D100': 'not told apart'} and o['e_lower'] == ['L100']
    # (f): a run outside the bands
    runs.loc[(runs['kind'] == 'B400') & (runs['wiring'] == 'M2'), 'in_bands'] = False
    assert X.outcomes(cl, comps, X.flags(cl, comps), runs)['f']
    # incomplete arms give None
    o = X.outcomes({'B25': c('equal')}, {}, {}, runs)
    assert o['a'] is None and o['b'] is None and o['c'] is None


def test_differ_is_e129s():
    E129 = pytest.importorskip('experiments.E129_block_survivor_static')
    for vals in ({'a': 0.050, 'b': 0.055}, {'a': 0.0551, 'b': 0.0550}, {'a': 0.056, 'b': 0.055}):
        rd, runs = _rd_runs(dict(vals, family=0.048, E110_A=0.063))
        assert X.differ(rd, runs, 'a', 'b')['verdict'] == E129.differ(rd, runs, 'a', 'b')['verdict']
    rd, runs = _rd_runs({'a': 0.050, 'b': 0.055})
    assert X.differ(rd, runs, 'a', 'b')['verdict'] == 'a lower'
