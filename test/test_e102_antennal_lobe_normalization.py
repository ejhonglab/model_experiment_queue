#!/usr/bin/env python3
"""Tests for E102 (the antennal lobe's normalization at the strength the flies' APL asks for;
experiments/E102_antennal_lobe_normalization.py), on made-up rows and runs: the step names and PN kwargs; the three level
rules, with a flagged step and with none existing; the refit classes (equal, better, price) and T's rule for better; the
low-concentration classes; the equal rule's inputs; the refusal to run past a failed (or missing) identity, in
`require_identity` and in part B's runner; part B's runner on made-up held and extended runs (its files in E98's layout,
which E98's loader reads, the check that the extended run's panel columns equal the held run's, E98's amended blank
identity, resuming); the chain of totals and its ratios; the results tables' caveat line; part D's refits: the init
file (read back exactly by E10's `read_params`), the options and the tag E10's own main gives each kind of wiring, the
start identity and the start check's stop on a difference, the refit's parameters from its history, which refits are
due, the per-wiring class counts and the comparison table."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
M = pytest.importorskip('experiments.E102_antennal_lobe_normalization')


# ---------------------------------------------------------------------------
# steps

def test_step_tags_round_trip():
    assert M.step_tag(M.STEP_IDENTITY) == 'identity'
    assert [M.step_tag(a) for a in M.LADDER] == ['inhadd100', 'inhadd314', 'inhadd1000', 'inhadd3141', 'inhadd10000',
                                                 'gain_fixed']
    for a in M.LADDER:
        assert M.step_of_tag(M.step_tag(a)) == a
    assert M.step_of_tag('identity') == M.STEP_IDENTITY
    assert all(c.isascii() and (c.isalnum() or c in '_-.') for a in M.LADDER for c in M.step_tag(a))
    with pytest.raises(ValueError):
        M.step_of_tag('S1_inhadd100')


def test_step_kws():
    k0 = M.step_kws(M.STEP_IDENTITY)
    assert k0 == dict(pn_inhadd=M.PN_INHADD_DEFAULT, pn_inhsc=M.PN_INHSC_DEFAULT)
    from experiments.E77_input_screen import inhsc_for
    for a in M.LADDER:
        k = M.step_kws(a)
        assert k['pn_inhadd'] == a and k['pn_inhsc'] == inhsc_for(a)
    # pn_inhsc keeps the gain at a typical inhibition (E3's 2000): inhsc / (inhadd + 2000) is the default's
    g0 = M.PN_INHSC_DEFAULT / (M.PN_INHADD_DEFAULT + 2000.0)
    for a in M.LADDER:
        k = M.step_kws(a)
        assert k['pn_inhsc'] / (k['pn_inhadd'] + 2000.0) == pytest.approx(g0)
    assert M.step_value(M.STEP_IDENTITY) == M.PN_INHADD_DEFAULT


def test_gain_spread_falls_along_the_ladder():
    g = [M._gain_spread(a) for a in (M.PN_INHADD_DEFAULT,) + M.LADDER]
    assert all(x > y for x, y in zip(g, g[1:]))
    assert g[0] == pytest.approx(np.log((31.4088 + 5600) / (31.4088 + 1300)))
    assert g[-1] < 0.005


# ---------------------------------------------------------------------------
# part C's levels

def _steps(v2, v2_lin, eq, flagged, silent=None):
    vals = (M.PN_INHADD_DEFAULT,) + M.LADDER
    return pd.DataFrame(dict(value=vals, v2=v2, v2_linear=v2_lin, equal_S=eq, flagged=flagged,
                             silent=[False] * len(vals) if silent is None else silent))


def test_levels_basic():
    s = _steps(v2=['flatter', 'flatter', 'flatter', 'inside', 'inside', 'inside', 'steeper'],
               v2_lin=['flatter'] * 4 + ['flatter', 'inside', 'inside'],
               eq=[True, True, True, False, False, False, False], flagged=[False] * 7)
    lv = M.levels(s)
    assert lv['L_apl'] == 1000.0 and lv['L_lin'] == 10000.0 and lv['L_free'] == 314.0 and lv['flagged'] == []
    assert M.refit_levels(lv, s) == [1000.0, 10000.0]


def test_levels_flagged_step_is_not_a_level():
    s = _steps(v2=['flatter', 'flatter', 'flatter', 'inside', 'inside', 'inside', 'inside'],
               v2_lin=['flatter'] * 4 + ['inside'] * 3, eq=[True, True, True, True, False, False, False],
               flagged=[False, False, False, True, True, False, False])
    lv = M.levels(s)
    assert lv['L_apl'] == 10000.0                    # 1000 and 3141 are inside but flagged
    assert lv['L_lin'] == 10000.0
    assert lv['L_free'] == 314.0                     # 1000 is an equal but flagged
    assert lv['flagged'] == [1000.0, 3141.0]
    assert M.refit_levels(lv, s) == [10000.0]        # L_lin equals L_apl: one level


def test_levels_none_exist():
    s = _steps(v2=['flatter'] * 7, v2_lin=['flatter'] * 7, eq=[True] + [False] * 6,
               flagged=[False] * 5 + [True, False])
    lv = M.levels(s)
    assert lv['L_apl'] is None and lv['L_lin'] is None and lv['L_free'] is None   # the identity row is no level
    # without L_apl the refit runs at the last unflagged step (the fixed gain here; 10000 is flagged)
    assert M.refit_levels(lv, s) == [M.FIXED_GAIN_INHADD]
    s2 = s.assign(flagged=[False] + [True] * 6)
    assert M.refit_levels(M.levels(s2), s2) == []
    with pytest.raises(ValueError):
        M.levels(s.drop(columns=['equal_S']))
    with pytest.raises(ValueError):
        M.levels(s.drop(columns=['silent']))


def test_levels_step_with_a_silent_odor_is_not_a_level():
    # the fixed gain inside on V2 and its linear reading, an equal, unflagged, but with an odor silent in one run
    s = _steps(v2=['flatter'] * 6 + ['inside'], v2_lin=['flatter'] * 6 + ['inside'], eq=[True] * 7, flagged=[False] * 7,
               silent=[False] * 6 + [True])
    lv = M.levels(s)
    assert lv['L_apl'] is None and lv['L_lin'] is None and lv['L_free'] == 10000.0
    assert lv['silent'] == [M.FIXED_GAIN_INHADD] and lv['flagged'] == []
    # without L_apl the refit's fallback, the last unflagged step, skips the silent one too
    assert M.refit_levels(lv, s) == [10000.0]
    s2 = s.assign(flagged=[False] * 5 + [True, False])
    assert M.refit_levels(M.levels(s2), s2) == [3141.0]


# ---------------------------------------------------------------------------
# part D's classes

def _dp(z8, a=None, b=None, c=None, d=None):
    return dict(z8=z8, halfA=z8 if a is None else a, halfB=z8 if b is None else b, halfC=z8 if c is None else c,
                halfD=z8 if d is None else d)


def test_refit_class():
    assert M.refit_class(_dp(0.002, -0.001, 0.0029), 0.004, -0.009, 0) == 'equal'
    assert M.refit_class(_dp(-0.0025), 0.0, 0.0, 0) == 'equal'                       # lower, but by less than 0.003
    assert M.refit_class(_dp(-0.004), 0.01, 0.05, 0) == 'better'                     # better is read on the pattern part
    assert M.refit_class(_dp(-0.002, d=0.0001), 0.0, 0.0, 0) == 'equal'              # within 0.003; half D decides nothing
    # lower by more than 0.003 but not lower against half D: neither better nor within 0.003, so the rule as written
    # (E94's, |change| within 0.003) gives a price
    assert M.refit_class(_dp(-0.004, d=0.0001), 0.0, 0.0, 0) == 'price'
    assert M.refit_class(_dp(0.0031), 0.0, 0.0, 0) == 'price'
    assert M.refit_class(_dp(0.001), 0.006, 0.0, 0) == 'price'                       # the raw-count loss
    assert M.refit_class(_dp(0.001), 0.0, 0.011, 0) == 'price'                       # the mean APL effect
    assert M.refit_class(_dp(0.001), 0.0, 0.0, 1) == 'price'                         # a flagged fit
    assert M.refit_class(_dp(-0.01), 0.0, 0.0, 1) == 'price'
    with pytest.raises(ValueError):
        M.refit_class(dict(z8=0.0, halfA=0.0), 0.0, 0.0, 0)


def test_t_better_holds():
    single = {w: -0.001 for w in M.ALL_WIRINGS}
    assert M.t_better_holds(_dp(-0.003), single, 0.0)
    assert not M.t_better_holds(_dp(-0.0029), single, 0.0)
    assert not M.t_better_holds(_dp(-0.004, c=0.0), single, 0.0)
    assert not M.t_better_holds(_dp(-0.004), single, 0.0001)
    three_up = dict(single, P1=0.001, P2=0.002, M1=0.0)
    assert not M.t_better_holds(_dp(-0.004), three_up, -0.01)                       # 5 of 8 lower
    assert M.t_better_holds(_dp(-0.004), dict(single, P1=0.001, P2=0.002), -0.01)    # 6 of 8


def test_equal_rule():
    d = dict(z8=0.003, halfA=-0.003, halfB=0.0)
    assert M.equal_rule(d, 0.005, 0.01, 0, 0)['equal']
    assert not M.equal_rule(dict(d, halfB=0.0031), 0.0, 0.0, 0, 0)['equal']
    assert not M.equal_rule(d, 0.0, 0.0, 0, 1)['equal']                              # a run at a bracket's bound
    assert M.equal_rule(dict(z8=-0.001, halfA=-0.002, halfB=-0.0001), 0.0, 0.0, 0, 0)['lower_on_all']
    with pytest.raises(ValueError):
        M.equal_rule(dict(z8=0.0, halfA=0.0), 0.0, 0.0, 0, 0)


# ---------------------------------------------------------------------------
# the KCs at low concentration

def test_low_conc_class():
    # the flies 0.114 with APL output blocked (TNT); C2 0.054 (a factor 2.11 below)
    assert M.low_conc_class(0.054, 0.114, 0.054) == 'about the same'
    assert M.low_conc_class(0.029, 0.114, 0.054) == 'further'                         # 3.93 below: 1.86 times further
    assert M.low_conc_class(0.045, 0.114, 0.054) == 'about the same'                  # 1.2 times further
    assert M.low_conc_class(0.0675, 0.114, 0.054) == 'nearer'                         # exactly 1.25 nearer
    assert M.low_conc_class(0.2, 0.114, 0.054) == 'about the same'                    # overshoots by as much as C2 falls short
    assert M.low_conc_class(0.0, 0.114, 0.054) == 'further'
    assert M.low_conc_class(0.0, 0.114, 0.0) == 'about the same'
    assert M.low_conc_class(0.054, 0.114, 0.0) == 'nearer'
    with pytest.raises(ValueError):
        M.low_conc_class(0.05, 0.0, 0.05)


# ---------------------------------------------------------------------------
# refusing to run past a failed identity

def _identity(tmp_path, points=None, readers=True, e77=None):
    ok = dict(passes=readers)
    d = dict(base_population=ok, e98_reading=ok, e100_reading=ok, points=points or {}, e77=e77 or {}, b4={})
    p = tmp_path / 'identity.json'
    p.write_text(json.dumps(d))
    return p


def test_require_identity(tmp_path):
    with pytest.raises(RuntimeError, match='part 0 has not run'):
        M.require_identity('point', ['real'], path=tmp_path / 'nothing.json')
    p = _identity(tmp_path, points={'real': dict(exact=True), 'P1': dict(exact=False)})
    M.require_identity('point', ['real'], path=p)
    with pytest.raises(RuntimeError, match='identity point P1: failed'):
        M.require_identity('point', ['real', 'P1'], path=p)
    with pytest.raises(RuntimeError, match='identity point P2: not run'):
        M.require_identity('point', ['P2'], path=p)
    M.require_identity('readers', path=p)
    with pytest.raises(RuntimeError, match='e98_reading: failed'):
        M.require_identity('readers', path=_identity(tmp_path, readers=False))
    p = _identity(tmp_path, points={'real': dict(exact=True)}, e77={'inhadd100': dict(exact=False)})
    with pytest.raises(RuntimeError, match='e77 inhadd100: failed'):
        M.require_identity('all', ['real'], path=p)


def test_run_point_refuses_past_a_failed_identity(tmp_path, monkeypatch):
    p = _identity(tmp_path, points={'P1': dict(exact=False)})
    monkeypatch.setattr(M, 'identity_path', lambda: p)

    def boom(*a, **k):
        raise AssertionError('simulated past a failed identity')
    monkeypatch.setattr(M, 'hold_point', boom)
    with pytest.raises(RuntimeError, match='refusing to go on'):
        M.run_point('P1', 100.0, root=tmp_path)
    with pytest.raises(RuntimeError, match='refusing to go on'):
        M.run_point('P2', 100.0, root=tmp_path)                                       # never run


# ---------------------------------------------------------------------------
# part B's runner on made-up runs

PANEL = list(H.APL16_ODORS)
EXTRA = ['ma [kiwi-control/diag] @ -7', '2h [kiwi-control/diag] @ -6', 'blank @ 0']
N_KC, N_CLAW = 6, 9


def _kc_index():
    return pd.MultiIndex.from_arrays([np.arange(N_KC), ['g'] * N_KC], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])


def _counts(cols, seed):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.integers(0, 4, size=(N_KC, len(cols))).astype(float), index=_kc_index(), columns=cols)


def _fake_held(seed=0):
    out = {}
    for i, c in enumerate(('on', 'off')):
        s = _counts(PANEL, seed + i)
        out[c] = ((s >= 2).astype(float), s)
    return dict(out, params=dict(fixed_thr=120.0, wAPLKC_scale=2.0, wKCAPL_scale=0.01, wAPLPN_scale=0.5, wPNAPL_scale=0.1,
                                 wKCKC_scale=40.0),
                params_off=dict(fixed_thr=120.0, wAPLKC_scale=0.0, wKCAPL_scale=0.01, wAPLPN_scale=0.0, wPNAPL_scale=0.1,
                                wKCKC_scale=40.0),
                thr_factor=1.1, apl_factor=0.9, n_on=10, n_off=20, target_on=10, target_off=20, at_bound_thr=False,
                at_bound_apl=False, n_sims=7, breadth_sha256='x', kws={}, history=[],
                pn_kws=dict(pn_inhadd=100.0, pn_inhsc=1.0), step='inhadd100', wiring='P1',
                breadth_series=pd.Series(np.ones(N_KC)))


def _fake_extended(held, differ=False, blank_off=0.0):
    labels = PANEL + EXTRA
    gl = list(H.apl16_orn_deltas().index)
    out = {}
    for i, (c, ct) in enumerate((('on', 'on'), ('off', 'acute_off'))):
        r16, s16 = held[c]
        s = pd.concat([s16, _counts(EXTRA, 10 + i)], axis=1)
        s.columns = pd.Index(labels, name='stimulus')
        if differ and c == 'on':
            s.iloc[0, 0] += 1.0
        r = (s >= 2).astype(float)
        win = pd.Series(np.linspace(9.0, 20.0, len(labels)), index=labels)
        win['blank @ 0'] = 8.0 + blank_off
        apl = pd.DataFrame(dict(win=win, pre=7.5, pre_settled=8.0))
        apl.index.name = 'stimulus'
        pn = pd.DataFrame(np.tile(np.linspace(10, 50, len(labels)), (len(gl), 1)), index=pd.Index(gl, name='glomerulus'),
                          columns=labels)
        claws = pd.DataFrame(np.tile(np.linspace(1, 5, len(labels)), (N_CLAW, 1)), columns=labels)
        trace = pd.DataFrame(np.ones((5, len(labels))), index=pd.Index(np.arange(5) * 0.001, name='time_s'), columns=labels)
        out[ct] = dict(responses=r, spike_counts=s, apl=apl, pn_win=pn, pn_pre=pn * 0.5, apl_trace=trace,
                       claw_win=claws, claw_pre=claws * 0.2, params=held['params' if c == 'on' else 'params_off'],
                       seconds=0.0)
    table = pd.DataFrame(dict(label=labels))
    meta = dict(claw_glomerulus=[gl[i % len(gl)] for i in range(N_CLAW)], kc_of_claw=[i % N_KC for i in range(N_CLAW)],
                bouton_glomerulus=gl[:3])
    return out, table, meta


@pytest.fixture
def runner(tmp_path, monkeypatch):
    p = _identity(tmp_path, points={'P1': dict(exact=True)})
    monkeypatch.setattr(M, 'identity_path', lambda: p)
    monkeypatch.setattr(M, 'run_commits', lambda: dict(queue_commit='abc'))
    calls = dict(hold=0)

    def hold(w, step):
        calls['hold'] += 1
        return _fake_held()
    monkeypatch.setattr(M, 'hold_point', hold)
    return tmp_path, calls


def test_run_point_on_made_up_runs(runner, monkeypatch):
    tmp_path, calls = runner
    monkeypatch.setattr(M, 'run_extended', lambda w, st, held, b: _fake_extended(held))
    rec = M.run_point('P1', 100.0, root=tmp_path)
    assert rec['complete'] and rec['extended_done'] and rec['extended_matches_held'] and rec['a0_identity_passes']
    assert rec['pn_inhadd'] == 100.0 and rec['params']['fixed_thr'] == 120.0 and rec['params_off']['wAPLKC_scale'] == 0.0
    d = tmp_path / 'B' / 'inhadd100'
    for f in ('P1_point.json', 'P1_family.csv', 'P1_runs.csv', 'b1_stimuli.csv', 'P1_units.json', 'P1_on_apl.csv',
              'P1_acute_off_spike_counts.parquet', 'P1_held_on_spike_counts.parquet', 'P1_on_pn_win.parquet'):
        assert (d / f).exists(), f
    held = M.load_held(100.0, 'P1', tmp_path)
    assert held['on'][1].equals(_fake_held()['on'][1])
    from experiments import E98_apl_response_strength as E98
    b1 = E98.load_b1(['P1'], roots=[d])                         # E98's loader reads the layout
    assert set(b1['cond']) == {'on', 'acute_off'} and len(b1) == 2 * (len(PANEL) + len(EXTRA))
    assert set(M.b1_rows(100.0, ['P1'], tmp_path)['stimulus']) == set(PANEL + EXTRA)
    fam = pd.read_csv(d / 'P1_family.csv')
    assert {'claw_input_00_win', 'apl_pre_settled', 'pn_total'} <= set(fam.columns)
    runs = pd.read_csv(d / 'P1_runs.csv')
    assert list(runs['cond']) == ['on', 'acute_off'] and runs['param_wAPLKC_scale'].tolist() == [2.0, 0.0]
    assert M.available_steps(tmp_path) == ['inhadd100']
    # resuming: a complete point is not run again
    rec2 = M.run_point('P1', 100.0, root=tmp_path)
    assert calls['hold'] == 1 and rec2['time'] == rec['time']


def test_run_point_stops_when_the_extended_run_differs(runner, monkeypatch):
    tmp_path, _ = runner
    monkeypatch.setattr(M, 'run_extended', lambda w, st, held, b: _fake_extended(held, differ=True))
    with pytest.raises(RuntimeError, match='a check failed'):
        M.run_point('P1', 100.0, root=tmp_path)
    rec = M.load_point(100.0, 'P1', tmp_path)
    assert not rec['complete'] and rec['extended_vs_held_on']['counts_differ'] == 1
    assert M.available_steps(tmp_path) == []


def test_run_point_stops_on_the_blank_identity(runner, monkeypatch):
    tmp_path, _ = runner
    monkeypatch.setattr(M, 'run_extended', lambda w, st, held, b: _fake_extended(held, blank_off=0.1))   # 1.25 % off
    with pytest.raises(RuntimeError, match='a check failed'):
        M.run_point('P1', 100.0, root=tmp_path)
    rec = M.load_point(100.0, 'P1', tmp_path)
    assert rec['extended_matches_held'] and not rec['a0_identity_passes']


# ---------------------------------------------------------------------------
# part A's chain on made-up runs

def test_chain_table_and_ratios(monkeypatch):
    gl = ['g1', 'g2', 'g3', 'g4']
    labels = ['a', 'b', 'blank @ 0']
    x = {'a': [40.0, 30.0, 20.0, 10.0], 'b': [20.0, 15.0, 10.0, 5.0], 'blank @ 0': [0.0] * 4}
    stim = pd.DataFrame([dict(label=k, **{f'orn_{g}': v for g, v in zip(gl, x[k])}) for k in labels])
    pn_w = pd.DataFrame({'a': [30.0, 25, 20, 15, 5], 'b': [22.0, 18, 15, 12, 5], 'blank @ 0': [10.0, 10, 10, 10, 5]},
                        index=pd.Index(gl + ['g5'], name='glomerulus'))
    apl = pd.DataFrame(dict(win=[18.0, 13.0, 8.0], pre=[7.0] * 3, pre_settled=[8.0] * 3), index=labels)
    resp = pd.DataFrame(np.zeros((3, 3)), columns=labels)
    runs = {'C2': {'on': dict(apl=apl, pn_win=pn_w, pn_pre=pn_w * 0 + 9.0, responses=resp)}}
    ch = M.chain_table(runs, stim)
    a = ch.set_index('stimulus').loc['a']
    assert a['input_total'] == 100.0 and a['input_median'] == 25.0 and a['input_broad'] == 100.0
    assert a['input_rest'] == pytest.approx(15.0 + 5.0)
    assert a['pn_total_blank'] == pytest.approx(50.0) and a['pn_total_blank_input'] == pytest.approx(50.0)
    assert a['apl_response'] == 10.0
    monkeypatch.setattr(M, 'flies_interval', lambda src: dict(flies_ratio=0.5, flies_lo90=0.4, flies_hi90=0.6))
    r = M.chain_ratios(ch, [dict(ratio_id='t', estimate='e', numerator='b', denominator='a', flies='x')]).iloc[0]
    assert r['input_total_ratio'] == pytest.approx(0.5)
    assert r["pn_total_blank_ratio"] == pytest.approx(27.0 / 50.0)          # (12 + 8 + 5 + 2 + 0) / 50
    assert r['apl_response_ratio'] == pytest.approx(0.5)
    assert r['exponent_model'] == pytest.approx(1.0) and r['exponent_flies'] == pytest.approx(1.0)
    assert r['al_share'] + r['mb_share'] == pytest.approx(1.0) or np.isnan(r['al_share'])
    assert r['model_class'] == 'inside'


def test_tables_carry_the_caveat(tmp_path):
    df = pd.DataFrame(dict(a=[1.0, 2.0], b=['x', 'y']))
    M.write_table(df, tmp_path / 't.csv', caveat=True)
    lines = (tmp_path / 't.csv').read_text().splitlines()
    assert lines[0] == f'# {M.YANG_QUALITY_CAVEAT}'
    assert M.read_table(tmp_path / 't.csv').equals(df)
    M.write_table(df, tmp_path / 'u.csv')
    assert M.read_table(tmp_path / 'u.csv').equals(df)


def test_update_identity_merges_points_under_its_lock(tmp_path, monkeypatch):
    p = tmp_path / 'identity.json'
    monkeypatch.setattr(M, 'identity_path', lambda: p)
    M.update_identity(points={'real': dict(exact=True)}, base_population=dict(passes=True))
    M.update_identity(points={'P1': dict(exact=False)})
    d = json.loads(p.read_text())
    assert d['points'] == {'real': dict(exact=True), 'P1': dict(exact=False)} and d['base_population'] == dict(passes=True)
    assert not (tmp_path / 'identity.json.lock').exists()
    (tmp_path / 'identity.json.lock').mkdir()                       # held by another lane
    with pytest.raises(TimeoutError):
        with M._Lock(p, timeout=0.5):
            pass


# ---------------------------------------------------------------------------
# odors silent in a run: their pairs handled explicitly, every other missing pair raising

def _silent_counts(rng, silent=()):
    x = pd.DataFrame(rng.poisson(1.0, size=(60, 16)).astype(float), columns=list(H.APL16_ODORS))
    x[list(silent)] = 0.0
    return x


def _sorted(s):
    """E83's `_cs` orientation (each pair's odors sorted), which differs from the panel's for many pairs."""
    s = s.copy()
    s.index = pd.MultiIndex.from_tuples([tuple(sorted(p)) for p in s.index], names=['odor_a', 'odor_b'])
    return s.sort_index()


def _data(rng):
    return H.canonical_pairs(pd.Series(rng.normal(0.18, 0.05, 120), index=pd.MultiIndex.from_tuples(H.PAIRS16)))


def test_silent_odors_and_label():
    rng = np.random.default_rng(0)
    o = {'on': (None, _silent_counts(rng, ['B-myr'])), 'off': (None, _silent_counts(rng))}
    assert M.silent_odors(o) == {'on': ['B-myr'], 'off': []}
    assert len(M.silent_pairs(['B-myr'])) == 15 and len(M.silent_pairs(['B-myr', 'ms'])) == 29
    assert M.silent_label({'M1': {'on': ['B-myr'], 'off': []}}) == 'B-myr (M1, APL on)'
    lab = M.silent_label({'M1': {'on': ['B-myr'], 'off': []}, 'P1': {'on': ['ms', 'B-myr'], 'off': ['ms']}})
    assert lab == 'ms (P1, APL on and off); B-myr (M1, APL on; P1, APL on)'
    assert M.silent_label({'M1': {'on': [], 'off': []}}) == ''


def test_unit_silence_and_mask_pairs():
    silent = {('c', 'M1'): {'on': ['B-myr'], 'off': []}, ('c', 'P1'): {'on': [], 'off': ['ms']},
              ('c', 'real'): {'on': [], 'off': []}}
    u = M.unit_silence(silent, 'c', ('real', 'P1', 'M1'))          # a population: every member's pairs
    assert u['n_pairs'] == 91 and u['odors'] == ['ms', 'B-myr'] and u['undefined'] == M.silent_pairs(['ms', 'B-myr'])
    assert u['silent_odors'] == 'ms (P1, APL off); B-myr (M1, APL on)'
    assert M.unit_silence(silent, 'c', ('real',)) == dict(undefined=frozenset(), n_pairs=120, silent_odors='', odors=[])
    rng = np.random.default_rng(1)
    d = {c: _sorted(H.readout_pair_distances(_silent_counts(rng), 'raw')) for c in ('on', 'off')}
    m = M.mask_pairs(d, u['undefined'])
    assert int(m['on'].isna().sum()) == 29 and int(m['off'].isna().sum()) == 29 and int(d['on'].isna().sum()) == 0
    assert 'all 120' == M.pairs_note(120, '')
    assert M.pairs_note(105, 'B-myr (M1, APL on)', 'population').startswith('105 of 120, the pairs every member')


def test_check_pairs_handles_a_silent_odor():
    rng = np.random.default_rng(2)
    on = H.readout_pair_distances(_silent_counts(rng, ['B-myr']), 'raw')            # B-myr's 15 pairs NaN
    off = H.readout_pair_distances(_silent_counts(rng), 'raw')
    assert int(on.isna().sum()) == 15
    data = _data(rng)
    und = M.silent_pairs(['B-myr'])
    assert M.check_pairs(on, off, data, und, 'x') == 105
    # the pairs absent instead of NaN (as in E83's per-pair rows), in E83's sorted orientation: the same
    on2, off2 = _sorted(on).dropna(), _sorted(off).reindex(_sorted(on).dropna().index)
    assert M.check_pairs(on2, off2, data, und, 'x') == 105
    # the statistics on the defined pairs are those of the 105 pairs alone
    st = H.apl_effect_stats(on, off, data)
    keep = [p for p in data.index if frozenset(p) not in und]
    st105 = H.apl_effect_stats(on.reindex(keep), off.reindex(keep), data.reindex(keep))
    assert st['apl_effect_n_pairs'] == 105 and st['apl_effect_rmse_centered'] == pytest.approx(st105['apl_effect_rmse_centered'])


def test_check_pairs_raises_on_an_ordering_fault():
    rng = np.random.default_rng(3)
    on, off = (H.readout_pair_distances(_silent_counts(rng), 'raw') for _ in range(2))
    data = _data(rng)
    assert M.check_pairs(on, off, data, frozenset(), 'x') == 120
    # the data in sorted orientation: pyolfsysm puts only the model's side in the panel's order, so pairs drop
    with pytest.raises(RuntimeError, match='ordering or alignment fault'):
        M.check_pairs(on, off, _sorted(data), frozenset(), 'x')
    # the same with a silent odor: the fault still raises
    on_s = H.readout_pair_distances(_silent_counts(rng, ['B-myr']), 'raw')
    with pytest.raises(RuntimeError, match='ordering or alignment fault'):
        M.check_pairs(on_s, off, _sorted(data), M.silent_pairs(['B-myr']), 'x')


def test_check_pairs_raises_on_pairs_no_silent_odor_explains():
    rng = np.random.default_rng(4)
    on, off = (H.readout_pair_distances(_silent_counts(rng), 'raw') for _ in range(2))
    data = _data(rng)
    with pytest.raises(RuntimeError, match='no silent odor behind them'):           # one pair dropped
        M.check_pairs(on.iloc[1:], off, data, frozenset(), 'x')
    on_nan = on.copy()
    on_nan.iloc[5] = np.nan
    with pytest.raises(RuntimeError, match='no silent odor behind them'):           # one pair NaN
        M.check_pairs(on_nan, off, data, frozenset(), 'x')
    with pytest.raises(RuntimeError, match="silent odors' pairs with a distance"):  # a silent odor's pairs present
        M.check_pairs(on, off, data, M.silent_pairs(['B-myr']), 'x')
    on_s = H.readout_pair_distances(_silent_counts(rng, ['B-myr']), 'raw')
    with pytest.raises(RuntimeError, match='no silent odor behind them'):           # silent, but not declared
        M.check_pairs(on_s, off, data, frozenset(), 'x')
    with pytest.raises(RuntimeError, match='repeated pairs'):
        M.check_pairs(pd.concat([on, on.iloc[:1]]), off, data, frozenset(), 'x')
    with pytest.raises(RuntimeError, match="data's per-pair effect"):
        M.check_pairs(on, off, data.iloc[1:], frozenset(), 'x')


def _pair_rows(change, unit, on, off):
    """E83's `score_runs` per-pair rows of one reading (sorted orientation, undefined pairs absent)."""
    me = _sorted(H.apl_effect(on, off))
    return pd.DataFrame(dict(change=change, unit=unit, readout='cal', odor_a=me.index.get_level_values(0),
                             odor_b=me.index.get_level_values(1), on=_sorted(on).reindex(me.index).to_numpy(),
                             off=_sorted(off).reindex(me.index).to_numpy(), effect=me.to_numpy()))


def test_effect_parts_with_a_silent_odor():
    rng = np.random.default_rng(5)
    off = H.readout_pair_distances(_silent_counts(rng), 'raw')
    on_s = H.readout_pair_distances(_silent_counts(rng, ['B-myr']), 'raw')
    on = H.readout_pair_distances(_silent_counts(rng), 'raw')
    pairs = pd.concat([_pair_rows('gain_fixed', 'M1', on_s, off), _pair_rows('gain_fixed', 'real', on, off)])
    tg = {'z8': _sorted(_data(rng)), 'halfA': _sorted(_data(rng))}         # E83's targets come sorted
    parts = M._effect_parts(pairs, tg, {('gain_fixed', 'M1'): M.silent_pairs(['B-myr'])}).set_index('unit')
    assert parts.loc['M1', 'n_pairs'] == 105 and parts.loc['real', 'n_pairs'] == 120
    sp = H.apl_effect_spread_shape(on_s, off, H.canonical_pairs(tg['z8']))
    assert parts.loc['M1', 'z8_shape'] == pytest.approx(sp['apl_effect_shape'])
    with pytest.raises(RuntimeError, match='gain_fixed M1 cal z8: 15 pairs missing'):    # not declared: raises
        M._effect_parts(pairs, tg)


def test_mark_silent():
    df = pd.DataFrame(dict(step=['inhadd100', 'gain_fixed'], x=[1, 2]))
    m = M.mark_silent(df, {'gain_fixed': 'B-myr (M1, APL on)'})
    assert m['silent_odors_16odor'].tolist() == ['', 'B-myr (M1, APL on)']
    with pytest.raises(ValueError):
        M.mark_silent(df.drop(columns=['step']), {})


# ---------------------------------------------------------------------------
# part D: the refits

PARAMS = dict(fixed_thr=118.95990186705896, wAPLKC_scale=0.1453349100781187, wKCAPL_scale=8.391161089960699e-05,
              wAPLPN_scale=0.3759560687748453, wPNAPL_scale=0.0010378099834196272, wKCKC_scale=26.43035318456299)
WIRING_ARGS = {'real': [], 'P1': ['--wiring', 'claw_permutation_within_type', '--wiring-seed', '1'],
               'M1': ['--wiring-file', 'WFILE', '--wiring-tag', 'curveball_within_type_matched-seed1']}


@pytest.fixture
def dfolders(tmp_path, monkeypatch):
    """local_dir and out_dir in tmp_path, the wiring options made up (M1's file a placeholder path)."""
    loc, res = tmp_path / 'outputs', tmp_path / 'results'
    res.mkdir()

    def local_dir(sub=''):
        d = loc / sub if sub else loc
        d.mkdir(parents=True, exist_ok=True)
        return d
    monkeypatch.setattr(M, 'local_dir', local_dir)
    monkeypatch.setattr(M, 'out_dir', lambda: res)
    monkeypatch.setattr(M, 'refit_wiring_args', lambda w: list(WIRING_ARGS[w]))
    return tmp_path


def _point(root, level, w, params=PARAMS, complete=True, held=None):
    d = M.step_dir(level, root)
    M._write_json(d / f'{w}_point.json', dict(wiring=w, params=dict(params), complete=complete))
    if held is not None:
        for c, (fr, fs) in M._held_paths(d, w).items():
            held[c][0].to_parquet(fr)
            held[c][1].to_parquet(fs)


def test_write_init_params_reads_back_exactly(dfolders):
    from experiments import E10_joint_fit as J
    root = dfolders / 'B_root'
    _point(root, 1000.0, 'real')
    _point(root, 1000.0, 'P1', dict(PARAMS, fixed_thr=108.57708630314373))
    _point(root, 1000.0, 'P2', complete=False)
    path = M.write_init_params(1000.0, root)
    assert path == M.init_params_path(1000.0) and path.name == 'init_inhadd1000.json'
    d = json.loads(path.read_text())
    assert set(d) == {'real', 'P1'}
    assert J.read_params(f'{path}:real') == PARAMS                                    # exactly, every parameter
    assert J.read_params(f'{path}:P1')['fixed_thr'] == 108.57708630314373
    with pytest.raises(RuntimeError, match='no complete part B point'):
        M.write_init_params(314.0, root)


def test_refit_args_and_tag_are_e10s(dfolders, monkeypatch):
    """The tag `refit_tag` predicts is the one E10's main gives the refit (its history file's name), for each kind of
    wiring; the options carry the level's PN kwargs exactly and the init file's entry for the wiring."""
    from experiments import E10_joint_fit as J
    wfile = dfolders / 'W.parquet'
    n = 4
    idx = pd.MultiIndex.from_arrays([[f'kc{i}' for i in range(n)]], names=[H.KC_ID_LEVEL])
    W = pd.DataFrame(np.ones((n, 2)), index=idx, columns=pd.MultiIndex.from_arrays([['g1', 'g2']], names=['glomerulus']))
    W.to_parquet(wfile)
    monkeypatch.setattr(M, 'refit_wiring_args', lambda w: [str(wfile) if a == 'WFILE' else a for a in WIRING_ARGS[w]])
    seen = {}

    class Stop(Exception):
        pass

    class Hist:
        def to_csv(self, path):
            seen['path'] = Path(path)
            raise Stop

    def fake_base(orn, out_dir, claw_thr=None, **kw):
        seen['base'] = kw
        initial = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.1, wPNAPL_scale=1e-3)
        wnorm = {k: pd.Series(np.ones(n), index=idx) for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')}
        return initial, wnorm, W, pd.Series(np.arange(n, dtype=float), index=idx.get_level_values(0))

    def fake_fit(k, orn, initial, wnorm, breadth, loss_fn, **kw):
        seen['fit'] = kw
        return dict(history=Hist()), {}, {}

    monkeypatch.setattr(J, 'bouton_base', fake_base)
    monkeypatch.setattr(J, 'fit_one', fake_fit)
    orig_r, orig_o = H.results_dir, H.outputs_dir
    monkeypatch.setattr(H, 'results_dir', lambda e: (dfolders / 'e10r') if e == 'E102_check' else orig_r(e))
    monkeypatch.setattr(H, 'outputs_dir', lambda e: (dfolders / 'e10o') if e == 'E102_check' else orig_o(e))
    for d in ('e10r', 'e10o'):
        (dfolders / d).mkdir(exist_ok=True)
    kws = M.step_kws(1000.0)
    for w in ('real', 'P1', 'M1'):
        args = M.refit_args(1000.0, w, experiment='E102_check')
        assert args[args.index('--pn-inhadd') + 1] == repr(kws['pn_inhadd'])
        assert float(args[args.index('--pn-inhsc') + 1]) == kws['pn_inhsc']
        assert args[args.index('--init-params') + 1] == f'{M.init_params_path(1000.0)}:{w}'
        assert args[args.index('--maxfev') + 1] == str(M.REFIT_MAXFEV) and args[args.index('--tag') + 1] == '_cap3'
        seen.clear()
        monkeypatch.setattr(sys, 'argv', ['E10_joint_fit.py', *args])
        with pytest.raises(Stop):
            J.main()
        assert seen['path'].name == f'history_{M.refit_tag(1000.0, w)}.csv'
        assert seen['base']['pn'] == kws and seen['fit']['pn'] == kws
        assert seen['fit']['init_params'] == f'{M.init_params_path(1000.0)}:{w}' and seen['fit']['aplpn_max'] == 3.0
    assert M.refit_tag(1000.0, 'P1') == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_pnadd1000sc544.444_k2_cap3'
    assert M.refit_tag(1000.0, 'real') == 'rate0.068_aplboth_pnadd1000sc544.444_k2_cap3'
    assert M.refit_tag(1000.0, 'M1') == 'wiring-curveball_within_type_matched-seed1_rate0.068_aplboth_pnadd1000sc544.444_k2_cap3'


def _held16(seed=0):
    out = {}
    for i, c in enumerate(('on', 'off')):
        s = _counts(PANEL, seed + i)
        out[c] = ((s >= 2).astype(float), s)
    return out


def test_compare_start():
    held = _held16()
    b = pd.Series(np.arange(N_KC, dtype=float), index=pd.Index(np.arange(N_KC), name=H.KC_ID_LEVEL))
    assert M.compare_start(held, held, b, b.copy())['exact']
    e10 = {c: (r.copy(), s.copy()) for c, (r, s) in held.items()}
    e10['off'][1].iloc[2, 3] += 1.0
    r = M.compare_start(e10, held, b, b)
    assert not r['exact'] and r['off']['counts_differ'] == 1 and r['on']['counts_differ'] == 0
    b2 = b.copy()
    b2.iloc[0] += 1.0
    r = M.compare_start(held, held, b2, b)
    assert not r['exact'] and r['breadth_max_abs_diff'] == 1.0
    r = M.compare_start(held, held, b.iloc[::-1], b)
    assert not r['exact'] and not r['breadth_index_equal']


def _fake_e10_run(held, tag, breadth_counts, differ=False):
    """A stand-in for E78's redirected E10 run: writes the outputs, the breadth run and the history an exact start
    writes."""
    import subprocess

    def run(cmd, cwd=None, stdout=None, stderr=None):
        root = Path(cmd[cmd.index('--run-e10') + 1])
        o, r = root / 'outputs', root / 'results'
        o.mkdir(parents=True, exist_ok=True)
        r.mkdir(parents=True, exist_ok=True)
        for c, ct in M.REFIT_CONDS:
            rr, ss = held[c][0].copy(), held[c][1].copy()
            if differ and c == 'on':
                ss.iloc[0, 0] += 1.0
            rr.to_parquet(o / f'{tag}_{ct}_responses.parquet')
            ss.to_parquet(o / f'{tag}_{ct}_spike_counts.parquet')
        breadth_counts.to_parquet(o / 'pnadd1000sc544.444_base_acute_off_spike_counts.parquet')
        pd.DataFrame([dict(eval=0, loss=1.5, **PARAMS)]).to_csv(r / f'history_{tag}.csv', index=False)
        assert '--maxfev' in cmd and cmd[cmd.index('--maxfev') + 1] == '1'
        return subprocess.CompletedProcess(cmd, 0)
    return run


@pytest.mark.parametrize('differ', [False, True])
def test_start_check(dfolders, monkeypatch, differ):
    """The start check passes on E10 outputs equal to part B's point and its breadth, and stops (recording why) on one
    count's difference."""
    import subprocess
    root = dfolders / 'B_root'
    held = _held16()
    _point(root, 1000.0, 'P1', held=held)
    off = _counts(PANEL, 7)
    monkeypatch.setattr(M, 'point_breadth', lambda w, kws: M._breadth_of(off))
    tag = M.refit_tag(1000.0, 'P1')
    monkeypatch.setattr(subprocess, 'run', _fake_e10_run(held, tag, off, differ))
    if differ:
        with pytest.raises(RuntimeError, match='does not reproduce'):
            M.start_check(1000.0, 'P1', root)
    else:
        rec = M.start_check(1000.0, 'P1', root)
        assert rec['exact'] and rec['params_max_rel_diff'] == 0.0 and rec['first_eval_loss'] == 1.5
    saved = json.loads(M.start_path(1000.0).read_text())['P1']
    assert saved['exact'] == (not differ)
    assert saved['on']['counts_differ'] == (1 if differ else 0)


def test_refit_params_from_history(dfolders):
    res = M.out_dir()
    tag = 'rate0.068_aplboth_pnadd1000sc544.444_k2_cap3'
    best = dict(PARAMS, fixed_thr=121.123456789012345)
    pd.DataFrame([dict(eval=0, loss=2.0, **PARAMS), dict(eval=1, loss=1.0, **best)]).to_csv(res / f'history_{tag}.csv', index=False)
    row = dict(tag=tag, **{k: v for k, v in best.items() if k != 'wKCKC_scale'}, kc_kc_inh=best['wKCKC_scale'])
    pd.DataFrame([row]).to_csv(res / 'fits.csv', index=False)
    assert M.refit_params(tag) == best
    pd.DataFrame([dict(row, fixed_thr=120.0)]).to_csv(res / 'fits.csv', index=False)
    with pytest.raises(RuntimeError, match='disagree'):
        M.refit_params(tag)


def test_refit_due(dfolders):
    res = M.out_dir()
    M._write_json(res / 'levels.json', dict(refit_levels=[1000.0]))
    assert M.refit_due(1000.0, M.S_WIRINGS)[0]
    assert not M.refit_due(3141.0, M.S_WIRINGS)[0]
    assert not M.refit_due(1000.0, M.T_WIRINGS)[0]                                     # S not scored
    for cls, ok in (('price', False), ('equal', True), ('better', True)):
        M._write_json(res / 'd_class_inhadd1000.json', dict(S={'class': cls}))
        assert M.refit_due(1000.0, ['P3', 'M2'])[0] == ok


def test_extended_counts():
    w98 = pd.DataFrame(dict(comparison=['V2, x'] * 4 + ['V2, x'], reading_set=['carrying'] * 4 + ['settled'],
                            wiring=['real', 'P1', 'P2', 'M1', 'real'], cls=['inside', 'inside', 'flatter', 'inside', 'steeper'],
                            cls_linear=['flatter'] * 5, ratio_linear=[0.7, 0.8, 0.9, 0.75, 0.1], flies_lo90=0.46,
                            flies_hi90=0.61))
    w100 = pd.DataFrame(dict(a0_level=['blank', 'blank', 'settled'], stimulus='farn', estimate='e1', wiring=['real', 'P1', 'real'],
                             cls=['inside', 'steeper', 'inside'], cls_linear=['inside'] * 3))
    w101 = pd.DataFrame(dict(a0_level=['blank'] * 2, ratio_id='r1', quantity='APL', estimate='e', wiring=['real', 'P1'],
                             cls=['flatter', 'flatter']))
    c = M.extended_counts({'refit': dict(e98_wirings=w98, e100_wirings=w100, e101_wirings=w101)})
    r98 = c[c['source'] == 'E98'].iloc[0]
    assert (r98['n_wirings'], r98['n_inside'], r98['n_flatter'], r98['n_steeper'], r98['n_flatter_linear']) == (4, 3, 1, 0, 4)
    assert r98['ratio_linear_min'] == 0.7 and r98['ratio_linear_max'] == 0.9 and r98['set'] == 'refit'
    r100 = c[c['source'] == 'E100'].iloc[0]
    assert (r100['n_wirings'], r100['n_inside'], r100['n_steeper']) == (2, 1, 1)
    r101 = c[c['source'] == 'E101'].iloc[0]
    assert r101['n_flatter'] == 2 and pd.isna(r101['n_inside_linear'])                 # E101's rows have no linear class


def test_compare_table():
    def rows(change, unit, base):
        out = []
        for ro in ('cal', 'raw'):
            d = dict(change=change, unit=unit, readout=ro, n_pairs=120, silent_odors='')
            d.update({c: base + i * 0.001 for i, c in enumerate(M.COMPARE_CAL)})
            out.append(d)
        return out
    rd = pd.DataFrame(rows('none', 'S', 0.05) + rows('refit_inhadd1000', 'S', 0.06) + rows('none', 'real', 0.04))
    fixed = pd.DataFrame(rows('inhadd1000', 'S', 0.055))
    t = M.compare_table(rd, fixed, 1000.0, ['S'])
    assert list(t['row']) == ['unchanged fits', 'fixed scales at inhadd1000 (part B)', 'refit at inhadd1000',
                              "the flies (floors, ceilings, the data's mean APL effect)"]
    assert t['cal_z8_pattern'].tolist()[:3] == [0.05, 0.055, 0.06]
    assert t['cal_mean_effect'].iloc[3] == pytest.approx(H.data_target('apl_effect_mean'))
    assert t['cal_z8_pattern'].iloc[3] == pytest.approx(H.data_target('apl_effect_pattern_floor'))
    t = M.compare_table(rd, None, 1000.0, ['S'])
    assert len(t) == 3
