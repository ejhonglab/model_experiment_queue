#!/usr/bin/env python3
"""Tests for E102 (the antennal lobe's normalization at the strength the flies' APL asks for;
experiments/E102_antennal_lobe_normalization.py), on made-up rows and runs: the step names and PN kwargs; the three level
rules, with a flagged step and with none existing; the refit classes (equal, better, price) and T's rule for better; the
low-concentration classes; the equal rule's inputs; the refusal to run past a failed (or missing) identity, in
`require_identity` and in part B's runner; part B's runner on made-up held and extended runs (its files in E98's layout,
which E98's loader reads, the check that the extended run's panel columns equal the held run's, E98's amended blank
identity, resuming); the chain of totals and its ratios; the results tables' caveat line."""

import json

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
