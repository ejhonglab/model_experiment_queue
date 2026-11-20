#!/usr/bin/env python3
"""Tests for E113 (the input of paa and B-myr without a search; experiments/E113_input_within_spread.py): the control
glomeruli rule, the arms' factors, the lanes, the gate of an arm and part 0's rule on made-up readings, the outcome
letters and the refit's choice; with the shared directory (no simulation), the sizes against E109's and E116's tables,
the control list on the lab's input, the per-fly deviations against E108's and the inputs."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E113_input_within_spread')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        H.data_targets_path().read_bytes()
        return True
    except (FileNotFoundError, OSError):
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='the shared directory is not available')


def test_interface():
    assert X.ARMS == ('D2', 'D1', 'D1x') and X.CLASS_ARM == 'D2'
    assert X.CONTROLS == ('C1', 'C2', 'C3', 'C4', 'D1inv')
    assert set(X.EXCLUDED) == {'DM1', 'DP1m', 'DL5', 'VM7d', 'DP1l', 'VC3', 'VM5d'}


def _lab():
    gl = ['DM1', 'DP1m', 'DL5', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6']
    cols = pd.MultiIndex.from_tuples([('validation2', 'paa @ -4'), ('validation2', 'B-myr @ -2'), ('megamat', '2h @ -6')],
                                     names=['panel', 'odor'])
    v = np.array([[100, 50, 10], [95, 49, 1], [99, 48, 1], [90, 10, 1], [80, 45, 1], [60, 30, 1], [101, 20, 1],
                  [70, 47, 1], [10, 5, 1]], dtype=float)
    return pd.DataFrame(v, index=gl, columns=cols)


def test_control_glomeruli_rule():
    cg = X.control_glomeruli(_lab())
    # paa: |x - 100| over G1..G6 = 10, 20, 40, 1, 30, 90 -> G4, G1, G2, G5; B-myr: |x - 50| = 40, 5, 20, 30, 3, 45
    assert cg['paa'] == ('G4', 'G1', 'G2', 'G5')
    assert cg['B-myr'] == ('G5', 'G2', 'G3', 'G4')


def test_entry_sizes_and_factors():
    s = X.entry_sizes(check=False)
    assert set(s) == set(X.ARMS + X.CONTROLS)
    assert s['D1x'] == {('paa', 'DM1'): -1.0, ('B-myr', 'DM1'): -1.0}
    assert s['D1inv'] == {k: -v for k, v in s['D1'].items()}
    assert set(s['D2']) == {('paa', 'DM1'), ('B-myr', 'DM1'), ('paa', 'DP1m')}
    for i in range(4):
        c = s[f'C{i + 1}']
        assert c == {('paa', X.CONTROL_GLOMERULI['paa'][i]): s['D1'][('paa', 'DM1')],
                     ('B-myr', X.CONTROL_GLOMERULI['B-myr'][i]): s['D1'][('B-myr', 'DM1')]}
    lab = pd.DataFrame(1.0, index=['DM1', 'DP1m', 'VL2a', 'DM6', 'DL2v', 'VA3', 'VM5v', 'DC1', 'DM4'],
                       columns=pd.MultiIndex.from_tuples([('v', 'paa @ -4'), ('v', 'B-myr @ -2')]))
    F = X.factors_of('D2', lab)
    assert list(F.columns) == ['paa', 'B-myr']
    assert F.loc['DP1m', 'paa'] == 2 ** X.D2_LOG2[('paa', 'DP1m')] and F.loc['DP1m', 'B-myr'] == 1.0
    assert (F.drop(['DM1', 'DP1m']) == 1.0).all().all()
    F = X.factors_of('C2', lab)
    assert F.loc['DM6', 'paa'] < 1 and F.loc['DL2v', 'B-myr'] < 1 and F.loc['DM1'].eq(1.0).all()


def test_lane_points_and_fly_names():
    pts = [dict(tag=str(i)) for i in range(10)]
    assert [p['tag'] for p in X.lane_points(pts, '1/3')] == ['1', '4', '7']
    assert X.lane_points(pts, None) == pts
    with pytest.raises(ValueError):
        X.lane_points(pts, '3/3')
    assert X.parse_fly('F3_pb') == (3, 'pb') and X.parse_fly('F0_bmyr') == (0, 'bmyr') and X.parse_fly('D1') is None
    names = X.input_names(bmyr=False)
    assert names[:3] == list(X.ARMS) and len(names) == 3 + 5 + 12
    assert len(X.input_names(bmyr=True)) == 26


# ---------------------------------------------------------------------------
# the gates on made-up readings

def _rows(change, units, d_on, d_off, per_target=None, extra=None):
    rows = []
    for unit in units:
        for ro in ('raw', 'cal'):
            r = dict(change=change, unit=unit, readout=ro, mean_effect=0.18, block=0.0, vcf_on=0.3, vcf_off=0.2,
                     struct_wt_on=0.3, struct_wt_off=0.2, paa_bmyr_on=0.5, paa_bmyr_off=0.4)
            for t in ('all', 'A', 'B'):
                a, b = (per_target or {}).get(t, (d_on, d_off))
                r[f'rho_on_{t}'] = 0.68 + a
                r[f'rho_off_{t}'] = 0.76 + b
            for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson_effect', 'rho_effect', 'loss'):
                r[f'all_{k}'] = 0.05
            r.update(extra or {})
            rows.append(r)
    return rows


def _runs(change, wirings, d_rate=0.0, flagged=False):
    return [dict(change=change, wiring=w, flagged=flagged, rate_rho_on=0.72 + d_rate, rate_rho_off=0.80) for w in wirings]


T = list(X.T_WIRINGS)
ALL = list(X.ALL_WIRINGS)


def _arm_tables(gain_on=0.03, gain_off=0.02, ctl_on=0.0, d_rate=0.0, flagged=False):
    rd = _rows('I0', ALL + ['T', 'all8', 'S'], 0, 0)
    runs = _runs('I0', ALL)
    for a in X.ARMS:
        rd += _rows(a, ALL + ['T', 'all8', 'S'], gain_on, gain_off)
        runs += _runs(a, ALL, d_rate, flagged)
    for c in X.CONTROLS:
        rd += _rows(c, ALL + ['T', 'all8', 'S'], ctl_on, 0.0)
        runs += _runs(c, ALL)
    return pd.DataFrame(rd), pd.DataFrame(runs)


def test_gate_arm_classes():
    g = X.gate_arm(*_arm_tables(), 'D2')
    assert g['complete'] and g['cls'] == 'passes' and g['gain_on'] and g['gain_off']
    assert g['all8']['n_wirings_gain_on_all'] == 8
    g = X.gate_arm(*_arm_tables(d_rate=-0.06), 'D2')
    assert g['cls'] == 'a price' and 'response rate ρ, APL on' in g['fails']
    g = X.gate_arm(*_arm_tables(flagged=True), 'D1')
    assert g['cls'] == 'a price'
    # one control at more than half the gain: the controls fail on that ρ
    rd, runs = _arm_tables(gain_off=0.0)
    rd.loc[(rd.change == 'C3') & rd.unit.isin(T + ['T']), 'rho_on_all'] = 0.68 + 0.02
    g = X.gate_arm(rd, runs, 'D2')
    assert not g['controls']['on']['ok'] and g['cls'] == 'does not carry'
    g = X.gate_arm(*_arm_tables(gain_on=0.008, gain_off=0.0), 'D2')
    assert g['cls'] == 'carries below the gate'
    rd, runs = _arm_tables()
    g = X.gate_arm(rd[rd.change != 'C4'], runs, 'D2')
    assert not g['complete'] and 'C4' in g['note']


def test_gate_arm_needs_each_target():
    rd, runs = _arm_tables(gain_off=0.0)
    rd.loc[(rd.change == 'D2') & (rd.unit == 'T'), 'rho_on_B'] = 0.68 + 0.005
    g = X.gate_arm(rd, runs, 'D2')
    assert not g['gain_on'] and g['cls'] != 'passes'


def _fly_tables(best_A=2, best_B=4, gain_held=0.03, kind='pb'):
    rd = _rows('I0', ALL + ['T', 'all8', 'S'], 0, 0)
    runs = _runs('I0', ALL)
    for f in range(X.N_FLIES):
        n = X.fly_name(f, kind)
        s = {'A': (0.02 if f == best_A else 0.0, 0.0), 'B': (0.02 if f == best_B else 0.0, 0.0), 'all': (0.0, 0.0)}
        rd += _rows(n, ['S'], 0, 0, per_target=s)
        # on T: fly best_A gains against B, fly best_B against A
        t = {'all': (0.0, 0.0), 'A': (gain_held if f == best_B else 0.0, 0.0), 'B': (gain_held if f == best_A else 0.0, 0.0)}
        rd += _rows(n, ALL + ['T', 'all8'], 0, 0, per_target=t)
        runs += _runs(n, ALL)
    return pd.DataFrame(rd), pd.DataFrame(runs)


def test_gate_flies_rule():
    rd, runs = _fly_tables()
    g = X.gate_flies(rd, runs, 'pb')
    assert g['complete'] and g['picks'] == {'A': 'F2_pb', 'B': 'F4_pb'}
    assert g['gain_on'] and g['cls'] == 'passes'
    assert g['directions']['AtoB']['pick'] == 'F2_pb' and g['directions']['AtoB']['read_on'] == 'B'
    rd, runs = _fly_tables(gain_held=0.005)
    assert X.gate_flies(rd, runs, 'pb')['cls'] == 'carries below the gate'
    # a fly that gains only on the half it was picked on does not pass
    rd, runs = _fly_tables()
    rd.loc[(rd.change == 'F2_pb') & (rd.unit.isin(T + ['T'])), 'rho_on_B'] = 0.68
    g = X.gate_flies(rd, runs, 'pb')
    assert not g['gain_on']
    assert not X.gate_flies(rd, runs, 'paa')['complete']


def test_outcomes_and_refit_choice():
    rd, runs = _arm_tables()
    arms = {a: X.gate_arm(rd, runs, a) for a in X.ARMS}
    fr, fruns = _fly_tables()
    flies = {'pb': X.gate_flies(fr, fruns, 'pb'), 'paa': {'complete': True, 'cls': 'does not carry'}}
    assert X.outcomes(arms, flies) == ['a', 'c']
    assert X.refit_choice(arms, flies, rd) == 'D2'
    rd2, runs2 = _arm_tables(gain_on=0.008, gain_off=0.0)
    arms2 = {a: X.gate_arm(rd2, runs2, a) for a in X.ARMS}
    assert X.outcomes(arms2, {'pb': {}, 'paa': {}}) == ['d']
    assert X.refit_choice(arms2, {'pb': {}, 'paa': {}}, rd2) is None
    both = pd.concat([rd2, fr[fr.change != 'I0']], ignore_index=True)
    assert X.refit_choice(arms2, flies, both) in ('F2_pb', 'F4_pb')
    rd3, runs3 = _arm_tables(gain_on=0.0, gain_off=0.0, ctl_on=0.0)
    arms3 = {a: X.gate_arm(rd3, runs3, a) for a in X.ARMS}
    assert X.outcomes(arms3, {'pb': {}, 'paa': {}}) == ['e']
    arms2['D1'] = X.gate_arm(*_arm_tables(), 'D1')
    assert X.outcomes(arms2, {})[0] == 'b'


# ---------------------------------------------------------------------------
# with the shared directory (no simulation)

@needs_shared
def test_sizes_and_controls_from_their_sources():
    s = X.entry_sizes(check=True)
    assert abs(s['D1'][('paa', 'DM1')] + 0.60) < 0.005 and abs(s['D2'][('paa', 'DP1m')] + 0.49) < 0.005
    assert X.control_glomeruli() == X.CONTROL_GLOMERULI


@needs_shared
def test_fly_inputs():
    D = X.fly_deviations()
    assert len(D['keys']) == 6
    for o in ('paa', 'B-myr'):
        dev, V = D['dev'][o], D['values'][o]
        cols = np.isfinite(V).any(axis=0)
        assert np.allclose(np.where(np.isfinite(V), dev, 0).sum(axis=0)[cols], 0, atol=1e-9)
    lab = X.lab_input()
    for f in range(6):
        pb, pa = X.fly_input(f, 'pb'), X.fly_input(f, 'paa')
        ab = [H.odor_abbrev(c) for c in lab.columns]
        changed = {ab[j] for j in range(len(ab)) if not np.array_equal(pb.iloc[:, j].to_numpy(), lab.iloc[:, j].to_numpy())}
        assert changed <= {'paa', 'B-myr'}
        assert np.array_equal(pa.iloc[:, ab.index('B-myr')].to_numpy(), lab.iloc[:, ab.index('B-myr')].to_numpy())
        assert np.array_equal(pa.iloc[:, ab.index('paa')].to_numpy(), pb.iloc[:, ab.index('paa')].to_numpy())
    D1 = X.input_of('D1')
    diff = (D1.to_numpy() != lab.to_numpy())
    assert diff.sum() == 2
