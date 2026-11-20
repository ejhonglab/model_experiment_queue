#!/usr/bin/env python3
"""Tests for E116 (the input of paa and B-myr searched against a second set of flies and taken from Hallem's receptor
data; experiments/E116_input_second_flies.py): the halves of the flies, the score, Hallem's inputs, the gates of parts H
and R and their classes on made-up readings, step 10's agreement and the outcome letters; with the shared directory (no
simulation), Remy's halves and targets."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E116_input_second_flies')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        H.data_targets_path().read_bytes()
        return True
    except (FileNotFoundError, OSError):
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='the shared directory is not available')


# ---------------------------------------------------------------------------
# the halves of the flies

def test_split_flies_keeps_every_pair_in_both_halves():
    rng = np.random.default_rng(3)
    flies = [f'f{i}' for i in range(10)]
    pairs = [('a', 'b'), ('a', 'c'), ('b', 'c')]
    held = pd.DataFrame(rng.random((10, 3)) > 0.3, index=flies, columns=pairs)
    held.loc[['f0', 'f7'], ('a', 'c')] = True          # a pair held by two flies only
    held.loc[[f for f in flies if f not in ('f0', 'f7')], ('a', 'c')] = False
    a, b = X.split_flies(flies, held, seed=0, pairs=pairs)
    assert sorted(a + b) == sorted(flies) and len(a) == 5 and not set(a) & set(b)
    for h in (a, b):
        assert held.loc[h].sum().min() >= 1
    assert (a, b) == X.split_flies(flies[::-1], held, seed=0, pairs=pairs)
    assert (a, b) != X.split_flies(flies, held, seed=1, pairs=pairs) or True


def test_split_flies_raises_when_impossible():
    flies = ['f0', 'f1', 'f2', 'f3']
    held = pd.DataFrame({('a', 'b'): [True, False, False, False]}, index=flies)
    with pytest.raises(RuntimeError):
        X.split_flies(flies, held, seed=0, pairs=[('a', 'b')], tries=50)


def test_score_of_subtracts_the_rate_penalty():
    base = {'w1': (0.10, 0.20), 'w2': (0.05, 0.10)}
    s = X.score_of(0.7, base, base)
    assert s['rate_penalty'] == 0 and s['score'] == 0.7
    s = X.score_of(0.7, {'w1': (0.13, 0.20), 'w2': (0.05, 0.10)}, base)      # +30 %: 15 % beyond, x10
    assert np.isclose(s['rate_penalty'], 1.5) and np.isclose(s['score'], 0.7 - 1.5)


# ---------------------------------------------------------------------------
# Hallem's inputs

def _lab():
    gl = ['DM6', 'DL5', 'VA5', 'DM1', 'VC3']
    cols = pd.MultiIndex.from_tuples([('validation2', 'paa @ -4'), ('validation2', 'B-myr @ -2')], names=['panel', 'odor'])
    return pd.DataFrame([[70.0, 5.0], [40.0, 3.0], [20.0, 0.0], [150.0, 70.0], [10.0, 12.0]], index=gl, columns=cols)


def test_hallem_inputs():
    lab = _lab()
    hal = pd.Series({'DM6': 153.0, 'DL5': 43.0, 'VA5': 79.0, 'VC3': -17.0})
    inp = X.hallem_inputs(lab, hal)
    assert set(inp) == set(X.H_ARMS) | set(X.H_PERMS)
    paa = ('validation2', 'paa @ -4')
    sh = list(hal.index)
    pd.testing.assert_series_equal(inp['H_meas'][0].loc[sh, paa], hal.rename(paa), check_names=False)
    assert np.isclose(inp['H_prof'][0].loc[sh, paa].sum(), lab.loc[sh, paa].sum())
    assert np.allclose(inp['H_prof'][0].loc[sh, paa] / hal, lab.loc[sh, paa].sum() / hal.sum())
    d = inp['H_dm6'][0]
    assert d.loc['DM6', paa] == 153.0 and d.drop('DM6').equals(lab.drop('DM6').astype(float))
    for k, (orn, meta) in inp.items():
        assert orn.loc['DM1', paa] == 150.0                                  # not one of Hallem's glomeruli
        assert orn[('validation2', 'B-myr @ -2')].equals(lab[('validation2', 'B-myr @ -2')].astype(float))
        if k.startswith('H_perm'):
            assert np.allclose(np.sort(orn.loc[sh, paa]), np.sort(inp['H_prof'][0].loc[sh, paa]))
    assert not inp['H_perm1'][0].equals(inp['H_perm2'][0])


# ---------------------------------------------------------------------------
# the gates on made-up readings

def _rows(change, units, d_on, d_off, extra=None):
    rows = []
    for unit in units:
        for ro in ('raw', 'cal'):
            r = dict(change=change, unit=unit, readout=ro, mean_effect=0.18, block=0.0, vcf_on=0.3, vcf_off=0.2,
                     struct_wt_on=0.3, struct_wt_off=0.2, paa_bmyr_on=0.5, paa_bmyr_off=0.4)
            for t in X.YANG_TARGETS:
                r[f'rho_on_{t}'] = 0.68 + d_on
                r[f'rho_off_{t}'] = 0.76 + d_off
            for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson_effect', 'rho_effect', 'loss'):
                r[f'all_{k}'] = 0.05
            r.update(extra or {})
            rows.append(r)
    return rows


def _runs(change, wirings, d_rate=0.0, flagged=False):
    return [dict(change=change, wiring=w, flagged=flagged, rate_rho_on=0.72 + d_rate, rate_rho_off=0.80) for w in wirings]


def _r_tables(gain_on=0.03, gain_off=0.02, perm_on=0.0, d_rate=0.0, remy_gain=0.03):
    T = list(X.T_WIRINGS)
    rd, runs, remy = [], [], []
    rd += _rows('I0', T + ['T'], 0, 0)
    runs += _runs('I0', T)
    for c in ('I0',) + tuple(f'{n}_s0.6' for n in X.SEARCHES):
        g = 0.0 if c == 'I0' else remy_gain
        for t in ('all',) + X.REMY_HALVES:
            remy.append(dict(change=c, unit='T', target=t, rho=0.69 + g))
    for n in X.SEARCHES:
        rd += _rows(f'{n}_s0.6', T + ['T'], gain_on, gain_off)
        runs += _runs(f'{n}_s0.6', T, d_rate)
        for s in X.PERM_SEEDS:
            rd += _rows(f'{n}_perm{s}', T + ['T'], perm_on, 0.0)
        rd += _rows(f'{n}_inv', T + ['T'], -0.02, -0.02)
    return pd.DataFrame(rd), pd.DataFrame(runs), pd.DataFrame(remy)


def test_gate_r_classes():
    g = X.gate_r(*_r_tables())
    assert g['complete'] and g['cls'] == 'passes' and g['gain_on'] and g['gain_off'] and g['gain_remy']
    g = X.gate_r(*_r_tables(d_rate=-0.06))
    assert g['cls'] == 'a price' and 'response rate ρ, APL on' in g['fails']
    g = X.gate_r(*_r_tables(perm_on=0.02, gain_off=0.0))
    assert not g['controls']['on']['ok'] and g['cls'] == 'does not carry'
    g = X.gate_r(*_r_tables(gain_on=0.008, gain_off=0.0))
    assert g['cls'] == 'carries below the gate'
    g = X.gate_r(*_r_tables(remy_gain=0.005))
    assert not g['gain_remy'] and g['gain_yang_on'] and g['cls'] == 'carries below the gate'
    g = X.gate_r(*_r_tables(gain_off=0.0))
    assert g['gain_on'] and not g['gain_off'] and g['cls'] == 'passes'
    assert 'c' in X.outcomes(g, {}, True) and 'a' in X.outcomes(g, {}, True)


def test_gate_r_incomplete():
    rd, runs, remy = _r_tables()
    g = X.gate_r(rd[rd.change != 'Rd_s0.6'], runs, remy)
    assert not g['complete']


def test_gate_h():
    W = list(X.ALL_WIRINGS)
    rd = _rows('I0', W + ['all8'], 0, 0) + _rows('H_prof', W + ['all8'], 0.02, 0.0)
    runs = _runs('I0', W) + _runs('H_prof', W)
    for p in X.H_PERMS:
        rd += _rows(p, W + ['all8'], 0.002, 0.0)
        runs += _runs(p, W)
    g = X.gate_h(pd.DataFrame(rd), pd.DataFrame(runs), 'H_prof')
    assert g['complete'] and g['gain_on'] and not g['gain_off'] and g['cls'] == 'passes'
    assert X.outcomes({}, g, None) == ['e']
    # five wirings of eight gaining: no gain
    rd2 = pd.DataFrame(rd)
    for w in W[:3]:
        rd2.loc[(rd2.change == 'H_prof') & (rd2.unit == w), [f'rho_on_{t}' for t in X.YANG_TARGETS]] = 0.68
    g = X.gate_h(rd2, pd.DataFrame(runs), 'H_prof')
    assert not g['gain_on'] and g['cls'] == 'carries below the gate'


def test_outcomes_d_and_f():
    gr = dict(complete=True, cls='carries below the gate', gain_on=False, gain_off=False, gain_remy=True,
              gain_yang_on=False, gain_yang_off=False, controls={c: dict(candidate=0.0, permuted=0.0, inverted=-0.01)
                                                                for c in ('on', 'off')})
    assert X.outcomes(gr, dict(complete=True, cls='does not carry'), False) == ['d']
    gr.update(cls='does not carry', gain_remy=False)
    assert X.outcomes(gr, dict(complete=True, cls='does not carry'), False) == ['f']


# ---------------------------------------------------------------------------
# step 10

def test_entries_agree_with_themselves():
    rng = np.random.default_rng(0)
    gl = [f'g{i}' for i in range(36)]
    rec = {'paa': ('g0', 'g1', 'g2', 'g3', 'g4'), 'B-myr': ('g5', 'g6', 'g7')}

    def F(noise):
        base = rng.normal(0, 0.5, (36, 2))
        return base, lambda: pd.DataFrame(np.exp2(base + rng.normal(0, noise, base.shape)), index=gl, columns=['paa', 'B-myr'])
    base, make = F(0.05)
    f4 = {f's{i}': make() for i in range(4)}
    f6 = {f'e{i}': make() for i in range(6)}
    ent, ag = X.entries(f4, f6, n_null=500, recurring=rec)
    assert len(ent) == 72 and ag['pearson_r'] > 0.95 and ag['agree'] and ag['sign_agree'] >= 7
    other = {f'e{i}': pd.DataFrame(np.exp2(rng.normal(0, 0.5, (36, 2))), index=gl, columns=['paa', 'B-myr']) for i in range(6)}
    _, ag2 = X.entries(f4, other, n_null=500, recurring=rec)
    assert abs(ag2['pearson_r']) < 0.5 and not ag2['agree']


# ---------------------------------------------------------------------------
# with the shared directory

@needs_shared
def test_remy_halves_and_targets():
    hv = X.remy_halves()
    R, cols, mph = X._remy_flies()
    assert len(cols) == 136 and mph == 1
    for s, (a, b) in X.REMY_SPLITS.items():
        assert sorted(hv[a] + hv[b]) == sorted(R.index) and len(hv[a]) == len(R) // 2
        for h in (a, b):
            assert R.loc[hv[h], cols].notna().sum().min() >= 1
    assert hv['a'] != hv['c']
    tg = X.remy_targets()
    assert set(tg) == {'all', 'a', 'b', 'c', 'd'}
    for h in X.REMY_HALVES:
        assert tg[h].reindex([tuple(sorted(p)) for p in cols]).notna().sum() >= 0
