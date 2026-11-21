#!/usr/bin/env python3
"""Tests for E114 (a compressive map from the ORN recordings to the input rates; experiments/E114_input_compression.py):
the two maps on made-up inputs, the gate, the classes, part C's choice and the outcome letters on made-up readings; with
the model's input files, the inputs (pow1 is the lab's input, E10_joint_fit.py's reader gives each input back)."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E114_input_compression')

QUICK = os.environ.get('QUICK') == '1'


def _orn(seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    x = rng.gamma(0.6, 40.0, size=(12, 5))
    x[0, 0], x[3, 2], x[5, 4] = -7.0, 0.0, -1.5
    return pd.DataFrame(x, index=[f'g{i}' for i in range(12)], columns=[f'o{j}' for j in range(5)])


# ---------------------------------------------------------------------------
# the maps

def test_power_map_identity_sum_and_compression():
    orn = _orn()
    p1 = X.power_map(orn, 1.0)
    assert np.array_equal(p1.to_numpy(), orn.to_numpy())
    pos = orn.to_numpy() > 0
    for e in (0.6, 0.9, 1.15):
        y = X.power_map(orn, e).to_numpy()
        np.testing.assert_allclose(np.where(pos, y, 0).sum(0), np.where(pos, orn.to_numpy(), 0).sum(0))
        assert np.array_equal(y[~pos], orn.to_numpy()[~pos])           # zero and negative deltas unchanged
    # compression: each odor's largest input falls below 1 and rises above it; the participation ratio moves the other way
    for j in range(orn.shape[1]):
        v = orn.iloc[:, j].to_numpy()
        lo, hi = X.power_map(orn, 0.6).iloc[:, j].to_numpy(), X.power_map(orn, 1.15).iloc[:, j].to_numpy()
        assert lo.max() < v.max() < hi.max()
        assert X.participation_ratio(lo) > X.participation_ratio(v) > X.participation_ratio(hi)


def test_saturating_map():
    orn = _orn(1)
    pos = orn.to_numpy() > 0
    y100, y400 = X.saturating_map(orn, 100.0).to_numpy(), X.saturating_map(orn, 400.0).to_numpy()
    for y in (y100, y400):
        np.testing.assert_allclose(np.where(pos, y, 0).sum(0), np.where(pos, orn.to_numpy(), 0).sum(0))
        assert np.array_equal(y[~pos], orn.to_numpy()[~pos])
    # order kept within an odor, and a smaller ceiling compresses more
    for j in range(orn.shape[1]):
        v = orn.iloc[:, j].to_numpy()
        m = v > 0
        assert np.array_equal(np.argsort(v[m]), np.argsort(y100[m, j]))
        assert y100[m, j].max() <= y400[m, j].max() <= v[m].max()
    # the formula on one odor: x / (1 + x / c), rescaled
    x = np.array([10.0, 100.0, 300.0])
    df = pd.DataFrame({'a': x})
    y = x / (1 + x / 200.0)
    np.testing.assert_allclose(X.saturating_map(df, 200.0)['a'].to_numpy(), y * x.sum() / y.sum())


def test_names_and_maps():
    assert X.pow_name(0.6) == 'pow0.6' and X.pow_name(1.0) == 'pow1' and X.sat_name(100.0) == 'sat100'
    orn = _orn(2)
    np.testing.assert_array_equal(X.map_of('pow0.7')(orn).to_numpy(), X.power_map(orn, 0.7).to_numpy())
    np.testing.assert_array_equal(X.map_of('sat200')(orn).to_numpy(), X.saturating_map(orn, 200.0).to_numpy())
    with pytest.raises(ValueError):
        X.map_of('I0')
    assert X.CONTROL == 'pow1.15' and X.CONTROL in X.A_INPUTS and X.CONTROL not in X.CANDIDATES


# ---------------------------------------------------------------------------
# the gate on made-up readings

UNITS = tuple(X.ALL_WIRINGS) + ('all8', 'S', 'T')


def _rd(spec: dict) -> pd.DataFrame:
    """Readings of made-up inputs: spec {name: dict of changes from a base}: d_on, d_off (ρ, every unit and target,
    unless 'wirings_on' lists per-wiring changes), pattern (calibrated, every target), raw_loss, mean (calibrated)."""
    rows = []
    for name, s in spec.items():
        for u in UNITS:
            for ro in ('raw', 'cal'):
                row = dict(input=name, unit=u, readout=ro)
                for c in ('on', 'off'):
                    d = s.get(f'd_{c}', 0.0)
                    if u in X.ALL_WIRINGS and f'wirings_{c}' in s:
                        d = s[f'wirings_{c}'][X.ALL_WIRINGS.index(u)]
                    for t in X.TARGETS:
                        row[f'rho_{c}_{t}'] = 0.70 + d + (s.get(f'half_{c}', 0.0) if t != 'all' else 0.0)
                for t in X.TARGETS:
                    row[f'{t}_pattern'] = 0.06 + (s.get('pattern', 0.0) if ro == 'cal' else 0.0)
                row['all_loss'] = 0.09 + (s.get('raw_loss', 0.0) if ro == 'raw' else 0.0)
                row['mean_effect'] = 0.12 + (s.get('mean', 0.0) if ro == 'cal' else 0.0)
                for k in ('amplitude', 'shape', 'sd_ratio', 'rho_effect', 'pearson_effect'):
                    row[f'all_{k}'] = 0.5
                for k in ('block', 'vcf_on', 'vcf_off', 'vcf_diff', 'struct_wt_on', 'struct_wt_off'):
                    row[k] = 0.3 + s.get(k, 0.0)
                rows.append(row)
    return pd.DataFrame(rows)


def _runs(spec: dict) -> pd.DataFrame:
    rows = []
    for name, s in spec.items():
        for i, w in enumerate(X.ALL_WIRINGS):
            rows.append(dict(input=name, unit=w, flagged=i < s.get('n_flagged', 0), realism_flag_zcrit='',
                             rate_rho_on=0.6 + s.get('rate_rho_on', 0.0), rate_rho_off=0.6 + s.get('rate_rho_off', 0.0)))
    return pd.DataFrame(rows)


BASE = {'I0': {}, 'pow1.15': dict(d_on=-0.01, d_off=-0.01)}


def _g(spec: dict, name: str, **kw) -> dict:
    sp = dict(BASE, **spec)
    return X.gate(_rd(sp), _runs(sp), name, **kw)


def test_gate_better():
    g = _g({'pow0.8': dict(d_on=0.012, d_off=0.0)}, 'pow0.8')
    assert g['gain_on'] and not g['gain_off'] and g['control_ok'] and X.class_of(g) == 'better'
    assert X.best_gain(g) == pytest.approx(0.012)


def test_gate_gain_needs_wirings_halves_and_the_other_rho():
    w = [0.02] * 5 + [0.0] * 3                                       # 5 of 8 wirings
    g = _g({'pow0.8': dict(d_on=0.012, wirings_on=w)}, 'pow0.8')
    assert g['n_wirings_gain_on'] == 5 and not g['gain']
    g = _g({'pow0.8': dict(d_on=0.012, half_on=-0.005)}, 'pow0.8')   # one half below GAIN
    assert not g['gain']
    g = _g({'pow0.8': dict(d_on=0.012, d_off=-0.006)}, 'pow0.8')     # the other ρ lower by more than 0.005
    assert not g['gain'] and X.class_of(g) == 'no gain'


def test_gate_price_named():
    g = _g({'pow0.8': dict(d_on=0.012, pattern=0.004)}, 'pow0.8')
    assert X.class_of(g) == 'a price' and any(k.startswith('pattern part') for k in g['fails'])
    g = _g({'pow0.8': dict(d_on=0.012, rate_rho_off=-0.06, n_flagged=1)}, 'pow0.8')
    assert X.class_of(g) == 'a price'
    assert 'response rate ρ, APL off' in g['fails'] and 'realism (runs outside the bands)' in g['fails']
    # the wrong-way control that gains too is a failure, named
    g = X.gate(_rd(dict(BASE, **{'pow1.15': dict(d_on=0.003), 'pow0.8': dict(d_on=0.012)})),
               _runs(dict(BASE, **{'pow1.15': {}, 'pow0.8': {}})), 'pow0.8')
    assert X.class_of(g) == 'a price' and g['control_ok'] is False
    # the block and the alignment are columns only
    g = _g({'pow0.8': dict(d_on=0.012, block=-0.2, vcf_on=-0.3)}, 'pow0.8')
    assert X.class_of(g) == 'better' and g['d_block'] == pytest.approx(-0.2)


def test_gate_equal_and_no_gain():
    g = _g({'sat200': dict(d_on=0.004, d_off=-0.004)}, 'sat200')
    assert X.class_of(g) == 'an equal'
    g = _g({'sat200': dict(d_on=0.004, raw_loss=0.006)}, 'sat200')
    assert X.class_of(g) == 'no gain'
    g = _g({}, 'pow0.6')
    assert not g['complete'] and X.class_of(g) == ''


def _gates(spec: dict) -> tuple:
    sp = dict(BASE, **spec)
    rd, runs = _rd(sp), _runs(sp)
    names = [n for n in sp if n not in ('I0',)]
    gates = {n: X.gate(rd, runs, n, control=None if n == X.CONTROL else X.CONTROL) for n in names}
    return gates, rd


def test_choose_c():
    gates, _ = _gates({'pow0.6': dict(d_on=0.015), 'pow0.7': dict(d_on=0.02), 'pow0.8': dict(d_on=0.03, pattern=0.01),
                       'pow0.9': dict(d_on=0.02)})
    ch = X.choose_c(gates)
    assert ch['exponent'] == 0.9                                     # tie with 0.7 -> nearer 1; 0.8 is 'a price'
    gates, _ = _gates({'pow0.6': {}, 'pow0.7': {}, 'pow0.8': {}, 'pow0.9': {}})
    assert X.choose_c(gates)['exponent'] == X.FALLBACK_EXPONENT


def test_outcomes():
    a = {'pow0.6': {}, 'pow0.7': {}, 'pow0.8': {}, 'pow0.9': {}, 'sat100': {}, 'sat200': {}, 'sat400': {}}
    # (a) the refit is better
    gates, rd = _gates(dict(a, refit=dict(d_on=0.012)))
    assert X.outcome_of(gates, rd, {'exponent': 0.75})['letter'] == 'a'
    # (b) a price on the pattern part at fixed scales and after the refit
    gates, rd = _gates(dict(a, **{'pow0.8': dict(d_on=0.012, pattern=0.004)}, refit=dict(d_on=0.012, pattern=0.004)))
    assert X.outcome_of(gates, rd, {'exponent': 0.8})['letter'] == 'b'
    # (b) with the fallback 0.75, read on its held points (part Cstart), while pow0.7 has a gain at fixed scales
    gates, rd = _gates(dict(a, **{'pow0.7': dict(d_off=0.012, pattern=0.004), 'pow0.75': dict(d_off=0.012, pattern=0.004)},
                            refit=dict(d_off=0.012, pattern=0.004)))
    assert X.outcome_of(gates, rd, {'exponent': 0.75})['letter'] == 'b'
    del gates['pow0.75']                                             # without its gate the letter cannot be (b)
    assert X.outcome_of(gates, rd, {'exponent': 0.75})['letter'] == '-'
    # (c) a gain on the hemibrain wiring only
    w = [0.02] + [0.0] * 7
    gates, rd = _gates(dict(a, **{'pow0.8': dict(wirings_on=w)}))
    assert X.outcome_of(gates, rd)['letter'] == 'c'
    # (d) no gain and the saturating form differs
    gates, rd = _gates(dict(a, sat100=dict(d_off=-0.02), sat200=dict(d_off=-0.02), sat400=dict(d_off=-0.03)))
    oc = X.outcome_of(gates, rd)
    assert oc['letter'] == 'd' and oc['provisional']
    # none of these
    gates, rd = _gates(a)
    assert X.outcome_of(gates, rd)['letter'] == '-'


# ---------------------------------------------------------------------------
# with the model's input files

def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


@pytest.mark.skipif(not _shared(), reason='needs the model input files')
def test_inputs_from_the_lab_input():
    lab = X.lab_input()
    assert np.array_equal(X.input_of('pow1').to_numpy(), lab.to_numpy())
    for n in X.A_INPUTS + X.B_INPUTS:
        x = X.input_of(n)
        assert x.index.equals(lab.index) and x.columns.equals(lab.columns)
        pos = lab.to_numpy() > 0
        np.testing.assert_allclose(np.where(pos, x.to_numpy(), 0).sum(0), np.where(pos, lab.to_numpy(), 0).sum(0))
    s = X.inputs_summary(('I0', 'pow0.6', 'pow1.15'))
    pr = s.groupby('input')['pr'].mean()
    assert pr['pow0.6'] > pr['I0'] > pr['pow1.15']
