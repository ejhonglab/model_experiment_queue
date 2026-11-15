#!/usr/bin/env python3
"""Tests for E98 (the strength of APL's response in the flies against the model's APL;
experiments/E98_apl_response_strength.py), on made-up maps and made-up APL values: the interpolation on log10
concentration and its range rule; the concentration matching to one log unit, inclusive; the fly bootstrap reproducing
with its seed and averaging a fly's recordings on the log scale; the two indicator forms at n 1 equal to the evoked
ratio; inside, flatter and steeper on made-up rows and the comparison's 6 of 8; V1's hull and its rule over the
estimates of ms's input; the flies' ratios with the response threshold and the 'missing' rule; the blank identity (the
first-in-block baseline with the blank subtracted is the order-matched one without it) on a made-up recording; the
power table's count of flies; and, with the shared directory, the reader against E79's saved tables."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
M = pytest.importorskip('experiments.E98_apl_response_strength')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return (H.shared_dir() / 'results' / 'E79_apl_calyx_constraints' / 'E79_maps_across_baselines.csv').exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# concentrations

def test_value_at_imaged_and_interpolated():
    concs, vals = [-6.0, -5.0, -3.0], [0.1, 0.2, 0.6]
    assert M.value_at(concs, vals, -5.0) == (0.2, 'imaged')
    v, how = M.value_at(concs, vals, -4.0)                    # halfway between -5 and -3 on log10 concentration
    assert how == 'interpolated' and v == pytest.approx(0.4)
    v, how = M.value_at(concs, vals, -5.75)
    assert how == 'interpolated' and v == pytest.approx(0.125)
    v, how = M.value_at(concs[::-1], vals[::-1], -4.0)         # order does not matter
    assert v == pytest.approx(0.4)


def test_value_at_range_rule():
    concs, vals = [-6.0, -5.0, -3.0], [0.1, 0.2, 0.6]
    # outside the range: the nearest end stands in within one log unit, inclusive
    assert M.value_at(concs, vals, -2.0) == (0.6, 'substituted')
    assert M.value_at(concs, vals, -7.0) == (0.1, 'substituted')
    v, how = M.value_at(concs, vals, -1.999)
    assert np.isnan(v) and how == 'outside'
    v, how = M.value_at(concs, vals, -7.01)
    assert np.isnan(v) and how == 'outside'
    # inside the series' design range but outside this recording's: left out, not substituted
    v, how = M.value_at(concs, vals, -2.5, design_range=(-6.0, -2.0))
    assert np.isnan(v) and how == 'missing'
    # outside the design range too: the recording's own end, within one log unit of it
    assert M.value_at(concs, vals, -2.0, design_range=(-6.0, -2.5)) == (0.6, 'substituted')
    with pytest.raises(ValueError):
        M.value_at([-5.0, -5.0], [0.1, 0.2], -5.0)


def test_concentration_match_one_log_unit_inclusive():
    assert M.concentration_match(-3.0, [-4.0, -6.0]) == -4.0                 # exactly one unit: allowed
    assert M.concentration_match(-3.0, [-4.0001]) is None
    assert M.concentration_match(-2.0, [-2.954, -3.431]) == -2.954
    assert M.concentration_match(-3.0, [-4.386]) is None                     # pa: the model's 10^-3, the maps' 10^-4.386
    assert M.concentration_match(-8.0, [-7.725, -6.294]) == -7.725
    assert M.concentration_match(-3.0, []) is None


# ---------------------------------------------------------------------------
# the flies' ratio

def _per_rec():
    rng = np.random.default_rng(3)
    rows = []
    for f in range(6):
        for r in range(1 + f % 2):
            rows.append(dict(fly=f'fly{f}', log_ratio=float(np.log(1.3) + 0.2 * rng.standard_normal())))
    return pd.DataFrame(rows)


def test_fly_ratio_reproduces_with_its_seed():
    d = _per_rec()
    a, b = M.fly_ratio(d, n_boot=2000, seed=7), M.fly_ratio(d, n_boot=2000, seed=7)
    assert a == b
    c = M.fly_ratio(d, n_boot=2000, seed=8)
    assert c['ratio'] == a['ratio'] and (c['lo90'], c['hi90']) != (a['lo90'], a['hi90'])
    assert a['lo90'] < a['ratio'] < a['hi90'] and a['n_flies'] == 6 and a['n_recordings'] == 9


def test_fly_ratio_averages_a_flys_recordings_on_the_log_scale():
    d = pd.DataFrame(dict(fly=['a', 'a', 'b'], log_ratio=np.log([2.0, 8.0, 1.0])))
    r = M.fly_ratio(d, n_boot=100)
    # fly a: geometric mean 4; then the geometric mean over flies of 4 and 1
    assert r['ratio'] == pytest.approx(2.0)
    assert r['sd_log2_flies'] == pytest.approx(np.std([2.0, 0.0], ddof=1))
    e = M.fly_ratio(pd.DataFrame(dict(fly=[], log_ratio=[])))
    assert np.isnan(e['ratio']) and e['n_flies'] == 0


def _resp(rows):
    df = pd.DataFrame(rows)
    df['baseline'], df['blank_subtracted'], df['pulse_s'], df['blank'] = 'order_matched', True, 1.0, 0.0
    df['stimulus'] = [f'{o} {c:g}' for o, c in zip(df['odor'], df['conc'])]
    return df


def test_flies_ratios_threshold_and_missing_rule():
    rows = []
    for f, (paa, ms) in enumerate([(0.3, 0.2), (0.4, 0.2), (0.25, 0.25)]):
        for o, cs, v in (('paa', (-5.0, -3.0), paa), ('ms', (-4.0, -2.0), ms)):
            for c in cs:
                rows.append(dict(recording=f'r{f}', fly=f'fly{f}', odor=o, conc=c, response=v))
    # a fourth recording without ms at 10^-2: ms at 10^-3 is inside the series' range but not its own
    rows += [dict(recording='r3', fly='fly3', odor='paa', conc=c, response=0.3) for c in (-5.0, -3.0)]
    rows += [dict(recording='r3', fly='fly3', odor='ms', conc=-4.0, response=0.2)]
    # a fifth whose ms is below the threshold
    rows += [dict(recording='r4', fly='fly4', odor='paa', conc=c, response=0.3) for c in (-5.0, -3.0)]
    rows += [dict(recording='r4', fly='fly4', odor='ms', conc=c, response=0.015) for c in (-4.0, -2.0)]
    resp = _resp(rows)
    ratios = (('x', ('paa', -3.0), ('ms', -3.0), 'V1'),)
    s, per = M.flies_ratios(resp, ('order_matched', True), 'all', ratios)
    assert s.iloc[0]['n_flies'] == 3 and s.iloc[0]['n_recordings_left_out'] == 2
    assert s.iloc[0]['ratio'] == pytest.approx(np.exp(np.mean(np.log([1.5, 2.0, 1.0]))))
    p = per.set_index('recording')
    assert p.loc['r3', 'den_how'] == 'missing' and not p.loc['r3', 'entered']
    assert not p.loc['r4', 'entered']
    # the per-recording rule substitutes ms at 10^-4 for 10^-3 in r3
    s2, per2 = M.flies_ratios(resp, ('order_matched', True), 'all_recording_ranges', ratios)
    assert per2.set_index('recording').loc['r3', 'den_how'] == 'substituted' and s2.iloc[0]['n_flies'] == 4


# ---------------------------------------------------------------------------
# the blank identity on a made-up recording

def _made_up_recording(seed: int = 0):
    from experiments import E79_apl_calyx_constraints as E79
    rng = np.random.default_rng(seed)
    stims = [('ms', -3.0), ('pa', -4.0), ('paa', -5.0)]
    rows = []
    for b, (o, c) in enumerate(stims):
        rows += [dict(block=b, odor='pfo', log10_conc=np.nan)] + [dict(block=b, odor=o, log10_conc=c)] * 2
    rows += [dict(block=len(stims), odor='pfo', log10_conc=np.nan)] * 3
    fr = pd.DataFrame(rows)
    fr.index.name = 'frame'
    fr = E79.add_positions(fr)
    stack = rng.normal(0.05, 0.02, size=(len(fr), 8, 9))
    for i, (o, c) in zip(fr.index, zip(fr['odor'], fr['log10_conc'])):
        if o != 'pfo':
            stack[i] += 0.1 * (c + 7)
    stack[fr.index[-2:]] += 0.08                    # the paraffin oil block's 2nd and 3rd respond: a blank
    mask = np.zeros((8, 9), bool)
    mask[1:7, 2:8] = True
    return E79.Recording(name='rec', side='R', fly='fly', stack=stack, frames=fr, mask=mask)


def test_blank_identity_on_a_made_up_recording():
    rec = _made_up_recording()
    resp = M.stimulus_responses([rec], pulses={'rec': 1.0})
    assert len(resp) == 3 * len(M.BASELINES) * 2
    assert resp['blank'].iloc[0] > 0.05
    assert M.blank_identity(resp) < 1e-12
    om = M._reading(resp, M.MAIN_READING).set_index('stimulus')
    raw = M._reading(resp, ('order_matched', False)).set_index('stimulus')
    assert np.allclose(om['response'], raw['response'] - raw['blank'])
    assert np.allclose(raw['response'], (raw['rep1'] + raw['rep2']) / 2)


# ---------------------------------------------------------------------------
# the model's ratio and the comparison

def test_indicator_forms_at_n1_equal_the_evoked_ratio():
    rng = np.random.default_rng(0)
    for _ in range(20):
        a0 = rng.uniform(1, 10)
        a1, a2 = a0 + rng.uniform(0.5, 20, size=2)
        ev = (a1 - a0) / (a2 - a0)
        assert M.indicator_ratio(a1, a2, a0, 1.0, 'baseline') == pytest.approx(ev)
        assert M.indicator_ratio(a1, a2, a0, 1.0, 'evoked') == pytest.approx(ev)
    assert M.indicator_ratio(3.0, 2.0, 1.0, 2.0, 'baseline') == pytest.approx((9 - 1) / (4 - 1))
    assert M.indicator_ratio(3.0, 2.0, 1.0, 2.0, 'evoked') == pytest.approx(4.0)
    rs = M.model_ratio_set(3.0, 2.0, 1.0)
    assert len(rs) == len(M.INDICATOR_NS) * 2 and rs['ratio'].iloc[0] == pytest.approx(2.0)
    assert M.INDICATOR_NS[0] == 1.0 and M.INDICATOR_NS[-1] == 2.5 and len(M.INDICATOR_NS) == 7
    with pytest.raises(ValueError):
        M.indicator_ratio(3.0, 2.0, 1.0, 1.0, 'other')


def test_classify_inside_flatter_steeper():
    # an interval above 1 (the flies' APL grows from stimulus 2 to 1)
    assert M.classify([1.1, 1.6, 2.5], 1.5, 2.0) == 'inside'
    assert M.classify([1.1, 1.3], 1.5, 2.0) == 'flatter'
    assert M.classify([0.9, 1.2], 1.5, 2.0) == 'flatter'                 # the other side of 1
    assert M.classify([2.1, 3.0], 1.5, 2.0) == 'steeper'
    assert M.classify([1.2, 2.2], 1.5, 2.0) == 'mixed'
    # an interval below 1 (V2: paa at 10^-5 over paa at 10^-4)
    assert M.classify([0.8, 0.9], 0.43, 0.6) == 'flatter'
    assert M.classify([1.1], 0.43, 0.6) == 'flatter'
    assert M.classify([0.3, 0.4], 0.43, 0.6) == 'steeper'
    assert M.classify([0.5], 0.43, 0.6) == 'inside'
    assert M.classify([0.43], 0.43, 0.6) == 'inside'                      # the ends are inside
    # an interval holding 1: no 'flatter'
    assert M.classify([1.3, 1.5], 0.8, 1.25) == 'steeper'
    assert M.classify([0.5, 0.7], 0.8, 1.25) == 'steeper'
    assert M.classify([0.5, 1.5], 0.8, 1.25) == 'mixed'
    with pytest.raises(ValueError):
        M.classify([0.0], 0.8, 1.25)


def test_comparison_verdict_six_of_eight():
    assert M.comparison_verdict(['inside'] * 6 + ['flatter'] * 2) == 'inside'
    assert M.comparison_verdict(['inside'] * 5 + ['flatter'] * 3) == 'not decided'
    assert M.comparison_verdict(['flatter'] * 7 + ['mixed']) == 'flatter'
    assert M.comparison_verdict(['steeper'] * 6 + ['inside'] * 2) == 'steeper'
    assert M.comparison_verdict(['inside'] * 4 + ['steeper'] * 4) == 'not decided'


def test_v1_hull_and_its_rule_over_the_estimates():
    assert M.hull([(0.67, 0.96), (1.0, 1.27)]) == (0.67, 1.27)
    assert M.hull([(1.0, 2.0)]) == (1.0, 2.0)
    with pytest.raises(ValueError):
        M.hull([(1.0, np.nan)])
    assert M.v1_verdict({'megamat': 'steeper', 'diag1': 'inside'}) == 'inside'
    assert M.v1_verdict({'megamat': 'flatter', 'diag1': 'flatter', 'control': 'flatter'}) == 'flatter'
    assert M.v1_verdict({'megamat': 'steeper', 'diag1': 'steeper'}) == 'steeper'
    assert M.v1_verdict({'megamat': 'steeper', 'diag1': 'not decided'}) == 'not decided'
    assert M.v1_verdict({'megamat': 'steeper', 'diag1': 'flatter'}) == 'not decided'


def test_wiring_classes_on_made_up_apl():
    rows = []
    for i, w in enumerate(M.ALL_WIRINGS):
        for o, a in (('paa', 20.0 + 0.05 * i), ('ms', 18.0)):
            rows.append(dict(wiring=w, odor=o, a=a, a0=8.0, evoked=a - 8.0))
    wc = M.wiring_classes(pd.DataFrame(rows), 'paa', 'ms', 1.15, 1.25)
    assert list(wc['wiring']) == list(M.ALL_WIRINGS)
    assert wc['ratio_linear'].iloc[0] == pytest.approx(12.0 / 10.0)
    assert (wc['cls'] == 'inside').all() and M.comparison_verdict(wc['cls']) == 'inside'
    assert (wc['ratio_max'] >= wc['ratio_linear']).all()
    wc = M.wiring_classes(pd.DataFrame(rows), 'paa', 'ms', 1.15, 1.2125)    # linear ratios 1.2 to 1.235 by 0.005
    assert list(wc['cls']) == ['inside'] * 3 + ['steeper'] * 5 and M.comparison_verdict(wc['cls']) == 'not decided'


def test_flies_needed():
    n = M.flies_needed(0.5, 0.4)
    assert n == float(np.ceil(((1.6449 + 0.8416) * 0.4 / 0.5) ** 2))
    assert M.flies_needed(-0.5, 0.4) == n
    assert M.flies_needed(5.0, 0.3) == 2.0                                 # at least two flies
    assert M.flies_needed(0.0, 0.3) == float('inf')
    assert M.flies_needed(0.25, 0.4) > M.flies_needed(0.5, 0.4)


# ---------------------------------------------------------------------------
# with the shared directory: the reader against E79's saved tables

@pytest.mark.skipif(not _shared(), reason='shared directory (E79 results, the calyx maps) not available')
def test_reader_reproduces_e79():
    res = M.reader_identity()
    assert res['exact'], res
    assert res['maps_order_matched']['n'] == 190
