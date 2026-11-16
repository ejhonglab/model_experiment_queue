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
    om = M._reading(resp, ('order_matched', True)).set_index('stimulus')
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


# ---------------------------------------------------------------------------
# part B1: the stimuli, their labels and the saved layout

def test_b1_specs_and_labels():
    from experiments import E100_apl_weak_odors as E100
    sp = M.b1_specs()
    labels = [x['label'] for x in sp]
    assert len(sp) == 52 and len(set(labels)) == 52
    assert labels[:16] == list(H.APL16_ODORS)
    roles = pd.Series([x['role'] for x in sp]).value_counts().to_dict()
    assert roles == {'panel': 16, 'V1': 4, 'V2': 1, 'beside': 9, 'E100': 12, 'hallem': 9, 'blank': 1}
    assert M.V2_NUM in labels and all(v in labels for v in M.V1_ESTIMATES.values())
    for x in sp:
        lab = E100.parse_b1_label(x['label'])
        if x['source'] not in ('panel', 'blank'):
            assert lab['odor'] == x['odor'] and lab['conc'] == pytest.approx(x['conc']) and lab['source'] == x['source']
        if x['role'] == 'E100' or x['role'] == 'V2':
            assert x['session'] in lab['source']                          # E100's match_b1 matches on the session's name
    assert M.session_label('paa', -5, 'validation2', M.DIAG) == E100.estimate_label(('paa', -5.0), 'validation2', M.DIAG)
    e100 = [x for x in sp if x['role'] == 'E100']
    assert {(x['odor'], x['conc']) for x in e100} == {('4h-ol', -6.0), ('farn', -2.0), ('2h', -6.0), ('2-but', -6.0), ('ma', -7.0)}


def _b1_rows(paa5=0.8, ms_scale=(1.0, 1.0, 1.0, 1.0, 1.0), pre=8.0, n=8):
    """Made-up B1 rows (`load_b1`'s columns) on n wirings: paa's evoked 12, each ms estimate's 10 x its scale, paa at
    10^-5 `paa5` x paa's evoked, the blank zero."""
    rows = []
    for i in range(n):
        w = M.ALL_WIRINGS[i]
        ev = {M.V1_NUM: 12.0, M.V2_NUM: 12.0 * paa5, **{lab: 10.0 * sc for lab, sc in zip(M.V1_ESTIMATES.values(), ms_scale)},
              **{lab: 12.0 * paa5 for lab in M.V2_BESIDE.values()}, M.BLANK_LABEL: 0.0}
        for st, e in ev.items():
            for cond, rate in (('on', 0.001), ('acute_off', 0.002)):
                rows.append(dict(wiring=w, cond=cond, stimulus=st, win=pre + e, pre=pre, evoked=e, response_rate=rate))
    return pd.DataFrame(rows)


def _flies(v1a=(0.67, 0.96), v1b=(1.0, 1.27), v2=(0.43, 0.59)):
    rows = [dict(ratio_id='V1_paa-4', lo90=v1a[0], hi90=v1a[1]), dict(ratio_id='V1_paa-2.954', lo90=v1b[0], hi90=v1b[1]),
            dict(ratio_id='V2', lo90=v2[0], hi90=v2[1])]
    return pd.DataFrame(rows)


def _readings(**over):
    """The flies' ratios under every reading of READING_SETS, the same unless `over` gives one ({key: _flies(...)})."""
    return {k: over.get(k, _flies()) for k in M.READING_SETS}


def test_held_against_reading_dependent():
    p = tuple(M.PARTNERS)
    assert M.held_category(dict(carrying='inside', first_in_block='inside', own_range='inside', without_three='inside',
                                literal='steeper'), p) == 'inside'                  # the literal reading decides nothing
    assert M.held_category(dict(carrying='inside', first_in_block='inside', own_range='not decided',
                                without_three='inside'), p) == 'reading-dependent'
    assert M.held_category(dict(carrying='not decided', first_in_block='not decided', own_range='not decided',
                                without_three='not decided'), p) == 'reading-dependent'
    # V2 has no V1-only partner
    assert M.held_category(dict(carrying='flatter', first_in_block='flatter', without_three='flatter'),
                           ('first_in_block', 'without_three')) == 'flatter'


def test_v1_rule_over_the_estimates():
    assert M.v1_from_estimates({'a': 'reading-dependent', 'b': 'inside'}) == 'inside'
    assert M.v1_from_estimates({'a': 'steeper', 'b': 'steeper'}) == 'steeper'
    assert M.v1_from_estimates({'a': 'steeper', 'b': 'flatter'}) == 'outside, mixed'
    assert M.v1_from_estimates({'a': 'steeper', 'b': 'reading-dependent'}) == 'not decided'


def test_verdict_core_held_and_outcomes():
    # paa / ms 1.2 (inside the hull 0.67 to 1.27 under every reading), paa at 10^-5 at 0.5 of paa: outcome (a)
    v, comps, rows = M.verdict_core(_readings(), _b1_rows(paa5=0.5))
    assert v['V1'] == 'inside' and v['V2'] == 'inside' and v['outcomes'] == ['a']
    assert v['V1_lines']['megamat/megamat (the panel)']['carrying'] == 'inside (8 of 8)'
    assert list(v['V1_lines'])[0] == 'megamat/megamat (the panel)'
    assert set(comps['reading_set']) == set(M.READING_SETS) and 'linear_minus_nearer_end' in rows.columns
    # a partner whose interval excludes the model: V1 reading-dependent under every estimate, not decided, outcome (e)
    v, _, _ = M.verdict_core(_readings(own_range=_flies(v1a=(0.6, 0.9), v1b=(0.9, 1.1))), _b1_rows(paa5=0.5))
    assert set(v['V1_held_per_estimate'].values()) == {'reading-dependent'} and v['V1'] == 'not decided'
    assert 'e' in v['outcomes'] and 'a' not in v['outcomes']
    # the V1-only partner does not touch V2; one estimate held inside is enough for V1
    v, _, _ = M.verdict_core(_readings(own_range=_flies(v2=(0.9, 1.0))), _b1_rows(paa5=0.5, ms_scale=(1.0, 0.4, 0.4, 0.4, 0.4)))
    assert v['V2'] == 'inside' and v['V1'] == 'inside'
    assert v['V1_held_per_estimate']['megamat/diag'] == 'steeper'
    # V2 flatter held: (b); V2 flatter under the carrying reading only: reading-dependent, (e)
    v, _, _ = M.verdict_core(_readings(), _b1_rows(paa5=0.85))
    assert v['V2'] == 'flatter' and 'b' in v['outcomes']
    v, _, _ = M.verdict_core(_readings(first_in_block=_flies(v2=(0.43, 0.90))), _b1_rows(paa5=0.85))
    assert v['V2_held'] == 'reading-dependent' and v['V2'] == 'not decided' and 'e' in v['outcomes']
    # every estimate held steeper: V1 outside, (d)
    v, _, _ = M.verdict_core(_readings(), _b1_rows(paa5=0.5, ms_scale=(0.5, 0.45, 0.4, 0.5, 0.42)))
    assert v['V1'] == 'steeper' and 'd' in v['outcomes']
    # the order factor: the nearer end over the model's nearest ratio
    assert v['V2_order_factor']['factor'] == 1.0


def test_verdict_core_pulse_and_blank():
    pulse = {'V1_paa-4': dict(disagree=False), 'V1_paa-2.954': dict(disagree=False), 'V2': dict(disagree=True)}
    v, _, _ = M.verdict_core(_readings(), _b1_rows(paa5=0.5), pulse, solvent=dict(ctrl=0.012, tnt=0.007))
    assert v['V2'] == 'not decided' and v['V1'] == 'inside' and 'e' in v['outcomes'] and 'h' not in v['outcomes']
    assert v['blank_response_rate']['acute_off']['median'] == pytest.approx(0.002)
    pulse['V1_paa-2.954'] = dict(disagree=True)
    v, _, _ = M.verdict_core(_readings(), _b1_rows(paa5=0.5), pulse)
    assert v['V1'] == 'not decided'
    v, _, _ = M.verdict_core(_readings(), _b1_rows(paa5=0.5), solvent=dict(ctrl=0.012, tnt=0.001))
    assert 'h' in v['outcomes']
    # fewer than six wirings agree: not decided under the carrying reading, so not held
    b = _b1_rows(paa5=0.5)
    b.loc[(b['stimulus'] == M.V2_NUM) & b['wiring'].isin(M.ALL_WIRINGS[:3]), 'win'] = 8.0 + 12.0 * 0.85
    v, _, _ = M.verdict_core(_readings(), b)
    assert v['V2'] == 'not decided' and 'e' in v['outcomes']


def test_pulse_permutation_count():
    from math import comb
    flies = [f'f{i}' for i in range(8)]
    x = pd.Series([0.0, 0.1, -0.1, 0.05, -0.05, 1.0, 1.1, 0.9], index=flies)   # the last three far from the first five
    r = M.pulse_permutation(x, flies[5:])
    assert r['n_splits'] == comb(8, 3) == 56 and r['n_short'] == 3
    assert r['n_as_large'] == 1 and r['share'] == pytest.approx(1 / 56) and r['disagree']
    y = pd.Series([0.0, 0.1, -0.1, 0.05, -0.05, 0.02, 0.08, -0.03], index=flies)
    r = M.pulse_permutation(y, flies[5:])
    assert r['n_as_large'] >= 6 and not r['disagree']
    # a fly without a value is left out: 7 flies, 2 short, C(7, 2) splits
    r = M.pulse_permutation(x.drop('f7'), flies[5:])
    assert r['n_splits'] == comb(7, 2) == 21


def test_stimulus_classes_refuses_different_pre_odor_levels():
    b = _b1_rows()
    b.loc[b['stimulus'] == M.V2_NUM, 'pre'] = 9.0
    with pytest.raises(RuntimeError):
        M.stimulus_classes(b[b['cond'] == 'on'], M.V2_NUM, M.V2_DEN, 0.4, 0.6)


def test_b1_save_and_load_round_trip(tmp_path):
    labels = list(H.APL16_ODORS) + [M.V2_NUM, M.BLANK_LABEL]
    kc = pd.Index(np.arange(5), name='kc_id')
    gl = pd.Index(['DL1', 'DM1'], name='glomerulus')
    for ct in ('on', 'acute_off'):
        o = dict(responses=pd.DataFrame(np.eye(5, len(labels)), index=kc, columns=labels),
                 spike_counts=pd.DataFrame(2 * np.eye(5, len(labels)), index=kc, columns=labels),
                 pn_win=pd.DataFrame(1.0, index=gl, columns=labels), pn_pre=pd.DataFrame(0.5, index=gl, columns=labels),
                 apl=pd.DataFrame(dict(win=np.arange(len(labels)) + 8.0, pre=8.0), index=pd.Index(labels, name='stimulus')))
        M.save_b1(tmp_path, 'real', ct, o)
    run = M.load_b1_run('real', [tmp_path])
    assert set(run) == {'on', 'acute_off'} and list(run['on']['apl'].index) == labels
    df = M.load_b1(['real', 'P1'], roots=[tmp_path])
    assert set(df['cond']) == {'on', 'acute_off'} and len(df) == 2 * len(labels)
    assert df.attrs['missing'] == [('P1', 'acute_off'), ('P1', 'on')]
    r = df[(df['stimulus'] == M.V2_NUM) & (df['cond'] == 'on')].iloc[0]
    assert (r['b1_odor'], r['b1_conc'], r['b1_source']) == ('paa', -5.0, 'validation2/diag')
    assert df[(df['stimulus'] == '2h') & (df['cond'] == 'on')]['response_rate'].iloc[0] == pytest.approx(0.2)
    # E100's own reader reads the same folder the same way
    from experiments import E100_apl_weak_odors as E100
    e = E100.read_e98_b1(['real'], dirs=[tmp_path])
    assert len(e) == 2 * len(labels) and set(e['stimulus']) == set(labels)


@pytest.mark.skipif(not _shared(), reason='shared directory not available')
def test_b1_input_structure():
    orn, t = M.b1_input()
    panel = H.apl16_orn_deltas()
    assert orn.shape == (36, 52) and orn.index.equals(panel.index)
    assert orn.iloc[:, :16].equals(panel.astype(float))
    assert list(t.loc[t['included'], 'label']) == M.input_labels(orn) and t['included'].all()
    assert (orn[(M.ADDED_PANEL, M.BLANK_LABEL)] == 0).all()
    o46, t46 = M.b1_input(include_hallem_lower=False)
    assert o46.shape == (36, 46)
    assert t46.loc[~t46['included'], 'label'].tolist() == [M.hallem_label(o, c) for o, c in M.HALLEM_STIMULI if c != -2.0]


def test_hallem_glomerulus_names_and_table():
    assert M.hallem_glomerulus('da4m') == 'DA4m' and M.hallem_glomerulus('vc3') == 'VC3' and M.hallem_glomerulus('vm5d') == 'VM5d'
    with pytest.raises(ValueError):
        M.hallem_glomerulus('dm3+dm5')
    t, info = M.hallem_hc_table()
    assert t.shape[1] == 23 and not t.columns.duplicated().any() and len(info['sha256']) == 64
    assert (t.loc['pentyl acetate -2'] == t.loc['pentyl acetate']).all()
    rows = M.hallem_rows()
    assert set(rows.columns) == set(M.HALLEM_STIMULI)


@pytest.mark.skipif(not _shared(), reason='shared directory not available')
def test_hallem_identity_against_e38():
    assert M.hallem_identity()['exact']


def test_with_a0_definitions():
    b = _b1_rows(paa5=0.5)
    b['pre_settled'] = b['pre'] + 0.5
    b.loc[b['stimulus'] == M.BLANK_LABEL, 'win'] = 8.0 + 0.5                 # the blank's window: the settled level
    e = M.with_a0(b, 'e90')
    assert np.allclose(e['pre'], b['pre'])
    st = M.with_a0(b, 'settled')
    assert np.allclose(st['pre'], 8.5) and np.allclose(st.loc[st['stimulus'] == M.BLANK_LABEL, 'evoked'], 0.0)
    bl = M.with_a0(b, 'blank')
    assert np.allclose(bl['pre'], 8.5) and np.allclose(bl.loc[bl['stimulus'] == M.BLANK_LABEL, 'evoked'], 0.0)
    with pytest.raises(ValueError):
        M.with_a0(b, 'other')
    with pytest.raises(RuntimeError):
        M.with_a0(b.assign(pre_settled=np.nan), 'settled')


def test_load_b1_all_takes_e101_where_e98_has_not_run(tmp_path):
    labels = list(H.APL16_ODORS) + [M.V2_NUM, M.BLANK_LABEL]
    kc = pd.Index(np.arange(5), name='kc_id')
    gl = pd.Index(['DL1'], name='glomerulus')
    own = tmp_path / 'own'
    own.mkdir()
    for ct in ('on', 'acute_off'):
        o = dict(responses=pd.DataFrame(np.eye(5, len(labels)), index=kc, columns=labels),
                 spike_counts=pd.DataFrame(np.eye(5, len(labels)), index=kc, columns=labels),
                 pn_win=pd.DataFrame(1.0, index=gl, columns=labels), pn_pre=pd.DataFrame(0.5, index=gl, columns=labels),
                 apl=pd.DataFrame(dict(win=np.arange(len(labels)) + 8.0, pre=8.0, pre_settled=8.5),
                                  index=pd.Index(labels, name='stimulus')))
        M.save_b1(own, 'real', ct, o)
    b = M.load_b1(['real'], roots=[own])
    e = pd.concat([b, b.assign(wiring='P1')], ignore_index=True)
    f = tmp_path / 'a3_apl.csv'
    e.to_csv(f, index=False)
    out = M.load_b1_all(['real', 'P1', 'P2'], roots=[own], e101_path=f)
    assert out.attrs['sources'] == {'real': 'E98', 'P1': 'E101'}
    assert out.attrs['missing'] == [('P2', 'acute_off'), ('P2', 'on')]
    assert out.attrs['e101']['overlap_check']['real']['max_abs_diff'] == 0.0
    e.loc[(e['wiring'] == 'real') & (e['stimulus'] == 'paa'), 'win'] += 1e-9
    e.to_csv(f, index=False)
    with pytest.raises(RuntimeError):
        M.load_b1_all(['real', 'P1'], roots=[own], e101_path=f)
