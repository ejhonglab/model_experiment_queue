"""Tests of E121 step 0's reader (no simulation)."""
import numpy as np
import pandas as pd
import pytest

from experiments import E121_onset_race as E


def test_apl_rise_exponential():
    dt = 0.0005
    t = np.arange(-0.5, 0.75, dt)
    tau = 0.05
    trace = np.where(t > 0, 1 - np.exp(-np.clip(t, 0, None) / tau), 0.0) + 0.2
    amp, t90 = E.apl_rise(trace, t, plateau='late')
    # (measured from the onset value, 0.2)
    assert amp == pytest.approx(1.0, abs=1e-3)
    # 90 % of a plateau reached at tau * ln(10)
    assert t90 == pytest.approx(tau * np.log(10), abs=2 * dt)
    amp_m, t90_m = E.apl_rise(trace, t, plateau='max')
    assert t90_m >= t90 - dt
    assert np.isnan(E.apl_rise(np.zeros_like(t), t)[1])


def test_apl_rise_transient_late_plateau_lower():
    t = np.arange(-0.1, 0.75, 0.001)
    trace = np.interp(t, [-0.1, 0, 0.1, 0.3, 0.5, 0.75], [0, 0, 2.0, 1.0, 1.0, 0])
    amp, t90 = E.apl_rise(trace, t, plateau='late')
    assert amp == pytest.approx(1.0, abs=1e-6)
    assert t90 == pytest.approx(0.045, abs=0.0015)
    amp_m, t90_m = E.apl_rise(trace, t, plateau='max')
    assert amp_m == pytest.approx(2.0, abs=1e-6) and t90_m == pytest.approx(0.09, abs=0.0015)


def test_second_spike_and_shares():
    spikes = pd.DataFrame({'kc': [0, 0, 0, 1, 2, 2, 0, 0], 'odor': [0, 0, 0, 0, 0, 0, 1, 1],
                           'time_s': [0.03, 0.01, 0.2, 0.05, 0.08, 0.30, 0.02, 0.04]})
    c = E.second_spike_times(spikes)
    got = {(r.kc, r.odor): (r.t_complete, r.n_spikes) for r in c.itertuples()}
    assert got == {(0, 0): (0.03, 3), (2, 0): (0.30, 2), (0, 1): (0.04, 2)}
    comp = E.completion_table(spikes, ['IaA', '1-8ol'], n_kcs=3)
    assert dict(zip(zip(comp.kc, comp.odor), comp.breadth)) == {(0, 0): 2, (2, 0): 1, (0, 1): 2}
    rise = pd.DataFrame({'odor_name': ['IaA', '1-8ol'], 'plateau': ['late', 'late'], 't90': [0.1, 0.03]})
    m = E.mark_before(comp, rise, 'late')
    sh = E.shares(m).set_index(['by', 'level'])
    assert sh.loc[('overall', 'all'), 'share'] == pytest.approx(1 / 3)
    assert sh.loc[('group', 'ester'), 'share'] == pytest.approx(0.5)
    assert sh.loc[('group', 'alcohol'), 'share'] == 0
    assert sh.loc[('breadth', '1'), 'n_responses'] == 1


def test_groups_and_wirings():
    assert E.odor_group('ms') == 'ester' and E.odor_group('6al') == 'alcohol' and E.odor_group('paa') == 'other'
    assert E.odor_name('2h @ -3') == '2h'
    assert set(E.WIRINGS) == set(E.S_SET) | set(E.T_SET) and len(E.WIRINGS) == 8


def test_apl_level_at():
    t = np.round(np.arange(-0.01, 0.6, 0.001), 6)
    tr = np.interp(t, [-0.01, 0, 0.1, 0.4, 0.5, 0.6], [1, 1, 3.0, 2.0, 2.0, 1])
    apl = pd.DataFrame([tr], index=['eb'], columns=t)
    comp = pd.DataFrame({'odor_name': ['eb', 'eb'], 't_complete': [0.05, 0.1]})
    c = E.apl_level_at(comp, apl)
    assert c['apl_frac_peak'].tolist() == pytest.approx([0.5, 1.0])
    assert c['apl_frac_late'].tolist() == pytest.approx([1.0, 2.0])


# --- part reads (the planner's reads, 2026-11-23) ---

def _trace():
    """A trace at rest 1 before onset, rising linearly to 4 at 0.08 s, falling to 2 by 0.3 s and staying there."""
    t = np.round(np.arange(-0.01, 0.75, 0.0005), 6)
    tr = np.interp(t, [-0.01, 0, 0.08, 0.3, 0.5, 0.75], [1, 1, 4.0, 2.0, 2.0, 1.0])
    return t, tr


def test_apl_levels_count_the_level_before_the_odor():
    t, tr = _trace()
    lv = E.apl_levels(tr, t)
    assert lv['rest'] == pytest.approx(1.0)
    assert lv['peak'] == pytest.approx(4.0) and lv['t_peak'] == pytest.approx(0.08)
    # (the mean over the last 0.1 s of the pulse, 0.4 < t <= 0.5)
    assert lv['sustained'] == pytest.approx(2.0)
    assert lv['rest_over_peak'] == pytest.approx(0.25) and lv['peak_over_sustained'] == pytest.approx(2.0)
    with pytest.raises(ValueError):
        E.apl_levels(tr, t[:-1])


def test_level_at_reads_the_last_step_at_or_before():
    t = np.array([-0.001, -0.0005, 0.0, 0.0005, 0.001])
    tr = np.array([10.0, 11.0, 12.0, 13.0, 14.0])
    # (a spike time a hair under a step's reads that step; between steps, the earlier one; before the first, NaN)
    got = E.level_at(tr, t, [0.0005 - 1e-12, 0.0007, 0.001, -0.002])
    assert got[:3].tolist() == [13.0, 13.0, 14.0] and np.isnan(got[3])


def test_auc():
    assert E.auc([3, 4], [1, 2]) == 1.0
    assert E.auc([1, 2], [3, 4]) == 0.0
    assert E.auc([1, 1], [1, 1]) == 0.5
    # (pairs: 2>1 yes, 2>3 no, 4>1 yes, 4>3 yes: 3 of 4)
    assert E.auc([2, 4], [1, 3]) == pytest.approx(0.75)
    assert np.isnan(E.auc([], [1]))


def test_window_reads_half_open_windows():
    spikes = pd.DataFrame({
        'kc':     [0, 0, 0, 1, 1, 2, 2, 2, 2, 3],
        'odor':   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        # KC 0: two spikes in the first 100 ms, one at 0.1 (outside it); KC 1: two at 0.25 and 0.3 (last 250 ms);
        # KC 2: one early, one at 0.49, two at and after the pulse's end; KC 3: one spike, not a responder
        'time_s': [0.03, 0.05, 0.1, 0.25, 0.3, 0.04, 0.49, 0.5, 0.6, 0.02]})
    w = E.window_reads(spikes)
    assert w['n_responders'] == 3 and w['n_spikes'] == 10
    assert w['n_first100'] == 1 and w['n_last250'] == 1 and w['n_after'] == 1
    assert w['share_spikes_first100'] == pytest.approx(4 / 10)
    assert w['share_spikes_after'] == pytest.approx(2 / 10)
    assert w['share_responders_first100'] == pytest.approx(1 / 3)
    assert w['share_exactly_two'] == pytest.approx(1 / 3)
    assert w['median_spikes'] == 3 and w['max_spikes'] == 4
    # (second spikes 0.05, 0.3, 0.49)
    assert w['t2_median'] == pytest.approx(0.3)


def test_levels_at_completion_and_summary():
    t, tr = _trace()
    apl = pd.DataFrame([tr, tr], index=['eb', '1-8ol'], columns=t)
    lv = pd.DataFrame({o: E.apl_levels(tr, t) for o in apl.index}).T
    comp = pd.DataFrame({'kc': [0, 1, 2], 'odor': [0, 0, 1], 'odor_name': ['eb', 'eb', '1-8ol'],
                         'group': ['ester', 'ester', 'alcohol'], 'breadth_bin': ['1', '2-3', '1'],
                         't_complete': [0.02, 0.08, 0.04]})
    c = E.levels_at_completion(comp, apl, lv)
    # (at 0.02 s the trace is 1.75: 0.4375 of the peak, 0.875 of the sustained level; at the peak's own time, 1)
    assert c['frac_peak'].tolist() == pytest.approx([1.75 / 4, 1.0, 2.5 / 4])
    assert c['frac_sustained'].tolist() == pytest.approx([1.75 / 2, 2.0, 1.25])
    assert c['below_sustained'].tolist() == [True, False, False]
    assert c['before_peak'].tolist() == [True, False, True]
    s = E.level_summary(c).set_index(['by', 'level'])
    assert s.loc[('overall', 'all'), 'n'] == 3
    assert s.loc[('overall', 'all'), 'frac_peak_median'] == pytest.approx(2.5 / 4)
    assert s.loc[('group', 'ester'), 'share_below_sustained'] == pytest.approx(0.5)
    assert s.loc[('breadth', '1'), 'share_before_peak'] == 1.0
    assert s.loc[('odor', '1-8ol'), 'frac_sustained_median'] == pytest.approx(1.25)


def test_survival_reads():
    # APL off: four responders; (0,0) and (1,0) also respond with APL on (survivors), (2,0) and (3,1) do not
    off = pd.DataFrame({'kc': [0, 1, 2, 3], 'odor': [0, 0, 0, 1], 't_complete': [0.030, 0.040, 0.044, 0.050],
                        'n_spikes': [6, 3, 4, 2]})
    on = pd.DataFrame({'kc': [0, 1, 5], 'odor': [0, 0, 0], 't_complete': [0.035, 0.045, 0.04], 'n_spikes': [2, 2, 3]})
    r = E.survival_reads(on, off)
    assert (r['n_off_responders'], r['n_survivors'], r['n_removed']) == (4, 2, 2)
    # (survivors' second spikes 30, 40 ms all earlier than the removed's 44, 50 ms)
    assert r['auc_t2'] == 1.0
    # (counts 6, 3 against 4, 2: 6>4, 6>2, 3>2, 3<4)
    assert r['auc_count'] == pytest.approx(0.75)
    assert r['median_t2_survivors'] == pytest.approx(0.035) and r['median_t2_removed'] == pytest.approx(0.047)
    assert r['share_early_survivors'] == 1.0 and r['share_early_removed'] == pytest.approx(0.5)
    assert (r['n_on_responders'], r['n_on_not_off']) == (3, 1) and r['share_on_not_off'] == pytest.approx(1 / 3)


def test_check_planner_within_rounding():
    checks = {'a': ((0.80,), 1, 2), 'b': ((69, 72), 100, 0), 'c': ((41,), 1000, 0)}
    rows = {r['name']: r for r in E.check_planner({'a': (0.8049,), 'b': (0.6849, 0.7214), 'c': (0.0405,)}, checks)}
    assert rows['a']['agrees'] and rows['c']['agrees']
    assert not rows['b']['agrees']
    with pytest.raises(ValueError):
        E.check_planner({'a': (0.8, 0.9), 'b': (0.7, 0.7), 'c': (0.04,)}, checks)
    # every check the reads part fills has a value per planner's value
    assert all(len(v[0]) >= 1 for v in E.PLANNER_CHECKS.values())
    assert len(E.PLACEMENT_CHECKS['family_n_on_not_off_archive'][0]) == len(E.WIRINGS)


def test_on_not_off_counts():
    idx = pd.Index([10, 11, 12], name='kc_id')
    on = pd.DataFrame([[1, 0], [1, 1], [0, 0]], index=idx, columns=['a', 'b'])
    off = pd.DataFrame([[1, 1], [0, 1], [0, 1]], index=idx, columns=['a', 'b'])
    r = E.on_not_off_counts(on, off)
    # (on: (10,a), (11,a), (11,b); of these only (11,a) is not a response with APL off)
    assert (r['n_on_responders'], r['n_off_responders'], r['n_on_not_off']) == (3, 4, 1)
    assert r['share_on_not_off'] == pytest.approx(1 / 3)
    with pytest.raises(ValueError):
        E.on_not_off_counts(on, off[['b', 'a']])
    assert E.e110_arm_a_variant('S') == 'A_real' and E.e110_arm_a_variant('M2') == 'A_M2'


def test_kc_apl_weight_sums_claws_in_kc_order():
    mi = pd.MultiIndex.from_tuples([(7, 0), (7, 1), (5, 0), (9, 0)], names=['kc_id', 'claw_id'])
    w = pd.DataFrame({'wAPLKC': [1.0, 2.0, 0.5, 4.0]}, index=mi)
    assert E.kc_apl_weight(w, [5, 7, 9]).tolist() == [0.5, 3.0, 4.0]
    with pytest.raises(ValueError):
        E.kc_apl_weight(w, [5, 8])


def test_survival_by_weight():
    # APL off: four responders; KCs 0 and 1 survive (respond with APL on too), KCs 2 and 3 are removed
    off = pd.DataFrame({'kc': [0, 1, 2, 3], 'odor': [0, 0, 0, 1], 't_complete': [0.03, 0.04, 0.044, 0.05],
                        'n_spikes': [6, 3, 4, 2]})
    on = pd.DataFrame({'kc': [0, 1, 5], 'odor': [0, 0, 0], 't_complete': [0.035, 0.045, 0.04], 'n_spikes': [2, 2, 3]})
    weight = np.array([0.1, 0.3, 0.2, 0.9, 1.0, 1.5])
    breadth = np.array([1.0, 3.0, 3.0, 2.0, 0.0, 0.0])
    r = E.survival_by_weight(on, off, weight, breadth)
    # (survivors' weights 0.1, 0.3 against 0.2, 0.9, lower favouring survival: 0.1<0.2, 0.1<0.9, 0.3>0.2, 0.3<0.9)
    assert r['auc_apl_weight'] == pytest.approx(0.75)
    # (breadth 1, 3 against 3, 2: 1<3, 1<2, 3=3 half, 3>2)
    assert r['auc_breadth'] == pytest.approx(0.625)
    assert r['mean_kc_weight'] == pytest.approx(4.0 / 6)
    assert r['weight_survivors_over_mean'] == pytest.approx(0.2 / (4.0 / 6))
    assert r['weight_removed_over_mean'] == pytest.approx(0.55 / (4.0 / 6))
