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
