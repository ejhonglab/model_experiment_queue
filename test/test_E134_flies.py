"""Tests of experiments/E134_flies.py (step 1, the flies alone) on planted populations: the populations' rows (floor
below one fly against the others), validation2's odor set, a step between two stages with no change, the no-step
yardstick and the megamat-against-control reading. The data module's populations are replaced by made-up ones."""
import numpy as np
import pandas as pd
import pytest

from data_prep import stage_population_data as D
from experiments import E134_flies as E
from kcq import population_fit as P

SPEC = dict(units=D.UNITS_BASELINE, criterion=1.0, window={}, notes=[], recording_site='somata')


def _pop(key, n_flies=6, n_units=60, seed=0, odors=None, gain=1.0, missing=None):
    rng = np.random.default_rng(seed)
    odors = list(D.PANEL_ODORS[key[1]] if odors is None else odors)
    tuning = rng.gamma(0.6, 1.0, size=(n_units, len(odors)))   # shared across flies: the flies agree on the pattern
    rows = []
    for f in range(n_flies):
        X = gain * tuning[None] + rng.normal(scale=0.3, size=(3, n_units, len(odors)))
        if missing and f in missing:
            X[:, :, [odors.index(o) for o in missing[f]]] = np.nan
        rows.append(D.long_from_trials(f'{key[0]}_f{f}', [str(u) for u in range(n_units)], odors, X))
    d = pd.concat(rows, ignore_index=True)
    return D.population_from_long(key, d, spec=dict(SPEC, recording_site=D.RECORDING_SITE[key[0]]))


@pytest.fixture
def planted(monkeypatch):
    pops = {}
    for i, key in enumerate(D.POPULATIONS):
        n = 1 if key[2] == 'same_fly' else 6
        miss = {0: ['paa', 'euc'], 1: ['paa'], 2: ['euc']} if key == ('kc_soma', 'validation2', 'wild_type') else None
        pops[key] = _pop(key, n_flies=n, seed=i, missing=miss)
    monkeypatch.setattr(E, 'recorded', lambda key: pops[key])
    return pops


def test_population_rows_floor_below_one_fly(planted):
    v, loo = E.population_rows(('kc_soma', 'apl16', 'control'), n_boot=100)
    err = v[(v.level == 'population') & (v.name == 'error')].iloc[0]
    assert err.floor < err.loo_mean
    assert 0.6 < err.loo_mean < 1.6                              # each fly against the others is about 1 by construction
    parts = v[v.level == 'part']
    assert set(parts.name) >= {'strength_rate', 'sharing_who', 'breadth_bins'}
    assert parts.one_fly_unit.notna().all()
    assert len(loo) == 6 and 'error' in loo


def test_v2_odor_set_largest_shared(planted):
    flies = E.flies_of(('kc_soma', 'validation2', 'wild_type'))
    odors, keep = E.v2_odor_set(flies, min_flies=4)
    assert 'paa' not in odors or 'euc' not in odors
    assert len(keep) >= 4
    assert all(set(odors) <= set(f.odors) for f in flies if f.name in keep)
    odors3, keep3 = E.v2_odor_set(flies, min_flies=6)
    assert len(odors3) <= len(odors) and len(keep3) == 6
    with pytest.raises(ValueError):
        E.v2_odor_set(flies, min_flies=7)


def test_no_change_between_identical_stages(planted, monkeypatch):
    a = _pop(('kc_claw', 'megamat', 'wild_type'), seed=5)
    b = _pop(('kc_soma', 'megamat', 'wild_type'), seed=5)       # the same tuning and noise: no step
    pops = {('kc_claw', 'megamat', 'wild_type'): a, ('kc_soma', 'megamat', 'wild_type'): b}
    monkeypatch.setattr(E, 'recorded', lambda key: pops[key])
    monkeypatch.setattr(E, 'STEPS', ((('kc_claw', 'megamat', 'wild_type'), ('kc_soma', 'megamat', 'wild_type')),))
    st, ss = E.recorded_steps(n_boot=50)
    assert np.allclose(ss.rms_change, 0.0)
    assert (ss.n_interval_excludes_0 == 0).all()
    ns = E.no_step_yardstick(n_boot=50)
    e = ns[ns.name == 'error'].iloc[0]
    assert e.no_step < 1e-9                                     # the upstream flies predict the downstream exactly


def test_step_with_a_gain_moves_the_mean_part_not_the_rate(planted, monkeypatch):
    a = _pop(('pn_bouton', 'megamat', 'wild_type'), seed=7)
    b = _pop(('kc_claw', 'megamat', 'wild_type'), seed=8)
    pops = {('pn_bouton', 'megamat', 'wild_type'): a, ('kc_claw', 'megamat', 'wild_type'): b}
    monkeypatch.setattr(E, 'recorded', lambda key: pops[key])
    monkeypatch.setattr(E, 'STEPS', ((('pn_bouton', 'megamat', 'wild_type'), ('kc_claw', 'megamat', 'wild_type')),))
    st, ss = E.recorded_steps(n_boot=50)
    assert set(ss.part) >= {'strength/rate', 'sharing/who', 'breadth/bins'}
    assert (ss.rms_change > 0).all()                            # different tunings: every part moves
    ns = E.no_step_yardstick(n_boot=50)
    assert ns[ns.name == 'error'].iloc[0].no_step > ns[ns.name == 'error'].iloc[0].down_floor


def test_same_fly_step_has_no_interval(planted):
    st, ss = E.recorded_steps(n_boot=20)
    s = ss[ss.step.str.contains('same_fly')]
    assert len(s) and s.n_interval_excludes_0.isna().all()


def test_megamat_vs_control_runs(planted):
    m = E.megamat_vs_control(n_boot=50)
    assert set(m.reference) == {'control', 'megamat'}
    assert set(m.criterion) == set(E.MEGAMAT_CRITERIA)
    e = m[(m.criterion == 'lock') & (m.reference == 'control') & (m.name == 'error')].iloc[0]
    assert np.isfinite(e.mean_views) and e.n_odors == 12
    assert e.megamat_criterion == 1.0                           # the planted population's criterion
    o = m[(m.criterion == 'apl16') & (m.reference == 'control') & (m.name == 'error')].iloc[0]
    assert o.megamat_criterion == D.CRITERIA['kc_soma']
    assert o.megamat_rate_17 < e.megamat_rate_17                # a higher criterion: fewer responses


def test_soma_criteria_lines():
    assert E.soma_criteria_lines({}) == []
    row = dict(cohort='validation2, the 13 tsub03 recordings', n_flies=13, criterion=5.9, units=D.UNITS_BASELINE,
               factor=1.03, rate_mean=0.05, rate_min=0.02, rate_max=0.08, apl16_criterion=1.9684,
               rate_mean_at_apl16_criterion=0.15)
    g = dict(noise_bridge_median=0.88, noise_bridge_min=0.77, noise_bridge_max=1.09)
    f = dict(noise_bridge_reference=dict(median=0.91), table=[row, dict(row, cohort='kmix (session z)', factor=None,
             units=D.UNITS_SESSION_Z)], validation2=dict(groups=dict(tsub03=g, other=g)),
             same_fly=dict(noise_bridge_per_fly={'x': 0.84}))
    L = E.soma_criteria_lines(dict(soma_criteria=f))
    assert any('tsub03' in x and '1.030' in x for x in L)
    assert any(x.startswith('| kmix (session z)') and ' z |' in x and '| — |' in x for x in L)


def test_energy_rows(planted, monkeypatch):
    keys = [('kc_soma', 'apl16', 'control')]
    monkeypatch.setattr(D, 'POPULATIONS', tuple(keys))
    en = E.flies_energy(max_units=50)
    assert set(en.reading) == {'raw', 'ranked'}
    assert len(en) == 12 and (en.ratio > 0).all()
