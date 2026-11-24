"""Tests of E134's stubs by atlas (kcq/stage_readings.py and experiments/E134_population_baseline.py): the modules
import, their main names exist, their names agree with viki's data module and hal's measures, and the trivial helpers
work. No shared-directory reads, no model run."""
import json

import numpy as np
import pandas as pd
import pytest

from data_prep import stage_population_data as D
from experiments import E134_population_baseline as B
from kcq import stage_readings as R

STAGE_READINGS_NAMES = (
    'STAGES', 'PANELS', 'CONDITIONS', 'Reading', 'DEFAULT_READINGS', 'ALTERNATIVE_READINGS', 'SavedRun', 'StageCurve',
    'NoiseModel', 'ModelReading', 'cohort_of', 'readings_of', 'reading_by_name', 'reading_definitions',
    'pooled_values', 'load_e90_run', 'load_e107_run', 'load_step3_run', 'run_panels', 'stage_changes',
    'driver_glomeruli', 'sample_driver_glomeruli', 'pool_boutons', 'pool_claws', 'sample_units', 'pooled_fly_values',
    'fit_stage_curve', 'apply_curve', 'fly_noise', 'add_measurement_noise', 'fly_side', 'read_stage', 'read_run',
    'read_family', 'SomaReading', 'soma_readings', 'reading_provenance',
)
RUNNER_NAMES = (
    'STEP2_SOURCES', 'BESIDE_MODELS', 'require_lock', 'step2_runs', 'score_family', 'gate_wirings', 'part_step2',
    'tag_of', 'parse_tag', 'planned_step3', 'claim', 'release', 'panel_input', 'held_targets', 'build_point',
    'run_stages', 'check_identity', 'save_run', 'load_run', 'run_point', 'part_step3', 'runs_table', 'stage_table',
    'step_table', 'classify_stages', 'classify_steps', 'classify_panels', 'realism_rows', 'part_score', 'provenance',
    'print_commands', 'main',
)


def test_names_exist():
    for n in STAGE_READINGS_NAMES:
        assert hasattr(R, n), n
    for n in RUNNER_NAMES:
        assert hasattr(B, n), n


def test_names_agree_with_the_data_module():
    assert R.STAGES == D.STAGES
    assert R.PANELS == D.PANELS
    assert R.NAMED_UNITS == D.NAMED_UNITS
    cohorts = {c for _, _, c in D.POPULATIONS}
    for panel in R.PANELS:
        for cond in R.CONDITIONS:
            c = R.cohort_of(panel, cond)
            assert c is None or c in cohorts
    assert R.cohort_of('apl16', 'off') == 'tnt'
    assert R.cohort_of('megamat', 'off') is None
    with pytest.raises(ValueError):
        R.cohort_of('megamat', 'acute_off')


def test_readings():
    for s in R.STAGES:
        rs = R.readings_of(s)
        assert rs[0] is R.DEFAULT_READINGS[s]
        assert len({r.name for r in rs}) == len(rs)
        for r in rs:
            assert r.stage == s
            assert r.curve in R.CURVES and r.noise in R.NOISE_MODELS and r.criterion in R.CRITERIA
            assert r.curve_set_on in R.CURVE_SET_ON and r.curve_scope in R.CURVE_SCOPES
            assert r.fly_weighting in R.FLY_WEIGHTING
            assert R.reading_by_name(s, r.name) is r
    assert R.DEFAULT_READINGS['kc_claw'].pool in R.CLAW_POOLS
    assert R.DEFAULT_READINGS['pn_bouton'].pool in R.BOUTON_POOLS
    with pytest.raises(KeyError):
        R.reading_by_name('kc_soma', 'no such reading')
    d = R.reading_definitions()
    assert set(d['stages']) == set(R.STAGES)
    assert json.loads(json.dumps(d)) == d


def test_trivial_helpers():
    a = pd.DataFrame({'2h': [1.0, np.nan], 'eb': [3.0, 4.0]}, index=['u1', 'u2'])
    b = pd.DataFrame({'2h': [5.0], 'eb': [6.0]})
    assert sorted(R.pooled_values([a, b])) == [1.0, 3.0, 4.0, 5.0, 6.0]
    assert sorted(R.pooled_values([a, b], odors=['eb'])) == [3.0, 4.0, 6.0]
    g = pd.DataFrame({'2h': [1.0, 2.0, 3.0]}, index=['DM1', 'DL5', 'VA2'])
    assert list(R.sample_driver_glomeruli(g, ['VA2', 'DM1']).index) == ['VA2', 'DM1']
    with pytest.raises(KeyError):
        R.sample_driver_glomeruli(g, ['DA1'])
    run = R.SavedRun('E107:B', 'megamat_fitted_real', 'C2', 'real', 'megamat',
                     conditions={'on': {'spike_counts': b, 'responses': b}, 'off': {'spike_counts': b, 'responses': b}})
    assert run.stages() == ('kc_soma',)
    assert set(R.soma_readings()) == {r.name for r in R.readings_of('kc_soma')}
    assert isinstance(R.ModelReading('C2', 'real').populations, dict)


def test_points_and_tags():
    pts = B.planned_step3()
    assert len(pts) == 47
    assert len(set(pts)) == len(pts)
    assert sum(pt.model in B.WIRINGS for pt in pts) == 32
    assert pts[0] == B.Point('megamat', 'fitted', 'real')
    assert B.planned_step3(reverse=True) == pts[::-1]
    for pt in pts:
        assert B.parse_tag(B.tag_of(pt)) == pt
    assert B.tag_of(B.Point('validation2', 'fitted', 'published_hemibrain')) == 'validation2_fitted_published_hemibrain'
    assert len(B.planned_step3(models=['real'], readings=['held'])) == 2
    assert all(pt.reading == 'fitted' for pt in B.planned_step3(models=list(B.BESIDE_MODELS)))
    with pytest.raises(ValueError):
        B.tag_of(B.Point('megamat', 'searched', 'real'))


def test_claims(tmp_path):
    assert B.claim('megamat_fitted_real', 'test', root=tmp_path)
    assert not B.claim('megamat_fitted_real', 'test2', root=tmp_path)
    assert B.claim('megamat_fitted_real', 'test2', reclaim=True, root=tmp_path)
    c = json.loads((tmp_path / 'claims' / 'megamat_fitted_real' / 'claim.json').read_text())
    assert c['lane'] == 'test2'
    B.release('megamat_fitted_real', root=tmp_path)
    assert not (tmp_path / 'claims' / 'megamat_fitted_real').exists()


def test_guard_refuses_before_the_lock():
    if B.LOCK_COPY.exists():
        pytest.skip('the lock exists; the guard is then viki\'s check')
    with pytest.raises(PermissionError):
        B.require_lock()
    for part in ('step2', 'step3', 'score'):
        with pytest.raises(PermissionError):
            B.main(['--part', part])


def test_commands(capsys):
    B.main(['--part', 'commands'])
    out = capsys.readouterr().out.splitlines()
    assert len(out) == sum(len(lanes) for _, lanes in B.LANES.values())
    assert all('--part step3 --lane ' in line for line in out)
