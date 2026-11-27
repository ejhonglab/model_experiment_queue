#!/usr/bin/env python3
"""Tests for archive_tools/archive_model_outputs.py's folder prefixes and excludes (E102: the step folders repeat their
files' names): the exclusion of paths and glob patterns, the archived name of a file under a prefixed folder, and an
archive of a made-up outputs folder (renamed copies, the earlier bare-named copy removed, the default excludes left
out, a variant list with every step's variants apart, params.csv from variants.csv). Needs the shared directory (the
archiver's module reads it at import); skipped without it."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.shared_dir()
    A = pytest.importorskip('archive_tools.archive_model_outputs')
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory is not available', allow_module_level=True)

E102 = 'E102_antennal_lobe_normalization'


def test_excluded():
    ex = ('D/start_check', 'D/runs/B/*/*_held_*', 'E/B/*/refit_D_*', 'tmp/')
    assert A.excluded('D/start_check/inhadd1000_real/outputs/x_on_spike_counts.parquet', ex)
    assert A.excluded('D/start_check', ex)
    assert not A.excluded('D/start_checks.json', ex)
    assert A.excluded('D/runs/B/inhadd1000/real_held_on_spike_counts.parquet', ex)
    assert not A.excluded('D/runs/B/inhadd1000/real_on_spike_counts.parquet', ex)
    assert A.excluded('E/B/inhadd1000/refit_D_on_apl.csv', ex) and not A.excluded('E/B/inhadd1000/refit_F_on_apl.csv', ex)
    assert A.excluded('tmp/a.txt', ex) and not A.excluded('tmpx/a.txt', ex)


@pytest.mark.parametrize('rel, out', [
    ('B/inhadd1000/real_held_on_spike_counts.parquet', 'B/inhadd1000/B_inhadd1000_real_held_on_spike_counts.parquet'),
    ('B/gain_fixed/M1_point.json', 'B/gain_fixed/B_gain_fixed_M1_point.json'),
    ('B/identity/b1_stimuli.csv', 'B/identity/B_identity_b1_stimuli.csv'),
    ('D/runs/B/inhadd1000/real_on_apl.csv', 'D/runs/B/inhadd1000/D_inhadd1000_real_on_apl.csv'),
    ('E/B/inhadd1000/refit_F_acute_off_spike_counts.parquet', 'E/B/inhadd1000/E_inhadd1000_refit_F_acute_off_spike_counts.parquet'),
    ('H/B/inhadd1000/H1_real_held_on_spike_counts.parquet', 'H/B/inhadd1000/H_inhadd1000_H1_real_held_on_spike_counts.parquet'),
    ('H/identity_inhadd1000.json', None),
    ('rate0.068_aplboth_pnadd1000sc544.444_k2_cap3_on_spike_counts.parquet', None),
    ('logs/lane_a.log', None), ('D/init_inhadd1000.json', None), ('B/x.csv', None)])
def test_archive_rel(rel, out):
    assert A.archive_rel(E102, Path(rel)) == Path(out or rel)
    assert A.archive_rel('E94_bouton_share', Path(rel)) == Path(rel)              # no prefixes for other experiments


def test_archive_experiment_with_folder_prefixes(tmp_path, monkeypatch):
    """A made-up E102 outputs folder archived: the step folders' files renamed with their step, a copy archived earlier
    under the bare name removed, D/start_check and the excluded patterns left out, the README's variants apart per
    step, params.csv from the results' variants.csv."""
    here, shared = tmp_path / 'repo', tmp_path / 'shared'
    src, res = here / 'outputs' / E102, shared / 'results' / E102
    arch = shared / 'data' / 'model_outputs'
    monkeypatch.setattr(H, 'HERE', here)
    monkeypatch.setattr(H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(H, 'results_dir', lambda e: shared / 'results' / e)
    monkeypatch.setattr(A, 'ARCHIVE', arch)
    s = pd.DataFrame(np.ones((3, 2)), columns=['a', 'b'])
    files = ['B/identity/real_held_on', 'B/identity/real_on', 'B/inhadd1000/real_held_on', 'B/inhadd1000/real_on',
             'D/runs/B/inhadd1000/real_on', 'D/runs/B/inhadd1000/real_held_on', 'D/start_check/inhadd1000_real/outputs/t_on',
             'E/B/inhadd1000/refit_F_held_on', 'E/B/inhadd1000/refit_D_on', 'tagx_on']
    for f in files:
        for kind in ('spike_counts', 'responses'):
            p = src / f'{f}_{kind}.parquet'
            p.parent.mkdir(parents=True, exist_ok=True)
            s.to_parquet(p)
    res.mkdir(parents=True)
    pd.DataFrame(dict(tag=['B_identity_real', 'B_inhadd1000_real', 'D_inhadd1000_real', 'E_inhadd1000_refit_F'],
                      fixed_thr=[1.0, 2.0, 3.0, 4.0], pn_inhadd=[31.4088, 1000.0, 1000.0, 1000.0], stray=1)).to_csv(
        res / 'variants.csv', index=False)
    old = arch / E102 / 'B' / 'inhadd1000' / 'real_held_on_spike_counts.parquet'
    old.parent.mkdir(parents=True)
    s.to_parquet(old)
    r = A.archive_experiment(E102, 'host')
    dst = arch / E102
    got = sorted(p.relative_to(dst).as_posix() for p in dst.rglob('*_spike_counts.parquet'))
    assert got == ['B/identity/B_identity_real_held_on_spike_counts.parquet', 'B/identity/B_identity_real_on_spike_counts.parquet',
                   'B/inhadd1000/B_inhadd1000_real_held_on_spike_counts.parquet', 'B/inhadd1000/B_inhadd1000_real_on_spike_counts.parquet',
                   'D/runs/B/inhadd1000/D_inhadd1000_real_on_spike_counts.parquet',
                   'E/B/inhadd1000/E_inhadd1000_refit_F_held_on_spike_counts.parquet', 'tagx_on_spike_counts.parquet']
    assert not old.exists()
    assert r['n_variants'] == 7
    readme = (dst / 'README.md').read_text()
    for v in ('B_identity_real', 'B_identity_real_held', 'B_inhadd1000_real', 'D_inhadd1000_real', 'E_inhadd1000_refit_F_held'):
        assert f'- `{v}`' in readme
    assert 'B_{1}_' in readme
    p = pd.read_csv(dst / 'params.csv')
    assert list(p['tag']) == ['B_identity_real', 'B_inhadd1000_real', 'D_inhadd1000_real', 'E_inhadd1000_refit_F']
    assert 'stray' not in p.columns and list(p['pn_inhadd']) == [31.4088, 1000.0, 1000.0, 1000.0]


def test_withheld_folders_are_never_copied():
    """A folder of withheld outputs (withheld_*, e.g. E97's of 2026-11-15) is skipped at any depth; a file whose own
    name starts so is not a folder and is kept; E97's claims and caches are left out by default."""
    A = pytest.importorskip('archive_tools.archive_model_outputs')
    assert A.excluded('withheld_2026-11-15/runs/B2_mean_real_on_responses.parquet', ())
    assert A.excluded('B/withheld_x/y.csv', ())
    assert not A.excluded('withheld_note.txt', ())
    ex = A.DEFAULT_EXCLUDE['E97_glomerulus_claw_shares']
    assert A.excluded('claims/B2_mean2_real/owner', ex) and A.excluded('distances_cache.pkl', ex)
    assert not A.excluded('runs_mcns/B2_mean2_real_on_responses.parquet', ex)


def test_params_keep_input_change_columns(tmp_path, monkeypatch):
    """E116's variants.csv names each variant's input by change, search, size and panel (which search's factors, at
    what fraction of their log size, on which odor panel); params.csv keeps those columns, not only the fit's."""
    name = 'E116_input_second_flies'
    here, shared = tmp_path / 'repo', tmp_path / 'shared'
    src, res = here / 'outputs' / name, shared / 'results' / name
    arch = shared / 'data' / 'model_outputs'
    monkeypatch.setattr(H, 'HERE', here)
    monkeypatch.setattr(H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(H, 'results_dir', lambda e: shared / 'results' / e)
    monkeypatch.setattr(A, 'ARCHIVE', arch)
    s = pd.DataFrame(np.ones((3, 2)), columns=['a', 'b'])
    for f in ('R/Ra_s0.6_P3_on', 'V/Ra_s1_P3_v2_on'):
        p = src / f'{f}_spike_counts.parquet'
        p.parent.mkdir(parents=True, exist_ok=True)
        s.to_parquet(p)
    res.mkdir(parents=True)
    pd.DataFrame(dict(tag=['Ra_s0.6_P3', 'Ra_s1_P3_v2'], change=['Ra_s0.6', 'Ra_s1'], search=['Ra', 'Ra'], size=[0.6, 1.0],
                      panel=['APL16', 'validation2'], fixed_thr=[1.0, 2.0], kind_text='x')).to_csv(res / 'variants.csv', index=False)
    A.archive_experiment(name, 'host')
    p = pd.read_csv(arch / name / 'params.csv')
    assert list(p['change']) == ['Ra_s0.6', 'Ra_s1'] and list(p['search']) == ['Ra', 'Ra']
    assert list(p['size']) == [0.6, 1.0] and list(p['panel']) == ['APL16', 'validation2']
    assert 'kind_text' not in p.columns


def test_params_tables_override(tmp_path, monkeypatch):
    """E120's params.csv is built from variants.csv only: its a1_ / a3_ / b_stage_summary tables are readings of other
    runs (rows with a form or a reading's variant but no saved model), which the migrate would otherwise keep or set
    aside as unidentified."""
    name = 'E120_bouton_stage_measured'
    here, shared = tmp_path / 'repo', tmp_path / 'shared'
    src, res = here / 'outputs' / name, shared / 'results' / name
    arch = shared / 'data' / 'model_outputs'
    monkeypatch.setattr(H, 'HERE', here)
    monkeypatch.setattr(H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(H, 'results_dir', lambda e: shared / 'results' / e)
    monkeypatch.setattr(A, 'ARCHIVE', arch)
    p = src / 'B' / 'inhadd10000' / 'real_on_spike_counts.parquet'
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(np.ones((3, 2)), columns=['a', 'b']).to_parquet(p)
    res.mkdir(parents=True)
    pd.DataFrame(dict(form=['baseline', 'evoked'], n=[1.0, 2.0])).to_csv(res / 'a1_summary.csv', index=False)
    pd.DataFrame(dict(variant=['ref', 'B_sub_k2_d2'], form=['sub', 'sub'])).to_csv(res / 'a3_summary.csv', index=False)
    pd.DataFrame(dict(tag=['real'], wiring=['real'], pn_inhadd=[1e4], fixed_thr=[132.4])).to_csv(res / 'variants.csv', index=False)
    A.archive_experiment(name, 'host')
    q = pd.read_csv(arch / name / 'params.csv')
    assert list(q['tag']) == ['real'] and list(q['source_table']) == ['variants.csv']
    assert list(q['pn_inhadd']) == [1e4] and 'form' not in q.columns


def test_in_place(tmp_path, monkeypatch):
    """E136's script writes its runs into the archive from every host: nothing of this host's outputs/ is copied, and
    params.csv (from results/<name>/measure_strength/variants.csv), the provenance files at their paths, the folder's
    README and the file count come from the archive."""
    name = 'E136_stage_population_fit'
    here, shared = tmp_path / 'repo', tmp_path / 'shared'
    src, res = here / 'outputs' / name, shared / 'results' / name
    arch = shared / 'data' / 'model_outputs'
    monkeypatch.setattr(H, 'HERE', here)
    monkeypatch.setattr(H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(H, 'results_dir', lambda e: shared / 'results' / e)
    monkeypatch.setattr(A, 'ARCHIVE', arch)
    local = src / 'measure_strength' / 'fit_A0' / 'A0_e5_real_on_spike_counts.parquet'
    local.parent.mkdir(parents=True)
    pd.DataFrame(np.ones((3, 2)), columns=['a', 'b']).to_parquet(local)
    for tag, sub in (('A1_megamat_P3', 'read'), ('A1_e122_real', 'fit_A1')):
        for cond in ('on', 'acute_off'):
            p = arch / name / 'measure_strength' / sub / f'{tag}_{cond}_spike_counts.parquet'
            p.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(np.ones((3, 2)), columns=['a', 'b']).to_parquet(p)
    (res / 'measure_strength').mkdir(parents=True)
    pd.DataFrame(dict(tag=['A1_megamat_P3', 'A1_e122_real'], subdir=['read', 'fit_A1'], pn_inhadd=[1017.3, 1017.3],
                      off_wKCKC_scale=[138.3, 40.1], reading=['read', 'fit'])).to_csv(
        res / 'measure_strength' / 'variants.csv', index=False)
    (res / 'measure_strength' / 'provenance_score.json').write_text('{}')
    r = A.archive_experiment(name, 'host')
    assert not (arch / name / 'measure_strength' / 'fit_A0').exists()
    assert (r['n_files'], r['n_variants']) == (4, 2) and r['conditions'] == 'acute_off, on'
    q = pd.read_csv(arch / name / 'params.csv')
    assert list(q['tag']) == ['A1_megamat_P3', 'A1_e122_real'] and list(q['off_wKCKC_scale']) == [138.3, 40.1]
    assert (arch / name / 'measure_strength' / 'provenance_score.json').exists()
    assert 'Written in place' in (arch / name / 'README.md').read_text()


def test_e142_in_place_from_step1_variants(tmp_path, monkeypatch):
    """E142's step-1 runs are written into the archive from every host: nothing of this host's outputs/ is copied, and
    params.csv comes from step1_variants.csv only (step 0's *_summary.csv tables are readings of the flies), keeping the
    table a run imposed, its route, the APL->bouton factor and the PN stage."""
    name = 'E142_dendrite_bouton_step'
    here, shared = tmp_path / 'repo', tmp_path / 'shared'
    src, res = here / 'outputs' / name, shared / 'results' / name
    arch = shared / 'data' / 'model_outputs'
    monkeypatch.setattr(H, 'HERE', here)
    monkeypatch.setattr(H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(H, 'results_dir', lambda e: shared / 'results' / e)
    monkeypatch.setattr(A, 'ARCHIVE', arch)
    local = src / 'step1' / 'megamat_imposed_a1_real_on_spike_counts.parquet'
    local.parent.mkdir(parents=True)
    pd.DataFrame(np.ones((3, 2)), columns=['a', 'b']).to_parquet(local)
    (src / 'logs').mkdir()
    (src / 'logs' / 'step1_hal-s1.log').write_text('x')
    tags = ['megamat_a1stage_imposed_P3', 'apl16_a1stage_imposed_P3']
    for tag in tags:
        for cond in ('on', 'acute_off'):
            p = arch / name / 'step1' / f'{tag}_{cond}_spike_counts.parquet'
            p.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(np.ones((3, 2)), columns=['a', 'b']).to_parquet(p)
    res.mkdir(parents=True)
    pd.DataFrame(dict(odor=['pa', 'va'], rate=[0.4, -0.4])).to_csv(res / 'flies_steps_summary.csv', index=False)
    pd.DataFrame(dict(tag=tags, subdir=['step1'] * 2, table_arm=['imposed'] * 2, route=[1, 1], aplpn_factor=[0.155] * 2,
                      pn_stage=['A1'] * 2, pn_inhadd=[1017.3] * 2)).to_csv(res / 'step1_variants.csv', index=False)
    (res / 'provenance_score_step1.json').write_text('{}')
    r = A.archive_experiment(name, 'host')
    assert not (arch / name / 'logs').exists() and not (arch / name / 'step1' / local.name).exists()
    assert (r['n_files'], r['n_variants']) == (4, 2)
    q = pd.read_csv(arch / name / 'params.csv')
    assert list(q['tag']) == tags and set(q['source_table']) == {'step1_variants.csv'}
    assert list(q['pn_stage']) == ['A1', 'A1'] and list(q['aplpn_factor']) == [0.155, 0.155]
    assert list(q['table_arm']) == ['imposed'] * 2 and list(q['route']) == [1, 1] and 'odor' not in q.columns
    assert (arch / name / 'provenance_score_step1.json').exists()
