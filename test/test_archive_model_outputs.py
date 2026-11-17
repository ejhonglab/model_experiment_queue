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
