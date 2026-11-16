#!/usr/bin/env python3
"""Tests for the rebuild archives of fitted models (archive_tools/rebuild_archive.py, pyolfsysm.rebuild).
Without simulations: E10's `breadth_weights` and E27's `gained` are pyolfsysm.rebuild's; the presets; every written
archive of the presets loads with its files' sha256 intact, holds its fit's absolute parameters (fits.csv) and the
fit's saved outputs (the experiment's archived files) exactly, and its recorded checks found 0 differing spike counts
and responses. Slow (skipped with QUICK=1): C2's archive rebuilt from its files alone reproduces the saved fit
(`check_archive`, both routes; a few minutes). Needs the shared directory; skipped without it."""

import json
import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
RA = pytest.importorskip('archive_tools.rebuild_archive')
from pyolfsysm import rebuild as R

QUICK = os.environ.get('QUICK') == '1'
VARIANTS = [v for vs in RA.presets().values() for v in vs]


def test_moved_functions_are_pyolfsysms():
    from experiments import E10_joint_fit as J
    from experiments import E27_glomerulus_gain as E27
    assert J.breadth_weights is R.breadth_weights
    assert E27.gained is R.glomerulus_gain


def test_presets():
    p = RA.presets()
    assert p['C2'] == [RA.C2] and p['E44'] == [RA.E44]
    assert len(p['C2-wirings']) == 7 and all(f == 'E72_wiring_population' for f, _ in p['C2-wirings'])
    assert RA.resolve(['C2', 'E44_operating_point:x']) == [RA.C2, ('E44_operating_point', 'x')]


@pytest.mark.parametrize('folder,tag', VARIANTS, ids=[t for _, t in VARIANTS])
def test_written_archive(folder, tag):
    d = RA.archive_dir(RA.ARCHIVE, folder, tag)
    if not (d / R.ARCHIVE_JSON).exists():
        pytest.skip(f'no archive for {folder}:{tag}')
    a = R.load_archive(d)
    row = RA.fits_row(folder, tag)
    expect = {c: float(row[c]) for c in RA.PARAM_COLS}
    expect['wKCKC_scale'] = float(row['kc_kc_inh'])
    assert a['params'] == expect
    assert a['model_kws']['n_spikes_for_response'] == 2 and a['meta']['response_criterion']['n_spikes_for_response'] == 2
    arch = RA.ARCHIVE / folder
    for c, prefix in R.CONDITIONS.items():
        for i, kind in enumerate(('responses', 'spike_counts')):
            saved = pd.read_parquet(arch / f'{tag}_{prefix}_{kind}.parquet')
            assert a['outputs'][c][i].equals(saved)
    assert a['orn_deltas'].equals(H.apl16_orn_deltas())
    checks = a['meta']['checks']
    assert checks and all(ch['ok'] for ch in checks)
    for route in ('fit_mb_model', 'simulator'):
        for c in ('on', 'off'):
            assert checks[-1][route][c]['spike_count_diffs'] == 0 and checks[-1][route][c]['response_diffs'] == 0


@pytest.mark.skipif(QUICK, reason="C2's rebuild from its archive (simulations, minutes)")
def test_c2_rebuild_from_its_archive():
    d = RA.archive_dir(RA.ARCHIVE, *RA.C2)
    if not (d / R.ARCHIVE_JSON).exists():
        pytest.skip('no archive for C2')
    res = R.check_archive(d)
    assert res['ok'], json.dumps(res)
