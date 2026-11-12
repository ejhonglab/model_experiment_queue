"""Tests of the repository's layout (kcq/layout.py: the flat module names from before the move into folders on
2026-11-12, and `package_imports`, which rewrites an old script's imports), and of E78_claw_knee.py's baseline E10
(`baseline_e10_source` at commits from the three layouts, and `_run_e10_redirected` loading a baseline as
experiments.E10_joint_fit). The E78 tests need the shared directory (E10 reads the data-target registry at import)."""
import ast
import subprocess
import sys
from pathlib import Path

import pytest

from kcq.layout import PACKAGES, ROOT, new_module, package_imports

REPO = Path(__file__).resolve().parents[1]  # the repository root
# commits of E10_joint_fit.py in each layout: this repository's root (before the move into folders) and al_analysis's
# scripts/claude_pair_remapping_experiments (before 2026-11-11)
QUEUE_ROOT_COMMIT = '6d7b597'
AL_ANALYSIS_COMMIT = '6b0c928'


def test_root_and_packages():
    assert ROOT == REPO
    for pkg in PACKAGES:
        assert (ROOT / pkg / '__init__.py').is_file(), pkg


def test_new_module():
    assert new_module('kc_remapping_harness') == 'kcq.harness'
    assert new_module('apl_effect_scores') == 'kcq.apl_effect'
    assert new_module('noise_floor_trial_distances') == 'kcq.noise_floor'
    assert new_module('E10_joint_fit') == 'experiments.E10_joint_fit'
    assert new_module('build_summary') == 'reports.build_summary'
    assert new_module('numpy') is None and new_module('json') is None
    # every module of the packages maps back to itself
    for pkg in PACKAGES:
        for p in (ROOT / pkg).glob('*.py'):
            if p.stem != '__init__' and p.stem not in ('harness', 'apl_effect', 'noise_floor'):
                assert new_module(p.stem) == f'{pkg}.{p.stem}'


def test_package_imports():
    src = '\n'.join([
        'import os, sys',
        'import numpy as np',
        'import kc_remapping_harness as H',
        'import apl_effect_scores',
        'import noise_proxy  # comment',
        'from E3_drive_dependent_gain import (BOUTON,',
        '    rate)',
        'def f():',
        '    import E72_wiring_population as E72',
        '    from apl_effect_scores import frac_data',
        '    import al_analysis.mb_model as mb',
        '    return H, E72, frac_data, mb, apl_effect_scores, noise_proxy',
        ''])
    out = package_imports(src)
    assert out.splitlines() == [
        'import os, sys',
        'import numpy as np',
        'from kcq import harness as H',
        'from kcq import apl_effect as apl_effect_scores',
        'from kcq import noise_proxy  # comment',
        'from experiments.E3_drive_dependent_gain import (BOUTON,',
        '    rate)',
        'def f():',
        '    from experiments import E72_wiring_population as E72',
        '    from kcq.apl_effect import frac_data',
        '    import al_analysis.mb_model as mb',
        '    return H, E72, frac_data, mb, apl_effect_scores, noise_proxy']
    ast.parse(out)
    # a file in the package layout comes back unchanged
    cur = (ROOT / 'experiments' / 'E78_claw_knee.py').read_text()
    assert package_imports(cur) == cur


def _flat_imports(src):
    bad = []
    for n in ast.walk(ast.parse(src)):
        if isinstance(n, ast.Import):
            bad += [a.name for a in n.names if new_module(a.name.split('.')[0])]
        elif isinstance(n, ast.ImportFrom) and n.module and new_module(n.module.split('.')[0]):
            bad.append(n.module)
    return bad


@pytest.fixture(scope='module')
def e78():
    H = pytest.importorskip('kcq.harness')
    try:
        H.data_targets_path().read_bytes()
    except (FileNotFoundError, OSError):
        pytest.skip('the shared directory (data-target registry) is not available')
    return pytest.importorskip('experiments.E78_claw_knee')


def test_baseline_e10_source(e78):
    """E10 at HEAD (experiments/), at a commit from before the move into folders (the root) and at an al_analysis
    commit (before 2026-11-11): each comes back with absolute imports only."""
    head = e78.baseline_e10_source('HEAD')
    committed = subprocess.run(['git', 'show', 'HEAD:experiments/E10_joint_fit.py'], cwd=REPO, capture_output=True,
        text=True).stdout
    assert head == package_imports(committed) and not _flat_imports(head)
    for commit in (QUEUE_ROOT_COMMIT, AL_ANALYSIS_COMMIT):
        try:
            src = e78.baseline_e10_source(commit)
        except ValueError:
            pytest.skip(f'{commit} is not in this clone')
        assert 'from kcq import harness as H' in src and not _flat_imports(src), commit
    with pytest.raises(ValueError):
        e78.baseline_e10_source('0' * 40)


def test_run_e10_redirected_baseline(e78, monkeypatch, tmp_path):
    """With E78_E10_DIR, `_run_e10_redirected` runs that folder's E10_joint_fit.py as experiments.E10_joint_fit, with
    the scratch experiment's outputs under the given root."""
    import experiments
    from kcq import harness as H
    E10 = pytest.importorskip('experiments.E10_joint_fit')
    # the function replaces these for the rest of its process (a separate one when E78 runs it): restored afterwards
    monkeypatch.setitem(sys.modules, 'experiments.E10_joint_fit', E10)
    monkeypatch.setattr(experiments, 'E10_joint_fit', E10)
    monkeypatch.setattr(H, 'results_dir', H.results_dir)
    monkeypatch.setattr(H, 'outputs_dir', H.outputs_dir)
    monkeypatch.setattr(sys, 'argv', list(sys.argv))
    old = tmp_path / 'old'
    old.mkdir()
    (old / 'E10_joint_fit.py').write_text(package_imports('\n'.join([
        'import sys',
        'import kc_remapping_harness as H',
        'CALLS = []',
        'def main():',
        '    import E10_joint_fit as me',
        f'    CALLS.append((me.__file__, sys.argv[1:], H.outputs_dir({e78.CHECK_SCRATCH!r})))',
        ''])))
    monkeypatch.setenv('E78_E10_DIR', str(old))
    e78._run_e10_redirected(tmp_path / 'root', ['--x', '1'])
    mod = sys.modules['experiments.E10_joint_fit']
    assert Path(mod.__file__) == old / 'E10_joint_fit.py'
    assert mod.CALLS == [(mod.__file__, ['--x', '1'], tmp_path / 'root' / 'outputs')]
