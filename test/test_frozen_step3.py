"""The queue's generally useful functions that moved to hong2p and pyolfsysm (2026-11-12) give the queue the same
outputs as before the move: every case of frozen_step3/cases.py, run through the queue's names (kcq.harness and the
experiments' re-exports, with the queue's defaults), against the outputs frozen from the code before the move
(frozen_step3/make_frozen.py). To the last bits (rtol 1e-12) with the library versions of the frozen files (the queue's
venv has them), looser with others (frozen_outputs.check)."""
from pathlib import Path

import pytest

from frozen_step3 import cases as C
from frozen_step3 import frozen_outputs as F

HERE = Path(__file__).resolve().parent / 'frozen_step3'
GROUPS = {'hong2p': C.hong2p_cases, 'pyolfsysm': C.pyolfsysm_cases, 'queue': C.queue_cases}
_FROZEN = {}


def _frozen(group):
    if group not in _FROZEN:
        fz = F.load(HERE / f'frozen_{group}.npz')
        _FROZEN[group] = (fz, F.objects(fz['in']), GROUPS[group]())
    return _FROZEN[group]


@pytest.mark.parametrize('group,case', [(g, c) for g, f in GROUPS.items() for c in f()])
def test_frozen_case(group, case):
    fz, O, cases = _frozen(group)
    fn, solver = cases[case]
    if group == 'queue':
        # the queue's own data (the shared directory) is read here
        pytest.importorskip('kcq.harness').shared_dir()
    F.check(fz, case, fn(O), solver=solver)


def test_versions_match():
    """The queue's venv has the library versions the files were frozen with, so every comparison above is to the last bits."""
    for g in GROUPS:
        assert F.same_versions(_frozen(g)[0]), (g, _frozen(g)[0]['version'], F.versions())
