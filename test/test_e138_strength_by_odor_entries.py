"""Tests of experiments/E138_strength_by_odor_entries.py (E138, the odors that carry the strength error and the glomeruli
that feed them): the module imports, its constants, and the held design item 6 is refused without Tom's or Yang's word."""
import pytest

from experiments import E138_strength_by_odor_entries as E


def test_imports_and_constants():
    assert E.EXPERIMENT == 'E138_strength_by_odor_entries'
    assert E.WIRINGS == E.S_WIRINGS + E.T_WIRINGS and len(E.WIRINGS) == 8
    assert E.STRENGTH_PARTS == ('strength/rate', 'strength/mean')
    assert E.CLASSES == ('ir', 'not_dendrite28', 'few_hemibrain')


def test_regression_refused_without_word():
    with pytest.raises(SystemExit):
        E.main(['--part', 'regression'])
