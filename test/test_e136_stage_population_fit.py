"""Tests of experiments/E136_stage_population_fit.py (E136, the best family fitted to all megamat stages at once). The
stub's tests: the arms, the wirings and the open definitions."""
import pytest

from experiments import E136_stage_population_fit as F


def test_arms_and_wirings():
    F.check_wirings()
    assert set(F.WIRINGS) == set(F.S_WIRINGS) | set(F.T_WIRINGS) and not set(F.S_WIRINGS) & set(F.T_WIRINGS)
    assert tuple(F.ARMS) == ('A0', 'A1', 'A2', 'A3', 'A4')
    for name, arm in F.ARMS.items():
        assert arm.name == name
        assert 0 < len(arm.free) <= F.MAX_FREE and len(set(arm.free)) == len(arm.free)
    assert 'kckc' not in F.ARMS['A4'].free and F.ARMS['A4'].kws == {'kckc_fixed': 0.0}
    assert F.CONTROL not in F.ARMS


def test_open_definitions_named():
    assert {'halves', 'operating_point', 'free_parameters', 'claw_reading'} <= set(F.OPEN_DEFINITIONS)


def test_stub_refuses():
    with pytest.raises(NotImplementedError):
        F.main(['--part', 'fit', '--arm', 'A0'])
