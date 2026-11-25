"""Tests of E137 (the held-out panels read like for like): the interface and the readings' rules on made-up tables."""
import pytest

from experiments import E137_heldout_readings as E


def test_interface():
    assert E.READINGS[0] == 'as_uploaded'
    assert len(E.READINGS) == 5
    assert E.POSITION_READINGS == E.READINGS[1:]
    assert E.V2_RULES == ('common_mean', 'fly_scale')
    for name in ('load_lock', 'fly_mixture_view', 'part_lock', 'part_position', 'part_readings', 'part_rescore',
                 'part_rise', 'part_validation2', 'part_order', 'part_classes'):
        assert callable(getattr(E, name))
