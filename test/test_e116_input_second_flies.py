#!/usr/bin/env python3
"""Tests for E116 (the input of paa and B-myr searched against a second set of flies and taken from Hallem's receptor
data; experiments/E116_input_second_flies.py)."""

import pytest

X = pytest.importorskip('experiments.E116_input_second_flies')


def test_interface():
    for name in ('split_flies', 'remy_halves', 'remy_targets', 'RemyObjective', 'search', 'hallem_inputs', 'part0',
                 'part_h', 'part_read', 'part_remy', 'gate_r', 'gate_h', 'entries', 'outcomes', 'score'):
        assert hasattr(X, name)
