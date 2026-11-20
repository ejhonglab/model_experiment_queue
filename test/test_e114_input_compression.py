#!/usr/bin/env python3
"""Tests for E114 (a compressive map from the ORN recordings to the input rates; experiments/E114_input_compression.py)."""

import pytest

X = pytest.importorskip('experiments.E114_input_compression')


def test_interface():
    for f in ('power_map', 'saturating_map', 'input_of', 'write_inputs', 'identity', 'fit_path_check', 'part_sims',
              'readings', 'gate', 'class_of', 'choose_c', 'outcome_of', 'part_c', 'part_d', 'part_beside', 'score',
              'write_variants', 'main'):
        assert callable(getattr(X, f))
