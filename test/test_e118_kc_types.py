#!/usr/bin/env python3
"""Tests for E118 (the KC types in the best family and in the recordings; experiments/E118_kc_types.py)."""

import pytest

X = pytest.importorskip('experiments.E118_kc_types')


def test_interface():
    for name in ('kc_weights', 'permuted_labels', 'weighted_pair_distances', 'weighted_rate', 'score_distances',
                 'reading', 'part0', 'part_a', 'rule_b', 'outcome_a'):
        assert callable(getattr(X, name))
