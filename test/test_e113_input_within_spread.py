#!/usr/bin/env python3
"""Tests for E113 (the input of paa and B-myr without a search; experiments/E113_input_within_spread.py): the interface."""

import pytest

X = pytest.importorskip('experiments.E113_input_within_spread')


def test_interface():
    for name in ('out_dir', 'local_dir', 'entry_sizes', 'control_glomeruli', 'factors_of', 'fly_deviations', 'fly_input',
                 'input_of', 'write_inputs', 'identity', 'plan', 'part_sims', 'readings', 'gate_arm', 'gate_flies',
                 'outcomes', 'refit_choice', 'part_refit', 'part_remy', 'score', 'write_variants', 'main'):
        assert callable(getattr(X, name))
    assert X.ARMS == ('D2', 'D1', 'D1x') and X.CLASS_ARM == 'D2'
    assert len(X.CONTROLS) == 5
