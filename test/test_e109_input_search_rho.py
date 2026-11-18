#!/usr/bin/env python3
"""Tests for E109 (the input of paa, B-myr and 2-but searched for the two ρ's on held-out wirings and recordings;
experiments/E109_input_search_rho.py). Stub: the interface imports."""

import pytest

X = pytest.importorskip('experiments.E109_input_search_rho')


def test_interface():
    for name in ('e108_class', 'require_e108_class', 'start_input', 'two_rhos_population', 'pattern_penalty', 'search',
                 'part0', 'gate', 'class_of', 'recurring_entries', 'step7'):
        assert callable(getattr(X, name))
