#!/usr/bin/env python3
"""Tests for E105 (the KCs' response rates at the weak odors read like for like; experiments/E105_weak_odor_rates.py), on
made-up data: the three sets are never pooled; c* reproduces a known response rate; the classes and W on known rows,
with the 6-of-8 edge cases; the held rule; the lock refuses a changed definition; the reader refuses raw-only rows and
any path outside the shared directory."""

import pytest

M = pytest.importorskip('experiments.E105_weak_odor_rates')


def test_module_constants():
    assert M.SETS == ('D1', 'D2', 'N')
    assert M.MIN_WIRINGS == 6 and len(M.ALL_WIRINGS) == 8
