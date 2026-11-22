#!/usr/bin/env python3
"""Tests for E129 (the inhibition that outlasts the block reduced to a static level;
experiments/E129_block_survivor_static.py). No model simulation."""

import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E129_block_survivor_static')


def test_import():
    assert X.STEPS['step1a'] == ('G', 'T', 'P', 'S')
