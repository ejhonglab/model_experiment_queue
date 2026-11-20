"""Tests of E120 (experiments/E120_bouton_stage_measured.py); synthetic data only."""
from experiments import E120_bouton_stage_measured as E


def test_imports():
    assert E.EXPERIMENT == 'E120_bouton_stage_measured'
