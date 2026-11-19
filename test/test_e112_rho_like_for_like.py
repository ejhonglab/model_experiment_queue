"""E112 (the two ρ's on odor-odor distances read as the recordings are made): the interface imports."""
from experiments import E112_rho_like_for_like as E


def test_interface():
    assert E.READINGS == ('R0', 'R1', 'R2', 'R3', 'R4')
    for f in ('reading_distances', 'pseudo_recording_distances', 'trial_noise_distances', 'reading_scores', 'bound_scores',
              'moves_verdict', 'reread_gate', 'e109_holds', 'part_a', 'part_b_prepare', 'part_b_score'):
        assert callable(getattr(E, f))
