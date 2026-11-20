"""E115 (the megamat gap against two sets of flies and the yardstick of the price in response rate ρ), on made-up
matrices: the interface."""
from experiments import E115_two_yardsticks as E


def test_interface():
    assert E.SHARED_MIN == 0.30 and E.LEVEL_TOL == 0.03 and E.POOLED_MAX_FALL == 0.02
    assert (E.CLASS_NEAR, E.CLASS_FAR) == (10, 15)
    for f in ('residual_correlation', 'two_set_rank', 'pair_classes', 'class_shares', 'reader_row', 'resample_part_a',
              'input_size', 'fractions', 'rate_rho', 'pooled_rate_rho', 'percentile_of', 'price_verdict', 'verdicts',
              'part0', 'part_a', 'part_b', 'registry_entries'):
        assert callable(getattr(E, f))
