#!/usr/bin/env python3
"""Tests for E118 (the KC types in the best family and in the recordings; experiments/E118_kc_types.py): the composition
weights, the permuted labels, the weighted distances and response rates on made-up counts, the rule for part B and the
outcome on made-up readings; with the shared directory, part 0's identities."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E118_kc_types')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        H.data_targets_path().read_bytes()
        return True
    except (FileNotFoundError, OSError):
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='the shared directory is not available')


def _counts(n_kc=60, seed=0):
    rng = np.random.default_rng(seed)
    types = np.array(['ab', "a'b'", 'g'] * (n_kc // 3), dtype=object)
    idx = pd.MultiIndex.from_arrays([np.arange(n_kc), types], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    c = pd.DataFrame(rng.poisson(1.5, (n_kc, len(H.APL16_ODORS))).astype(float), index=idx, columns=list(H.APL16_ODORS))
    return c, types


def test_kc_weights():
    labels = np.array(['ab', "a'b'", 'g', 'g'], dtype=object)
    assert X.kc_weights(labels, g=4, apbp=0).tolist() == [1, 0, 4, 4]
    with pytest.raises(ValueError):
        X.kc_weights(np.array(['x'], dtype=object))


def test_permuted_labels_keep_each_claw_counts_labels():
    rng = np.random.default_rng(1)
    labels = rng.choice(['ab', "a'b'", 'g'], 300).astype(object)
    claws = rng.integers(3, 9, 300)
    p = X.permuted_labels(labels, claws, seed=0)
    assert not np.array_equal(p, labels)
    for c in np.unique(claws):
        m = claws == c
        assert sorted(p[m]) == sorted(labels[m])
    assert np.array_equal(p, X.permuted_labels(labels, claws, seed=0))


def test_weighted_distances_unit_weights_are_the_plain_reading():
    c, _ = _counts()
    a = X.weighted_pair_distances(c, np.ones(len(c)))
    b = H.readout_pair_distances(c, 'raw')
    assert np.allclose(H.canonical_pairs(a), H.canonical_pairs(b), atol=1e-14)
    a = X.weighted_pair_distances(c, np.ones(len(c)), 'saturating', 3.0)
    b = H.readout_pair_distances(c, 'saturating', 3.0)
    assert np.allclose(H.canonical_pairs(a), H.canonical_pairs(b), atol=1e-14)


def test_integer_weight_is_duplication():
    c, types = _counts()
    w = X.kc_weights(types, g=3, apbp=0)
    dup = pd.concat([c[types == 'ab']] + [c[types == 'g']] * 3)
    a = H.canonical_pairs(X.weighted_pair_distances(c, w))
    b = H.canonical_pairs(H.readout_pair_distances(dup, 'raw'))
    assert np.allclose(a, b, atol=1e-12)


def test_weighted_rate():
    r = pd.DataFrame([[1, 0], [0, 0], [1, 1]], dtype=float)
    assert X.weighted_rate(r, np.ones(3)) == pytest.approx(0.5)
    assert X.weighted_rate(r, np.array([0, 0, 1.0])) == pytest.approx(1.0)
    assert X.weighted_rate(r, np.array([2, 1, 1.0])) == pytest.approx((2 * 0.5 + 1) / 4)


def _comp(gain_on=0.0, gain_off=0.0, n_up=8, half_on=None):
    """A made-up composition table: base ρ's 0.69 / 0.76; at γ x4 (α′β′ x1) the APL-on ρ higher by `gain_on` on `n_up`
    wirings, the APL-off by `gain_off`."""
    wir = ['real', 'P1', 'M1', 'P2', 'M2', 'P3', 'M3', 'P4']
    rows = []
    for apbp in X.APBP_WEIGHTS:
        for g in X.G_WEIGHTS:
            on = gain_on if (g, apbp) == (4.0, 1.0) else 0.0
            off = gain_off if (g, apbp) == (4.0, 1.0) else 0.0
            for i, w in enumerate(wir):
                rows.append(dict(w_g=g, w_apbp=apbp, wiring=w, rho_on_120=0.69 + (on if i < n_up else 0.0),
                                 rho_off_120=0.76 + off, rho_on_120_halfA=0.69, rho_off_120_halfA=0.76,
                                 rho_on_120_halfB=0.69, rho_off_120_halfB=0.76))
            h = on if half_on is None else half_on
            rows.append(dict(w_g=g, w_apbp=apbp, wiring='population', rho_on_120=0.69 + on, rho_off_120=0.76 + off,
                             rho_on_120_halfA=0.69 + h, rho_off_120_halfA=0.76 + off, rho_on_120_halfB=0.69 + h,
                             rho_off_120_halfB=0.76 + off))
    return pd.DataFrame(rows)


def _a1(g_rate=0.034):
    return pd.DataFrame([dict(kc_type='ab', wiring='population', rate_on=0.089),
                         dict(kc_type='g', wiring='population', rate_on=g_rate)])


def _perm(gain_on=0.0):
    p = _comp(gain_on=gain_on)
    base = (p.w_g == 1.0) & (p.w_apbp == 1.0)
    # (as part A3 writes it: the all-KC setting once, for the first permutation only)
    return pd.concat([p.assign(perm_seed=s) if s == 0 else p[~base].assign(perm_seed=s) for s in range(3)], ignore_index=True)


def test_rule_composition_passes_and_fails():
    pub = X.published_type_rates()
    r = X.rule_b(_comp(gain_on=0.015), _perm(0.0), _a1(), pub)
    assert r['composition'] and r['run_b'] and r['outside'] is None
    assert X.outcome_a(r)['letter'] == 'b'
    # too few wirings
    r = X.rule_b(_comp(gain_on=0.015, n_up=5), _perm(0.0), _a1(), pub)
    assert not r['composition']
    # the other ρ lower by more than 0.005
    r = X.rule_b(_comp(gain_on=0.015, gain_off=-0.006), _perm(0.0), _a1(), pub)
    assert not r['composition']
    # a half without the gain
    r = X.rule_b(_comp(gain_on=0.015, half_on=0.005), _perm(0.0), _a1(), pub)
    assert not r['composition']
    # the permuted control gaining as much
    r = X.rule_b(_comp(gain_on=0.015), _perm(0.015), _a1(), pub)
    assert not r['composition']
    assert X.outcome_a(r)['letter'].startswith('a')


def test_rule_outside_with_published_ratios():
    pub = pd.DataFrame([dict(kc_type='g', rate=0.1, ratio_to_ab=0.8, source='made up'),
                        dict(kc_type='g', rate=0.12, ratio_to_ab=1.2, source='made up')])
    r = X.rule_b(_comp(), _perm(), _a1(0.034), pub)
    assert r['outside'] is True and r['run_b']
    r = X.rule_b(_comp(), _perm(), _a1(0.089), pub)
    assert r['outside'] is False and not r['run_b']
    assert X.outcome_a(r)['letter'] == 'a'


def test_no_published_rates_on_record():
    assert len(X.published_type_rates()) == 0


@needs_shared
@pytest.mark.skipif(QUICK, reason='reads the eight archived fits (about 1 minute)')
def test_part0(tmp_path):
    res = X.part0(tmp_path)
    assert res['passes']
