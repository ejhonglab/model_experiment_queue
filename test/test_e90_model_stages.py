#!/usr/bin/env python3
"""Tests for E90 (the model's stages against the flies' stages; experiments/E90_model_stages.py): Remy's per-fly window
means reproduce data/remy_rdm_timeseries/windows/; the block-average estimate recovers a known correlation on synthetic
units at the flies' reliability (and the as-measured reading is attenuated); the claw inputs formed from the bouton
rates and APL's activity reproduce olfsysm's on toy arrays, and its claw_sims on the real wiring (slow, skipped with
QUICK=1); the 28 shared glomeruli are those of the note of 2026-11-07; the groups have their sizes; the verdict function
on made-up tables. The synthetic tests need no shared data; the others are skipped without the shared directory."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E90_model_stages')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='needs the shared directory')


def test_groups_have_their_sizes():
    m = X.group_masks(X.megamat_pairs())
    assert {g: int(v.sum()) for g, v in m.items()} == dict(ester=15, alcohol=6, ester_alcohol=24, other=21, near=6, all66=66)
    assert not (m['ester'] & m['alcohol']).any() and not (m['other'] & m['ester_alcohol']).any()


@needs_shared
@pytest.mark.parametrize('stage,name', [('bouton', 'pn_boutons'), ('claw', 'kc_claws')])
@pytest.mark.parametrize('window', ['odor_0.5-3', 'early_0.5-1', 'late_2-3'])
def test_window_means_reproduce_windows(stage, name, window):
    lo, hi = (float(x) for x in window.split('_')[1].split('-'))
    wm = X.remy_window_means(stage, (lo, hi))
    _, _, odors = X.remy_matrices(stage, 'stim')
    for i in range(wm.shape[0]):
        c = pd.read_csv(X.remy_dir() / 'windows' / f'{name}_fly{i + 1}_corrdist_{window}.csv', index_col=0)
        np.testing.assert_allclose(c.loc[odors, odors].to_numpy(), wm[i], atol=1e-12)


def _units_with_corr(rng, rho, n_units=3000, n_odors=4):
    """Odor vectors (odors x units) whose population correlation between odors 0 and 1 is rho (the others independent)."""
    base = rng.standard_normal((n_odors, n_units))
    base[1] = rho * base[0] + np.sqrt(1 - rho ** 2) * base[1]
    return base


@pytest.mark.parametrize('rel', [0.65, 0.72])
@pytest.mark.parametrize('rho', [0.2, 0.6])
def test_block_average_estimate_recovers_known_correlation(rel, rho):
    rng = np.random.default_rng(1)
    est, meas = [], []
    for _ in range(8):
        x = _units_with_corr(rng, rho)
        sig = np.sqrt(x.var(axis=1) * (1 - rel) / rel)
        trials = x[None] + sig[None, :, None] * rng.standard_normal((3,) + x.shape)
        m, nr = X.block_average_distances(trials)
        est.append(nr[0, 1])
        meas.append(m[0, 1])
        # the diagonal of the estimate is 0 (each odor fully correlated with itself once the noise is removed)
        assert abs(nr[0, 0]) < 1e-12
    truth = 1 - np.corrcoef(x)[0, 1]
    assert abs(np.mean(est) - (1 - rho)) < 0.03
    # the trial means are attenuated (farther apart than the truth) by about the three-trial reliability
    r3 = 3 * rel / (1 + 2 * rel)
    assert abs(np.mean(meas) - (1 - rho * r3)) < 0.03
    assert np.mean(meas) > np.mean(est)


def test_noisy_reading_at_reliability():
    rng = np.random.default_rng(2)
    odors = list(X.ODORS)
    Z = pd.DataFrame(rng.standard_normal((2000, len(odors))), columns=odors)
    Z.iloc[:, 1] = 0.7 * Z.iloc[:, 0] + np.sqrt(1 - 0.49) * Z.iloc[:, 1]
    rel = pd.Series(0.65, index=odors)
    nz = X.noisy_reading(Z, rel, 4, rng, n_reps=3)
    truth = X.odor_distances(Z)
    assert float((nz['noise_removed'] - truth).abs().max()) < 0.03
    pair = X._canon(odors[0], odors[1])
    # the correlated pair is farther apart as measured (r 0.7 attenuated by the three-trial reliability)
    assert nz['measured'][pair] - truth[pair] > 0.05


def test_claw_inputs_toy():
    rng = np.random.default_rng(3)
    W = (rng.random((50, 7)) < 0.2).astype(float)
    b = rng.uniform(0, 40, (7, 30))
    inh = rng.uniform(0, 5, 30)
    w = rng.uniform(0, 3, 50)
    x = X.claw_inputs(W, b, inh, w)
    expect = np.maximum(W @ b - w[:, None] * inh[None, :], 0)
    assert np.array_equal(x, expect)
    assert (x >= 0).all() and (x == 0).any()


@needs_shared
def test_shared_glomeruli_are_the_notes():
    note = pd.read_csv(H.shared_dir() / 'notes' / 'gh146_pn_vs_orn_gain_2026-11-07.csv')
    assert X.shared_glomeruli() == sorted(note.glomerulus)
    assert len(X.shared_glomeruli()) == 28


@needs_shared
def test_orn_flies_reproduce_consensus():
    m = X.flies_stage_mean('orn')
    c = X._cs(H.orn16_consensus_pairs()).reindex(X.megamat_pairs())
    assert float((m - c).abs().max()) < 1e-9


def test_differs_rule():
    # far, outside the interval, 7 of 8 wirings on the side: differs
    r = X.differs(0.70, [0.70, 0.72, 0.68, 0.71, 0.69, 0.73, 0.74, 0.55], 0.60, 0.57, 0.63)
    assert r['differs'] and r['n_same_side'] == 7
    # 0.05 from the flies but inside a wide interval: does not
    assert not X.differs(0.65, [0.65] * 8, 0.60, 0.50, 0.70)['differs']
    # outside the interval but less than 0.05 away: does not
    assert not X.differs(0.64, [0.64] * 8, 0.60, 0.58, 0.62)['differs']
    # five of eight on the side: does not
    assert not X.differs(0.70, [0.8, 0.8, 0.8, 0.8, 0.8, 0.5, 0.5, 0.5], 0.60, 0.58, 0.62)['differs']


def _res(differs=False, mixed=False):
    return dict(differs=differs, mixed=mixed, n_differ=0, n_readings=3, signs=[])


def test_group_verdict_on_made_up_tables():
    v = X.group_verdict
    base = dict(bouton=_res(), claw=_res(), step_on=_res(), step_off=_res())
    assert v(dict(base, claw=_res(True)), 0.0) == 'upstream'
    assert v(dict(base, bouton=_res(True), step_on=_res(True)), 0.0) == 'upstream'
    assert v(dict(base, claw=_res(mixed=True), step_on=_res(True)), 0.0) == 'not resolved'
    assert v(dict(base, step_on=_res(True)), 0.0) == 'the step with APL on'
    assert v(dict(base, step_off=_res(True)), 0.0) == 'the step without APL'
    assert v(dict(base, step_on=_res(True), step_off=_res(True)), 0.0) == 'both steps'
    assert v(dict(base, step_on=_res(True), step_off=_res(mixed=True)), 0.0) == 'not resolved'
    assert v(base, 0.05) == 'not resolved'
    assert v(base, 0.01) == 'none'


def test_element_result():
    rows = pd.DataFrame(dict(differs=[True, True, False], diff=[0.1, 0.08, 0.02]))
    r = X.element_result(rows)
    assert r['mixed'] and not r['differs'] and r['n_differ'] == 2
    r = X.element_result(rows.iloc[:2])
    assert r['differs'] and not r['mixed'] and r['signs'] == [1]


@needs_shared
@pytest.mark.skipif(QUICK, reason='slow (builds a simulator on the real wiring); skipped with QUICK=1')
def test_claw_inputs_reproduce_claw_sims():
    res = X.claw_identity()
    assert res['counts_differ'] == 0
    assert res['claw_values_differ'] == 0 and res['exact']


def test_cross_half_recovers_the_mean_pattern():
    # flies = a shared pattern plus fly-specific deviations; the cross-half estimate recovers the shared pattern's
    # distances, while the per-fly distances are larger
    rng = np.random.default_rng(4)
    n_odors, n_gl, n_flies = len(X.ODORS), 400, 8
    base = rng.standard_normal((n_odors, n_gl))
    base[1] = 0.8 * base[0] + 0.6 * base[1]
    flies = base[None] + 0.7 * rng.standard_normal((n_flies, n_odors, n_gl))
    est = X.cross_half_distances(flies)
    truth = X._pairs_from_square(1 - np.corrcoef(base), X.ODORS)
    assert float((est - truth).abs().mean()) < 0.03
    pair = X._canon(X.ODORS[0], X.ODORS[1])
    per_fly = np.mean([X._pairs_from_square(1 - np.corrcoef(f), X.ODORS)[pair] for f in flies])
    assert abs(est[pair] - truth[pair]) < 0.03
    # the correlated pair is farther apart per fly (the fly-specific deviations attenuate it)
    assert per_fly - truth[pair] > 0.1
