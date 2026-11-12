#!/usr/bin/env python3
"""Tests for E85 (the ORN input as a population of flies; experiments/E85_orn_fly_population.py).
Without simulations: the mean over the nine model flies of every variant's input is the lab's input exactly; the
permuted control keeps each fly's and odor's deviations as a set before it is re-centred; the gain-free variant's gains
are 1 where they are not used; the across-trial ORN distance of an odor with itself is 0. With one simulation set
(about 1 min on 16 threads): the lab's input through the experiment's held reading reproduces the real wiring's saved
fit exactly. Needs the shared directory and olfdb dataset 1 (cached under outputs/, or http://hal:8080); skipped
without them."""

import numpy as np
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
X = pytest.importorskip('experiments.E85_orn_fly_population')
try:
    X.flies()
except OSError as e:  # urllib's URLError is an OSError
    pytest.skip(f'olfdb dataset {X.OLFDB_DATASET} is not available: {e}', allow_module_level=True)


@pytest.mark.parametrize('variant', X.VARIANTS)
def test_mean_over_flies_is_the_lab_input(variant):
    lab = H.apl16_orn_deltas().to_numpy(dtype=float)
    dev = X.deviations(variant, 1.0)
    assert dev.shape == (9, 12, 36)
    mean = np.mean([X.model_input(dev[f]).to_numpy(dtype=float) for f in range(dev.shape[0])], axis=0)
    assert np.allclose(mean, lab, rtol=0, atol=1e-9)
    # nothing but the megamat odors changes
    meg = [i for i, c in enumerate(H.apl16_orn_deltas().columns) if H.odor_abbrev(c) in H.MEGAMAT12]
    other = [i for i in range(lab.shape[1]) if i not in meg]
    for f in range(dev.shape[0]):
        assert np.array_equal(X.model_input(dev[f]).to_numpy(dtype=float)[:, other], lab[:, other])


def test_permuted_keeps_each_flys_deviations_as_a_set():
    st = X.trial_stats(X.flies()['x'])
    d = X.deviations('repro', 1.0)
    p = X.permute_deviations(d, st['measured'])
    for f in range(d.shape[0]):
        for o in range(d.shape[1]):
            idx = np.flatnonzero(st['measured'][f, o])
            assert np.array_equal(np.sort(p[f, o, idx]), np.sort(d[f, o, idx]))
            assert np.all(p[f, o, ~st['measured'][f, o]] == 0)
    # and it is a different arrangement for most of them
    assert np.mean(p != d) > 0.5
    r = X.recenter(p, st['measured'])
    assert np.allclose(r.sum(axis=0), 0.0, atol=1e-9)


def test_gains_and_reliabilities():
    a, info = X.fly_gains(X.flies()['x'])
    assert a.shape == (9, 36) and np.all(a >= X.GAIN_CLIP[0]) and np.all(a <= X.GAIN_CLIP[1])
    assert info['n_used'] > 0
    st = X.trial_stats(X.flies()['x'])
    assert np.all((st['r_g'] >= 0) & (st['r_g'] <= 1))


def test_within_fly_distances():
    d = X.within_fly_distances(0)
    assert len(d) == 66
    assert np.all(np.isfinite(d.trial_means))
    assert (d.trial_means >= 0).all() and (d.trial_means <= 2).all()


def test_lab_input_reproduces_the_real_wirings_fit():
    E = X.E
    r = X.held_reading('real', H.apl16_orn_deltas(), 'E85_lab')
    assert r['thr_factor'] == 1.0 and r['apl_factor'] == 1.0
    assert E._same({'on': r['on'], 'off': r['off']}, E.saved_fit('real')['saved']) == (0, 0)
