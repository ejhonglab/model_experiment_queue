#!/usr/bin/env python3
"""Tests for E87 (the form of APL's inhibition in the flies' own responses; experiments/E87_inhibition_form.py), on the
flies' recordings only (no simulations): with no inhibition the data's effect is the registry's thresholded value and a
form's effect within a fly is zero; the uniform divisive form reproduces the density-matched row of the data-only
confound controls (E7 part C: criterion 1.85, mean distance 0.845 with APL output blocked); every form's constant gives
the control flies' response rate; the divisive kind does not depend on a rescaling of a recording (with its
criterion); the permuted controls keep the strengths and the breadths as sets; the halves are disjoint, cover the
recordings and hold both genotypes. Needs the shared directory; skipped without it."""

import numpy as np
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
X = pytest.importorskip('experiments.E87_inhibition_form')


@pytest.fixture(scope='module')
def rd():
    return X.Reading('z', 'trial_means')


def test_no_inhibition_is_the_registry_value(rd):
    assert np.isclose(np.nanmean(rd.data), H.data_target('apl_effect_thresholded_mean'), atol=1e-12)
    pred, c = rd.predict(X.UNIFORM, c=0.0)
    assert c == 0.0 and np.allclose(pred, 0.0, atol=1e-12)


def test_uniform_divisive_is_e7_density_matched(rd):
    pred, c = rd.predict(X.UNIFORM_DIV)
    assert abs((1 + c) - 1.85) < 0.01
    assert abs(float(np.nanmean(rd.d_off + pred)) - 0.845) < 0.001


@pytest.mark.parametrize('form', [X.UNIFORM, X.Form('subtractive', 1, 0.25), X.Form('divisive', 2, 1, 'summed'),
                                  ('similar', 2)])
def test_constant_gives_the_control_rate(rd, form):
    pred, c = rd.predict(form)
    red = rd.reduced(form, c)
    rate = np.mean([(m > 0).mean() for m in red])
    assert abs(rate - rd.rate) < 2e-4


def test_constant_gives_the_control_rate_single_trials():
    rn = X.Reading('z', 'noise_free')
    form = X.Form('subtractive', 1, 0.5)
    inhs = rn.inhibitions(form)
    c = X.fit_c(rn.tr_tnt, [np.broadcast_to(i, x.shape) for i, x in zip(inhs, rn.tr_tnt)], rn.crit1, rn.rate1)
    red = [X.apply_form(x, c, i[None], 'subtractive', rn.crit1) for x, i in zip(rn.tr_tnt, inhs)]
    assert abs(np.mean([(r > 0).mean() for r in red]) - rn.rate1) < 2e-4


def test_divisive_kind_does_not_depend_on_scale(rd):
    m = rd.tm_tnt[0]
    f = X.Form('divisive', 1, 0.5, 'summed')
    inh = X.inhibition(m, f, rd.crit)
    c = X.fit_c([m], [inh], rd.crit, 0.07)
    s = 3.7
    inh_s = X.inhibition(s * m, f, s * rd.crit)
    c_s = X.fit_c([s * m], [inh_s], s * rd.crit, 0.07)
    a = X.cos_dist(X.apply_form(m, c, inh, 'divisive', rd.crit))
    b = X.cos_dist(X.apply_form(s * m, c_s, inh_s, 'divisive', s * rd.crit))
    assert np.allclose(a, b, equal_nan=True)


def test_permuted_controls_keep_sets(rd):
    f = X.Form('subtractive', 1, 1, 'count')
    A = X._perm_strengths(rd, 3)
    B = X._perm_breadths(rd, f, 3)
    for n, m in enumerate(rd.tm_tnt):
        assert np.allclose(np.sort(A[n]), np.sort(X.strengths(m, rd.crit)))
        b = X.breadths(m, rd.crit, 'count')
        assert np.allclose(np.sort(B[n]), np.sort(b))
        assert np.array_equal(B[n] > 0, b > 0)


@pytest.mark.parametrize('seed', [0, 1, 5])
def test_halves(seed):
    h = X._halves(seed)
    tm = X.trialmean_mats('z')
    for g, key in (('ctrl', 'TNT3c'), ('tnt', 'TNT')):
        a, b = set(h['A'][g]), set(h['B'][g])
        assert a and b and not (a & b)
        assert a | b == set(tm[key])
