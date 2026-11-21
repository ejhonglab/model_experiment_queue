#!/usr/bin/env python3
"""Tests for E123 (the form of the inhibition among KCs; experiments/E123_kckc_form.py) on made-up matrices and
readings: the arms' matrices (total weight, what each keeps), the planned points, the classes, the draws' verdict, the
step-2 gate, and E77_input_screen.hold_rates's off_extra. No model simulation."""

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E123_kckc_form')


def _toy(seed=0):
    rng = np.random.default_rng(seed)
    types = np.array(['g'] * 7 + ["a'b'"] * 4 + ['ab'] * 9, dtype=object)
    n = len(types)
    ref = rng.poisson(0.6, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.4)
    np.fill_diagonal(ref, 0.0)
    ref *= n / ref.sum()
    return ref, types


@pytest.mark.parametrize('arm', ['SH1', 'SH2', 'SH3', 'SH0', 'AK', 'K4', 'K16', 'HBab', 'SH0ab', 'U'])
def test_total_weight_matched(arm):
    ref, types = _toy()
    M = X.arm_array(arm, ref, types)
    assert M.shape == ref.shape and (M >= 0).all()
    assert np.isclose(M.sum(), ref.sum())


def test_seeds_and_e122():
    E122 = pytest.importorskip('experiments.E122_feedback_neuron')
    ref, types = _toy()
    assert np.array_equal(X.arm_array('SH0', ref, types), E122.arm_array('SH', ref, types))
    assert np.array_equal(X.arm_array('SH0ab', ref, types), X.arm_array('SH0', ref, types))
    assert np.array_equal(X.arm_array('HBab', ref, types), ref)
    assert np.array_equal(X.arm_array('U', ref, types), E122.arm_array('U', ref, types))
    a, b = X.arm_array('SH1', ref, types), X.arm_array('SH2', ref, types)
    assert not np.array_equal(a, b)


def test_amount_kept():
    ref, types = _toy()
    M = X.amount_kept_shuffle(ref, types, 0)
    assert np.allclose(np.diag(M), 0.0)
    assert not np.array_equal(M, ref)
    for t in X.TYPES:
        m = types == t
        assert np.allclose(M[:, m].sum(axis=1), ref[:, m].sum(axis=1))
        # each row's values from type t are a permutation of the reference's
        for i in range(len(types)):
            c = np.flatnonzero(m & (np.arange(len(types)) != i))
            assert np.allclose(np.sort(M[i, c]), np.sort(ref[i, c]))


@pytest.mark.parametrize('k', [1, 4, 16])
def test_random_k_partners(k):
    n, total = 20, 20.0
    M = X.random_k_partners(n, k, total, 0)
    assert np.allclose(np.diag(M), 0.0)
    assert ((M > 0).sum(axis=1) == k).all()
    assert np.allclose(M.sum(axis=1), total / n)
    assert np.isclose(M.sum(), total)
    assert len(np.unique(M[M > 0])) == 1
    with pytest.raises(ValueError):
        X.random_k_partners(n, n, total)


def test_planned_order_and_reads():
    p = X.planned()
    assert p[0] == ('SH1', 'real')
    assert ('U', 'real') not in p and ('SH0', 'P1') not in p
    assert ('U', 'P3') in p and ('SH0', 'M3') in p
    assert len(p) == len(X.ARMS) * 8 - 2 * 4
    first40 = {a for a, _ in p[:40]}
    assert first40 == {'SH1', 'SH2', 'SH3', 'HBab', 'SH0ab'}
    assert X.planned(wirings=['P3', 'P4']) == [(a, w) for a in X.RUN_ORDER for w in ('P3', 'P4')]


def _rd(arm_pattern, fam=0.050, zero=0.061):
    rows = []
    for kind, val in (('family', fam), ('E110_A', zero), ('X', arm_pattern)):
        for unit in list(X.WIRINGS) + ['S4', 'T4', 'all8']:
            v = val[unit] if isinstance(val, dict) else val
            rows.append(dict(kind=kind, unit=unit, readout='cal', mean_effect=0.17,
                             **{f'{t}_pattern': v for t in X.TARGETS}, all_loss=0.08))
            rows.append(dict(kind=kind, unit=unit, readout='raw', mean_effect=0.12,
                             **{f'{t}_pattern': v for t in X.TARGETS}, all_loss=0.08))
    return pd.DataFrame(rows)


def _runs(in_bands=True, rho=0.7):
    rows = [dict(kind=k, wiring=w, in_bands=(in_bands if k == 'X' else True), rate_rho_on=rho)
            for k in ('family', 'E110_A', 'X') for w in X.WIRINGS]
    return pd.DataFrame(rows)


def test_classes():
    c = X.classify(_rd(0.051), _runs(), 'X')
    assert c['cls'] == 'equal' and c['keeps_gain'] and c['equal_two_sided']
    c = X.classify(_rd(0.051), _runs(in_bands=False), 'X')
    assert c['cls'] == 'keeps the gain' and not c['equal']
    c = X.classify(_rd(0.054), _runs(), 'X')
    assert c['cls'] == 'keeps the gain'            # 0.054 <= half-way 0.0555, but 0.004 above the family
    c = X.classify(_rd(0.058), _runs(), 'X')
    assert c['cls'] == 'neither'
    c = X.classify(_rd(0.040), _runs(), 'X')
    assert c['cls'] == 'equal' and not c['equal_two_sided']
    # 3 of 8 wirings above the family by more than 0.003: not equal
    v = {w: (0.060 if w in ('real', 'P1', 'P3') else 0.050) for w in X.WIRINGS}
    v.update(S4=0.051, T4=0.051, all8=0.051)
    assert X.classify(_rd(v), _runs(), 'X')['cls'] == 'keeps the gain'
    assert not X.classify(_rd(0.05), _runs()[lambda d: d.wiring != 'M3'], 'X')['complete']


def test_draws_and_step2():
    runs = pd.concat([_runs(rho=0.7)] + [_runs(rho=r).assign(kind=s) for s, r in
                                          zip(X.DRAWS, (0.6, 0.62, 0.71, 0.6))], ignore_index=True)
    cls = {s: dict(complete=True, cls='equal') for s in X.DRAWS}
    d = X.draws_verdict(cls, runs)
    assert d['agree'] and d['seeds_below_family'] == 3 and d['costs_rate_rho']
    cls['SH2'] = dict(complete=True, cls='neither')
    assert not X.draws_verdict(cls, runs)['agree']
    classes = {a: dict(complete=True, cls='neither') for a in X.ARMS}
    classes.update(K16=dict(complete=True, cls='equal'), K64=dict(complete=True, cls='equal'),
                   AK=dict(complete=True, cls='equal'), HBab=dict(complete=True, cls='equal'))
    p = X.step2_plan(classes)
    assert p['complete'] and p['refit_sets'] == ['K16', 'HBab']
    classes.update(K16=dict(complete=True, cls='neither'), K64=dict(complete=True, cls='neither'))
    assert X.step2_plan(classes)['refit_sets'] == ['AK', 'HBab']


def test_off_extra(monkeypatch):
    E77 = pytest.importorskip('experiments.E77_input_screen')
    from pyolfsysm import tuning
    seen = {}
    monkeypatch.setattr(E77, 'c2_parts', lambda: dict(p=dict(fixed_thr=1.0), target_off=1, target_on=1, size=1))

    def fake(sim, p, **kw):
        seen.update(kw)
        return dict(params=p, params_off=p)
    monkeypatch.setattr(tuning, 'hold_rates', fake)
    E77.hold_rates(lambda: (lambda q: None), params=dict(fixed_thr=1.0))
    assert seen['off_params'] == dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    E77.hold_rates(lambda: (lambda q: None), params=dict(fixed_thr=1.0), off_extra=dict(wKCKC_scale=0.0))
    assert seen['off_params'] == dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0, wKCKC_scale=0.0)
