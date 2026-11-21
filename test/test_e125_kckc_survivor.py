#!/usr/bin/env python3
"""Tests for E125 (the part of the inhibition among KCs that outlasts the block; experiments/E125_kckc_survivor.py): the
survivor forms' matrices and totals, the two conditions' matrices, the dispatch of the APL-off condition to its own
simulator (on pyolfsysm's hold_rates with a made-up simulator), the planned points, the "does what the matrix does" rule,
the corner rule and the outcomes' flags. No model simulation."""

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
X = pytest.importorskip('experiments.E125_kckc_survivor')
E122 = pytest.importorskip('experiments.E122_feedback_neuron')


def _toy(seed=0, n_g=60, n_apbp=30, n_ab=80):
    rng = np.random.default_rng(seed)
    types = np.array(['g'] * n_g + ["a'b'"] * n_apbp + ['ab'] * n_ab, dtype=object)
    n = len(types)
    ref = rng.poisson(0.5, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.15)
    ref[types[:, None] == types[None, :]] *= 3.0
    np.fill_diagonal(ref, 0.0)
    ref *= n / ref.sum()
    return ref, types


def test_survivor_forms():
    ref, types = _toy()
    for f in ('PO', 'DG', 'TY'):
        assert np.isclose(X.survivor_array(f, ref, types).sum(), ref.sum())
    assert np.array_equal(X.survivor_array('PO', ref, types), E122.arm_array('U', ref, types))
    assert np.array_equal(X.survivor_array('TY', ref, types), E122.arm_array('TY', ref, types))
    D = X.survivor_array('DG', ref, types)
    assert np.count_nonzero(D - np.diag(np.diag(D))) == 0
    for t in X.TYPES:
        m = types == t
        # each KC its type's mean amount (the hemibrain's mean row sum over the type)
        assert np.allclose(np.diag(D)[m], ref[m].sum(axis=1).mean())
        assert np.isclose(D[m].sum(), ref[m].sum())
    assert not X.survivor_array('Z', ref, types).any()
    with pytest.raises(ValueError):
        X.survivor_array('nope', ref, types)


@pytest.mark.parametrize('arm', list(X.ARMS))
def test_condition_arrays(arm):
    ref, types = _toy()
    d = X.dose_of(arm)
    c = X.condition_arrays(arm, ref, types)
    Xs = X.survivor_array(X.form_of(arm), ref, types)
    assert np.allclose(c['on'], (1 - d) * ref + d * Xs)
    assert np.allclose(c['off'], d * Xs)
    assert np.isclose(c['on'].sum(), ref.sum()) and np.isclose(c['off'].sum(), d * ref.sum())
    assert (c['on'] >= 0).all() and (c['off'] >= 0).all()


def test_z0_and_doses():
    ref, types = _toy()
    c = X.condition_arrays('Z0', ref, types)
    assert np.array_equal(c['on'], ref) and not c['off'].any()
    assert [X.dose_of(a) for a in X.ARMS] == [0.25, 0.5, 0.75] * 2 + [0.5]
    assert [X.form_of(a) for a in X.ARMS] == ['PO'] * 3 + ['DG'] * 3 + ['TY']


def test_plan():
    pl = X.planned()
    assert len(pl) == 56 and len(set(pl)) == 56
    assert [a for a, _ in pl[:24]] == [a for a in X.PARTS['arm1'] for _ in X.WIRINGS]
    assert X.planned(wirings=['P3']) == [(a, 'P3') for a in X.RUN_ORDER]
    assert X.runnable_arms(dict(arm2_testable=False)) == X.PARTS['arm1'] + X.PARTS['arm3']
    assert X.runnable_arms(dict(arm2_testable=True)) == X.ARMS


def test_dispatch_on_hold_rates():
    """The APL-off condition (marked by off_extra) goes to its own simulator, the APL-on calls to the other, and the
    mark never reaches a simulator."""
    from pyolfsysm import tuning
    calls = {'on': [], 'off': []}
    kcs, odors = pd.RangeIndex(200), pd.RangeIndex(4)

    def fake(cond, strength):
        def sim(p):
            assert X.OFF_MARK not in p
            calls[cond].append(dict(p))
            drive = np.linspace(0, 400, len(kcs))[:, None] * np.ones((1, len(odors)))
            counts = np.floor(np.clip(drive / p['fixed_thr'] - strength * p['wAPLKC_scale'], 0, None))
            c = pd.DataFrame(counts, index=kcs, columns=odors)
            return (c >= 2).astype(float), c
        return sim

    sim = X.dispatch(fake('on', 1.0), fake('off', 1.0))
    p = dict(fixed_thr=100.0, wAPLKC_scale=1.0, wKCAPL_scale=1.0, wPNAPL_scale=1.0, wAPLPN_scale=0.5, wKCKC_scale=2.0)
    res = tuning.hold_rates(sim, p, target_off=300, target_on=150, size=800,
                            off_params=dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0, **{X.OFF_MARK: 1.0}),
                            tied_on=('wAPLKC_scale', 'wKCAPL_scale', 'wPNAPL_scale'))
    assert calls['off'] and calls['on']
    assert all(c['wAPLKC_scale'] == 0.0 for c in calls['off'])
    assert all(c['wAPLKC_scale'] > 0.0 for c in calls['on'])
    assert all(c['wKCKC_scale'] == 2.0 for c in calls['on'] + calls['off'])
    assert res['n_off'] > 0
    with pytest.raises(RuntimeError):
        sim(dict(p, wAPLKC_scale=0.0))
    with pytest.raises(RuntimeError):
        sim(dict(p, **{X.OFF_MARK: 1.0}))


def _rd(vals):
    """readings rows: vals[kind] = pattern part (all targets, both populations, each wiring)."""
    rows = []
    for k, v in vals.items():
        for u in X.POPS + ('all8',) + tuple(X.WIRINGS):
            r = dict(kind=k, unit=u, readout='cal', mean_effect=0.17)
            for t in X.TARGETS:
                r[f'{t}_pattern'] = v
            rows.append(r)
    return pd.DataFrame(rows)


def test_matrix_like():
    rd = _rd({'HBoff50': 0.0537, 'HBoff75': 0.0529, 'PO50': 0.0560, 'PO75': 0.0535, 'PO25': 0.05})
    assert X.matrix_like(rd, 'PO25', True)['matrix_like'] is None
    r = X.matrix_like(rd, 'PO50', True)
    assert r['matrix_like'] and r['rung'] == 'HBoff50'
    assert not X.matrix_like(rd, 'PO50', False)['matrix_like']
    rd2 = _rd({'HBoff75': 0.0529, 'PO75': 0.0560})
    assert not X.matrix_like(rd2, 'PO75', True)['matrix_like']


def test_corners():
    runs = pd.DataFrame(dict(kind='DG50', wiring=list(X.WIRINGS), at_bound_thr=[True, False, True] + [False] * 5,
                             at_bound_apl=[False, False, False, True] + [False] * 4))
    c = X.corners(runs, 'DG50')
    assert c['corner_wirings'] == sorted(X.WIRINGS[i] for i in (0, 2, 3))
    assert c['corner']
    runs['at_bound_thr'] = False
    c = X.corners(runs, 'DG50')
    assert len(c['corner_wirings']) == 1 and not c['corner']


def _c(kg, ml=None, corner=False):
    return dict(complete=True, keeps_gain=kg, matrix_like=ml, corner=corner, cls='x')


def test_outcomes():
    ladder = {'HBoff25': dict(keeps_gain=False), 'HBoff50': dict(keeps_gain=True), 'HBoff75': dict(keeps_gain=True)}
    # (a) and (c): the pooled form does what the matrix does, the own feedback does not
    cl = {'PO25': _c(False), 'PO50': _c(True, True), 'PO75': _c(True, True), 'DG25': _c(False), 'DG50': _c(False, False),
          'DG75': _c(False, False), 'TY50': _c(True, False)}
    o = X.outcomes(cl, ladder)
    assert o['a'] and o['c'] and not o['b'] and not o['d']
    # (b): neither keeps the gain where the ladder does (a keep at 0.25, where the ladder does not, does not count)
    cl.update(PO25=_c(True), PO50=_c(False, False), PO75=_c(False, False))
    o = X.outcomes(cl, ladder)
    assert o['b'] and not o['a'] and not o['c']
    # (d): a corner arm; it counts for nothing else
    cl.update(DG75=_c(True, True, corner=True))
    o = X.outcomes(cl, ladder)
    assert o['d'] and o['forms']['DG']['corner_arms'] == ['DG75'] and not o['a']
    # arm 2 not run: only the pooled form is read
    o = X.outcomes({k: v for k, v in cl.items() if not k.startswith('DG')}, ladder)
    assert o['complete_forms'] == ['PO', 'TY'] and not o['c']
