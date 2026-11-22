#!/usr/bin/env python3
"""Tests for E126 (the ladder of the pooled survivor; experiments/E126_kckc_survivor_ladder.py): the arms' names, the
random groups (the types' sizes, seeded, the survivor on the types E122's TY), the two conditions' matrices (E125's for
TY), the planned points, the "differ" rule, the gate and the outcomes' flags. No model simulation."""

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
X = pytest.importorskip('experiments.E126_kckc_survivor_ladder')
E122 = pytest.importorskip('experiments.E122_feedback_neuron')
E125 = pytest.importorskip('experiments.E125_kckc_survivor')


def _toy(seed=0, n_g=60, n_apbp=30, n_ab=80):
    rng = np.random.default_rng(seed)
    types = np.array(['g'] * n_g + ["a'b'"] * n_apbp + ['ab'] * n_ab, dtype=object)
    n = len(types)
    ref = rng.poisson(0.5, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.15)
    ref[types[:, None] == types[None, :]] *= 3.0
    np.fill_diagonal(ref, 0.0)
    ref *= n / ref.sum()
    return ref, types


def test_names():
    assert X.STEP1_ARMS == ('TY75', 'TY25') and X.STEP2_ARMS == ('RG75s0', 'RG75s1')
    assert [X.dose_of(a) for a in X.ARMS] == [0.75, 0.25, 0.75, 0.75]
    assert [X.form_of(a) for a in X.ARMS] == ['TY', 'TY', 'RG', 'RG']
    assert [X.seed_of(a) for a in X.ARMS] == [None, None, 0, 1]
    assert [X.step_of(a) for a in X.ARMS] == [1, 1, 2, 2]
    assert X.dose_of('Z0') == 0 and X.form_of('Z0') == 'Z'
    # E125's arms read with this experiment's parsing
    assert [X.dose_of(a) for a in E125.ARMS] == [E125.dose_of(a) for a in E125.ARMS]
    assert X.tag_of('TY25', 'P3') == 'S1_TY25_P3' and X.tag_of('RG75s1', 'real') == 'S2_RG75s1_real'


def test_random_groups():
    ref, types = _toy()
    for s in X.GROUP_SEEDS:
        g = X.random_groups(types, s)
        assert sorted(g) == sorted(types)
        assert np.array_equal(g, X.random_groups(types, s))             # seeded
        assert (g != types).any()
    assert not np.array_equal(X.random_groups(types, 0), X.random_groups(types, 1))
    # the survivor on the types is E122's TY exactly; on the groups it keeps the total
    assert np.array_equal(X.group_survivor(ref, types), E122.arm_array('TY', ref, types))
    G = X.group_survivor(ref, X.random_groups(types, 0))
    assert np.isclose(G.sum(), ref.sum())
    g = X.random_groups(types, 0)
    for t in X.TYPES:
        m = g == t
        assert np.allclose(G[np.ix_(m, m)], ref[m].sum() / m.sum() ** 2)
        assert not G[np.ix_(m, ~m)].any()
    info = X.group_info('RG75s0', types)
    assert info['group_seed'] == 0 and 0 < info['group_same_as_type'] < 1
    assert X.group_info('TY75', types) == {}


@pytest.mark.parametrize('arm', list(X.ARMS) + ['TY50', 'Z0'])
def test_condition_arrays(arm):
    ref, types = _toy()
    d = X.dose_of(arm)
    c = X.condition_arrays(arm, ref, types)
    if X.form_of(arm) == 'TY' or arm == 'Z0':
        e = E125.condition_arrays(arm, ref, types)
        assert all(np.array_equal(c[k], e[k]) for k in ('on', 'off'))
    if arm == 'Z0':
        return
    Xs = (X.group_survivor(ref, X.random_groups(types, X.seed_of(arm))) if X.form_of(arm) == 'RG'
          else E125.survivor_array('TY', ref, types))
    assert np.allclose(c['on'], (1 - d) * ref + d * Xs)
    assert np.allclose(c['off'], d * Xs)
    assert np.isclose(c['on'].sum(), ref.sum()) and np.isclose(c['off'].sum(), d * ref.sum())


def test_plan():
    pl = X.planned()
    assert len(pl) == 32 and len(set(pl)) == 32
    assert [a for a, _ in pl[:8]] == ['TY75'] * 8                        # the gate's arm first
    assert X.planned(X.STEP1_ARMS, ['P3', 'real']) == [('TY75', 'real'), ('TY75', 'P3'), ('TY25', 'real'),
                                                         ('TY25', 'P3')]
    assert len(X.planned(X.STEP2_ARMS)) == 16


def _rd(vals):
    """readings rows: vals[kind] = pattern part, either a number (every unit and target) or a dict per unit."""
    rows = []
    for k, v in vals.items():
        for u in X.POPS + ('all8',) + tuple(X.WIRINGS):
            x = v.get(u, v.get('default')) if isinstance(v, dict) else v
            r = dict(kind=k, unit=u, readout='cal', mean_effect=0.17)
            for t in X.TARGETS:
                r[f'{t}_pattern'] = x
            rows.append(r)
    return pd.DataFrame(rows)


def _runs(kinds, bound=()):
    return pd.DataFrame([dict(kind=k, wiring=w, at_bound_thr=(k, w) in bound, at_bound_apl=False)
                         for k in kinds for w in X.WIRINGS])


def test_differ():
    runs = _runs(['TY75', 'PO75'])
    rd = _rd({'TY75': 0.050, 'PO75': 0.053})
    c = X.differ(rd, runs, 'TY75', 'PO75')
    assert c['verdict'] == 'TY75 lower' and c['n_wirings_a_lower'] == 8
    assert X.differ(rd, runs, 'PO75', 'TY75')['verdict'] == 'TY75 lower'
    # lower on 6 of 8 wirings only: not told apart
    v = {'default': 0.050, **{w: 0.060 for w in X.WIRINGS[:2]}}
    c = X.differ(_rd({'TY75': v, 'PO75': 0.053}), runs, 'TY75', 'PO75')
    assert c['n_wirings_a_lower'] == 6 and c['verdict'] == 'not told apart'
    # 7 of 8 wirings: differ, if lower on both populations
    v = {'default': 0.050, X.WIRINGS[0]: 0.060}
    assert X.differ(_rd({'TY75': v, 'PO75': 0.053}), runs, 'TY75', 'PO75')['verdict'] == 'TY75 lower'
    # higher on T (on every target): not told apart
    v = {'default': 0.050, 'T4': 0.054}
    assert X.differ(_rd({'TY75': v, 'PO75': 0.053}), runs, 'TY75', 'PO75')['verdict'] == 'not told apart'
    # a corner on three wirings: no class
    runs_c = _runs(['TY75', 'PO75'], bound={('PO75', w) for w in X.WIRINGS[:3]})
    assert X.differ(rd, runs_c, 'TY75', 'PO75')['verdict'] == 'no class (corner)'
    # incomplete
    assert X.differ(rd, runs[runs['wiring'] != 'P3'], 'TY75', 'PO75')['verdict'] == 'incomplete'


def _c(kg, ml=None, corner=False, low=()):
    return dict(complete=True, keeps_gain=kg, matrix_like=ml, corner=corner, cls='x', thr_below_60=list(low),
                fixed_thr_min=80.0, corner_wirings=[])


def test_gate():
    assert not X.gate({})['complete']
    assert X.gate({'TY75': _c(True, True)})['open']
    assert X.gate({'TY75': _c(True, False)})['open']                     # between (a) and (b): step 2 still runs
    assert not X.gate({'TY75': _c(False, False)})['open']
    assert not X.gate({'TY75': _c(True, True, low=('P1',))})['open']
    assert not X.gate({'TY75': _c(True, True, corner=True)})['open']


def test_outcomes():
    cl = {'TY75': _c(True, True), 'TY25': _c(True)}
    g = X.gate(cl)
    o = X.outcomes(cl, {}, g)
    assert o['a'] and not o['b'] and not o['between_a_b'] and o['types_needed'] is None and o['c'] is None
    # (a) and the draws fail: the types are needed
    cl2 = dict(cl, RG75s0=_c(False, False), RG75s1=_c(True, True))
    o = X.outcomes(cl2, {}, g)
    assert o['types_needed'] and not o['c']
    # (c): both draws keep the gain, TY75 not lower
    cl3 = dict(cl, RG75s0=_c(True, True), RG75s1=_c(True, False))
    o = X.outcomes(cl3, {'RG75s0|TY75': dict(verdict='not told apart')}, g)
    assert o['c'] and not o['types_needed']
    # TY75 lower than a draw: types needed, not (c)
    o = X.outcomes(cl3, {'RG75s1|TY75': dict(verdict='TY75 lower')}, g)
    assert o['types_needed'] and not o['c'] and o['ty75_lower_than'] == ['RG75s1']
    # (b): a threshold under 60 Hz
    cl4 = {'TY75': _c(True, True, low=('P3',))}
    o = X.outcomes(cl4, {}, X.gate(cl4))
    assert o['b'] and not o['a'] and not o['gate_open']
    # between (a) and (b)
    cl5 = {'TY75': _c(True, False)}
    o = X.outcomes(cl5, {}, X.gate(cl5))
    assert o['between_a_b'] and not o['a'] and not o['b'] and o['gate_open']
    # (d)
    o = X.outcomes({'TY75': _c(True, True), 'TY25': _c(True, corner=True)}, {}, g)
    assert o['d']
