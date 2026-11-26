"""Tests of E142's lock, stub and step 0's pure functions (experiments/E142_dendrite_bouton_step.py and
experiments/E142_step0.py): the weighted strength parts against repeated rows, the yardsticks against hal's no-step
score, the null's summary, the glomerulus units, the flies' group units, the classes read again against E136's
`classify`, and the lock's step 1 against hal's runner. No model is run and no recording is read."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments import E142_dendrite_bouton_step as E142
from experiments import E142_step0 as S0
from kcq import population_fit as P

ODORS = [f'o{i}' for i in range(6)]


def _named_pop(rng, glomeruli, name='f', crit=2.0, nan_frac=0.0):
    V = pd.DataFrame(rng.normal(1.5, 2.0, size=(len(glomeruli), len(ODORS))), index=glomeruli, columns=ODORS)
    if nan_frac:
        V = V.mask(rng.random(V.shape) < nan_frac)
    return P.Population.from_criterion(V, crit, name=name, named_units=True)


def test_weighted_strength_equals_repeated_rows():
    rng = np.random.default_rng(1)
    gl = [f'g{i}' for i in range(7)]
    pop = _named_pop(rng, gl, nan_frac=0.1)
    w = {g: int(k) for g, k in zip(gl, rng.integers(0, 5, len(gl)))}
    ws = S0.weighted_strength(pop, w, ODORS)
    rep = S0.repeated_population(pop, w)
    pv = P.population_views(rep, ODORS, parts=list(E142.STRENGTH_PARTS))
    for p in E142.STRENGTH_PARTS:
        np.testing.assert_allclose(ws[p].to_numpy(), pv[p].reindex(ODORS).to_numpy(), rtol=0, atol=1e-12)


def test_weighted_arrays_ones_are_population_views():
    rng = np.random.default_rng(2)
    gl = [f'g{i}' for i in range(5)]
    pop = _named_pop(rng, gl)
    a = S0.weighted_arrays(np.ones((1, 5)), pop.values.to_numpy(), pop.responds.to_numpy())
    pv = P.population_views(pop, ODORS, parts=list(E142.STRENGTH_PARTS))
    for p in E142.STRENGTH_PARTS:
        np.testing.assert_allclose(a[p][0], pv[p].to_numpy(), atol=1e-12)


def test_yardstick_scores_ones_equal_no_step_score():
    rng = np.random.default_rng(3)
    gl = [f'g{i}' for i in range(8)]
    up = [_named_pop(rng, [g for g in gl if rng.random() > 0.15] or gl, name=f'u{i}') for i in range(4)]
    down = [P.Population.from_criterion(pd.DataFrame(rng.normal(1.0, 2.0, size=(40, len(ODORS))), columns=ODORS), 2.0,
                                        name=f'd{i}') for i in range(3)]
    ref = P.fly_reference(down, ODORS, parts=list(E142.STRENGTH_PARTS))
    mats = S0.fly_matrices(up, gl, ODORS)
    got = S0.yardstick_scores(np.ones((1, len(gl))), mats, ref).iloc[0].view_strength
    want = P.no_step_score(up, ref).views['strength']
    assert got == pytest.approx(want, abs=1e-12)
    # many weightings at once give the same as one at a time
    W = rng.integers(0, 4, size=(5, len(gl))).astype(float)
    many = S0.yardstick_scores(W, mats, ref).view_strength.to_numpy()
    one = [S0.yardstick_scores(W[i:i + 1], mats, ref).iloc[0].view_strength for i in range(5)]
    np.testing.assert_allclose(many, one, atol=1e-12)


def test_fly_matrices_extra_rows():
    rng = np.random.default_rng(4)
    pop = _named_pop(rng, ['a', 'b'])
    extra = pd.DataFrame([[3.0, 1.0, np.nan, 5.0, 0.0, 2.5]], index=['c'], columns=ODORS)
    (V, Rr), = S0.fly_matrices([pop], ['a', 'b', 'c'], ODORS, extra, 2.0)
    np.testing.assert_allclose(V[2], extra.loc['c'].to_numpy(), equal_nan=True)
    np.testing.assert_allclose(Rr[2], [1, 0, np.nan, 1, 0, 1], equal_nan=True)


def test_permutations_and_null_summary():
    p = S0.permutations(5, 10, 0)
    assert p.shape == (10, 5) and all(sorted(r) == list(range(5)) for r in p)
    np.testing.assert_array_equal(p, S0.permutations(5, 10, 0))
    s = S0.null_summary(0.5, np.arange(1, 101) / 100.0)
    assert s['null_percentile_of_actual'] == pytest.approx(49.0)
    assert not s['below_p05'] and S0.null_summary(0.01, np.arange(1, 101) / 100.0)['below_p05']


def test_outside_zero():
    assert S0.outside_zero(0.1, 0.3) == 1 and S0.outside_zero(-0.3, -0.1) == -1
    assert S0.outside_zero(-0.1, 0.2) == 0 and S0.outside_zero(np.nan, 1.0) == 0


def test_glomerulus_units():
    V = pd.DataFrame([[1.0, 3.0], [3.0, 1.0], [5.0, 0.0]], columns=['x', 'y'])
    pop = P.Population.from_criterion(V, 2.0)
    eq = S0.glomerulus_units(pop, ['A', 'A', 'B'], 2.0)
    assert list(eq.values.index) == ['A', 'B']
    np.testing.assert_allclose(eq.values.to_numpy(), [[2.0, 2.0], [5.0, 0.0]])
    np.testing.assert_allclose(eq.responds.to_numpy(), [[0, 0], [1, 0]])
    with pytest.raises(ValueError):
        S0.glomerulus_units(pop, ['A'], 2.0)


def test_fly_group_units_reads_the_held_out_repeat():
    rng = np.random.default_rng(5)
    n = 30
    base = np.zeros((n, len(ODORS)))
    base[:10, :2] = 6.0     # one tuning
    base[10:20, 3:5] = 6.0  # another
    cols = pd.MultiIndex.from_product([ODORS, [0, 1, 2]], names=['odor', 'repeat'])
    X = np.concatenate([base[:, [i]] + rng.normal(0, 0.3, size=(n, 3)) for i in range(len(ODORS))], axis=1)
    tr = pd.DataFrame(X, columns=cols)
    pop, lab = S0.fly_group_units(tr, ODORS, ODORS, 2)
    assert set(np.unique(lab)) >= {0}
    assert len(pop.values) == len(np.unique(lab))
    held = tr.xs(2, axis=1, level='repeat')[ODORS]
    for i, g in enumerate(sorted(np.unique(lab))):
        np.testing.assert_allclose(pop.values.iloc[i].to_numpy(), held[lab == g].mean().to_numpy())


def test_reference_from_views_units():
    rng = np.random.default_rng(6)
    pops = [P.Population.from_criterion(pd.DataFrame(rng.normal(1, 2, size=(20, len(ODORS))), columns=ODORS), 2.0,
                                        name=f'f{i}') for i in range(4)]
    want = P.fly_reference(pops, ODORS, parts=list(E142.STRENGTH_PARTS))
    views = {p.name: P.population_views(p, ODORS, parts=list(E142.STRENGTH_PARTS)) for p in pops}
    got = S0.reference_from_views(views, ODORS, P.APL16_BREADTH_EDGES)
    for p in E142.STRENGTH_PARTS:
        assert got.one_fly[p] == pytest.approx(want.one_fly[p])
        np.testing.assert_allclose(got.mean[p].to_numpy(), want.mean[p].to_numpy())


def _planted_mega(rng, shift_arm, stages):
    w8 = list(E142.WIRINGS)
    fam = {w: {s: 1.0 + 0.02 * rng.normal() for s in stages} for w in w8}
    arm = {w: {s: v + shift_arm.get(s, 0.0) for s, v in fam[w].items()} for w in E142.T_WIRINGS}
    return S0.with_stage(fam, stages), S0.with_stage(arm, stages)


@pytest.mark.parametrize('shift', [{'pn_bouton': -0.3}, {'pn_bouton': -0.3, 'kc_claw': 0.2}, {'kc_soma': 0.01}, {}])
def test_reclass_matches_e136_classify(shift):
    from experiments import E136_stage_population_fit as E136
    rng = np.random.default_rng(7)
    stages = list(E142.LOSS_STAGES)
    fam, arm = _planted_mega(rng, shift, stages)
    spreads = dict(stage=S0._sd(v['stage'] for v in fam.values()), **{s: S0._sd(v[s] for v in fam.values())
                                                                       for s in stages})
    panels = {a: {k: {w: 1.0 for w in E142.WIRINGS} for k in E136.SCORED_PANELS} for a in ('A9', E142.FAMILY)}
    for k in E136.SCORED_PANELS:
        spreads['/'.join(k)] = 0.1
    orig = E136.classify('A9', stages, {'A9': arm, E142.FAMILY: fam}, panels, spreads, {}, {}, 0.05)
    got = S0.reclass(orig, arm, fam, spreads, stages, 0.05)
    assert got['cls'] == orig['cls']
    assert got['megamat']['change'] == pytest.approx(orig['megamat']['change'])


def test_reclass_drops_the_guard_stage():
    rng = np.random.default_rng(8)
    stages = [s for s in E142.LOSS_STAGES if s != 'pn_bouton']
    fam, arm = _planted_mega(rng, {'kc_claw': -0.4}, stages)
    spreads = dict(stage=0.01, **{s: 0.01 for s in stages})
    orig = dict(cls='x', reasons=[], guard_worse=['pn_bouton'], panel_worse=[], panel_missing=[])
    assert S0.reclass(orig, arm, fam, spreads, stages, None)['cls'] == "strength at another view's cost"
    assert S0.reclass(orig, arm, fam, spreads, stages, None, drop_guard=('pn_bouton',))['cls'] == 'holds'


def test_readings_class():
    assert S0.readings_class(-0.01, 0.02)['holds'] and S0.readings_class(0.05, 0.01)['holds']
    assert not S0.readings_class(-0.19, 0.015)['holds']


def test_counts_class():
    rows = []
    for s, c, n in (('all', 0.8, 1.0), ('fit', 0.85, 1.0), ('held', 0.95, 1.0)):
        rows.append(dict(calyx='hemibrain-first', stage='pn_bouton', glomeruli='dendrites', fly_set=s, counts=c,
                         no_step=n, null_p05=0.9, null_percentile_of_actual=1.0, below_p05=c < 0.9))
    r = S0.counts_class(pd.DataFrame(rows), 'hemibrain-first')
    assert r['margin_met'] == {'fit': True, 'held': False} and not r['holds']
    rows[2]['counts'] = 0.85
    assert S0.counts_class(pd.DataFrame(rows), 'hemibrain-first')['holds']


def test_lock_constants_and_step1_names():
    from experiments import E142_imposed_calyx_runs as S1
    assert E142.N_STEP1_RUNS == len(S1.STEP1_ARMS) * 8 == 56
    tags = {S1.tag_of(a, w) for a, w in S1.planned_step1()}
    for name in (*E142.STEP1_APL_ARMS, E142.STEP1_CARRY_ARM, *E142.STEP1_CONTROLS):
        assert all(f'{name}_{w}' in tags for w in E142.WIRINGS), name
    assert set(E142.CLASSES) >= {'counts_make_step', 'apl_hides_step', 'dendrites_carry', 'readings', 'no_class'}
    assert E142.COUNTS_MARGIN == 0.10 and E142.NULL_PERCENTILE == 5.0 and E142.N_PERM == 1000
    assert S0.PIECES == ('flies', 'models', 'yardsticks', 'sampling', 'equal', 'probe', 'classes')


def test_committed_lock_copy_is_consistent():
    lock = json.loads(E142.LOCK_COPY.read_text())
    assert lock['experiment'] == E142.EXPERIMENT
    assert lock['classes']['counts_make_step']['margin'] == E142.COUNTS_MARGIN
    assert len(lock['step1']['tags']) == E142.N_STEP1_RUNS
    assert lock['glomeruli']['carried'] and 'VA4' in lock['glomeruli']['dendrites']
    assert lock['narrow_readings'] == list(E142.NARROW_READINGS)
    assert lock['amendments'] == list(E142.AMENDMENTS)


def test_cli_parts():
    with pytest.raises(SystemExit):
        E142.main(['--help'])
    with pytest.raises(NotImplementedError):
        E142.score_step1()
