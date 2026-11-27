"""Tests of E146's pure functions (experiments/E146_bouton_step_forms.py) on planted data: the weighted entries against
E142's weighted arrays and hal's no-step score, the Jacobian and the fit's gradient against finite differences, the fit
recovering planted weights on left-out odors, the permuted predictions, the linear left-out prediction, the features,
E128's deals, the uneven draws, the power rule and the penalty's choice, the intervals between sets, the classes and the
outcome letters, and the committed lock against the module. No model is run and no recording is read."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments import E142_step0 as S0
from experiments import E146_bouton_step_forms as E
from kcq import harness as H
from kcq import population_fit as P

N = 9
ODORS = [f'o{i}' for i in range(N)]
PARTS = list(E.STRENGTH_PARTS)


def _named_pop(rng, glomeruli, name='f', crit=2.0):
    V = pd.DataFrame(rng.normal(2.0, 2.0, size=(len(glomeruli), N)), index=glomeruli, columns=ODORS)
    return P.Population.from_criterion(V, crit, name=name, named_units=True)


def _flies(rng, n_flies=4, n_gl=7, drop=True):
    gl = [f'g{i}' for i in range(n_gl)]
    pops = []
    for i in range(n_flies):
        keep = [g for g in gl if not (drop and i % 2 and g == gl[i % n_gl])]
        pops.append(_named_pop(rng, keep, name=f'd{i}'))
    return gl, pops


def test_weighted_entries_ones_are_the_no_step():
    rng = np.random.default_rng(1)
    gl, pops = _flies(rng)
    mats = S0.fly_matrices(pops, gl, ODORS)
    down = [P.Population.from_criterion(pd.DataFrame(rng.normal(1.0, 2.0, size=(30, N)), columns=ODORS), 2.0,
                                        name=f'b{i}') for i in range(3)]
    ref = P.fly_reference(down, ODORS, parts=PARTS)
    ent = E.weighted_entries(np.ones(len(gl)), mats, np.arange(N))
    got = E.score_entries({p: v for p, (v, _) in ent.items()}, ref)['view_strength']
    assert got == pytest.approx(P.no_step_score(pops, ref).views['strength'], abs=1e-12)


def test_weighted_entries_equal_e142_weighted_arrays():
    rng = np.random.default_rng(2)
    gl, pops = _flies(rng)
    mats = S0.fly_matrices(pops, gl, ODORS)
    w = rng.integers(0, 5, size=len(gl)).astype(float) + 0.5
    ent = E.weighted_entries(w, mats, np.arange(N))
    for p in PARTS:
        want = np.mean([S0.weighted_arrays(w[None, :], V, Rr)[p][0] for V, Rr in mats], axis=0)
        np.testing.assert_allclose(ent[p][0], want, atol=1e-12)


@pytest.mark.parametrize('norm', [np.arange(N), np.array([0, 2, 3, 5, 6, 8])])
def test_jacobian_against_finite_differences(norm):
    rng = np.random.default_rng(3)
    gl, pops = _flies(rng)
    prep = E.prepare(S0.fly_matrices(pops, gl, ODORS))
    w = rng.uniform(0.5, 2.0, size=len(gl))
    ent = E.weighted_entries(w, prep, norm, jacobian=True)
    h = 1e-6
    for g in range(len(gl)):
        dw = np.zeros_like(w)
        dw[g] = h
        up = E.weighted_entries(w + dw, prep, norm)
        dn = E.weighted_entries(w - dw, prep, norm)
        for p in PARTS:
            np.testing.assert_allclose(ent[p][1][g], (up[p][0] - dn[p][0]) / (2 * h), rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize('penalty', [0.0, 0.3])
def test_fit_loss_gradient(penalty):
    rng = np.random.default_rng(4)
    gl, pops = _flies(rng)
    prep = E.prepare(S0.fly_matrices(pops, gl, ODORS))
    fit = np.array([i for i in range(N) if i != 4])
    targets = {p: rng.uniform(0.5, 1.5, size=N) for p in PARTS}
    units = {'strength/rate': 0.3, 'strength/mean': 0.4}
    w = rng.uniform(0.5, 2.0, size=len(gl))
    f0, g0 = E.fit_loss(w, prep, targets, units, fit, penalty)
    h = 1e-6
    for g in range(len(gl)):
        dw = np.zeros_like(w)
        dw[g] = h
        num = (E.fit_loss(w + dw, prep, targets, units, fit, penalty)[0]
               - E.fit_loss(w - dw, prep, targets, units, fit, penalty)[0]) / (2 * h)
        assert g0[g] == pytest.approx(num, rel=1e-5, abs=1e-8)
    # the penalty is the squared coefficient of variation
    f_eq, _ = E.fit_loss(np.ones(len(gl)), prep, targets, units, fit, penalty)
    f_eq0, _ = E.fit_loss(np.ones(len(gl)), prep, targets, units, fit, 0.0)
    assert f_eq == pytest.approx(f_eq0)


def _planted(rng, n_flies=4, n_gl=7):
    gl, pops = _flies(rng, n_flies, n_gl)
    prep = E.prepare(S0.fly_matrices(pops, gl, ODORS))
    w = rng.gamma(1.5, 1.0, size=len(gl)) + 0.05
    truth = E.plant_truth(w, prep)
    parts = {p: np.repeat(truth[p][None], 3, axis=0) for p in PARTS}
    return gl, pops, prep, w, truth, parts


@pytest.mark.parametrize('fit', [np.arange(N), np.array([0, 2, 3, 5, 6, 8])])
def test_fit_entries_jacobian(fit):
    rng = np.random.default_rng(14)
    gl, pops = _flies(rng)
    prep = E.prepare(S0.fly_matrices(pops, gl, ODORS))
    w = rng.uniform(0.5, 2.0, size=len(gl))
    ent = E.fit_entries(w, prep, fit, jacobian=True)
    full = E.weighted_entries(w, prep, np.arange(N))
    h = 1e-6
    for p in PARTS:
        np.testing.assert_allclose(ent[p][0], full[p][0] / full[p][0][fit].mean(), atol=1e-12)
        for g in range(len(gl)):
            dw = np.zeros_like(w)
            dw[g] = h
            num = (E.fit_entries(w + dw, prep, fit)[p][0] - E.fit_entries(w - dw, prep, fit)[p][0]) / (2 * h)
            np.testing.assert_allclose(ent[p][1][g], num, rtol=1e-5, atol=1e-7)


def test_planted_weights_reproduce_the_truth_on_the_fit_odors():
    # the lock's amendment of 2026-11-26 20:28: a plant's own weights give a zero fit loss on every fold
    rng = np.random.default_rng(15)
    gl, pops, prep, w, truth, parts = _planted(rng)
    units = {p: 0.3 for p in PARTS}
    for o in range(N):
        fit = np.array([i for i in range(N) if i != o])
        assert E.fit_loss(w, prep, E.relative_parts(parts, fit), units, fit, 0.0)[0] < 1e-24


def test_fit_recovers_planted_weights_on_left_out_odors():
    rng = np.random.default_rng(5)
    gl, pops, prep, w, truth, parts = _planted(rng, n_gl=4)
    units = {p: 0.3 for p in PARTS}
    pred, W, recs = E.left_out(prep, parts, units, 0.0)
    ref = E.synthetic_reference(parts, units, ODORS)
    lo = E.score_entries(pred, ref)['view_strength']
    ns = E.score_entries({p: v for p, (v, _) in E.weighted_entries(np.ones(len(gl)), prep, np.arange(N)).items()},
                         ref)['view_strength']
    assert ns > 0.05 and lo < 0.01 * ns
    np.testing.assert_allclose(W, np.repeat((w / w.mean())[None], N, axis=0), rtol=1e-3)
    assert W.shape == (N, len(gl)) and np.allclose(W.mean(axis=1), 1.0)
    # a strong penalty pulls the weights toward equal and the reading toward the no-step
    pred_h, W_h, _ = E.left_out(prep, parts, units, 100.0)
    assert np.max(np.abs(W_h - 1.0)) < 0.2
    assert E.score_entries(pred_h, ref)['view_strength'] > lo


def test_predictions_of_matches_left_out_and_batches():
    rng = np.random.default_rng(6)
    gl, pops, prep, w, truth, parts = _planted(rng)
    units = {p: 0.3 for p in PARTS}
    pred, W, _ = E.left_out(prep, parts, units, 0.1)
    again = E.predictions_of(W, prep)
    for p in PARTS:
        np.testing.assert_allclose(again[p], pred[p], atol=1e-12)
    perms = S0.permutations(len(gl), 5, 0)
    batch = E.predictions_of(np.transpose(W[:, perms], (1, 0, 2)), prep)
    for i, pm in enumerate(perms):
        one = E.predictions_of(W[:, pm], prep)
        for p in PARTS:
            np.testing.assert_allclose(batch[p][i], one[p], atol=1e-12)


def test_relative_parts():
    A = {'x': np.array([[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])}
    r = E.relative_parts(A, [0, 1])['x']
    np.testing.assert_allclose(r, np.mean([[1 / 1.5, 2 / 1.5, 3 / 1.5], [1.0, 1.0, 1.0]], axis=0))


def test_linear_left_out_recovers_an_exact_linear_step():
    rng = np.random.default_rng(7)
    x = rng.uniform(0.5, 1.5, size=N)
    f = rng.uniform(0.0, 1.0, size=N)
    raw = 0.4 + 0.8 * x + 0.5 * f
    T = raw / raw.mean()
    flies = np.repeat(T[None], 3, axis=0)
    np.testing.assert_allclose(E.linear_left_out(x, f, flies), T, atol=1e-10)
    # without the feature the fit is not exact
    assert np.max(np.abs(E.linear_left_out(x, None, flies) - T)) > 1e-3


def test_feature_prediction_both_ways_is_the_mean_of_the_two():
    rng = np.random.default_rng(8)
    x = {s: {p: rng.uniform(0.5, 1.5, N) for p in PARTS} for s in E.FLY_SETS}
    halves = {'other_part': {h: {p: x[h][E.OTHER_PART[p]] for p in PARTS} for h in E.HALVES}}
    single = {'orn_summed': {p: rng.uniform(1, 5, N) for p in PARTS}}
    I = E.FeatureInputs(x, single, halves)
    assert I.features() == ['orn_summed', 'other_part']
    fp = {p: rng.uniform(0.5, 1.5, (3, N)) for p in PARTS}
    got = E.feature_prediction(I, 'other_part', fp)
    for p in PARTS:
        a = E.linear_left_out(x['fit'][p], x['held'][E.OTHER_PART[p]], fp[p])
        b = E.linear_left_out(x['held'][p], x['fit'][E.OTHER_PART[p]], fp[p])
        np.testing.assert_allclose(got[p], (a + b) / 2)
    pm = np.roll(np.arange(N), 2)
    got = E.feature_prediction(I, 'orn_summed', fp, perm=pm)
    np.testing.assert_allclose(got['strength/rate'], E.linear_left_out(x['all']['strength/rate'],
                                                                       single['orn_summed']['strength/rate'][pm],
                                                                       fp['strength/rate']))


def test_top_share_and_summed_positive():
    V = pd.DataFrame({'a': [5.0, 3.0, 1.0, 1.0, -2.0], 'b': [0.0, 0.0, 0.0, 0.0, 0.0]})
    t = E.top_share(V, 3)
    assert t['a'] == pytest.approx(9.0 / 10.0) and np.isnan(t['b'])
    s = E.summed_positive(V)
    assert s['a'] == pytest.approx(10.0) and s['b'] == 0.0


def test_flies_mean_table():
    a = P.Population.from_criterion(pd.DataFrame([[1.0, 2.0]], index=['g1'], columns=['x', 'y']), 2.0, name='a',
                                    named_units=True)
    b = P.Population.from_criterion(pd.DataFrame([[3.0, 4.0], [5.0, 6.0]], index=['g1', 'g2'], columns=['x', 'y']),
                                    2.0, name='b', named_units=True)
    t = E.flies_mean_table([a, b], ['x', 'y'])
    np.testing.assert_allclose(t.loc['g1'].to_numpy(), [2.0, 3.0])
    np.testing.assert_allclose(t.loc['g2'].to_numpy(), [5.0, 6.0])


def test_deal_weights_are_e128s():
    from experiments import E128_bouton_groups as E128
    sizes = np.array([5, 1, 12, 3, 7, 2])
    w = E.deal_weights(sizes, 10, np.random.default_rng(9))
    assert w.sum() == sizes.sum() and len(w) == 10
    U = np.arange(30, dtype=float)[:, None] * np.ones((1, 3))
    gl = np.repeat([f'g{i}' for i in range(10)], 3)
    gloms = np.unique(gl)
    _, wts = E128.deal_population(U, gl, sizes, np.random.default_rng(9), gloms)
    got = {str(g): int(k) for g, k in zip(gloms, w) if k}
    assert got == wts


@pytest.mark.parametrize('cv', [0.0, 0.3, 1.2])
def test_uneven_weights(cv):
    z = np.random.default_rng(10).standard_normal(29)
    w = E.uneven_weights(cv, z)
    assert E.cv_of(w) == pytest.approx(cv, abs=1e-9)
    assert np.all(w > 0)


def test_noise_fly_parts_and_synthetic_reference():
    rng = np.random.default_rng(11)
    truth = {p: rng.uniform(0.5, 1.5, N) for p in PARTS}
    bp = {p: rng.uniform(0.5, 1.5, (5, N)) for p in PARTS}
    out = E.noise_fly_parts(truth, bp, np.random.default_rng(12))
    idx = np.random.default_rng(12).integers(0, 5, size=5)
    for p in PARTS:
        np.testing.assert_allclose(out[p], truth[p][None] + (bp[p] - bp[p].mean(axis=0))[idx])
    ref = E.synthetic_reference(out, {p: 0.25 for p in PARTS}, ODORS)
    assert ref.one_fly == {p: 0.25 for p in PARTS}
    np.testing.assert_allclose(ref.mean['strength/rate'].to_numpy(), out['strength/rate'].mean(axis=0))


def _rows(kind, no_step, left_out):
    return pd.DataFrame(dict(kind=kind, no_step=no_step, left_out=left_out))


def test_power_rule_and_choose_penalty():
    ok = pd.concat([_rows('hemibrain', [0.5] * 20, [0.1] * 20),                       # vacuous
                    _rows('deals', [1.0] * 10 + [0.5] * 10, [0.85] * 10 + [0.4] * 10),  # 10 eligible, all pass
                    _rows('uneven', [1.2] * 10, [1.0] * 9 + [1.15])])                   # 9 of 10
    r = E.power_rule(ok)
    assert r['vacuous'] == ['hemibrain'] and r['n_eligible'] == 20 and r['holds']
    assert r['kinds']['uneven']['share'] == pytest.approx(0.9)
    bad = pd.concat([_rows('deals', [1.0] * 10, [0.95] * 10), _rows('uneven', [1.0] * 10, [0.5] * 10)])
    rb = E.power_rule(bad)
    assert not rb['holds'] and rb['min_share'] == 0.0
    few = E.power_rule(_rows('deals', [1.0] * 5, [0.5] * 5))
    assert not few['holds']   # fewer than POWER_MIN_ELIGIBLE eligible plants
    by = {0.0: rb, 0.1: r, 1.0: dict(r), 3.0: few}
    assert E.choose_penalty(by) == 1.0   # ties to the larger penalty; too few eligible never chosen


def test_resampled_and_measured_errors():
    rng = np.random.default_rng(13)
    a = {p: rng.uniform(0.5, 1.5, (4, N)) for p in PARTS}
    b = {p: rng.uniform(0.5, 1.5, (6, N)) for p in PARTS}
    pops = [P.Population.from_criterion(pd.DataFrame(rng.normal(1, 2, size=(20, N)), columns=ODORS), 2.0,
                                        name=f'r{i}') for i in range(6)]
    ref = P.fly_reference(pops, ODORS, parts=PARTS)
    units = {p: float(ref.one_fly[p]) for p in PARTS}
    m = E.measured_error(a, b, units)
    ref_b = E.synthetic_reference(b, units, ODORS)
    want = P.score_views({p: pd.Series(a[p].mean(axis=0), index=ODORS) for p in PARTS}, ref_b).views['strength']
    assert m['view_strength'] == pytest.approx(want, abs=1e-12)
    d1 = E.resampled_errors(a, b, units, 50, 0)
    d2 = E.resampled_errors(a, b, units, 50, 0)
    assert d1.shape == (50,) and np.array_equal(d1, d2)
    iv = E.interval(d1)
    assert iv[0] <= np.median(d1) <= iv[1]


def test_overlap_and_outcome():
    assert E.overlap([0.9, 1.2], [1.1, 1.5]) and not E.overlap([0.9, 1.0], [1.1, 1.5])
    assert E.outcome(True, True, False, True) == ['a']
    assert E.outcome(True, True, True, False) == ['a', 'b']
    assert E.outcome(False, True, True, False) == ['b']
    assert E.outcome(True, False, False, True) == ['c']
    assert E.outcome(True, False, False, False) == ['d']
    assert E.outcome(False, False, False, True) == ['e']


def test_weight_and_feature_classes():
    bound = dict(left_out=dict(all=0.8, fit=1.0, held=0.95), no_step=dict(all=1.052, fit=1.184, held=1.084),
                 null=dict(all=dict(null_p05=0.9)), penalty=0.1)
    c = E.weight_class(bound, True)
    assert c['holds'] and c['cls'] == 'makes the step'
    assert E.weight_class(bound, False)['cls'] == 'not readable' and not E.weight_class(bound, False)['holds']
    bound['left_out']['held'] = 1.0   # 0.084 below: the margin fails on the held half
    assert not E.weight_class(bound, True)['holds']
    forms = {f: dict(left_out=dict(all=1.0, fit=1.15, held=1.05), no_step=dict(all=1.052, fit=1.184, held=1.084))
             for f in E.FEATURES}
    forms['orn_summed'] = dict(left_out=dict(all=0.7, fit=1.0, held=0.9), no_step=forms['other_part']['no_step'])
    fc = E.feature_class(dict(forms=forms, null_best_p05=0.8))
    assert fc['holds'] and fc['passing'] == ['orn_summed'] and fc['n_forms'] == 4
    assert not E.feature_class(dict(forms=forms, null_best_p05=0.6))['holds']


def test_remap_path():
    assert E.remap_path('/home/yang/src/kc_queue/model_experiment_queue/experiments/E134_lock.json') == \
        H.HERE / 'experiments' / 'E134_lock.json'
    assert E.remap_path('repo:experiments/x.json') == H.HERE / 'experiments' / 'x.json'
    try:
        sd = H.shared_dir()
    except FileNotFoundError:
        pytest.skip('no shared directory')
    assert E.remap_path('/local/tensor/Tom/claude_kc_remapping/results/a.csv') == sd / 'results' / 'a.csv'
    assert E.remap_path('shared:results/a.csv') == sd / 'results' / 'a.csv'


def test_committed_lock_matches_the_module():
    if not E.LOCK_COPY.exists():
        pytest.skip('no lock yet')
    lock = json.loads(E.LOCK_COPY.read_text())
    c = lock['constants']
    assert c['penalties'] == list(E.PENALTIES) and c['n_perm'] == E.N_PERM and c['margin'] == E.MARGIN
    assert c['power_eligible'] == E.POWER_ELIGIBLE and c['power_share'] == E.POWER_SHARE
    assert c['n_hemibrain_plants'] == E.N_HEMIBRAIN_PLANTS and c['n_deals_per_fly'] == E.N_DEALS_PER_FLY
    assert c['n_uneven'] == E.N_UNEVEN and c['features'] == list(E.FEATURES)
    assert c['seven_odors'] == list(E.SEVEN_ODORS) and c['sets_q'] == list(E.SETS_Q)
    assert lock['classes'] == json.loads(json.dumps(E.CLASSES)) and lock['outcomes'] == E.OUTCOMES
    assert lock['choices'] == list(E.CHOICES) and lock.get('amendments', []) == list(E.AMENDMENTS)
    assert len(lock['glomeruli']['dendrites']) == 29 and len(lock['odors']) == E.N_ODORS
    assert not any(f.startswith('2022-10-21__fly02') for f in lock['recordings']['pn_bouton']['flies'])


def test_cli_parts(monkeypatch):
    called = []
    for p in E.PARTS:
        monkeypatch.setattr(E, f'part_{p}', (lambda q: (lambda *a, **k: called.append(q)))(p))
    E.main(['--part', 'sets', 'identities'])
    assert called == ['identities', 'sets']
