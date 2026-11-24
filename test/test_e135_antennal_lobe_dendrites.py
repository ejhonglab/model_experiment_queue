"""Tests of experiments/E135_antennal_lobe_dendrites.py on planted data: the measured step (gains, r, residuals,
candidate partners, intervals, halves), the lock's folds, the gate's refusal, the numpy port of olfsysm's ORN, LN and PN
layers against olfsysm itself (on planted parameters and input, no fit_mb_model), the glomerulus terms, the bounds,
route 1 (olfsysm's PN rates given directly) carrying the port's terms, route 2 (the solved input) and its identity, the
imposed table and its bounds policy, the rung fits and the cross-fit's assembly on planted dendrites, and the classes.
No recording and no model of the experiment is read."""
import json

import numpy as np
import pandas as pd
import pytest

import olfsysm as osm
from data_prep import stage_population_data as D
from experiments import E135_antennal_lobe_dendrites as E

pytestmark = pytest.mark.filterwarnings('ignore')

SHORT = dict(time_pre_start=-0.6, time_start=-0.3, time_end=0.7)


def _mp(G=10, O=4, seed=1, zero_rows=2, n_phys=12):
    rng = np.random.default_rng(seed)
    mp = osm.ModelParams()
    for k, v in SHORT.items():
        setattr(mp, k, v)
    spont = rng.uniform(5, 30, (G, 1))
    delta = rng.uniform(-20, 150, (G, O))
    if zero_rows:
        delta[-zero_rows:] = 0
    mp.orn.data.spont = spont
    mp.orn.data.delta = delta
    mp.orn.n_physical_gloms = n_phys
    mp.kc.N = 10
    rv = osm.RunVars(mp)
    return mp, rv


def _stage(mp):
    G, O = np.asarray(mp.orn.data.delta).shape
    return E.PNStage.from_mp(mp, [f'g{i}' for i in range(G)], [f'o{j}' for j in range(O)])


@pytest.fixture(scope='module')
def planted():
    mp, rv = _mp()
    return mp, rv, _stage(mp)


# ---------------------------------------------------------------------------
# step 0's readings

def _flies(gains, n_orn=6, n_pn=5, n_odors=12, seed=0, noise=0.02):
    rng = np.random.default_rng(seed)
    G = len(gains)
    x = rng.gamma(1.5, 1.0, (G, n_odors))
    orn = x[None] + rng.normal(scale=noise, size=(n_orn, G, n_odors))
    pn = (np.asarray(gains)[:, None] * x)[None] + rng.normal(scale=noise, size=(n_pn, G, n_odors))
    return orn, pn, x


def test_measured_step_recovers_planted_gains():
    gains = np.array([0.5, 1.0, 1.0, 2.0, 4.0])
    orn, pn, _ = _flies(gains)
    gl, od = [f'g{i}' for i in range(5)], [f'o{j}' for j in range(12)]
    st = E.measured_step(orn, pn, gl, od, total=np.nansum(orn.mean(axis=0), axis=0))
    np.testing.assert_allclose(st.gain.to_numpy(), gains / np.median(gains), rtol=0.03)
    assert (st.r > 0.98).all()
    assert {'own_coef', 'total_coef', 'r2_own_total'} <= set(st.columns)
    iv, ev = E.step_intervals(orn, pn, gl, od, n_boot=200)
    assert ((iv.gain_lo <= st.gain) & (st.gain <= iv.gain_hi)).all()
    assert ((iv.r_lo <= st.r) & (st.r <= iv.r_hi)).all()
    assert len(ev) == 5 * 12


def test_candidate_partner_found_and_weighted():
    rng = np.random.default_rng(3)
    O = 15
    X = pd.DataFrame(rng.gamma(1.5, 1.0, (6, O)), index=[f'p{i}' for i in range(6)])
    own = X.loc[['p0']].rename(index={'p0': 'g'})
    den = 1.5 * own + 0.6 * X.loc[['p3']].to_numpy()
    st = E.measured_step(own.to_numpy()[None], den.to_numpy()[None], ['g'], list(range(O)))
    res = E.residuals(own.to_numpy()[None], den.to_numpy()[None], st, ['g'], list(range(O)))
    res.columns = X.columns
    cp = E.candidate_partners(res, X.drop(index='p0'), own.set_axis(X.columns, axis=1),
                              den.set_axis(X.columns, axis=1), n=2, n_perm=200)
    top = cp[cp['rank'] == 1].iloc[0]
    assert top.partner == 'p3'
    assert top.r2 > top.chance_r2_95
    assert top.weight_rel == pytest.approx(0.6 / 1.5, rel=1e-6)


def test_halves_and_folds_and_follows_rule():
    flies = [f'f{i}' for i in range(7)]
    a, b = E.halves(flies)
    assert len(a) == 3 and len(b) == 4 and set(a) | set(b) == set(flies) and E.halves(flies) == (a, b)
    folds = E.odor_folds(5)
    assert sorted(o for f in folds for o in f) == sorted(E.MEGAMAT_ODORS) and len(folds) == 5
    assert len(E.odor_folds(17)) == 17
    with pytest.raises(ValueError):
        E.write_lock(None)
    st = pd.DataFrame(dict(r_lo=[0.2, -0.1], gain=[1.5, 0.3]), index=['a', 'b'])
    assert E.follows_rule(st).tolist() == [True, False]


def test_gate_refuses_until_e134_records_usable(tmp_path):
    with pytest.raises(PermissionError, match='measure_gate.json'):
        E.require_e134_gate(tmp_path)
    (tmp_path / E.E134_GATE_FILE).write_text(json.dumps(dict(verdict='not usable')))
    (tmp_path / E.E134_SCORED_FILE).write_text('{}')
    with pytest.raises(PermissionError, match="'not usable'"):
        E.require_e134_gate(tmp_path)
    (tmp_path / E.E134_GATE_FILE).write_text(json.dumps(dict(verdict='usable')))
    assert E.require_e134_gate(tmp_path)['verdict'] == 'usable'
    (tmp_path / E.E134_SCORED_FILE).unlink()
    with pytest.raises(PermissionError, match='absent'):
        E.require_e134_gate(tmp_path)


# ---------------------------------------------------------------------------
# the port against olfsysm

def test_port_matches_olfsysm(planted):
    mp, rv, st = planted
    r0 = st.check_identity(mp, rv)
    assert r0['ok'] and r0['max_abs_diff']['pn'] < 1e-9 and r0['max_abs_diff']['evoked_preset'] < 1e-9
    r1 = st.check_identity(mp, rv, params=dict(ln_inhsc=900.0, ln_inhadd=120.0, pn_inhadd=300.0, pn_tanhsc=3.0,
                                               pn_inhsc=500.0))
    assert r1['ok']


def test_bounds_match_olfsysm(planted):
    mp, rv, st = planted
    osm.run_ORN_LN_sims(mp, rv)
    osm.run_PN_sims(mp, rv)
    lo, hi = osm.pn_rate_bounds(mp, rv)
    plo, phi = st.rate_bounds()
    np.testing.assert_allclose(plo.to_numpy(), lo, atol=1e-9)
    np.testing.assert_allclose(phi.to_numpy(), hi, atol=1e-9)


def test_terms_act_on_the_pn_drive_only(planted):
    mp, rv, st = planted
    base = st.run()
    t = E.Terms(gains={'g1': 2.0}, laterals=[('g2', 'g4', 0.5)])
    x = st.run(terms=t)
    other = [g for g in st.glomeruli if g not in ('g1', 'g2')]
    pd.testing.assert_frame_equal(x.loc[other], base.loc[other])     # the LN keeps the recorded input
    assert (x.loc['g1'] > base.loc['g1']).all()
    d2 = st.delta.copy()
    d2.loc['g1'] *= 2
    y = st.run(delta=d2)                                              # the input changed: the LN reads it too
    assert not np.allclose(y.loc[other].to_numpy(), base.loc[other].to_numpy())
    M = t.matrix(st.glomeruli)
    v = st.run_virtual(M[:, [1, 2]], np.array([1, 2]))
    np.testing.assert_allclose(v, x.loc[['g1', 'g2']].to_numpy(), atol=1e-12)
    with pytest.raises(ValueError):
        E.Terms(laterals=[('g1', 'g2', 0.1), ('g1', 'g3', 0.1), ('g1', 'g4', 0.1)]).matrix(st.glomeruli)


def test_route1_carries_the_ports_terms(planted):
    """olfsysm's PN rates given directly keep the LN and shift a PN's drive by a constant times the unit ORN kernel:
    a gain or a lateral term of the port is such a shift, so route 1 given the port's rates returns the port's
    traces."""
    mp, rv, st = planted
    t = E.Terms(gains={'g0': 1.8, 'g3': 0.5}, laterals=[('g5', 'g6', -0.4)])
    target = st.run(terms=t, reading='preset')
    tab = E.route1_pn_rates(target, ['g0', 'g3', 'g5'])
    assert np.isnan(tab.loc['g1']).all() and np.isfinite(tab.loc['g0']).all()
    try:
        rv.pn.preset_rates = tab.to_numpy()
        osm.run_ORN_LN_sims(mp, rv)
        osm.run_PN_sims(mp, rv)
        got = np.asarray(osm.pn_evoked_rates(mp, rv))
        traces = np.stack([np.asarray(x) for x in rv.pn.pn_sims])
    finally:
        rv.pn.preset_rates = np.zeros((0, 0))
    m = np.isfinite(tab.to_numpy())
    np.testing.assert_allclose(got[m], tab.to_numpy()[m], atol=1e-8)
    sp = st.sp
    port = E.pn_layer(sp, st.spont, st.orn(), *st.ln(sp), t.matrix(st.glomeruli), np.arange(len(st.glomeruli)),
                      full=True)
    np.testing.assert_allclose(traces, port, atol=1e-6)


def test_solve_identity_and_a_new_target(planted):
    mp, rv, st = planted
    ident = E.solve_identity(st, glomeruli=st.glomeruli[:8])
    assert ident['ok'] and ident['max_delta_diff'] < 1e-3
    d2 = st.delta.copy()
    d2.loc[['g0', 'g1', 'g2']] = d2.loc[['g0', 'g1', 'g2']].to_numpy() * np.array([[0.5], [1.7], [1.2]])
    st2 = E.PNStage(d2, pd.Series(st.spont, index=st.glomeruli), st.sp)
    target = st2.run(reading='preset')
    solved, rec = E.solve_input(st, target.loc[st.glomeruli[:8]], tol=0.01)
    assert rec['converged'] and rec['fixed_point'] and rec['max_abs_diff'] <= 0.01
    np.testing.assert_allclose(solved.loc[st.glomeruli[:8]].to_numpy(), d2.loc[st.glomeruli[:8]].to_numpy(), atol=0.5)
    pd.testing.assert_frame_equal(solved.loc[st.glomeruli[8:]], st.delta.loc[st.glomeruli[8:]])


def test_imposed_table_rank_map_and_bounds(planted):
    mp, rv, st = planted
    own = st.run(reading='preset')
    gl = st.glomeruli[:8]
    rng = np.random.default_rng(0)
    den = pd.DataFrame(rng.normal(size=(8, len(st.odors))), index=gl, columns=st.odors)
    tab, rec = E.imposed_rates(own, den, gl)
    np.testing.assert_allclose(np.sort(tab.loc[gl].to_numpy().ravel()), np.sort(own.loc[gl].to_numpy().ravel()),
                               atol=1e-9)
    assert (np.argsort(tab.loc[gl].to_numpy().ravel()) == np.argsort(den.to_numpy().ravel())).all()
    pd.testing.assert_frame_equal(tab.loc[st.glomeruli[8:]], own.loc[st.glomeruli[8:]])
    lo, hi = st.rate_bounds()
    tight = (lo, hi.clip(upper=own.max().max() * 0.5))
    with pytest.raises(ValueError, match='outside'):
        E.imposed_rates(own, den, gl, bounds=tight)
    tab2, rec2 = E.imposed_rates(own, den, gl, bounds=tight, out_of_bounds='clip')
    assert rec2['n_outside'] > 0
    assert (tab2.loc[gl].to_numpy() < tight[1].loc[gl].to_numpy()).all()


# ---------------------------------------------------------------------------
# the fits on planted dendrites

def _pn_population(rates: pd.DataFrame, n_flies=4, seed=0, scale=0.05, offset=0.2, noise=0.15):
    """Planted PN flies: the rates through an affine map, plus trial noise, 3 trials, in viki's layout."""
    rng = np.random.default_rng(seed)
    gl, od = list(rates.index), list(rates.columns)
    rows = []
    for f in range(n_flies):
        X = (scale * rates.to_numpy() + offset)[None] + rng.normal(scale=noise, size=(3, len(gl), len(od)))
        rows.append(D.long_from_trials(f'f{f}', gl, od, X))
    d = pd.concat(rows, ignore_index=True)
    spec = dict(units=D.UNITS_TRIAL, criterion=2.0, window={}, notes=[], recording_site='dendrites')
    pop = D.population_from_long(('pn_dendrite', 'megamat', 'wild_type'), d.assign(
        odor=d.odor.map(dict(zip(od, D.MEGAMAT_ODORS)))), spec=spec, trials=True)
    return pop, dict(zip(od, D.MEGAMAT_ODORS))


@pytest.fixture(scope='module')
def fit_case():
    mp, rv = _mp(G=8, O=10, seed=5, zero_rows=0, n_phys=10)
    st = _stage(mp)
    true = E.Terms(gains={'g0': 2.5, 'g2': 0.5, 'g5': 1.6})
    rates = st.run(terms=true)
    pop, names = _pn_population(rates)
    st = E.PNStage(st.delta.rename(columns=names), pd.Series(st.spont, index=st.glomeruli), st.sp)
    return st, pop, true, list(names.values())


def test_fit_gains_recovers_planted_gains(fit_case):
    st, pop, true, odors = fit_case
    tgt = E.Target.make(pop, odors, st.glomeruli)
    start = {g: 1.0 for g in st.glomeruli}
    loss0 = E.matrix_loss(st.run(), tgt)
    f = E.fit_gains(st, tgt, {}, start, n_grid=31)
    assert f['loss'] < loss0
    g = f['terms'].gains
    assert g['g0'] > 1.5 and g['g2'] < 0.8
    sc = E.read_prediction(st.run(terms=f['terms']), tgt)
    sc0 = E.read_prediction(st.run(), tgt)
    assert sc.parts['matrix/values'] < sc0.parts['matrix/values']
    lat = E.fit_lateral(st, tgt, {}, f['terms'], 'g1', 'g3', n_grid=21)
    assert E.LATERAL_BOUNDS[0] <= lat['weight'] <= E.LATERAL_BOUNDS[1]
    assert lat['loss'] <= f['loss'] + 1e-6


def test_fit_global_moves_from_the_defaults(fit_case):
    st, pop, true, odors = fit_case
    tgt = E.Target.make(pop, odors, st.glomeruli)
    r = E.fit_global(st, tgt, names=('pn_inhadd', 'pn_tanhsc'), maxiter=40)
    assert r['loss'] <= E.matrix_loss(st.run(), tgt) + 1e-12
    assert set(r['params']) == {'pn_inhadd', 'pn_tanhsc'}


def test_cross_fit_assembles_held_out_columns(fit_case):
    st, pop, true, odors = fit_case
    A, B = E.halves([str(f) for f in pop.flies])
    tg = {h: (lambda train, fl=fl: E.Target.make(E.flies_subset(pop, fl), train, st.glomeruli))
          for h, fl in (('half_a', A), ('half_b', B))}
    rt = {h: E.Target.make(E.flies_subset(pop, fl), odors, st.glomeruli) for h, fl in (('half_a', A), ('half_b', B))}
    seen = []

    def fit_one(target, train):
        seen.append(tuple(train))
        return st.run(terms=true), dict(n_train=len(train))
    folds = [odors[:5], odors[5:]]
    out = E.cross_fit(fit_one, folds, tg, rt, n_boot=30)
    assert len(seen) == 4 and all(len(s) == 5 for s in seen)
    k = 'half_a->half_b'
    assert list(out[k]['table'].columns) == odors and len(out[k]['draws']) == 30
    assert np.isfinite(out['point'])
    base = E.held_out_draws(st.run(), rt['half_b'], 30)
    kept, se = E.keep_rung(out[k]['draws'], base, out[k]['point'], float(np.mean(base)))
    assert se >= 0 and isinstance(kept, bool)


def test_matrix_loss_is_the_affine_reading(fit_case):
    from kcq import stage_readings as R
    st, pop, true, odors = fit_case
    tgt = E.Target.make(pop, odors, st.glomeruli)
    U = st.run().loc[st.glomeruli, odors]
    fv, fw = R.pooled_fly_values(pop, odors)
    c = R.fit_stage_curve(U.to_numpy().ravel(), fv, fw, 'affine')
    m = tgt.ref.mean['matrix/values'].reindex([f'{g}|{o}' for g in st.glomeruli for o in odors]).to_numpy()
    want = np.sqrt(np.nanmean((c(U.to_numpy()).ravel() - m) ** 2)) / tgt.ref.one_fly['matrix/values']
    assert E.matrix_loss(U, tgt) == pytest.approx(want, rel=1e-12)


def test_run_rungs_end_to_end(fit_case):
    st, pop, true, odors = fit_case
    rng = np.random.default_rng(2)
    orn = (st.delta.loc[st.glomeruli, odors].to_numpy()[None] * 0.02
           + rng.normal(scale=0.01, size=(5, len(st.glomeruli), len(odors))))
    flies = [str(f) for f in pop.flies]
    lock = dict(glomeruli=st.glomeruli, folds=[odors[:5], odors[5:]], halves=dict(half_a=flies[:2], half_b=flies[2:]))
    out = E.run_rungs(st, pop, orn, orn, st.glomeruli, lock, n_boot=10, n_perm=1, max_lateral_steps=1,
                      global_maxiter=12, verbose=False)
    r = out['rungs'].set_index('rung')
    assert list(r.index) == ['R0', 'R1', 'R2m', 'R2f', 'R3']
    assert r.loc['R2f', 'point'] < r.loc['R0', 'point']          # the planted gains are found out of sample
    assert out['r2'] in ('R2m', 'R2f') and out['highest_kept'] in ('R0', 'R1', out['r2'], 'R3')
    assert {'R1', 'R2m', 'R2f'} <= set(out['controls'].columns)
    assert set(out['tables']['R2f']) == {'half_a->half_b', 'half_b->half_a'}
    assert all(len(v) == 2 for v in out['folds']['R1'].values())


def test_permutations_and_measured_start():
    rng = np.random.default_rng(0)
    x = pd.DataFrame(np.arange(12.0).reshape(4, 3), index=list('abcd'))
    pop, _ = _pn_population(pd.DataFrame(rng.normal(size=(3, 4)), index=['a', 'b', 'c'],
                                         columns=['o0', 'o1', 'o2', 'o3']))
    q = E.permute_odor_labels(pop, np.random.default_rng(1), ['2h', 'IaA', 'pa'])
    f = pop.flies[0]
    assert q.tables[f]['2-but'].equals(pop.tables[f]['2-but'])
    assert sorted(q.tables[f][['2h', 'IaA', 'pa']].to_numpy().ravel()) == sorted(
        pop.tables[f][['2h', 'IaA', 'pa']].to_numpy().ravel())
    np.testing.assert_allclose(q.trials[f].T.groupby(level=0).mean().T[['2h', 'IaA', 'pa']].to_numpy(),
                               q.tables[f][['2h', 'IaA', 'pa']].to_numpy())
    p = E.permute_glomerulus_labels(x, ['a', 'b', 'c'], rng)
    assert sorted(p.loc[['a', 'b', 'c']].to_numpy().ravel()) == sorted(x.loc[['a', 'b', 'c']].to_numpy().ravel())
    assert p.loc['d'].equals(x.loc['d'])
    t = E.permute_terms(E.Terms(gains={'a': 2.0, 'b': 0.5}), ['a', 'b', 'c'], rng)
    assert sorted(t.gains.values()) == [0.5, 2.0]
    meas = pd.DataFrame(dict(gain=[2.0, 0.1, 1.0], follows_provisional=[True, True, False]), index=['a', 'b', 'c'])
    s = E.measured_gain_start(meas, pd.Series([0.5, 1.0, 1.0], index=['a', 'b', 'c']))
    assert s == {'a': 4.0, 'b': 0.25, 'c': 1.0}


# ---------------------------------------------------------------------------
# step 2's points and the classes

def test_planned_points():
    pts = E.planned_step2()
    assert len(pts) == 73 and len(E.planned_step2(with_rung=False)) == 49
    assert len({E.tag_of(a, w) for a, w in pts}) == 73
    assert E.planned_step2(reverse=True)[0] == pts[-1]


def _scores(soma_drop, bouton_ratio, perm_drop=0.0, r1_share=0.0, rung_share=0.0):
    rows = []
    for w in E.WIRINGS:
        for arm, sd, br in (('imposed', soma_drop, bouton_ratio), ('imposed_glom_perm', perm_drop, 1.0),
                            ('imposed_odor_perm', perm_drop, 1.0), ('R1', r1_share * soma_drop, 1.0),
                            ('rung', rung_share * soma_drop, 1.0), ('rung_glom_perm', 0.0, 1.0)):
            rows.append(dict(arm=arm, wiring=w, stage='kc_soma', error=1.0 - sd, error_base=1.0, strength=1.0,
                             strength_base=1.0))
            bdrop = 0.6 if arm == 'imposed' else (r1_share * 0.6 if arm == 'R1' else
                                                  rung_share * 0.6 if arm == 'rung' else 0.0)
            rows.append(dict(arm=arm, wiring=w, stage='pn_bouton', error=1.0 - bdrop, error_base=1.0,
                             strength=br * 2.0, strength_base=2.0))
    return pd.DataFrame(rows)


def test_classes():
    c = E.classify(_scores(0.3, 0.4))
    assert c['upstream'] and not c['boutons_only'] and not c['global_carry'] and not c['terms_carry']
    c = E.classify(_scores(0.3, 0.4, r1_share=0.6))
    assert c['upstream'] and c['global_carry']
    c = E.classify(_scores(0.3, 0.4, rung_share=0.7))
    assert c['terms_carry'] and not c['global_carry']
    c = E.classify(_scores(0.02, 0.4))
    assert c['boutons_only'] and not c['upstream']
    c = E.classify(_scores(0.3, 0.9))
    assert c['not_upstream']
    c = E.classify(_scores(0.3, 0.4, perm_drop=0.35))
    assert not c['upstream']
    thr = E.class_thresholds(pd.DataFrame(dict(stage='kc_soma', wiring=list(E.WIRINGS),
                                               error=np.linspace(1.0, 1.07, 8))))
    assert thr['kc_soma']['threshold'] == pytest.approx(3 * np.std(np.linspace(1.0, 1.07, 8), ddof=1) / np.sqrt(8))
