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
    assert len(E.odor_folds(17)) == 17 and all(len(f) == 1 for f in E.odor_folds(E.N_FOLDS))
    assert E.odor_folds(6) == E.odor_folds(6) and len(E.odor_folds(E.N_FOLDS_FALLBACK)) == 6
    # ruling 2: the slope's interval, not r's
    st = pd.DataFrame(dict(r_lo=[0.2, 0.3, -0.1], slope_lo=[0.1, -0.05, 0.2], gain=[1.5, 0.3, 1.0]),
                      index=['a', 'b', 'c'])
    assert E.follows_rule(st).tolist() == [True, False, True]


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
    out = E.run_rungs(st, pop, orn, orn, st.glomeruli, lock, n_boot=10, n_perm_lateral=2, n_perm_control=1,
                      max_lateral_terms=1, check_fifth=True, global_maxiter=12, verbose=False)
    r = out['rungs'].set_index('rung')
    assert list(r.index) == ['R0', 'R1', 'R2m', 'R2f', 'R3']
    assert r.loc['R2f', 'point'] < r.loc['R0', 'point']          # the planted gains are found out of sample
    assert out['r2'] in ('R2m', 'R2f') and out['highest_kept'] in ('R0', 'R1', out['r2'], 'R3')
    # each rung against the highest kept rung below it, in both readings
    assert r.loc['R1', 'below'] == 'R0'
    low = 'R1' if r.loc['R1', 'kept'] else 'R0'
    assert r.loc['R2m', 'below'] == low and r.loc['R2f', 'below'] == low
    assert {'kept_affine_noise', 'reading_dependent', 'point_affine_noise', 'kept_vs_R2m'} <= set(r.columns)
    assert r.loc['R1', 'reading_dependent'] == (r.loc['R1', 'kept'] != r.loc['R1', 'kept_affine_noise'])
    assert set(out['highest_kept_by_reading']) == {E.KEEP_READING, E.BESIDE_READING}
    assert {'R1', 'R2m', 'R2f', 'R1_affine_noise'} <= set(out['controls'].columns)
    assert out['controls'].seed.iloc[0] == E.control_seed(0)
    assert set(out['tables']['R2f']) == {'half_a->half_b', 'half_b->half_a'}
    assert all(len(v) == 2 for v in out['folds']['R1'].values())
    assert all('not_following' in f for v in out['folds']['R2m'].values() for f in v)
    assert len(out['selected']) <= 1
    for s_ in out['lateral_steps']:
        assert s_['needed'] == E.lateral_min_exceeded(2) and len(s_['perm_gains']) == 2
        assert s_['passes'] == (s_['n_exceeded'] >= s_['needed'])
    if out['selected']:
        assert out['fifth'] is not None and 'would_pass' in out['fifth']
    assert E.step2_rung(out['highest_kept']) in (None, 'R2', 'R3')


def test_rung_table_and_arms(fit_case):
    st, pop, true, odors = fit_case
    tA, tB = E.Terms(gains={'g0': 2.0}), E.Terms(gains={'g0': 0.5})
    fits = {'0': dict(held_out=odors[:5], R1={}, R2=tA.as_dict()), '1': dict(held_out=odors[5:8], R1={},
            R2=tB.as_dict()), 'all': dict(held_out=[], R1={'pn_inhadd': 60.0}, R2=E.Terms().as_dict())}
    tab = E.rung_table(st, fits, 'R2', scaled=False)
    pd.testing.assert_series_equal(tab[odors[0]], st.run(terms=tA, reading='preset')[odors[0]])
    pd.testing.assert_series_equal(tab[odors[6]], st.run(terms=tB, reading='preset')[odors[6]])
    pd.testing.assert_series_equal(tab[odors[9]], st.run({'pn_inhadd': 60.0}, E.Terms(), reading='preset')[odors[9]])
    with pytest.raises(KeyError, match='scale'):
        E.rung_table(st, fits, 'R2')
    for k, f in fits.items():                                        # ruling 5: one factor per fit
        f['scale'] = {'R2': dict(zip(fits, (2.0, 0.5, 1.5)))[k], 'R1': dict(zip(fits, (1.0, 1.0, 3.0)))[k]}
    tab_s = E.rung_table(st, fits, 'R2')
    np.testing.assert_allclose(tab_s[odors[0]], 2.0 * tab[odors[0]])
    np.testing.assert_allclose(tab_s[odors[6]], 0.5 * tab[odors[6]])
    np.testing.assert_allclose(tab_s[odors[9]], 1.5 * tab[odors[9]])
    gl = st.glomeruli[:6]
    den = pd.DataFrame(np.random.default_rng(0).normal(size=(6, len(odors))), index=gl, columns=odors)
    imp, set_gl, rec = E.arm_table(E.Arm('megamat', 'imposed'), st, den, gl)
    assert set_gl == gl and rec['n_mapped'] == 6 * len(odors)
    assert rec['clip']['policy'] == 'clip' and rec['clip']['n_set'] == 6 * len(odors)      # ruling 7, by default
    lo, hi = st.rate_bounds()
    assert (imp.loc[gl].to_numpy() >= lo.loc[gl].to_numpy() + E.BOUNDS_MARGIN_HZ - 1e-9).all()
    perm, _, _ = E.arm_table(E.Arm('megamat', 'imposed_odor_perm'), st, den, gl, seed=E.arm_seed('real'))
    assert not np.allclose(perm.loc[gl].to_numpy(), imp.loc[gl].to_numpy())
    perm2, _, r2_ = E.arm_table(E.Arm('megamat', 'imposed_odor_perm'), st, den, gl, seed=E.arm_seed('P1'))
    assert r2_['seed'] == E.SEED + 1 and not np.allclose(perm2.loc[gl].to_numpy(), perm.loc[gl].to_numpy())
    again, _, _ = E.arm_table(E.Arm('megamat', 'imposed_odor_perm'), st, den, gl, seed=E.arm_seed('P1'))
    pd.testing.assert_frame_equal(again, perm2)
    dff, _, rd = E.arm_table(E.Arm('megamat', 'imposed_dff'), st, den * 3.0 + 1.0, gl)
    assert rd['units'] == 'dff' and E.arm_units(E.Arm('megamat', 'imposed')) == E.FIT_UNITS
    pd.testing.assert_frame_equal(dff, imp)                           # rank for rank: a monotone map changes nothing
    r1, set_r1, rr = E.arm_table(E.Arm('megamat', 'R1'), st, None, gl, fits=fits)
    assert set_r1 == st.glomeruli and rr['rung'] == 'R1' and rr['clip']['n_set'] == r1.size
    assert rr['scale']['factors']['R1'] == {'0': 1.0, '1': 1.0, 'all': 3.0} and rr['scale']['named'] == {'R1': ['all']}
    kw, r1 = E.arm_input(st, imp, gl, route=1)
    assert np.isnan(kw['extra_kws']['pn_rates'].loc[st.glomeruli[6:]]).all().all()
    d2 = st.delta.copy()
    d2.loc[gl[:2]] = d2.loc[gl[:2]].to_numpy() * 1.3                  # a reachable table: another input's rates
    tab2 = E.PNStage(d2, pd.Series(st.spont, index=st.glomeruli), st.sp).run(reading='preset')
    kw2, r2 = E.arm_input(st, tab2, st.glomeruli, route=2)
    back = E.PNStage(kw2['orn'], pd.Series(st.spont, index=st.glomeruli), st.sp).run(reading='preset')
    assert np.abs(back.to_numpy() - tab2.to_numpy()).max() <= E.SOLVE_TOL_HZ
    np.testing.assert_allclose(kw2['orn'].to_numpy(), d2.to_numpy(), atol=0.5)


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
    meas = pd.DataFrame(dict(gain=[2.0, 0.1, 1.0], follows=[True, True, False]), index=['a', 'b', 'c'])
    s = E.measured_gain_start(meas, pd.Series([0.5, 1.0, 1.0], index=['a', 'b', 'c']))
    assert s == {'a': 4.0, 'b': 0.25, 'c': 1.0}


# ---------------------------------------------------------------------------
# step 2's points and the classes

def test_planned_points():
    pts = E.planned_step2()
    assert len(pts) == 73 and len(E.planned_step2(with_rung=False)) == 49
    assert len({E.tag_of(a, w) for a, w in pts}) == 73
    assert E.planned_step2(reverse=True)[0] == pts[-1]
    dff = E.planned_step2(with_dff_repeat=True)
    assert len(dff) == 97 and {a.arm for a, _ in dff[73:]} == set(E.DFF_REPEAT_ARMS)
    assert [E.arm_seed(w) for w in E.WIRINGS] == list(range(E.SEED, E.SEED + 8))


def _scores(soma_drop, bouton_ratio, perm_drop=0.0, r1_share=0.0, rung_share=0.0, perm_bouton=1.0):
    """Planted step-2 scores: the imposed arm lowers the somata's error by `soma_drop` and brings the boutons'
    strength-by-odor error to `bouton_ratio` of E134's; R1 and the rung give their shares of both gains."""
    rows = []
    for w in E.WIRINGS:
        for arm, sd, br in (('imposed', soma_drop, bouton_ratio), ('imposed_glom_perm', perm_drop, perm_bouton),
                            ('imposed_odor_perm', perm_drop, perm_bouton),
                            ('R1', r1_share * soma_drop, 1.0 - r1_share * (1.0 - bouton_ratio)),
                            ('rung', rung_share * soma_drop, 1.0 - rung_share * (1.0 - bouton_ratio)),
                            ('rung_glom_perm', 0.0, 1.0)):
            rows.append(dict(arm=arm, wiring=w, stage='kc_soma', error=1.0 - sd, error_base=1.0, strength=1.0,
                             strength_base=1.0))
            # the boutons' population error barely moves: the carry classes read their strength by odor
            rows.append(dict(arm=arm, wiring=w, stage='pn_bouton', error=0.99, error_base=1.0,
                             strength=br * 2.0, strength_base=2.0))
    return pd.DataFrame(rows)


def test_classes():
    c = E.classify(_scores(0.3, 0.4))
    assert c['upstream'] and not c['boutons_only'] and not c['global_carry'] and not c['terms_carry']
    assert E.outcome(c)['outcome'] == 'c'
    c = E.classify(_scores(0.3, 0.4, r1_share=0.6))
    assert c['upstream'] and c['global_carry'] and E.outcome(c)['outcome'] == 'a'
    c = E.classify(_scores(0.3, 0.4, r1_share=0.4))
    assert c['upstream'] and not c['global_carry']
    c = E.classify(_scores(0.3, 0.4, rung_share=0.7))
    assert c['terms_carry'] and not c['global_carry'] and E.outcome(c)['outcome'] == 'b'
    c = E.classify(_scores(0.02, 0.4))
    assert c['boutons_only'] and not c['upstream'] and E.outcome(c)['outcome'] == 'd'
    c = E.classify(_scores(0.3, 0.9, r1_share=0.9))
    assert c['not_upstream'] and not c['permuted_as_much']
    assert not c['global_carry']                                     # carry is read only when upstream holds
    assert E.outcome(c) == dict(outcome='needs the dF/F repeat', pending='e', cls='not_upstream')
    c = E.classify(_scores(0.3, 0.4, perm_drop=0.35))
    assert not c['upstream'] and c['permuted_as_much'] and c['not_upstream']
    assert E.outcome(c)['pending'] == 'f'
    thr = E.class_thresholds(pd.DataFrame(dict(stage='kc_soma', wiring=list(E.WIRINGS),
                                               error=np.linspace(1.0, 1.07, 8))))
    assert thr['kc_soma']['threshold'] == pytest.approx(3 * np.std(np.linspace(1.0, 1.07, 8), ddof=1) / np.sqrt(8))
    c = E.classify(_scores(0.3, 0.4), thr)
    assert c['thresholds']['soma_drop'] == c['thresholds']['soma_within'] == thr['kc_soma']['threshold']


def test_outcome_dff_guard():
    base = E.classify(_scores(0.3, 0.9))
    same = E.classify(_scores(0.25, 0.95))
    other = E.classify(_scores(0.3, 0.4))
    assert E.outcome(base, same) == dict(outcome='e', cls='not_upstream', cls_dff='not_upstream')
    assert E.outcome(base, other)['outcome'] == 'unit-dependent' and E.outcome(base, other)['pending'] == 'e'
    sc = _scores(0.3, 0.9)
    d = sc[sc.arm.isin(E.IMPOSED_ARMS)].assign(arm=lambda x: x.arm + E.DFF_SUFFIX, error=0.6)
    both = pd.concat([sc, d], ignore_index=True)
    ds = E.dff_scores(both)
    assert set(ds.arm) == set(sc.arm) and (ds[ds.arm == 'imposed'].error == 0.6).all()
    assert len(ds) == len(sc)


def test_class_thresholds_read_the_held_runs():
    """Ruling 10: the somata's thresholds from E134's runs held at the megamat flies' response rate (group
    'megamat_held'), the megamat panel, APL on, the default reading, the eight wirings only."""
    from kcq import stage_readings as R
    rows = []
    rng = np.random.default_rng(0)
    for group, scale in (('megamat_held', 0.02), ('megamat_fitted', 0.2)):
        for st in E.CLASS_STAGES:
            for w in list(E.WIRINGS) + ['mean', 'common']:
                for cond in ('on', 'off'):
                    for rd in (R.DEFAULT_READINGS[st].name, 'affine_noise'):
                        rows.append(dict(group=group, stage=st, panel='megamat', condition=cond, reading=rd, wiring=w,
                                         error=1.0 + scale * rng.normal()))
    t = pd.DataFrame(rows)
    thr = E.class_thresholds(t)
    held = t[(t.group == 'megamat_held') & (t.stage == 'kc_soma') & (t.condition == 'on')
             & (t.reading == R.DEFAULT_READINGS['kc_soma'].name) & t.wiring.isin(E.WIRINGS)].error
    assert thr['kc_soma']['n_wirings'] == 8
    assert thr['kc_soma']['threshold'] == pytest.approx(3 * held.std(ddof=1) / np.sqrt(8))
    assert set(thr) == set(E.CLASS_STAGES)
    with pytest.raises(KeyError):
        E.class_thresholds(t[t.group == 'megamat_fitted'])
    agree = E.routes_agree({'pn_bouton': 1.0, 'kc_claw': 1.0, 'kc_soma': 1.0},
                           {'pn_bouton': 1.0 + 0.5 * thr['pn_bouton']['threshold'], 'kc_claw': 1.0,
                            'kc_soma': 1.0 + 2 * thr['kc_soma']['threshold']}, thr)
    assert agree['stages']['pn_bouton']['agree'] and not agree['stages']['kc_soma']['agree'] and not agree['agree']


# ---------------------------------------------------------------------------
# the rulings' other pieces

def test_step0_readings_and_unit_dependence():
    """Step 0 in a reading of planted recordings (`step0_units` on given data) and the gains whose side of 1 differs
    between readings (ruling 1)."""
    gains = np.array([0.5, 1.0, 1.2, 2.0, 0.8, 1.5])
    orn, pn, _ = _flies(gains, n_orn=6, n_pn=7, n_odors=17)
    gl, od = [f'g{i}' for i in range(6)], list(E.MEGAMAT_ODORS)
    data = dict(orn=orn, pn=pn, orn_all=orn, partners=gl, pn_flies=[f'f{i}' for i in range(7)], orn_units='dff',
                pn_units='noise')
    r = E.step0_units('orn_dff_pn_noise', gl, od, n_boot=50, data=data)
    st = r['step']
    assert set(st.flies) == {'all', 'half_a', 'half_b'} and {'follows', 'slope_lo', 'orn_units'} <= set(st.columns)
    assert (st.orn_units == 'dff').all() and (st.pn_units == 'noise').all()
    assert st.follows.all()
    # a second reading in which the PN side of each glomerulus is rescaled (as by a per-glomerulus noise SD)
    sc = np.array([1.0, 1.0, 0.7, 1.0, 1.6, 1.0])[None, :, None]
    r2 = E.step0_units('noise', gl, od, n_boot=50, data=dict(data, pn=pn * sc, orn_units='noise'))
    dep = E.unit_dependence(pd.concat([st, r2['step']], ignore_index=True), ['orn_dff_pn_noise', 'noise']
                            ).set_index('glomerulus')
    side = lambda g, u: dep.loc[g, f'side_{u}']
    for g in gl:
        assert dep.loc[g, 'unit_dependent'] == (side(g, 'orn_dff_pn_noise') != side(g, 'noise'))
    assert dep.unit_dependent.any() and not dep.loc['g3', 'unit_dependent']
    assert set(E.STEP_READINGS) == {'dff', 'noise', 'orn_dff_pn_noise'} and E.RUNG_STEP_READING == 'orn_dff_pn_noise'


def test_keep_rule_against_each_direction():
    rng = np.random.default_rng(0)
    below = 1.0 + 0.01 * rng.normal(size=200)
    rung = below - 0.2
    assert E.keep_rung(rung, below, 0.8, 1.0)[0]
    assert E.keep_rung(rung, below, 0.8, 1.0, dir_rung=[0.7, 0.9], dir_below=[0.95, 1.05])[0]
    assert not E.keep_rung(rung, below, 0.8, 1.0, dir_rung=[0.5, 1.1], dir_below=[0.95, 1.05])[0]
    assert not E.keep_rung(below - 0.001 + 0.01 * rng.normal(size=200), below, 0.999, 1.0)[0]
    o = dict(by_reading={'r': dict(point=0.8, draws=rung, directions={'a->b': 0.7, 'b->a': 0.9})}, **{'a->b': {},
             'b->a': {}})
    b = dict(by_reading={'r': dict(point=1.0, draws=below, directions={'a->b': 0.95, 'b->a': 1.05})})
    k = E.keep_against(o, b, 'r')
    assert k['kept'] and k['no_higher_each_direction'] and k['diff'] == pytest.approx(0.2)
    assert [E.lateral_min_exceeded(n) for n in (20, 10, 5, 2, 1, 0)] == [19, 10, 5, 2, 1, 0]
    assert len({E.control_seed(i) for i in range(5)}) == 5 and E.control_seed(3) == E.control_seed(3)


def test_scale_factor_holds_the_models_rate_scale(fit_case):
    st, pop, true, odors = fit_case
    gl, train = st.glomeruli[:6], odors[:7]
    fit = dict(R1={'pn_inhadd': 80.0, 'pn_tanhsc': 3.0}, R2=E.Terms(gains={'g0': 2.0}).as_dict())
    for rung in ('R1', 'R2'):
        f = E.rung_scale_factor(st, fit, rung, train, gl)
        x = st.run(fit['R1'], None if rung == 'R1' else E.Terms.from_dict(fit['R2']), reading='preset')
        base = st.run(reading='preset')
        assert np.std((f * x).loc[gl, train].to_numpy()) == pytest.approx(np.std(base.loc[gl, train].to_numpy()))
    rec = E.scale_record({'0': dict(scale={'R1': 0.4, 'R2': 1.0}), 'all': dict(scale={'R1': 1.1, 'R2': 2.5})},
                         ('R1', 'R2'))
    assert rec['named'] == {'R1': ['0'], 'R2': ['all']}


def test_clip_table_both_sides(planted):
    mp, rv, st = planted
    lo, hi = st.rate_bounds()
    t = st.run(reading='preset')
    t.loc['g0', 'o0'] = lo.loc['g0', 'o0'] - 5.0
    t.loc['g1', 'o1'] = hi.loc['g1', 'o1'] + 3.0
    t.loc['g2', 'o2'] = np.nan
    out, rec = E.clip_table(t, (lo, hi), ['g0', 'g1', 'g2'])
    assert rec['n_outside'] == 2 and rec['largest_excursion_hz'] == pytest.approx(5.0)
    assert out.loc['g0', 'o0'] == pytest.approx(lo.loc['g0', 'o0'] + E.BOUNDS_MARGIN_HZ)
    assert out.loc['g1', 'o1'] == pytest.approx(hi.loc['g1', 'o1'] - E.BOUNDS_MARGIN_HZ)
    assert np.isnan(out.loc['g2', 'o2']) and rec['n_set'] == 3 * len(st.odors) - 1
    assert rec['share_clipped'] == rec['n_clipped'] / rec['n_set']
    assert rec['named'] == (rec['share_clipped'] > E.CLIP_NAME_SHARE)
    pd.testing.assert_frame_equal(out.loc[st.glomeruli[3:]], t.loc[st.glomeruli[3:]])
    with pytest.raises(ValueError, match='outside'):
        E.clip_table(t, (lo, hi), ['g0'], policy=None)
    many = t.copy()
    many.loc['g0'] = lo.loc['g0'] - 1.0
    assert E.clip_table(many, (lo, hi), ['g0'])[1]['named']


def test_route2_leaves_glomeruli_without_input_and_reports(planted):
    """Ruling 9: the glomeruli without input are not solved (their drift and the LN's change are reported), and a
    target out of reach is reported as not converged, not forced."""
    mp, rv, st = planted
    noinput = [g for g in st.glomeruli if not np.any(st.delta.loc[g].to_numpy() != 0)]
    assert noinput == ['g8', 'g9']
    d2 = st.delta.copy()
    d2.loc[['g0', 'g1']] *= 1.4
    tab = E.PNStage(d2, pd.Series(st.spont, index=st.glomeruli), st.sp).run(reading='preset')
    kw, rec = E.arm_input(st, tab, st.glomeruli, route=2)                 # a reachable table on every glomerulus
    assert rec['converged'] and rec['n_solved'] == 8 and rec['n_without_input'] == 2
    assert rec['set_without_input'] == noinput
    pd.testing.assert_frame_equal(kw['orn'].loc[noinput], st.delta.loc[noinput])
    assert rec['no_input_max_drift_hz'] > 0 and rec['ln_change']['inhA_max_abs'] > 0
    assert rec['ln_change']['pn_inh_factor_max_rel'] > 0
    lo, hi = st.rate_bounds()
    bad = tab.copy()
    bad.loc['g0'] = hi.loc['g0'] + 50.0                               # out of reach
    kw3, rec3 = E.arm_input(st, bad, ['g0'], route=2)
    assert not rec3['converged'] and 'orn' in kw3                       # reported, not forced, and not raised


def test_imposed_check_reads_the_table_in_e134s_reading(planted):
    """Ruling 8: a route-1 table's PN rates in E134's reading (the port's prediction: the drive shifted by a constant
    times the unit kernel) against the table; the terms of the port are such shifts, so the prediction is exact."""
    mp, rv, st = planted
    t = E.Terms(gains={'g0': 1.6, 'g3': 0.6}, laterals=[('g5', 'g6', 0.3)])
    table = st.run(terms=t, reading='preset')
    pred = E.predicted_stage_rates(st, table, ['g0', 'g3', 'g5'])
    want = st.run(terms=t, reading='stage')
    np.testing.assert_allclose(pred.loc[['g0', 'g3', 'g5']].to_numpy(), want.loc[['g0', 'g3', 'g5']].to_numpy(),
                               atol=1e-5)
    pd.testing.assert_frame_equal(pred.loc[['g1', 'g2']], st.run(reading='stage').loc[['g1', 'g2']])
    chk = E.imposed_run_check(table, pred, ['g0', 'g3', 'g5'], pn_matrix_error=0.33)
    d = np.abs(table.loc[['g0', 'g3', 'g5']].to_numpy() - pred.loc[['g0', 'g3', 'g5']].to_numpy()).max()
    assert chk['max_abs_diff_hz'] == pytest.approx(d) and chk['named'] == (d > E.IMPOSED_CHECK_HZ)
    assert chk['pn_matrix_error'] == 0.33 and chk['flies_floor'] == E.PN_MATRIX_FLOOR_NOTE
    assert E.imposed_run_check(table, table + 0.6, ['g0'])['named']
    assert not E.imposed_run_check(table, table + 0.4, ['g0'])['named']


def test_timing_projection_and_budget():
    """The timing helper on a small planted case, the projection's arithmetic, and the fallback on cost only: 5
    controls and 17 folds, then 3 and 17, then 3 and 6 (flagged if still above the budget)."""
    small = dict(n_glomeruli=12, n_input=10, n_shared=8, n_odors=17, n_pn_flies=4, n_trials=3)
    c = E.planted_step1(small)
    assert c['orn'].shape == (9, 8, 17) and c['orn_partners'].shape == (9, 10, 17)
    assert len(c['lock']['folds']) == 17 and sorted(sum(c['lock']['halves'].values(), [])) == sorted(c['pop'].flies)
    tm = E.time_step1(case=c, global_maxiter=3)
    assert set(tm['unit_seconds']) == {'target', 'run', 'R1', 'measured', 'gains', 'lateral', 'read'}
    assert all(v >= 0 for v in tm['unit_seconds'].values())
    u = dict(target=1.0, run=1.0, R1=10.0, measured=1.0, gains=5.0, lateral=1.0, read=1.0)
    p = E.project_step1(u, n_folds=2, n_perm_control=1, n_perm_lateral=1, max_terms=1, check_fifth=False,
                        n_candidates=2, n_directions=2)
    F = 4
    want_r0 = F * 2 + 2
    want_r1 = F * 12 + 2
    want_r2m = F * 3 + 2
    want_r2f = F * 7 + 2
    want_r3 = 2 * (F * 3 + 2) + (F * 4 + 2)
    want_all = 3 * (1 + 10 + 1 + 5 + 1)
    assert p['breakdown_hours']['R3'] * 3600 == pytest.approx(want_r3)
    assert p['real_lane_hours'] * 3600 == pytest.approx(want_r0 + want_r1 + want_r2m + want_r2f + want_r3 + want_all)
    assert p['control_lane_hours'] * 3600 == pytest.approx(want_r1 + want_r2m + want_r2f)
    cheap = E.choose_budget({k: 1e-3 for k in u})
    assert (cheap['n_perm_control'], cheap['n_folds'], cheap['within_budget']) == (5, 17, True)
    p53 = [E.project_step1(u, 17, nc)['total_lane_hours'] for nc in (5, 3)]
    mid = E.choose_budget(u, budget=(p53[0] + p53[1]) / 2)
    assert (mid['n_perm_control'], mid['n_folds']) == (3, 17) and len(mid['tried']) == 2
    dear = E.choose_budget(u, budget=1e-6)
    assert (dear['n_perm_control'], dear['n_folds'], dear['within_budget']) == (3, 6, False)


def test_lock_records_the_rulings_and_is_written_once(tmp_path, monkeypatch):
    r = E.lock_rules(17, 5)
    assert r['directions'] == [['half_a', 'half_b'], ['half_b', 'half_a']] and r['n_folds'] == 17
    assert r['units']['measured_step'] == dict(orn='dff', pn='noise') and r['units']['fit_pn'] == 'noise'
    assert r['loss']['keep_reading'] == E.KEEP_READING and r['loss']['beside_reading'] == 'affine_noise'
    assert r['permutations'] == dict(r['permutations'], lateral=20, lateral_min_exceeded=19, control=5)
    assert r['clip']['policy'] == 'clip' and r['clip']['margin_hz'] == E.BOUNDS_MARGIN_HZ
    assert r['r3']['max_terms'] == 4
    lock = tmp_path / 'E135_lock.json'
    lock.write_text('{}')
    monkeypatch.setattr(E, 'LOCK_COPY', lock)
    with pytest.raises(FileExistsError):
        E.write_lock()


def test_measured_for_reads_the_slope_rule(fit_case):
    st, pop, true, odors = fit_case
    tgt = E.Target.make(pop, odors, st.glomeruli)
    rng = np.random.default_rng(4)
    orn = (st.delta.loc[st.glomeruli, odors].to_numpy()[None] * 0.02
           + rng.normal(scale=0.01, size=(5, len(st.glomeruli), len(odors))))
    orn[:, 3] = rng.normal(scale=0.5, size=(5, len(odors)))           # g3's terminals carry nothing of its dendrites
    m = E.measured_for(tgt, orn, orn, st.glomeruli, n_boot=100, odors_all=odors)
    s = m['step']
    assert (s.follows == (s.slope_lo > 0)).all()
    start = E.measured_gain_start(s, E.model_gains(st, {}, st.glomeruli, odors))
    for g in s.index[~s.follows]:
        assert start[g] == 1.0


def test_kept_stage_table_on_another_input(fit_case):
    """E136's request: the kept rung (or R1) on any ORN input, in route 1's reading, at each fit's scale factor and
    clipped to the bounds on that input; on a held-out panel every odor from the fit to all 17 odors."""
    st, pop, true, odors = fit_case
    fits = {'0': dict(held_out=odors[:5], R1={}, R2=E.Terms(gains={'g0': 2.0}).as_dict(),
                      scale={'R1': 1.0, 'R2': 0.8}),
            'all': dict(held_out=[], R1={'pn_inhadd': 60.0}, R2=E.Terms(gains={'g1': 0.5}).as_dict(),
                        scale={'R1': 1.2, 'R2': 1.5})}
    rng = np.random.default_rng(7)
    new = pd.DataFrame(rng.uniform(-10, 120, (6, 3)), index=st.glomeruli[:6], columns=['x @ -2', 'y @ -2', 'y @ -3'])
    tab, rec = E.kept_stage_table(fits, 'R2f', st, orn=new, panel='validation2')
    assert list(tab.columns) == ['x @ -2', 'y @ -2', 'y @ -3']                  # names kept where abbreviations collide
    on = E.stage_on_input(st, new)
    assert (on.delta.loc[st.glomeruli[6:]] == 0).all().all()
    want = on.run({'pn_inhadd': 60.0}, E.Terms(gains={'g1': 0.5}), reading='preset') * 1.5
    want, _ = E.clip_table(want, on.rate_bounds(), on.glomeruli)
    pd.testing.assert_frame_equal(tab, want[tab.columns])
    assert rec['rung'] == 'R2' and rec['fold_odors'] == [] and rec['clip']['policy'] == 'clip'
    # on the megamat input, the megamat odors from the fold fit that left them out
    mm, r2 = E.kept_stage_table(fits, 'R2', st, panel='megamat')
    raw = st.run({}, E.Terms(gains={'g0': 2.0}), reading='preset') * 0.8
    raw, _ = E.clip_table(raw, st.rate_bounds(), st.glomeruli)
    pd.testing.assert_series_equal(mm[odors[0]], raw[odors[0]])
    assert r2['fold_odors'] == sorted(odors[:5])
    with pytest.raises(ValueError):
        E.kept_stage_table(fits, 'R0', st)
    assert E.odor_labels(['2h @ -3', 'pa @ -3']) == ['2h', 'pa']
