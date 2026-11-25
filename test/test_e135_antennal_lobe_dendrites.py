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


def _auth_env(tmp_path, monkeypatch):
    """E134's results (no gate) and E135's results in tmp folders; the lock refused with a planted message, so that a
    part that passes its gate stops there (before any model build)."""
    e134, res = tmp_path / 'e134', tmp_path / 'res'
    e134.mkdir()
    res.mkdir()
    monkeypatch.setattr(E, 'e134_results_dir', lambda: e134)
    monkeypatch.setattr(E, 'out_dir', lambda: res)

    def no_lock(path=None):
        raise PermissionError('planted: no lock')

    monkeypatch.setattr(E, 'read_lock', no_lock)
    return e134, res


def test_step1_refused_without_the_gate_or_toms_authorization(tmp_path, monkeypatch):
    e134, res = _auth_env(tmp_path, monkeypatch)
    with pytest.raises(PermissionError, match='measure_gate.json') as ex:
        E.require_step1_allowed()
    assert E.STEP1_AUTH_FILE in str(ex.value)
    for bad in (dict(authorized_by='viki', when='2026-11-24 15:21', words='run it'),
                dict(authorized_by='Tom', when='2026-11-24 15:21')):
        (res / E.STEP1_AUTH_FILE).write_text(json.dumps(bad))
        with pytest.raises(PermissionError, match='absent or incomplete'):
            E.require_step1_allowed()
    for part in (E.part_step1, lambda: E.part_step1_control([0])):
        with pytest.raises(PermissionError, match='measure_gate.json'):
            part()
    with pytest.raises(PermissionError, match='measure_gate.json'):
        E.family_stage('megamat', purpose='step1')
    (e134 / E.E134_GATE_FILE).write_text(json.dumps(dict(verdict='usable')))
    (e134 / E.E134_SCORED_FILE).write_text('{}')
    assert E.require_step1_allowed()['how'] == 'gate'


def test_step1_runs_on_toms_authorization(tmp_path, monkeypatch):
    e134, res = _auth_env(tmp_path, monkeypatch)
    (e134 / E.E134_GATE_FILE).write_text(json.dumps(dict(verdict='not usable')))
    (e134 / E.E134_SCORED_FILE).write_text('{}')
    auth = dict(authorized_by='Tom', when='2026-11-24 15:21', words='Run step 1 now')
    (res / E.STEP1_AUTH_FILE).write_text(json.dumps(auth))
    a = E.require_step1_allowed()
    assert a['how'] == 'authorized' and a['authorization'] == auth and a['gate']['verdict'] == 'not usable'
    for part in (E.part_step1, lambda: E.part_step1_control([0])):           # past the gate, stopped at the lock
        with pytest.raises(PermissionError, match='planted: no lock'):
            part()


def test_step2_still_waits_on_the_gate_under_the_authorization(tmp_path, monkeypatch):
    e134, res = _auth_env(tmp_path, monkeypatch)
    (res / E.STEP1_AUTH_FILE).write_text(json.dumps(dict(authorized_by='Tom', when='2026-11-24 15:21', words='step 1')))
    assert E.require_step1_allowed()['how'] == 'authorized'
    for part in (E.part_step2, E.part_score, lambda: E.family_stage('megamat'),
                 lambda: E.family_stage('apl16', 'P1', purpose='step2')):
        with pytest.raises(PermissionError, match='measure_gate.json'):
            part()


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
    """Planted step-2 scores: the imposed arm lowers the somata's strength-by-odor error by `soma_drop` and brings the
    boutons' to `bouton_ratio` of E134's; R1 and the rung give their shares of both gains."""
    rows = []
    for w in E.WIRINGS:
        for arm, sd, br in (('imposed', soma_drop, bouton_ratio), ('imposed_glom_perm', perm_drop, perm_bouton),
                            ('imposed_odor_perm', perm_drop, perm_bouton),
                            ('R1', r1_share * soma_drop, 1.0 - r1_share * (1.0 - bouton_ratio)),
                            ('rung', rung_share * soma_drop, 1.0 - rung_share * (1.0 - bouton_ratio)),
                            ('rung_glom_perm', 0.0, 1.0)):
            # the somata read on strength by odor (the rulings of 15:41); their population error follows here
            rows.append(dict(arm=arm, wiring=w, stage='kc_soma', error=1.0 - sd, error_base=1.0, strength=1.0 - sd,
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


# ---------------------------------------------------------------------------
# step 2's runner on a planted model path (the PN stage through olfsysm's own layers, the KC model faked), and the
# score on planted runs

def _planted_orn(st):
    """The planted stage's input in a panel input's layout: its glomeruli with input, (panel, 'odor @ -3') columns."""
    gl = [g for g in st.glomeruli if np.any(st.delta.loc[g].to_numpy() != 0)]
    cols = pd.MultiIndex.from_tuples([('megamat', f'{o} @ -3') for o in st.odors], names=['panel', 'odor'])
    return pd.DataFrame(st.delta.loc[gl].to_numpy(), index=pd.Index(gl, name='glomerulus'), columns=cols)


def _planted_wiring(st, n_kc=3):
    """A claw x bouton wiring with one bouton per glomerulus (in the stage's order) and two claws per KC, and its units'
    metadata."""
    cols = pd.MultiIndex.from_tuples([(g, i) for i, g in enumerate(st.glomeruli)], names=['glomerulus', 'bouton_id'])
    idx = pd.MultiIndex.from_tuples([(100 + k, c) for k in range(n_kc) for c in range(2)], names=['kc_id', 'claw_id'])
    W = pd.DataFrame(0.0, index=idx, columns=cols)
    bout = [(2 * k + c) % len(cols) for k in range(n_kc) for c in range(2)]
    for i, b in enumerate(bout):
        W.iloc[i, b] = 1.0
    meta = dict(bouton_glomerulus=np.asarray(st.glomeruli), claw_glomerulus=np.asarray([st.glomeruli[b] for b in bout]),
                kc_of_claw=np.asarray([100 + k for k in range(n_kc) for _ in range(2)]), claw_bouton=np.asarray(bout))
    return W, meta


def _kc_outputs(rng, odors, n_kc=6):
    idx = pd.MultiIndex.from_tuples([(100 + i, 'ab') for i in range(n_kc)], names=['kc_id', 'kc_type'])
    s = pd.DataFrame(rng.integers(0, 4, (n_kc, len(odors))).astype(float), index=idx, columns=list(odors))
    return (s >= 2).astype(float), s


class _FakeBuilt:
    """A built point whose PN stage is olfsysm's own (a fresh copy of the planted layers: route 1's pn_rates set as
    rv.pn.preset_rates, or route 2's solved input as the delta), the KC model faked."""

    def __init__(self, st, orn, extra, W, meta):
        mp, rv = _mp()
        if extra:
            t = extra['pn_rates']
            assert list(t.columns) == list(orn.columns)
            rv.pn.preset_rates = t.reindex(index=st.glomeruli).to_numpy(dtype=float)
        else:
            mp.orn.data.delta = orn.reindex(index=st.glomeruli).fillna(0.0).to_numpy(dtype=float)
        osm.run_ORN_LN_sims(mp, rv)
        osm.run_PN_sims(mp, rv)
        self.mp, self.rv, self.orn, self.W, self.meta = mp, rv, orn, W, meta
        self.stages, self.glomeruli, self.info = ('pn', 'bouton', 'claw'), list(st.glomeruli), dict(planted=True)

    def vars(self):
        return self.mp, self.rv, np.arange(np.asarray(self.mp.orn.data.delta).shape[1])


def _fake_out(st, W, hold, bouton_scale=1.0):
    rng = np.random.default_rng(3)
    out = {}
    for c in ('on', 'off'):
        d = {}
        for k, idx in (('pn', pd.Index(st.glomeruli, name='glomerulus')), ('bouton', W.columns), ('claw', W.index)):
            for x in ('win', 'pre', 'pre_settled'):
                v = rng.uniform(0, 50, (len(idx), len(st.odors))) * (bouton_scale if k == 'bouton' else 1.0)
                d[f'{k}_{x}'] = pd.DataFrame(v, index=idx, columns=st.odors)
        d['apl'] = pd.DataFrame({x: rng.uniform(size=len(st.odors)) for x in ('win', 'pre', 'pre_settled')},
                                index=st.odors)
        d['responses'], d['spike_counts'] = hold[c]
        out[c] = d
    return out


@pytest.fixture
def step2_env(tmp_path, monkeypatch, planted):
    import types
    from experiments import E134_population_baseline as B
    mp, rv, st = planted
    orn = _planted_orn(st)
    res = tmp_path / 'results'
    res.mkdir()
    dirs = [tmp_path / 'local', tmp_path / 'archive']
    monkeypatch.setattr(E, 'out_dir', lambda: res)
    monkeypatch.setattr(E, 'step2_dirs', lambda: dirs)
    gl = st.glomeruli[:6]
    den = pd.DataFrame(np.random.default_rng(11).normal(size=(6, len(st.odors))), index=gl, columns=st.odors)
    monkeypatch.setattr(E, 'dendrite_means', lambda units='noise', glomeruli=None: den * (3.0 if units == 'dff' else 1))
    calls = dict(family=0, identity=0, hold=[], build=[])

    def family_stage(panel, wiring='real'):
        calls['family'] += 1
        return st, dict(sim=None, identity=dict(ok=True))

    def family_identity(panel, wiring, sim):
        calls['identity'] += 1
        return dict(exact=True, checks=[], panel=panel, wiring=wiring)

    p = dict(fixed_thr=100.0, wAPLKC_scale=1.0, wKCAPL_scale=1.0, wAPLPN_scale=1.0, wPNAPL_scale=1.0, wKCKC_scale=0.5)
    fam = dict(p=p, breadth=None, saved=None, n_kc=6, target_on=3, target_off=6, size=6 * len(st.odors), folder='f',
               tag='t')
    W, meta = _planted_wiring(st)
    state = {}
    kc_odors = list(orn.columns.get_level_values('odor'))

    def hold(wiring, orn_, fam_, extra, tg, label):
        rng = np.random.default_rng(len(calls['hold']))
        h = dict(params=dict(p, fixed_thr=110.0), params_off=dict(p, fixed_thr=110.0, wAPLKC_scale=0.0, wAPLPN_scale=0.0),
                 on=_kc_outputs(rng, kc_odors), off=_kc_outputs(rng, kc_odors), thr_factor=1.1, apl_factor=0.9,
                 at_bound_thr=False, at_bound_apl=False, n_sims=4)
        calls['hold'].append(dict(wiring=wiring, orn=orn_, extra=extra, targets=tg))
        state['hold'] = h
        return h

    def build(panel, wiring, params, orn_in, extra):
        calls['build'].append(dict(orn=orn_in, extra=extra))
        return _FakeBuilt(st, orn if orn_in is None else orn_in, extra, W, meta)

    monkeypatch.setattr(E, 'family_stage', family_stage)
    monkeypatch.setattr(E, 'family_identity', family_identity)
    monkeypatch.setattr(E, 'panel_orn', lambda panel: orn)
    monkeypatch.setattr(E, 'family_parts', lambda w: dict(fam))
    monkeypatch.setattr(B, 'held_targets', lambda pt, size, lock=None: (8, 4, 8 / size, 4 / size))
    monkeypatch.setattr(E, '_hold', hold)
    monkeypatch.setattr(E, '_build', build)
    monkeypatch.setattr(E, '_run', lambda panel, wiring, built, params, params_off: _fake_out(st, W, state['hold']))
    monkeypatch.setattr(E, '_provenance', lambda orn_: dict(planted=True))
    fits = {'0': dict(held_out=st.odors[:2], R1={}, R2=E.Terms(gains={'g0': 1.5}).as_dict(), scale={'R1': 1.0, 'R2': 1.0}),
            'all': dict(held_out=[], R1={'pn_inhadd': 60.0}, R2=E.Terms(gains={'g1': 0.7}).as_dict(),
                        scale={'R1': 1.2, 'R2': 0.9})}
    fits_rec = json.loads(json.dumps(dict(all_flies=fits, highest_kept='R2m', step2_rung='R2')))
    lock = dict(glomeruli=list(gl))
    ctx = E.Step2(lock, fits_rec, 'lane-t', dict(sha256='x'))
    return types.SimpleNamespace(st=st, orn=orn, ctx=ctx, calls=calls, dirs=dirs, res=res, den=den, lock=lock,
                                 fits_rec=fits_rec, W=W)


def test_step2_point_route1_runs_checks_and_saves(step2_env):
    """A planted imposed point by route 1: the table enters as pn_rates on the input's own columns, set on the shared
    glomeruli only; olfsysm's PN rates are the table's there and the model's elsewhere; the run reproduces the hold;
    the run, its table, row, identity and provenance are saved; the lane's identity is checked once per wiring."""
    from kcq import stage_readings as R
    v = step2_env
    tab_path, prov_path = v.res / 'runs_2_h_lane-t.csv', v.res / 'provenance_2_h_lane-t.json'
    arm = E.Arm('megamat', 'imposed')
    row = E.run_point2(arm, 'real', v.ctx, tab_path, prov_path)
    assert row['route'] == 1 and row['identity_exact'] and row['lane_identity_exact']
    assert row['pn_check_max_set_hz'] <= E.ROUTE1_TOL_HZ and row['pn_check_max_unset_hz'] <= E.ROUTE1_TOL_HZ
    assert row['n_set_glomeruli'] == 6 and row['n_set_odors'] == len(v.st.odors) and row['imposed_check_max_hz'] >= 0
    pr = v.calls['hold'][0]['extra']['pn_rates']
    assert list(pr.columns) == list(v.orn.columns)
    assert pr.loc[v.st.glomeruli[:6]].notna().all().all() and pr.loc[v.st.glomeruli[6:]].isna().all().all()
    assert v.calls['hold'][0]['targets']['target_on'] == 4 and v.calls['build'][0]['orn'] is None
    run = R.load_step3_run(E.tag_of(arm, 'real'), dirs=v.dirs)
    assert run.info['arm'] == 'imposed' and run.info['route'] == 1 and run.input_panel == 'megamat'
    assert set(run.stages()) >= {'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma'}
    for d in v.dirs:
        t = pd.read_parquet(d / 'megamat_imposed_real_pn_table.parquet')
        assert t.shape == (len(v.st.glomeruli), len(v.st.odors))
    assert E.runs2_table(v.res).tag.tolist() == ['megamat_imposed_real']
    prov = json.loads(prov_path.read_text())['points']['megamat_imposed_real']
    assert prov['model'] == dict(planted=True) and prov['table_sha256'] == row['table_sha256']
    ids = json.loads(E.identity_path('lane-t').read_text())
    assert ids['family:real']['exact'] and ids['megamat_imposed_real']['exact']
    # the next point on the same wiring reuses the stage and the lane's identity; a new process re-reads the identity
    row2 = E.run_point2(E.Arm('megamat', 'R1'), 'real', v.ctx, tab_path, prov_path)
    assert row2['n_set_glomeruli'] == len(v.st.glomeruli) and row2['scale_min'] == 1.0 and row2['scale_max'] == 1.2
    assert v.calls['family'] == 1 and v.calls['identity'] == 1
    ctx2 = E.Step2(v.lock, v.fits_rec, 'lane-t', {})
    E.run_point2(E.Arm('megamat', 'rung'), 'real', ctx2, tab_path, prov_path)
    assert v.calls['family'] == 2 and v.calls['identity'] == 1
    assert set(E.runs2_table(v.res).tag) == {'megamat_imposed_real', 'megamat_R1_real', 'megamat_rung_real'}


def test_step2_point_route2_repeat(step2_env):
    """The repeat by route 2: the solved input in the input's layout and no pn_rates; olfsysm's PN rates on the solved
    input are the solve's (the port's) on the set glomeruli; the solve's record, converged or not, is in the row."""
    v = step2_env
    row = E.run_point2(E.REPEAT_ARM, 'real', v.ctx, v.res / 'runs_2_h_l.csv', v.res / 'prov.json')
    h, b = v.calls['hold'][0], v.calls['build'][0]
    assert h['extra'] is None and b['extra'] is None
    assert b['orn'].index.equals(v.orn.index) and b['orn'].columns.equals(v.orn.columns)
    assert row['route'] == 2 and 'solve_converged' in row and row['imposed_check_max_hz'] >= 0
    # olfsysm on the solved input gives the port's solve exactly; here the planted dendrites are out of the solve's
    # reach, which is reported (ruling 9), not forced
    assert row['pn_check_max_set_hz'] == pytest.approx(row['solve_max_abs_diff_hz'], abs=1e-6)
    assert row['solve_converged'] == (row['solve_max_abs_diff_hz'] <= E.SOLVE_TOL_HZ)
    assert row['ln_inhA_max_abs'] > 0 and row['solve_no_input_max_drift_hz'] >= 0


def test_step2_claims_and_release(step2_env, monkeypatch):
    v = step2_env
    arm = E.Arm('megamat', 'imposed_odor_perm')
    tag = E.tag_of(arm, 'P1')
    assert E.claim(tag, 'other')
    assert E.run_point2(arm, 'P1', v.ctx, v.res / 'r.csv', v.res / 'p.json') is None and not v.calls['hold']
    assert E.run_point2(arm, 'P1', v.ctx, v.res / 'r.csv', v.res / 'p.json', done={tag}) is None
    E.release(tag)

    def boom(*a, **k):
        raise RuntimeError('planted failure')

    monkeypatch.setattr(E, '_hold', boom)
    with pytest.raises(RuntimeError, match='planted failure'):
        E.run_point2(arm, 'P1', v.ctx, v.res / 'r.csv', v.res / 'p.json')
    assert not (v.res / 'claims' / tag).exists()
    assert E.claim(tag, 'again')


def test_step2_and_score_refuse_before_the_gate_lock_and_step1(tmp_path, monkeypatch):
    e134 = tmp_path / 'e134'
    e134.mkdir()
    monkeypatch.setattr(E, 'e134_results_dir', lambda: e134)
    monkeypatch.setattr(E, 'out_dir', lambda: tmp_path)
    for part in (E.part_step2, E.part_score):
        with pytest.raises(PermissionError, match='measure_gate.json'):
            part()
    (e134 / E.E134_GATE_FILE).write_text(json.dumps(dict(verdict='usable')))
    (e134 / E.E134_SCORED_FILE).write_text('{}')

    def no_lock(path=None):
        raise PermissionError('no lock yet (planted)')

    monkeypatch.setattr(E, 'read_lock', no_lock)
    for part in (E.part_step2, E.part_score):
        with pytest.raises(PermissionError, match='no lock yet'):
            part()
    monkeypatch.setattr(E, 'read_lock', lambda path=None: dict(glomeruli=['g0']))
    for part in (E.part_step2, E.part_score):
        with pytest.raises(PermissionError, match="step 1's results"):
            part()


def test_step2_points_rung_and_dff_repeat(tmp_path, monkeypatch):
    monkeypatch.setattr(E, 'out_dir', lambda: tmp_path)
    kept = dict(all_flies={'all': {}}, highest_kept='R3', step2_rung='R3')
    none = dict(all_flies={'all': {}}, highest_kept='R1', step2_rung=None)
    assert len(E.step2_points(kept)) == 73 and len(E.step2_points(none)) == 49
    assert not E.dff_repeat_asked()
    (tmp_path / E.SCORE_VERDICT_FILE).write_text(json.dumps(dict(outcome=dict(outcome='needs the dF/F repeat'))))
    assert E.dff_repeat_asked() and len(E.step2_points(kept)) == 97 and len(E.step2_points(none)) == 73
    assert E.step2_points(kept, reverse=True)[0] == E.step2_points(kept)[-1]
    with pytest.raises(PermissionError, match="step 1's results"):
        E.read_step1_fits(tmp_path)


def test_input_columns_and_layout():
    orn = pd.DataFrame(np.arange(6.0).reshape(2, 3), index=pd.Index(['a', 'b'], name='glomerulus'),
                       columns=pd.MultiIndex.from_tuples([('megamat', f'{o} @ -3') for o in ('2h', 'IaA', 'pa')],
                                                         names=['panel', 'odor']))
    t = pd.DataFrame([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], index=['a', 'b', 'c'], columns=['pa', '2h'])
    x = E.input_columns(t, orn)
    assert x.columns.equals(orn.columns) and np.isnan(x[('megamat', 'IaA @ -3')]).all()
    assert x.loc['c', ('megamat', 'pa @ -3')] == 5.0
    d = pd.DataFrame(np.ones((3, 3)), index=['a', 'b', 'c'], columns=['2h', 'IaA', 'pa'])
    with pytest.raises(RuntimeError, match='outside the input'):
        E.input_layout(d, orn)
    d.loc['c'] = 0.0
    y = E.input_layout(d, orn)
    assert y.index.equals(orn.index) and y.columns.equals(orn.columns) and (y.to_numpy() == 1.0).all()


def _score_runs(drops, bouton_hot=()):
    """Planted saved runs of every planned point (in memory): the KC outputs, the PN and bouton stages on four megamat
    odors, the arm in the info; `bouton_hot` (arm, wiring) of the megamat input get boutons above 300 Hz."""
    from kcq import stage_readings as R
    groups = {}
    rng = np.random.default_rng(0)
    od = list(E.MEGAMAT_ODORS[:4])
    gl = ['g0', 'g1', 'g2']
    for arm, w in E.planned_step2(with_rung=True, with_dff_repeat=True):
        conds = {}
        for c in ('on', 'off'):
            r, s = _kc_outputs(rng, od, n_kc=20)
            hot = 400.0 if (arm.panel == 'megamat' and (arm.arm, w) in bouton_hot) else 50.0
            conds[c] = dict(responses=r, spike_counts=s,
                            pn_win=pd.DataFrame(rng.uniform(10, 50, (3, 4)), index=gl, columns=od),
                            pn_pre=pd.DataFrame(rng.uniform(0, 5, (3, 4)), index=gl, columns=od),
                            bouton_win=pd.DataFrame(rng.uniform(10, 20, (3, 4)) + hot, columns=od),
                            bouton_pre=pd.DataFrame(rng.uniform(0, 5, (3, 4)), columns=od))
        run = R.SavedRun(source='E134:step3', tag=E.tag_of(arm, w), model='C2', wiring=w, input_panel=arm.panel,
                         conditions=conds, meta=dict(bouton_glomerulus=np.asarray(gl)),
                         info=dict(arm=arm.arm, route=2 if arm == E.REPEAT_ARM else 1, params=dict(fixed_thr=100.0)))
        groups.setdefault((arm.panel, arm.arm), []).append(run)
    return groups


def _score_rows(groups, drops, second=None, rise=None):
    """Planted score rows (what `score_group` gives): for each run the four megamat stages, APL on, in their default
    readings and the somata and boutons in their second readings (SECOND_READINGS); the somata's errors and the
    boutons' strength by odor fall by the arm's planted (soma drop, bouton ratio) (`second` for the second reading,
    default the same); `rise` (arm -> (stage, rise)) adds to a stage's strength by odor; an APL16 run's somata on the
    APL16 panel with APL on and off besides."""
    from kcq import stage_readings as R
    rows = []
    second = drops if second is None else second
    for (panel, arm), runs in groups.items():
        for run in runs:
            for st in ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma'):
                for rd_set, dr in (('default', drops), ('second', second)):
                    if rd_set == 'second' and st not in E.SECOND_READINGS:
                        continue
                    rd = R.DEFAULT_READINGS[st].name if rd_set == 'default' else E.SECOND_READINGS[st]
                    sd, br = dr.get(arm, (0.0, 1.0))
                    err = 1.0 - sd if st == 'kc_soma' else 0.99
                    strength = br * 2.0 if st == 'pn_bouton' else err
                    if rise and arm in rise and rise[arm][0] == st:
                        strength += rise[arm][1]
                    rows.append(dict(input_panel=panel, arm=arm, group=f'{panel}_{arm}', model='C2',
                                     wiring=run.wiring, stage=st, panel='megamat', condition='on', reading=rd,
                                     error=err, error_with_matrix=err, view_strength=strength, view_sharing=err,
                                     view_breadth=0.5, view_matrix=0.4))
            if panel == 'apl16':
                for cond in ('on', 'off'):
                    rows.append(dict(input_panel=panel, arm=arm, group=f'{panel}_{arm}', model='C2', wiring=run.wiring,
                                     stage='kc_soma', panel='apl16', condition=cond, reading='curve_noise', error=0.95,
                                     error_with_matrix=0.95, view_strength=0.95, view_sharing=0.95))
    return pd.DataFrame(rows)


def _e134_tables(path):
    from kcq import stage_readings as R
    rng = np.random.default_rng(1)
    rows = []
    for st in ('pn_dendrite',) + E.CLASS_STAGES:
        for rd in (R.DEFAULT_READINGS[st].name,) + ((E.SECOND_READINGS[st],) if st in E.SECOND_READINGS else ()):
            for w in E.WIRINGS:
                e = 1.0 + 0.01 * rng.normal()
                rows.append(dict(group='megamat_held', stage=st, panel='megamat', condition='on', reading=rd, wiring=w,
                                 error=e, error_with_matrix=e, view_strength=2.0 * e if st == 'pn_bouton' else e,
                                 view_sharing=e, view_breadth=0.5, view_matrix=0.4))
    pd.DataFrame(rows).to_csv(path / E.E134_STAGE_TABLE, index=False)
    pd.DataFrame([dict(source='E90:B', group='B', wiring=w, stage='kc_soma', panel='apl16', condition=c,
                       reading='curve_noise', error=1.0 + 0.01 * rng.normal(), view_strength=1.0 + 0.01 * rng.normal())
                  for c in ('on', 'off') for w in list(E.WIRINGS) + ['mean']]).to_csv(path / E.E134_STEP2_SCORES,
                                                                                   index=False)


@pytest.fixture
def score_env(tmp_path, monkeypatch):
    e134, res = tmp_path / 'e134', tmp_path / 'res'
    e134.mkdir()
    res.mkdir()
    _e134_tables(e134)
    monkeypatch.setattr(E, 'e134_results_dir', lambda: e134)
    monkeypatch.setattr(E, 'out_dir', lambda: res)
    fits_rec = dict(all_flies={'all': {}}, highest_kept='R2f', step2_rung='R2')
    monkeypatch.setattr(E, 'require_step2', lambda: (dict(glomeruli=['g0', 'g1']), fits_rec))
    monkeypatch.setattr(E, 'require_e134_lock', lambda: dict(remy_megamat_criterion_rate=dict(mean=0.0915)))
    monkeypatch.setattr(E, '_flies', lambda: None)
    monkeypatch.setattr(E, 'apl16_effect', lambda o: dict(raw_apl_effect_rmse=0.1))
    monkeypatch.setattr(E, 'family_saved', lambda w: {c: _kc_outputs(np.random.default_rng(2), ['2h']) for c in
                                                       ('on', 'off')})
    flag = {('R1', 'P1')}

    def bands(o):
        return dict(rate_rho_on=0.5, realism_flag=bands.next)

    bands.next = ''
    real_realism = E.run_realism

    def realism(run, gl, lock134=None):
        bands.next = 'planted band' if (run.input_panel == 'apl16' and (run.info['arm'], run.wiring) in flag) else ''
        return real_realism(run, gl, lock134)

    monkeypatch.setattr(E, 'run_bands', bands)
    monkeypatch.setattr(E, 'run_realism', realism)
    E.freeze_thresholds()
    return res, monkeypatch


def _set_scores(mp, drops, second=None, rise=None, groups=None):
    mp.setattr(E, 'score_group', lambda runs, flies, lock134=None, readings=None, n_floor=None: _score_rows(
        {(runs[0].input_panel, runs[0].info['arm']): runs}, drops, second, rise).drop(
        columns=['input_panel', 'arm', 'group']))
    if groups is not None:
        mp.setattr(E, 'step2_groups', lambda dirs=None: groups)


UPSTREAM_A = {'imposed': (0.3, 0.4), 'imposed_route2_repeat': (0.3, 0.4), 'R1': (0.18, 0.64), 'imposed_dff': (0.3, 0.4)}


def test_score_on_planted_runs_gives_the_classes(score_env):
    """Planted runs scored on strength by odor: the imposed dendrites lower the somata's strength-by-odor error well
    beyond the frozen threshold and halve the boutons', the permuted controls do not, R1 carries more than half:
    "upstream" and (a), with no guard raised and the same class in the second reading. One APL16 run outside a band
    and one megamat run with hot boutons lose their class (excluded, named), and the classes still hold on the
    remaining wirings."""
    res, mp = score_env
    _set_scores(mp, UPSTREAM_A, groups=_score_runs(UPSTREAM_A, bouton_hot={('imposed', 'M3')}))
    v = E.part_score()
    oc = v['outcome']
    assert v['complete'] and v['cls'] == 'upstream' and oc['outcome'] == 'a' and oc['improvement']
    assert not oc['reading_dependent'] and oc['outcome_second'] == 'a' and oc['guard_raised'] == dict(R1=False,
                                                                                                      rung=False)
    assert v['classes']['global_carry'] and not v['classes']['terms_carry']
    assert v['classes']['thresholds']['soma_col'] == 'strength' and v['measure']['column'] == 'view_strength'
    # R1 carries on the six wirings where neither R1 (P1, its APL16 run outside a band) nor the imposed arm (M3) is
    # excluded
    assert v['classes']['counts']['soma_upstream'] == 7 and v['classes']['counts']['r1_carry'] == 6
    assert set(v['excluded']) == {'R1|P1', 'imposed|M3'}
    thr = json.loads((res / E.THRESHOLDS_FILE).read_text())
    assert v['thresholds']['kc_soma']['threshold'] == thr['classes']['kc_soma']['threshold'] < 0.1
    assert v['routes']['available'] and v['routes']['agree'] and v['routes']['column'] == 'view_strength'
    assert v['routes']['population_error']['kc_soma']['diff'] == 0
    se = v['megamat_stage_error']['megamat_imposed']
    assert se['stage_error_strength_change'] < 0 and se['stage_error_strength_3_change'] < 0
    assert se['stage_error_three_view_change'] < 0 and se['stage_error_strength_sharing'] is not None
    assert v['megamat_stage_error']['megamat_R1']['stage_error_strength_3'] is None
    for f in ('score_scores.csv', 'score_class_rows.csv', 'score_class_rows_second.csv', 'score_guard.csv',
              'score_drops.csv', 'score_stage_error.csv', 'score_realism.csv', 'score_apl_effect.csv',
              'score_table.md', 'provenance_score.json'):
        assert (res / f).exists(), f
    rows = pd.read_csv(res / 'score_class_rows.csv')
    assert rows[(rows.arm == 'imposed') & (rows.wiring == 'M3')].excluded.all()
    assert {'view_sharing_change', 'view_breadth_change', 'error_strength_sharing_change'} <= set(rows.columns)
    g = pd.read_csv(res / 'score_guard.csv')
    assert set(g.panel) == {'megamat', 'apl16'} and len(g) == 2 * (3 * 2 + 2 * 2) and not g.raised.any()
    eff = pd.read_csv(res / 'score_apl_effect.csv')
    assert set(eff.arm) == {'R1', 'rung', 'imposed', 'family'} and len(eff) == 4 * 8
    assert not E.dff_repeat_asked()


def test_score_in_parallel_is_the_serial_score(score_env):
    """`part_score(jobs=3)` (`score_groups` with a pool; a thread pool here, so that the planted `score_group` is the
    one the workers call, and `_score_tags` finding the planted runs by tag) writes the serial score: the same
    score_scores.csv, row for row and in the groups' order, and the same verdict; `jobs` is in the provenance."""
    import concurrent.futures as cf
    res, mp = score_env
    groups = _score_runs(UPSTREAM_A)
    _set_scores(mp, UPSTREAM_A, groups=groups)
    serial = E.part_score()
    s_scores = pd.read_csv(res / 'score_scores.csv')
    by_tag = {r.tag: r for runs in groups.values() for r in runs}
    calls = []

    def tags(ts, lock134, readings, n_floor):
        calls.append(list(ts))
        return E.score_group([by_tag[t] for t in ts], None, lock134, readings, n_floor)

    mp.setattr(E, '_score_tags', tags)
    mp.setattr(E, '_score_pool', lambda jobs: cf.ThreadPoolExecutor(max_workers=jobs))
    par = E.part_score(jobs=3)
    pd.testing.assert_frame_equal(pd.read_csv(res / 'score_scores.csv'), s_scores)
    assert len(calls) == len(groups) and sorted(map(tuple, calls)) == sorted(
        tuple(r.tag for r in runs) for runs in groups.values())
    for k in ('outcome', 'cls', 'classes', 'classes_second', 'guard', 'megamat_stage_error', 'excluded'):
        assert par[k] == serial[k], k
    assert json.loads((res / 'provenance_score.json').read_text())['jobs'] == 3


def test_score_guard_second_reading_and_step1_on_strength(score_env):
    """Rule 3's guard: R1's claws rise on strength by odor on every wiring, and step 1's R1 is above R0 on held-out
    strength in both directions: R1's guard is raised at both, (a) is named at a cost beside and not an improvement.
    The second reading disagrees on the somata: the outcome is reading-dependent and the default decides. Step 1's
    kept R2f above R1 on strength in both directions is named."""
    res, mp = score_env
    d = {'half_a->half_b': 0, 'half_b->half_a': 0}
    s1 = pd.DataFrame([dict(rung='R0', kept=True, **{f'{k}_view_strength': 1.0 for k in d}),
                       dict(rung='R1', kept=True, below='R0', **{f'{k}_view_strength': 1.1 for k in d}),
                       dict(rung='R2m', kept=False, below='R1', **{f'{k}_view_strength': 0.9 for k in d}),
                       dict(rung='R2f', kept=True, below='R1', **{f'{k}_view_strength': 1.2 for k in d})])
    s1.to_csv(res / 'step1_rungs.csv', index=False)
    second = dict(UPSTREAM_A, imposed=(0.0, 0.4))
    _set_scores(mp, UPSTREAM_A, second=second, rise={'R1': ('kc_claw', 0.2)}, groups=_score_runs(UPSTREAM_A))
    v = E.part_score()
    oc = v['outcome']
    assert oc['outcome'] == 'a' and oc['reading_dependent'] and oc['outcome_second'] != 'a'
    assert oc['guard_raised']['R1'] and not oc['improvement']
    assert 'view_strength at kc_claw (megamat, APL on)' in oc['at_a_cost_beside']
    assert any('pn_dendrite' in x for x in oc['at_a_cost_beside'])
    assert v['guard']['R1']['dendrites']['raised'] and v['guard']['rung']['dendrites']['rung'] == 'R2f'
    assert [x['rung'] for x in v['step1_strength']] == ['R1', 'R2f']
    assert v['step1_strength'][0]['name'] == 'kept on the matrix, higher on strength'


def test_score_asks_for_the_dff_repeat_and_waits_for_every_run(score_env):
    """Not upstream: the outcome waits on the dF/F repeat, which the runner then plans; with the repeat's runs
    present and the same class the outcome is (e), unless the second reading differs ('reading-dependent'); with a
    run missing it is 'incomplete'."""
    res, mp = score_env
    drops = {'imposed': (0.0, 0.9)}
    groups = _score_runs(drops)
    dff = {k: v for k, v in groups.items() if k[1] in E.DFF_REPEAT_ARMS}
    first = {k: v for k, v in groups.items() if k[1] not in E.DFF_REPEAT_ARMS}
    _set_scores(mp, drops, groups=first)
    v = E.part_score()
    assert v['complete'] and v['outcome']['outcome'] == 'needs the dF/F repeat' and v['outcome']['pending'] == 'e'
    assert E.dff_repeat_asked()
    mp.setattr(E, 'step2_groups', lambda dirs=None: {**first, **dff})
    v = E.part_score()
    assert {k: v['outcome'][k] for k in ('outcome', 'cls', 'cls_dff', 'reading_dependent')} == dict(
        outcome='e', cls='not_upstream', cls_dff='not_upstream', reading_dependent=False)
    _set_scores(mp, drops, second={'imposed': (0.3, 0.4)})
    v = E.part_score()
    assert v['outcome']['outcome'] == 'reading-dependent' and v['outcome']['pending'] == 'e'
    _set_scores(mp, drops)
    short = dict(first)
    short[('megamat', 'imposed')] = short[('megamat', 'imposed')][:-1]
    mp.setattr(E, 'step2_groups', lambda dirs=None: short)
    v = E.part_score()
    assert not v['complete'] and v['outcome']['outcome'] == 'incomplete' and v['missing'] == ['megamat_imposed_M3']


def test_thresholds_frozen_once_on_strength(score_env):
    """The thresholds are 3 SE (ddof 1) over the eight wirings of view_strength in E134's held runs, written with the
    tables' sha256 before the first arm, and never recomputed once written; the score refuses without them."""
    from kcq import stage_readings as R
    res, mp = score_env
    thr = json.loads((res / E.THRESHOLDS_FILE).read_text())
    t = pd.read_csv(E.e134_results_dir() / E.E134_STAGE_TABLE)
    x = t[(t.stage == 'kc_soma') & (t.reading == R.DEFAULT_READINGS['kc_soma'].name)].view_strength
    assert thr['classes']['kc_soma']['threshold'] == pytest.approx(3 * x.std(ddof=1) / np.sqrt(8))
    assert thr['classes']['kc_soma']['column'] == 'view_strength' and thr['error']['kc_soma']['column'] == 'error'
    assert set(thr['classes']) == {'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma'}
    assert thr['second_reading']['kc_soma']['reading'] == 'raw'
    assert set(thr['apl16_somata']) == {'on', 'off'} and thr['apl16_somata']['on']['view_strength']['n_wirings'] == 8
    assert thr['e134_stage_table_sha256'] == E._sha256(E.e134_results_dir() / E.E134_STAGE_TABLE)
    t.assign(view_strength=t.view_strength * 5).to_csv(E.e134_results_dir() / E.E134_STAGE_TABLE, index=False)
    assert E.freeze_thresholds() == thr
    (res / E.THRESHOLDS_FILE).unlink()
    with pytest.raises(PermissionError, match='no thresholds frozen'):
        E.read_thresholds()
    with pytest.raises(PermissionError, match='no thresholds frozen'):
        E.part_score()


def test_commands_are_module_lanes(capsys, monkeypatch):
    monkeypatch.setattr(E, 'read_lock', lambda path=None: dict(permutations=dict(control=3)))
    E.print_commands()
    out = capsys.readouterr().out
    step2 = [l for l in out.splitlines() if '--part step2' in l]
    assert [l.split(':')[0] for l in step2] == ['viki', 'hal', 'hal', 'atlas', 'disco']
    assert all('-m experiments.E135_antennal_lobe_dendrites' in l for l in step2)
    assert sum('--reverse' in l for l in step2) == 2 and 'OMP_THREAD_LIMIT=8' in step2[3]
    ctl = [l for l in out.splitlines() if '--part step1_control' in l]
    assert [l.split(':')[0] for l in ctl] == ['hal', 'disco', 'hal'] and '--part step1 --lane hal-1' in out


def test_write_variants_one_row_per_run(tmp_path):
    """`write_variants` writes one row per step-2 run of every lane's runs_2_*.csv (the last row per tag), with the
    arm's kind (a dF/F repeat named), both conditions' parameters and the family's build from the lane's provenance."""
    base = dict(panel='megamat', wiring='real', route=1, units='noise', seed=0, fixed_thr=100.0, wAPLKC_scale=0.2,
                off_wAPLKC_scale=0.0, rate_on=0.0915, rate_off=0.2)
    pd.DataFrame([dict(base, tag='megamat_imposed_real', arm='imposed'),
                  dict(base, tag='megamat_imposed_dff_real', arm='imposed_dff', units='dff')]).to_csv(
        tmp_path / 'runs_2_hal_hal-1.csv', index=False)
    pd.DataFrame([dict(base, tag='megamat_imposed_real', arm='imposed', fixed_thr=101.0)]).to_csv(
        tmp_path / 'runs_2_viki_viki.csv', index=False)
    built = dict(model='C2', wiring='real', fit_folder='E67_apl_effect_fits', fit_tag='t', kws=dict(kc_kc_inh=1.0),
                 breadth='b')
    (tmp_path / 'provenance_2_hal_hal-1.json').write_text(json.dumps(dict(points={
        'megamat_imposed_real': dict(built_info=built), 'megamat_imposed_dff_real': dict(built_info=built)})))
    v = E.write_variants(tmp_path)
    assert (tmp_path / 'variants.csv').exists() and len(v) == 2 and set(v.tag) == set(v.variant)
    r = v.set_index('tag')
    assert r.loc['megamat_imposed_real', 'fixed_thr'] == 101.0
    assert 'in dF/F' in r.loc['megamat_imposed_dff_real', 'kind'] and 'in dF/F' not in r.loc['megamat_imposed_real', 'kind']
    assert json.loads(r.loc['megamat_imposed_real', 'kws']) == dict(kc_kc_inh=1.0)
    assert (v.subdir == 'step2').all() and r.loc['megamat_imposed_real', 'fit_folder'] == 'E67_apl_effect_fits'
