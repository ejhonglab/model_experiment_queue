"""Tests of kcq/stage_readings.py (E134, the model read as each stage is recorded) on planted runs and planted flies;
the loader tests read E90's and E107's saved runs where the shared directory has them (shapes, names, identities), and
nothing is scored against recorded flies."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from data_prep import stage_population_data as D
from kcq import harness as H
from kcq import population_fit as P
from kcq import stage_readings as R

ODORS = tuple(D.MEGAMAT_ODORS[:6])        # (a planted megamat run: the panel's names, so run_panels reads it)


def planted_flies(stage, panel='megamat', cohort='wild_type', odors=ODORS, n_flies=4, n_units=60, criterion=1.0,
                  seed=0, named=False, noise_sd=0.5, missing=None):
    """Planted recorded flies (viki's Population): gamma-distributed trial means plus Gaussian trial noise of SD
    `noise_sd`, three repeats; `missing` maps a fly to units it lacks (NaN rows)."""
    rng = np.random.default_rng(seed)
    tables, trials = {}, {}
    for i in range(n_flies):
        idx = [f'G{j}' for j in range(n_units)] if named else [str(j) for j in range(n_units)]
        base = rng.gamma(0.6, 1.5, size=(n_units, len(odors)))
        reps = base[:, :, None] + noise_sd * rng.standard_normal((n_units, len(odors), 3))
        for u in (missing or {}).get(f'f{i}', ()):
            reps[idx.index(u)] = np.nan
        tables[f'f{i}'] = pd.DataFrame(reps.mean(axis=2), index=idx, columns=list(odors))
        cols = pd.MultiIndex.from_tuples([(o, r) for o in odors for r in range(3)], names=['odor', 'repeat'])
        trials[f'f{i}'] = pd.DataFrame(reps.reshape(n_units, -1), index=idx, columns=cols)
    return D.Population(stage=stage, panel=panel, cohort=cohort, odors=tuple(odors), flies=tuple(tables),
                        tables=tables, criterion=criterion, units='baseline-noise SD', trials=trials,
                        named_units=named)


def planted_run(seed=0, odors=ODORS, n_gl=6, boutons_per_gl=3, claws_per_bouton=10, n_kc=40, tag='real'):
    """A planted saved run in E90's layout, both conditions (APL off: everything larger)."""
    rng = np.random.default_rng(seed)
    gl = [f'G{i}' for i in range(n_gl)]
    b_gl = np.repeat(gl, boutons_per_gl)
    b_idx = pd.MultiIndex.from_arrays([b_gl, np.repeat(np.arange(n_gl), boutons_per_gl),
                                       np.arange(len(b_gl))], names=['glomerulus', 'pn_id', 'bouton_id'])
    c_bouton = np.repeat(np.arange(len(b_gl)), claws_per_bouton)
    n_claw = len(c_bouton)
    kc = np.sort(rng.integers(0, n_kc, n_claw))
    c_idx = pd.MultiIndex.from_arrays([kc, np.arange(n_claw), np.array(['ab'] * n_claw)],
                                      names=[H.KC_ID_LEVEL, 'claw_id', 'kc_type'])
    pn = rng.gamma(1.0, 20.0, size=(n_gl, len(odors)))
    conds = {}
    for c, f in (('on', 1.0), ('off', 1.5)):
        st = {}
        st['pn_win'] = pd.DataFrame(pn + 5.0, index=pd.Index(gl, name='glomerulus'), columns=list(odors))
        st['pn_pre'] = pd.DataFrame(np.full(pn.shape, 5.0), index=st['pn_win'].index, columns=list(odors))
        bw = pn[np.repeat(np.arange(n_gl), boutons_per_gl)] * f * rng.uniform(0.5, 1.5, (len(b_gl), 1))
        st['bouton_win'] = pd.DataFrame(bw + 2.0, index=b_idx, columns=list(odors))
        st['bouton_pre'] = pd.DataFrame(np.full(bw.shape, 2.0), index=b_idx, columns=list(odors))
        cw = bw[c_bouton] * rng.uniform(0.0, 1.0, (n_claw, 1))
        st['claw_win'] = pd.DataFrame(cw, index=c_idx, columns=list(odors))
        st['claw_pre'] = pd.DataFrame(np.zeros(cw.shape), index=c_idx, columns=list(odors))
        counts = np.floor(rng.gamma(0.3, 4.0 * f, size=(n_kc, len(odors))))
        kidx = pd.MultiIndex.from_arrays([np.arange(n_kc), ['ab'] * n_kc], names=[H.KC_ID_LEVEL, 'kc_type'])
        st['spike_counts'] = pd.DataFrame(counts, index=kidx, columns=[f'{o} @ -3' for o in odors])
        st['responses'] = (st['spike_counts'] >= 2).astype(float)
        conds[c] = st
    meta = dict(bouton_glomerulus=b_gl, claw_glomerulus=b_gl[c_bouton], kc_of_claw=kc,
                claw_bouton=np.array([f'b{i}' for i in c_bouton]))
    orn = pd.DataFrame(rng.gamma(1.0, 30.0, size=(n_gl, len(odors))), index=gl, columns=list(odors))
    return R.SavedRun(source='planted', tag=tag, model='C2', wiring=tag, input_panel='megamat', conditions=conds,
                      meta=meta, orn=orn)


def planted_fly_cache(odors=ODORS):
    fixed = {}
    for stage in R.STAGES:
        named = R.NAMED_UNITS[stage]
        fixed[(stage, 'megamat', 'on')] = planted_flies(stage, odors=odors, named=named, n_units=6 if named else 60,
                                                        seed=hash(stage) % 1000)
        fixed[(stage, 'megamat', 'off')] = None
    return R.FlyCache(fixed, fallback=False)


# ---------------------------------------------------------------------------

def test_names_agree_with_the_data_module():
    assert R.STAGES == D.STAGES and R.PANELS == D.PANELS and R.NAMED_UNITS == D.NAMED_UNITS
    for s in R.STAGES:
        for r in R.readings_of(s):
            assert r.stage == s
    assert json.loads(json.dumps(R.reading_definitions())) == R.reading_definitions()


def test_reading_validation():
    with pytest.raises(ValueError):
        R.Reading('x', 'kc_claw', pool='g99')
    with pytest.raises(ValueError):
        R.Reading('x', 'kc_claw', pool='g14', subset='all', criterion='model_spikes')
    with pytest.raises(ValueError):
        R.Reading('x', 'pn_dendrite', glomeruli='every_fly', curve='none')
    assert not R.reading_by_name('kc_soma', 'raw').needs_flies


def test_pool_claws_regions_of_14():
    # glomerulus A: 30 claws (14, 14, and a remainder of 2, dropped); B: 21 claws (14 and a remainder of 7 = 14 / 2, kept)
    g = np.array(['A'] * 30 + ['B'] * 21)
    vals = pd.DataFrame({'a': np.where(g == 'A', 1.0, 10.0) + np.arange(51) * 1e-3, 'b': 1.0},
                        index=pd.MultiIndex.from_arrays([np.arange(51), np.arange(51)], names=[H.KC_ID_LEVEL, 'claw_id']))
    meta = dict(claw_glomerulus=g, claw_bouton=np.array([f'b{i // 7}' for i in range(51)]))
    V = R.pool_claws(vals, meta, 'g14', 'all', seed=0)
    assert len(V) == 4 and V.index.name == 'region' and list(V.index) == list(range(4))
    assert sorted(np.round(V['a'].to_numpy()).astype(int).tolist()) == [1, 1, 10, 10]   # one glomerulus per region
    # each region is the mean of its claws: the column of ones stays one
    assert np.allclose(V['b'], 1.0)
    V2 = R.pool_claws(vals, meta, 'g14', 'all', seed=0)
    pd.testing.assert_frame_equal(V, V2)
    assert len(R.pool_claws(vals, meta, 'r14', 'all', seed=0)) == 51 // 14
    assert len(R.pool_claws(vals, meta, '1', 'all', seed=0)) == 51
    # by bouton: 7 claws a bouton (the last 2), groups of 14 within a bouton: 7 = 14 / 2 is kept, 2 is not
    assert len(R.pool_claws(vals, meta, 'b14', 'all', seed=0)) == len(set(meta['claw_bouton'])) - 1
    # silent regions are dropped
    z = vals.copy()
    z.loc[g == 'B', :] = 0.0
    assert len(R.pool_claws(z, meta, 'g14', 'all', seed=0)) == 2
    with pytest.raises(ValueError):
        R.pool_claws(vals, dict(claw_glomerulus=g), 'b14', 'all')


def test_pool_boutons_regions_without_labels():
    run = planted_run()
    X = R.stage_changes(run, 'pn_bouton', 'on', ODORS)
    V = R.pool_boutons(X, run.meta['bouton_glomerulus'], '1', 'shared', subset_glomeruli=['G0', 'G1'])
    assert len(V) == 6 and V.index.name == 'region'
    assert 'glomerulus' not in (V.index.names or [])
    np.testing.assert_allclose(V.to_numpy(), X.iloc[:6].to_numpy())
    assert len(R.pool_boutons(X, run.meta['bouton_glomerulus'], 'g3', 'all')) == 6
    with pytest.raises(ValueError):
        R.pool_boutons(X, run.meta['bouton_glomerulus'], '1', 'driver')


def test_driver_line_sampling():
    fl = planted_flies('pn_dendrite', named=True, n_units=6, missing={'f1': ['G2']})
    assert R.driver_glomeruli('pn_dendrite', 'megamat', 'every_fly', flies=fl) == ['G0', 'G1', 'G3', 'G4', 'G5']
    assert R.driver_glomeruli('pn_dendrite', 'megamat', 'any_fly', flies=fl) == [f'G{i}' for i in range(6)]
    assert R.driver_glomeruli('pn_dendrite', 'megamat', 'any_fly', flies=fl, model_glomeruli=['G1', 'G9']) == ['G1']
    run = planted_run()
    rd = R.DEFAULT_READINGS['pn_dendrite']
    U = R.sample_units(run, 'pn_dendrite', 'on', ODORS, rd, glomeruli=['G4', 'G0', 'G99'])
    assert list(U.index) == ['G4', 'G0']
    np.testing.assert_allclose(U.loc['G4'].to_numpy(), (run.conditions['on']['pn_win'] -
                                                        run.conditions['on']['pn_pre']).loc['G4'].to_numpy())
    with pytest.raises(ValueError):
        R.sample_units(run, 'pn_dendrite', 'on', ODORS, rd)


def test_curve_monotone_and_set_on_pooled_values_only():
    rng = np.random.default_rng(1)
    m = np.concatenate([np.zeros(900), rng.gamma(1.0, 5.0, 100)])       # many ties at 0, as the somata
    f = rng.normal(0.0, 1.0, 3000) + np.concatenate([np.zeros(2700), rng.gamma(2.0, 2.0, 300)])
    c = R.fit_stage_curve(m, f, kind='quantile')
    grid = np.linspace(-1, 80, 500)
    assert np.all(np.diff(c(grid)) >= -1e-12)
    # equal model values stay equal; the zeros go to the flies' mean over the quantile range they hold
    assert np.ptp(c(np.zeros(5))) == 0
    assert c(np.zeros(1))[0] == pytest.approx(np.sort(f)[:2700].mean(), abs=0.02)
    # a permutation of either array gives the same curve
    c2 = R.fit_stage_curve(rng.permutation(m), rng.permutation(f), kind='quantile')
    np.testing.assert_array_equal(c.x, c2.x)
    np.testing.assert_allclose(c.y, c2.y)
    # the map takes distinct model values onto the flies' distribution
    m2 = rng.gamma(2.0, 3.0, 5000)
    c3 = R.fit_stage_curve(m2, f, kind='quantile')
    q = np.linspace(0.05, 0.95, 19)
    np.testing.assert_allclose(np.quantile(c3(m2), q), np.quantile(f, q), atol=0.05)
    for kind in ('affine', 'saturating'):
        ck = R.fit_stage_curve(m2, f, kind=kind)
        assert np.all(np.diff(ck(np.sort(m2))) >= -1e-12)
    ca = R.fit_stage_curve(m2, f, kind='affine')
    assert ca(m2).mean() == pytest.approx(f.mean()) and ca(m2).std() == pytest.approx(f.std(), rel=1e-6)
    assert np.array_equal(R.fit_stage_curve(m, f, kind='none')(m), m)


def test_curve_unchanged_by_an_odor_permutation():
    run = planted_run()
    fl = planted_flies('kc_soma')
    rd = R.DEFAULT_READINGS['kc_soma']
    perm = [ODORS[i] for i in (3, 0, 5, 1, 4, 2)]
    pa, ca, _, _ = R.read_stage(run, 'kc_soma', 'megamat', ODORS, rd, {'on': fl, 'off': None})
    pb, cb, _, _ = R.read_stage(run, 'kc_soma', 'megamat', perm, rd, {'on': fl, 'off': None})
    np.testing.assert_array_equal(ca['on'].x, cb['on'].x)
    np.testing.assert_allclose(ca['on'].y, cb['on'].y)
    # (the noise draws follow the table's layout, so compare the values before the noise through the curve)
    np.testing.assert_allclose(ca['on'](R.stage_changes(run, 'kc_soma', 'on', ODORS).to_numpy())[:, [3, 0, 5, 1, 4, 2]],
                               cb['on'](R.stage_changes(run, 'kc_soma', 'on', perm).to_numpy()))
    # relabelling the odors (of the run and of the flies alike) leaves the curve unchanged too
    new = [f'x{i}' for i in range(len(ODORS))]
    run2 = planted_run()
    for st in run2.conditions.values():
        st['spike_counts'].columns = [f'{o} @ -3' for o in new]
    fl2 = planted_flies('kc_soma', odors=tuple(new))
    pc, cc, _, _ = R.read_stage(run2, 'kc_soma', 'megamat', new, rd, {'on': fl2, 'off': None})
    np.testing.assert_allclose(ca['on'].y, cc['on'].y)
    np.testing.assert_allclose(pa['on'].values.to_numpy(), pc['on'].values.to_numpy())


def test_noise():
    fl = planted_flies('kc_soma', n_units=200, noise_sd=0.5)
    nm = R.fly_noise(fl, 'residual')
    assert nm.n_trials == 3 and nm.sd == pytest.approx(0.5, rel=0.05)
    nb = R.fly_noise(fl, 'residual_binned')
    assert len(nb.bin_sd) == len(nb.bin_edges) - 1 and np.allclose(nb.bin_sd, 0.5, rtol=0.15)
    nr = R.fly_noise(fl, 'reliability')
    assert set(nr.reliability.index) == set(ODORS) and (nr.reliability > 0).all()
    X = pd.DataFrame(np.ones((4000, 3)), columns=list(ODORS[:3]))
    V, T = R.add_measurement_noise(X, nm, np.random.default_rng(0), return_trials=True)
    assert T.shape == (3, 4000, 3)
    assert (V - 1).to_numpy().std() == pytest.approx(0.5 / np.sqrt(3), rel=0.05)
    V2 = R.add_measurement_noise(X, nm, np.random.default_rng(0))
    pd.testing.assert_frame_equal(V, V2)
    pd.testing.assert_frame_equal(R.add_measurement_noise(X, R.NoiseModel('none', 1), np.random.default_rng(0)), X)
    fl.trials = None
    with pytest.raises(ValueError):
        R.fly_noise(fl, 'residual')


def test_noise_aware_curve_matches_the_flies_after_the_noise():
    rng = np.random.default_rng(2)
    fl = planted_flies('kc_soma', n_units=300, noise_sd=0.8, seed=3)
    nm = R.fly_noise(fl, 'residual')
    fv, fw = R.pooled_fly_values(fl)
    q = np.linspace(0.02, 0.98, 49)
    target = np.quantile(fv, q)
    for m, gain in ((rng.gamma(1.0, 2.0, 3000), 0.5), (np.concatenate([np.zeros(2500), rng.gamma(1, 3, 500)]), 1.0)):
        err = {}
        for aware in (False, True):
            c = R.fit_stage_curve(m, fv, fw, 'quantile', noise=nm if aware else None)
            assert np.all(np.diff(c.y) >= -1e-12)
            noisy = R.add_measurement_noise(pd.DataFrame(c(m)[:, None]), nm, np.random.default_rng(5)).to_numpy().ravel()
            err[aware] = float(np.sqrt(np.mean((np.quantile(noisy, q) - target) ** 2)))
        # nearer the flies after the noise (much nearer for distinct values; no worse with the somata's many zeros)
        assert err[True] <= gain * err[False] + 1e-9


def test_read_run_end_to_end():
    run = planted_run()
    flies = planted_fly_cache()
    mr = R.read_run(run, flies=flies, subset_glomeruli=['G0', 'G1', 'G2', 'G3'])
    assert not mr.skipped
    for stage in R.STAGES:
        for c in ('on', 'off'):
            k = (stage, 'megamat', c)
            pop = mr.populations[k]
            assert isinstance(pop, P.Population) and pop.named_units == R.NAMED_UNITS[stage]
            assert list(pop.values.columns) == list(ODORS)
            crit = flies.get((stage, 'megamat', 'on')).criterion
            pd.testing.assert_frame_equal(pop.responds, (pop.values > crit).astype(float), check_names=False)
            pv = mr.provenance[k]
            assert pv['stage'] == stage and pv['condition'] == c and pv['curve']['kind'] == 'quantile'
            assert pv['noise']['kind'] == 'residual' and pv['flies']['cohort'] == 'wild_type'
        # the curve is set on APL on and applied unchanged to APL off
        assert mr.curves[(stage, 'megamat', 'off')] is mr.curves[(stage, 'megamat', 'on')]
    assert mr.populations[('kc_claw', 'megamat', 'on')].values.index.name == 'region'
    assert len(mr.populations[('pn_bouton', 'megamat', 'on')].values) == 12
    # the same claws form the same regions in both conditions
    n_on = len(mr.before_curve[('kc_claw', 'megamat', 'on')])
    assert n_on == len(mr.before_curve[('kc_claw', 'megamat', 'off')])
    # hal's measures read the populations
    fly_pops = P.fly_populations(flies.get(('kc_soma', 'megamat', 'on')))
    ref = P.fly_reference(fly_pops, list(ODORS), P.APL16_BREADTH_EDGES)
    assert np.isfinite(P.score(mr.populations[('kc_soma', 'megamat', 'on')], ref).error)
    # a family: a wiring's seed index changes its pooling and noise, not its stages' units
    fam = R.read_family([planted_run(tag='real'), planted_run(seed=1, tag='P1')], flies=flies,
                        subset_glomeruli=['G0', 'G1'])
    assert set(fam) == {'real', 'P1'}
    assert fam['P1'].provenance[('kc_claw', 'megamat', 'on')]['pool_seed'] == R.POOL_SEED + 1
    # a stage without flies is skipped and named
    fixed = dict(flies.fixed)
    fixed[('kc_claw', 'megamat', 'on')] = None
    mr2 = R.read_run(run, flies=R.FlyCache(fixed, fallback=False), subset_glomeruli=['G0'])
    assert ('kc_claw', 'megamat') in mr2.skipped and ('kc_claw', 'megamat', 'on') not in mr2.populations


def test_family_scope_curve_and_raw_reading():
    run = planted_run()
    flies = planted_fly_cache()
    rd = R.Reading('fam', 'kc_soma', curve_scope='family')
    fam = R.read_family([planted_run(tag='real'), planted_run(seed=1, tag='P1')], stages=['kc_soma'],
                        readings={'kc_soma': rd}, flies=flies)
    assert fam['real'].curves[('kc_soma', 'megamat', 'on')] is fam['P1'].curves[('kc_soma', 'megamat', 'on')]
    raw = R.read_run(run, stages=['kc_soma'], readings={'kc_soma': R.reading_by_name('kc_soma', 'raw')},
                     flies=R.FlyCache({}, fallback=False))
    pop = raw.populations[('kc_soma', 'megamat', 'on')]
    np.testing.assert_array_equal(pop.values.to_numpy(), run.conditions['on']['spike_counts'].to_numpy())
    np.testing.assert_array_equal(pop.responds.to_numpy(), (run.conditions['on']['spike_counts'] >= 2).to_numpy())


def test_soma_reading():
    run = planted_run()
    on, off = run.conditions['on']['spike_counts'], run.conditions['off']['spike_counts']
    raw = R.SomaReading('raw', 'apl16')(on)
    np.testing.assert_array_equal(raw.responds.to_numpy(), (on.to_numpy() >= 2).astype(float))
    fl = {'on': planted_flies('kc_soma', 'apl16', 'control'), 'off': planted_flies('kc_soma', 'apl16', 'tnt', seed=9)}
    sr = R.SomaReading('curve_noise', 'apl16', flies=fl)
    p_on = sr(on)
    assert list(p_on.values.columns) == list(ODORS)
    with pytest.raises(ValueError):
        sr(off, condition='off')
    p_off = sr(off, condition='off', curve_counts=on)
    assert p_off.values.shape == off.shape
    # the curve set on APL on maps the larger APL-off counts higher on average
    assert p_off.values.to_numpy().mean() > p_on.values.to_numpy().mean()
    own = R.SomaReading('curve_noise_own', 'apl16', flies=fl)(off, condition='off')
    assert own.values.shape == off.shape


def test_step3_layout_roundtrip(tmp_path):
    run = planted_run()
    conds = {}
    for c, st in run.conditions.items():
        x = {k: v for k, v in st.items()}
        for k in R.STAGE_FILES:
            x[f'{k}_pre_settled'] = x[f'{k}_pre'] + 0.5
        x['apl'] = pd.DataFrame({'win': np.ones(len(ODORS)), 'pre': 0.5, 'pre_settled': 0.6}, index=list(ODORS))
        conds[c] = x
    info = dict(model='C2', wiring='real', panel='megamat', reading='fitted', tag='megamat_fitted_real')
    R.save_step3_run('megamat_fitted_real', conds, run.meta, run.orn, info, [tmp_path])
    got = R.load_step3_run('megamat_fitted_real', dirs=[tmp_path])
    assert got.stages() == R.STAGES and got.input_panel == 'megamat' and got.wiring == 'real'
    for c in ('on', 'off'):
        pd.testing.assert_frame_equal(got.conditions[c]['spike_counts'], run.conditions[c]['spike_counts'],
                                      check_names=False)
        np.testing.assert_allclose(got.conditions[c]['claw_pre'].to_numpy(),
                                   run.conditions[c]['claw_pre'].to_numpy() + 0.5)
        np.testing.assert_allclose(got.conditions[c]['claw_pre_e90'].to_numpy(),
                                   run.conditions[c]['claw_pre'].to_numpy())
    assert list(got.meta['claw_bouton']) == list(run.meta['claw_bouton'])
    with pytest.raises(FileNotFoundError):
        R.load_step3_run('no_such_tag', dirs=[tmp_path])


# ---------------------------------------------------------------------------
# the saved runs (read only; nothing scored)

def _e90_available():
    try:
        R._find_dir(R.E90_EXP, 'B', 'real_on_claw_win.parquet')
        H.shared_dir()
        return True
    except FileNotFoundError:
        return False


@pytest.mark.skipif(not _e90_available(), reason="E90's saved stage runs are not on this host")
def test_load_e90_run_saved_matrices():
    run = R.load_e90_run('real', 'B')
    assert run.stages() == R.STAGES
    on = run.conditions['on']
    assert on['pn_win'].shape[0] == 54 and on['bouton_win'].shape[0] == len(run.meta['bouton_glomerulus'])
    assert on['claw_win'].shape[0] == len(run.meta['claw_glomerulus']) and list(on['pn_win'].columns) == H.APL16_ODORS
    d = R._find_dir(R.E90_EXP, 'B', 'real_on_claw_win.parquet')
    pd.testing.assert_frame_equal(on['spike_counts'], pd.read_parquet(d / 'real_on_spike_counts.parquet'))
    pd.testing.assert_frame_equal(on['claw_pre_e90'], R._abbrev(pd.read_parquet(d / 'real_on_claw_pre.parquet')))
    # the settled level: one value per unit for every odor
    assert (on['claw_pre'].nunique(axis=1) == 1).all()
    assert all(f['sha256'] and f['bytes'] > 0 for f in run.files)
    assert R.run_panels(run, 'kc_claw') == [('megamat', tuple(H.MEGAMAT12))]
    assert [p for p, _ in R.run_panels(run, 'kc_soma')] == ['apl16', 'megamat']
    unsettled = R.load_e90_run('real', 'B', settled=False)
    pd.testing.assert_frame_equal(unsettled.conditions['on']['claw_pre'], on['claw_pre_e90'])


@pytest.mark.skipif(not _e90_available(), reason="E107's saved runs are not on this host")
def test_load_e107_run():
    run = R.load_e107_run('megamat', 'fitted', 'real')
    assert run.stages() == ('orn', 'kc_soma')
    panel, odors = R.run_panels(run, 'kc_soma')[0]
    assert panel == 'megamat' and len(odors) == 17 and set(H.MEGAMAT12) <= set(odors)
    assert run.conditions['on']['spike_counts'].shape == run.conditions['off']['spike_counts'].shape
