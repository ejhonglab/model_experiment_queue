"""Tests of experiments/E134_population_baseline.py (E134's runner of step 2's saved runs and step 3's short runs):
the points, tags and claims, the held rates from the lock, the identity checks, the scoring and the classes on planted
readings and planted flies. No step-2 or step-3 run; one small real held run (two odors of the lab's input, the
hemibrain wiring) checks that a rebuilt held run reproduces the hold's outputs exactly, only with E134_REAL=1."""
import json
import os
import time

import numpy as np
import pandas as pd
import pytest

from data_prep import stage_population_data as D
from experiments import E134_population_baseline as B
from kcq import harness as H
from kcq import population_fit as P
from kcq import stage_readings as R

from test_e134_stage_readings import ODORS, planted_fly_cache, planted_flies, planted_run


def test_points_tags_and_panels():
    pts = B.planned_step3()
    assert len(pts) == 47 and len(set(pts)) == 47
    assert sum(B.is_family(pt) for pt in pts) == 32
    assert {pt.panel for pt in pts} == {'megamat', 'validation2', 'mixtures'}
    assert all(B.parse_tag(B.tag_of(pt)) == pt for pt in pts)
    assert all(pt.reading == 'fitted' for pt in pts if not B.is_family(pt))
    assert B.planned_step3(models=['P3'], panels=['mixtures']) == [B.Point('mixtures', 'fitted', 'P3'),
                                                                    B.Point('mixtures', 'held', 'P3')]
    with pytest.raises(ValueError):
        B.tag_of(B.Point('kmix', 'fitted', 'real'))
    B.check_wirings()


def test_held_rates_from_the_lock():
    lock = {'remy_megamat_criterion_rate': {'mean': 0.0915},
            'populations': {'kc_soma__kmix__wild_type': {'odors': list(range(8)), 'response_rate': {'mean_over_flies': 0.1}},
                            'kc_soma__cmix__wild_type': {'odors': list(range(8)), 'response_rate': {'mean_over_flies': 0.2}}}}
    assert B.held_rate('megamat', lock) == pytest.approx(0.0915)
    assert B.held_rate('mixtures', lock) == pytest.approx(0.15)
    with pytest.raises(RuntimeError):
        B.held_rate('mixtures', {'populations': {}})
    if B.LOCK_COPY.exists():
        real = json.loads(B.LOCK_COPY.read_text())
        pops = real['populations']
        assert B.held_rate('megamat', real) == pytest.approx(real['remy_megamat_criterion_rate']['mean'])
        assert B.held_rate('megamat', real) == pytest.approx(0.0915, abs=5e-4)
        rates = [pops[f'kc_soma__{p}__wild_type']['response_rate']['mean_over_flies'] for p in R.MIXTURE_PANELS]
        assert min(rates) <= B.held_rate('mixtures', real) <= max(rates)
        try:
            B.HELD_RATE_SOURCE, old = dict(B.HELD_RATE_SOURCE, megamat='lock_criterion'), B.HELD_RATE_SOURCE
            assert B.held_rate('megamat', real) == pytest.approx(
                pops['kc_soma__megamat__wild_type']['response_rate']['mean_over_flies'])
        finally:
            B.HELD_RATE_SOURCE = old


def test_held_targets(monkeypatch):
    from experiments import E83_cross_wiring_fit as E83
    monkeypatch.setattr(E83, 'saved_fit', lambda w: dict(target_on=100, target_off=220, size=1000))
    lock = {'remy_megamat_criterion_rate': {'mean': 0.09}}
    t_off, t_on, f_off, f_on = B.held_targets(B.Point('megamat', 'held', 'real'), 2000, lock)
    assert (t_on, f_on) == (180, 0.09) and f_off == pytest.approx(0.198) and t_off == 396
    monkeypatch.setattr(B, 'HELD_OFF_RULE', 'e107')
    assert B.held_targets(B.Point('megamat', 'held', 'real'), 2000, lock)[2] == pytest.approx(0.22)


def test_check_identity_of_a_held_point():
    run = planted_run()
    out = {c: dict(responses=run.conditions[c]['responses'], spike_counts=run.conditions[c]['spike_counts'])
           for c in ('on', 'off')}
    hold = {c: (run.conditions[c]['responses'].copy(), run.conditions[c]['spike_counts'].copy()) for c in ('on', 'off')}
    ok = B.check_identity(B.Point('megamat', 'held', 'real'), out, hold)
    assert ok['exact'] and len(ok['checks']) == 4
    hold['off'][1].iloc[0, 0] += 1
    bad = B.check_identity(B.Point('megamat', 'held', 'real'), out, hold)
    assert not bad['exact'] and [c['differ'] for c in bad['checks']] == [0, 0, 1, 0]


def test_pn_glomeruli():
    b = B.Built(sim=None, params={}, params_off={}, vars=None, W=None, glomeruli=[], stages=('pn',),
                orn=pd.DataFrame(np.zeros((3, 2)), index=['DA1', 'DM1', 'VA2']))
    assert B.pn_glomeruli(b, 3) == ['DA1', 'DM1', 'VA2']
    b.info['_last'] = dict(wPNKC=pd.DataFrame(np.zeros((2, 4)), columns=pd.Index(['DA1', 'DL5', 'DM1', 'VA2'],
                                                                                    name='glomerulus')))
    assert B.pn_glomeruli(b, 4) == ['DA1', 'DL5', 'DM1', 'VA2']
    with pytest.raises(RuntimeError):
        B.pn_glomeruli(b, 7)


def test_mixture_input():
    try:
        n = B.mixture_input()
    except FileNotFoundError:
        pytest.skip("pyolfsysm's mixture input is not on this host")
    names = [H.odor_abbrev(c[-1] if isinstance(c, tuple) else c) for c in n.columns]
    assert sorted(names) == sorted(D.KMIX_ODORS + D.CMIX_ODORS)
    pi = B.panel_input('mixtures')
    assert list(pi.index) == list(H.apl16_orn_deltas().index) and not pi.isna().any().any()


def _planted_family(n=8):
    flies = planted_fly_cache()
    runs = [planted_run(seed=i, tag=w) for i, w in enumerate(B.WIRINGS[:n])]
    return R.read_family(runs, flies=flies, subset_glomeruli=['G0', 'G1', 'G2', 'G3']), flies


def test_scoring_on_planted_readings():
    fam, flies = _planted_family()
    sc = B.score_family(fam, 'curve_noise', flies, n_floor=50)
    keys = set(zip(sc.stage, sc.panel, sc.condition))
    assert keys == {(s, 'megamat', 'on') for s in R.STAGES}         # APL off has no flies outside the APL16 panel
    for s in R.STAGES:
        x = sc[sc.stage == s]
        assert list(x.wiring) == list(B.WIRINGS) + ['mean', 'common']
        assert np.isfinite(x.error).all() and (x.flies_floor > 0).all() and (x.flies_loo_max >= x.flies_loo_mean).all()
        assert np.isfinite(x[x.wiring.isin(B.WIRINGS)].energy_ratio).all()
    assert np.isfinite(sc[sc.stage == 'orn'].error_with_matrix).all()
    g = B.gate_wirings(fam, 'curve_noise')
    assert set(g.condition) == {'on', 'off'} and set(g.passes) <= {True, False}
    assert (g.groupby(['stage', 'condition']).size() == len(B.WIRINGS)).all()
    # each wiring against the other seven in their own units; the gate is the mean over the wirings
    one = g[(g.stage == 'kc_soma') & (g.condition == 'on')]
    pops = [fam[w].populations[('kc_soma', 'megamat', 'on')] for w in B.WIRINGS]
    od = list(pops[0].values.columns)
    ref = P.fly_reference(pops[1:], od, B.breadth_edges(len(od)))
    assert one.error.iloc[0] == pytest.approx(P.score(pops[0], ref).error)
    assert one.mean_error.iloc[0] == pytest.approx(one.error.mean())
    assert one.n_above.iloc[0] == int((one.error > B.GATE_WIRING_MAX).sum()) and one.max_error.iloc[0] >= one.min_error.iloc[0]
    assert bool(one.passes.iloc[0]) == bool(one.error.mean() <= B.GATE_WIRING_MAX)
    st, inputs = B.step_table(fam, flies)
    steps = {f'{a}->{b}' for a, b in zip(R.STAGES[:-1], R.STAGES[1:])}
    assert set(st.step) == steps and set(inputs) == steps and set(st.level) == {'entry', 'part_distance', 'population'}
    pop = st[st.level == 'population']
    assert np.isfinite(pop.no_step_error).all() and np.isfinite(pop.model_downstream_error).all()
    dist = st[st.level == 'part_distance']
    assert (dist.groupby(['step', 'part']).size() == len(B.WIRINGS)).all() and (dist.distance >= 0).all()
    nulls = B.step_nulls(inputs)
    assert set(nulls) == steps
    assert all(np.isfinite(v) and v > 0 for thr in nulls.values() for v in thr.values())
    assert all(set(thr) <= set(dist.part) for thr in nulls.values())
    classes = B.classify_steps(st, nulls)
    assert {c['cls'] for c in classes.values()} <= {'explained', 'in part', 'no better than no step'}
    # the class from planted thresholds: everything inside, or nothing
    big = {s_: {p: np.inf for p in dist.part.unique()} for s_ in steps}
    zero = {s_: {p: -1.0 for p in dist.part.unique()} for s_ in steps}
    for s_, c in B.classify_steps(st, big).items():
        assert set(c['views']) <= {'strength', 'sharing', 'breadth', 'matrix'} and c['entries_inside']
        no_step = pop[pop.step == s_]
        expect = ('no better than no step' if float(no_step.model_downstream_error.iloc[0]) >=
                  float(no_step.no_step_error.iloc[0]) else 'explained')
        assert c['cls'] == expect
    assert all(c['cls'] in ('in part', 'no better than no step') for c in B.classify_steps(st, zero).values())
    sc.insert(0, 'group', 'megamat_fitted')
    t = B.stage_table(sc)
    assert list(dict.fromkeys(t.stage)) == list(R.STAGES)
    se = B.megamat_stage_error(sc)
    assert list(se.wiring) == list(B.WIRINGS) + ['model', 'common']
    w0 = sc[(sc.wiring == 'real') & (sc.reading == 'curve_noise')].set_index('stage')
    e4 = [w0.loc['pn_dendrite', 'error_with_matrix']] + [w0.loc[x, 'error'] for x in ('pn_bouton', 'kc_claw', 'kc_soma')]
    assert se.stage_error.iloc[0] == pytest.approx(np.sqrt(np.mean(np.square(e4))))
    assert se.stage_error_three.iloc[0] == pytest.approx(np.sqrt(np.mean(np.square(
        [w0.loc[x, 'error'] for x in B.MEGAMAT_ERROR_STAGES]))))
    assert se.stage_error.iloc[-2] == pytest.approx(se.stage_error.iloc[:-2].mean())
    cm = sc[(sc.wiring == 'common') & (sc.reading == 'curve_noise')].set_index('stage')
    assert se.stage_error.iloc[-1] == pytest.approx(np.sqrt(np.mean(np.square(
        [cm.loc['pn_dendrite', 'error_with_matrix']] + [cm.loc[x, 'error'] for x in ('pn_bouton', 'kc_claw', 'kc_soma')]))))


def test_classify_stages_and_panels():
    rows = []
    for stage in R.STAGES:
        for k, rd in enumerate(R.readings_of(stage)[:2]):
            for i, w in enumerate(B.WIRINGS):
                # orn inside under both readings; kc_soma inside under one and far under the other
                err = 0.5 if stage == 'orn' or (stage == 'kc_soma' and k == 0) else 3.0
                rows.append(dict(group='megamat_fitted', stage=stage, panel='megamat', condition='on', reading=rd.name,
                                 wiring=w, error=err, flies_loo_max=1.2))
    sc = pd.DataFrame(rows)
    c = B.classify_stages(sc)
    assert c['orn']['cls'] == 'inside' and c['kc_soma']['cls'] == 'not readable' and c['kc_claw']['cls'] == 'far'
    rd = R.DEFAULT_READINGS['kc_soma'].name
    s2 = pd.DataFrame([dict(source='E90:B', stage='kc_soma', panel='apl16', condition='on', reading=rd, wiring='mean',
                            error=1.0)])
    s3 = pd.DataFrame([dict(group=g, stage='kc_soma', panel=p, condition='on', reading=rd, wiring='mean', error=e)
                       for g, p, e in (('megamat_fitted', 'megamat', 1.15), ('validation2_fitted', 'validation2', 1.5),
                                       ('mixtures_fitted', 'kmix', 0.9), ('mixtures_fitted', 'cmix', 1.3))])
    cp = B.classify_panels(s3, s2)
    assert {k: v['cls'] for k, v in cp.items()} == {'megamat': 'holds', 'validation2': 'panel-specific',
                                                     'kmix': 'holds', 'cmix': 'panel-specific'}


def test_consistency_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(B, 'out_dir', lambda: tmp_path)
    assert B.consistency_gate() == {}
    rows = [dict(condition=c, measure=m, variant_set=v, split_half_rho=r)
            for c, r in (('on', 0.85), ('off', 0.82)) for m in ('error', 'view_strength') for v in ('all', 'inside_bands')]
    pd.DataFrame(rows).to_csv(tmp_path / 'archive_consistency_curve_noise.csv', index=False)
    raw = [dict(r, split_half_rho=0.5 if r['condition'] == 'off' else r['split_half_rho']) for r in rows]
    pd.DataFrame(raw).to_csv(tmp_path / 'archive_consistency_raw.csv', index=False)
    g = B.consistency_gate()
    assert g['gate'] is True and g['curve_noise']['split_half_rho'] == {'on': 0.85, 'off': 0.82}
    assert g['raw']['passes'] is False


def test_commands_and_guard(capsys):
    B.main(['--part', 'commands'])
    out = capsys.readouterr().out.splitlines()
    assert len(out) == 4 and all('--part step3 --lane ' in x for x in out)
    if not B.LOCK_COPY.exists():
        with pytest.raises(PermissionError):
            B.require_lock()


@pytest.mark.skipif(os.environ.get('E134_REAL') != '1', reason='a real held run (3 min at 8 threads on atlas, '
                    '2026-11-24): E134_REAL=1 to run it')
def test_held_rebuild_reproduces_the_hold(tmp_path):
    """The hemibrain wiring of the family on two odors of the lab's input (the breadth held at the full input's): the
    as-fitted run equals the saved fit on those odors; held at 10 % more responses in each condition, then built again
    at the held parameters and run with its stages read, it reproduces the hold's outputs exactly; its stages save and
    load back in step 3's layout."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments.E77_input_screen import hold_rates
    t0 = time.time()
    lab = H.apl16_orn_deltas()
    orn = lab.iloc[:, [4, 8]]                                           # eb, 1-8ol (E90's claw-identity odors)
    Pp = E83.wiring_parts('real')
    b0 = B.build_family('real', 'apl16', orn=orn)
    fitted = B.run_stages(B.Point('apl16', 'fitted', 'real'), b0, b0.params, b0.params_off)
    sv = E83.saved_fit('real')['saved']
    for c in ('on', 'off'):
        s = R._abbrev(fitted[c]['spike_counts'])
        ss = R._abbrev(sv[c][1])[list(s.columns)]
        assert int((s.to_numpy() != ss.to_numpy()).sum()) == 0
    n = {c: int((fitted[c]['responses'].to_numpy() > 0.5).sum()) for c in ('on', 'off')}
    hold = hold_rates(lambda: E83.build_simulator('real', orn=orn, breadth=Pp['breadth'], params=dict(Pp['p'])),
                      target_off=int(round(1.1 * n['off'])), target_on=int(round(1.1 * n['on'])), params=dict(Pp['p']),
                      label='E134 test')
    built = B.build_family('real', 'apl16', hold['params'], orn=orn)
    out = B.run_stages(B.Point('apl16', 'held', 'real'), built, hold['params'], hold['params_off'])
    ident = B.check_identity(B.Point('apl16', 'held', 'real'), out, {c: hold[c] for c in ('on', 'off')})
    assert ident['exact'], ident
    for c in ('on', 'off'):
        st = out[c]
        assert st['claw_win'].shape == (built.W.shape[0], 2) and st['bouton_win'].shape == (built.W.shape[1], 2)
        assert st['pn_win'].shape[0] == len(built.glomeruli) and set(st['apl'].columns) == {'win', 'pre', 'pre_settled'}
    R.save_step3_run('apl16_held_real', out, built.meta, built.orn, dict(model='C2', wiring='real', panel='apl16',
                     reading='held', params=hold['params']), [tmp_path])
    run = R.load_step3_run('apl16_held_real', dirs=[tmp_path])
    assert run.stages() == R.STAGES
    np.testing.assert_array_equal(run.conditions['on']['claw_pre'].to_numpy(),
                                  R._abbrev(out['on']['claw_pre_settled']).to_numpy())
    print(f'held rebuild identity: {time.time() - t0:.0f} s, {hold["n_sims"]} simulations in the hold')


def test_runs_csv_union_and_old_files(tmp_path, monkeypatch):
    """A held row after fitted rows widens the lane's table to the union of the columns (rewritten in place, read back
    whole); a file written before that fix, with a held row under a fitted row's header, is read with the held columns
    inserted before 'stages'; any other row length is refused."""
    base = dict(tag='a', panel='megamat', reading='fitted', fixed_thr=1.0, rate_on=0.1, stages='pn+kc', host='h',
                queue_commit='x')
    held = dict(base, tag='b', reading='held')
    held = {**{k: v for k, v in held.items() if k not in ('stages', 'host', 'queue_commit')},
            **{c: i + 0.5 for i, c in enumerate(B.HELD_EXTRA_COLUMNS)}, 'stages': 'pn+kc', 'host': 'h', 'queue_commit': 'x'}
    monkeypatch.setattr(B, 'out_dir', lambda: tmp_path)
    p = tmp_path / 'runs_3_h_l.csv'
    B._append_csv(p, base)
    B._append_csv(p, dict(base, tag='c'))
    B._append_csv(p, held)
    df = B.read_runs_csv(p)
    assert list(df.tag) == ['a', 'c', 'b'] and df.columns[-1] == B.HELD_EXTRA_COLUMNS[-1]
    assert df.loc[df.tag == 'b', 'thr_factor'].item() == 0.5 and df.loc[df.tag == 'a', 'thr_factor'].isna().item()
    B._append_csv(p, dict(base, tag='d'))
    assert list(B.read_runs_csv(p).tag) == ['a', 'c', 'b', 'd']
    # the old layout: the held row appended under the fitted header
    q = tmp_path / 'runs_3_old_l.csv'
    pd.DataFrame([base]).to_csv(q, index=False)
    pd.DataFrame([held]).to_csv(q, mode='a', header=False, index=False)
    with pytest.raises(pd.errors.ParserError):
        pd.read_csv(q)
    old = B.read_runs_csv(q)
    assert list(old.tag) == ['a', 'b'] and old.loc[1, 'fit_fixed_thr'] == len(B.HELD_EXTRA_COLUMNS) - 0.5
    assert old.loc[1, 'host'] == 'h' and old.loc[1, 'stages'] == 'pn+kc'
    assert set(B.runs_table().tag) == {'a', 'b', 'c', 'd'}
    with open(q, 'a') as f:
        f.write('z,1,2\n')
    with pytest.raises(ValueError):
        B.read_runs_csv(q)


def test_build_overrides_reach_the_simulator(monkeypatch):
    """build_point passes an ORN-input override and extra kwargs (E135's pn_rates) through build_family to E83's
    build_simulator, records them in the Built's info (a table by its sha256), and refuses them for a beside model."""
    from experiments import E83_cross_wiring_fit as E83
    seen = {}
    orn = pd.DataFrame(np.arange(6.0).reshape(3, 2), index=['DM1', 'DL5', 'VA2'], columns=['2h', 'pa'])
    rates = pd.DataFrame(np.ones((1, 2)), index=['DM1'], columns=['2h', 'pa'])

    def fake_build(label, orn=None, breadth=None, params=None, input_tag=None, extra_kws=None):
        seen.update(label=label, orn=orn, extra_kws=extra_kws)
        return 'sim', dict(cells='cells', W='W', kws=dict(a=1))

    monkeypatch.setattr(E83, 'wiring_parts', lambda label: dict(p=dict(fixed_thr=1.0), breadth=None, folder='f', tag='t'))
    monkeypatch.setattr(E83, 'build_simulator', fake_build)
    monkeypatch.setattr(B, '_simulator_built', lambda sim, cells, W, p, orn, stages, info: info)
    pt = B.Point(model='real', panel='megamat', reading='fitted')
    info = B.build_point(pt, orn=orn, extra_kws=dict(pn_rates=rates, scale=2.0))
    assert seen['orn'] is orn and seen['extra_kws']['pn_rates'] is rates and seen['extra_kws']['scale'] == 2.0
    from experiments.E77_input_screen import frame_sha256
    assert info['extra_kws.pn_rates'] == 'sha256:' + frame_sha256(rates) and info['extra_kws.scale'] == 2.0
    assert info['orn_override_sha256'] == frame_sha256(orn)
    plain = B.build_point(pt)
    assert seen['extra_kws'] is None and not any(k.startswith(('extra_kws.', 'orn_override')) for k in plain)
    with pytest.raises(ValueError):
        B.build_point(B.Point(model='uniform_draw', panel='megamat', reading='fitted'), extra_kws=dict(pn_rates=rates))


def test_step2_loads_claw_boutons_for_bouton_pools(monkeypatch):
    """Step 2 loads E90's runs with each claw's bouton when a claw reading pools by bouton (viki's step 2 stopped on
    'curve_noise_b14' at 13:14), and not when no reading to be read needs it; E107's runs are loaded as before."""
    calls = []
    monkeypatch.setattr(R, 'load_e90_run', lambda tag, part='B', settled=True, with_claw_bouton=False:
                        calls.append((tag, part, with_claw_bouton)) or tag)
    assert B.needs_claw_bouton() and B.needs_claw_bouton(['curve_noise_b14'])
    assert not B.needs_claw_bouton(['curve_noise']) and not B.needs_claw_bouton(['curve_noise_g3'])
    e90 = [s for s in B.STEP2_SOURCES if s.source.startswith('E90') and 'kc_claw' in s.stages][0]
    B.step2_runs(e90, with_claw_bouton=True)
    assert calls and all(cb for _, _, cb in calls)
    calls.clear()
    B.step2_runs(e90)
    assert calls and not any(cb for _, _, cb in calls)


def test_save_populations_with_a_multiindex(tmp_path, monkeypatch):
    """save_populations writes a population whose units have a MultiIndex (viki's step 2 stopped on the claws at
    13:25), keeping the levels and names as strings, as it writes a flat index."""
    monkeypatch.setattr(B, 'local_dir', lambda sub: tmp_path)
    mi = pd.MultiIndex.from_tuples([(1, 10), (1, 11), (2, 12)], names=['kc_id', 'claw_id'])
    flat = pd.Index(['DM1', 'DL5', 'VA2'], name='glomerulus')

    class Pop:
        def __init__(self, idx):
            self.values = pd.DataFrame(np.arange(6.0).reshape(3, 2), index=idx, columns=['2h', 'pa'])

    class Reading:
        model, wiring = 'C2', 'real'
        populations = {('kc_claw', 'apl16', 'on'): Pop(mi), ('pn_dendrite', 'apl16', 'on'): Pop(flat)}

    B.save_populations({'r': Reading()}, 'step2', 'E90:B')
    claw = pd.read_parquet(tmp_path / 'E90-B__C2__real__kc_claw__apl16__on.parquet')
    assert list(claw.index.names) == ['kc_id', 'claw_id'] and claw.index[0] == ('1', '10')
    assert claw.to_numpy().tolist() == [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]]
    pn = pd.read_parquet(tmp_path / 'E90-B__C2__real__pn_dendrite__apl16__on.parquet')
    assert list(pn.index) == ['DM1', 'DL5', 'VA2']


def test_write_variants(tmp_path, monkeypatch):
    """variants.csv has one row per run of the runs tables, the tag as the files' variant name, both conditions'
    absolute parameters, the held point's targets, and the model's kwargs from the lane's provenance."""
    monkeypatch.setattr(B, 'out_dir', lambda: tmp_path)
    common = dict(model='C2', fixed_thr=100.0, wAPLKC_scale=0.1, off_fixed_thr=100.0, off_wAPLKC_scale=0.0,
                  rate_on=0.09, rate_off=0.18, stages='pn+bouton+claw', host='h', lane='l')
    B._append_csv(tmp_path / 'runs_3_h_l.csv', dict(tag='megamat_fitted_real', panel='megamat', reading='fitted',
                                                     wiring='real', **common))
    B._append_csv(tmp_path / 'runs_3_h_l.csv', {**dict(tag='megamat_held_P1', panel='megamat', reading='held',
                                                       wiring='P1', **common),
                                                **{c: 1.0 for c in B.HELD_EXTRA_COLUMNS}})
    B._write_json(tmp_path / 'provenance_3_h_l.json', dict(points={'megamat_held_P1': dict(info=dict(
        kws="{'kc_kc_inh': 26.4}", fit_folder='E72_wiring_population', fit_tag='t', breadth='held'))}))
    df = B.write_variants()
    back = pd.read_csv(tmp_path / 'variants.csv')
    assert list(back.tag) == ['megamat_fitted_real', 'megamat_held_P1'] and list(back.variant) == list(back.tag)
    held = back.set_index('tag').loc['megamat_held_P1']
    assert held.thr_factor == 1.0 and held.fixed_thr == 100.0 and held.kws == "{'kc_kc_inh': 26.4}"
    assert held.fit_folder == 'E72_wiring_population' and 'kc_claw:' in held.readings
    assert back.set_index('tag').loc['megamat_fitted_real', 'kws'] != back.set_index('tag').loc['megamat_fitted_real', 'kws']
    assert (back.subdir == 'step3').all() and len(df) == 2


def test_claw_reader_hook(monkeypatch):
    """read_model_stages reads the claws through `claw_reader` when one is given (E136's divisive APL), with E90's
    claw_inputs' arguments, and through E90's claw_inputs otherwise."""
    from experiments import E90_model_stages as E90
    import pyolfsysm.mb_model as M
    n_t, glom = 30, ['DM1', 'DL5']
    W = pd.DataFrame(np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]),
                     columns=pd.MultiIndex.from_tuples([('DM1', 1), ('DL5', 2)], names=['glomerulus', 'bouton_id']))
    sims = dict(pn_sims=[np.ones((2, n_t))], bouton_sims=np.ones((1, 2, n_t)), inh_sims=[np.linspace(0, 1, n_t)])

    class V:
        def __init__(self, **kw):
            self.__dict__.update(kw)

    rv = V(pn=V(), kc=V(wAPLKC=np.array([0.1, 0.2, 0.3])))
    monkeypatch.setattr(M, 'get_sim_var', lambda obj, name: sims[name])
    monkeypatch.setattr(E90, 'step_times', lambda mp, rv: np.linspace(-1.0, 1.0, n_t))
    seen = []

    def reader(Wm, b, inh, wk):
        seen.append((Wm.shape, b.shape, inh.shape, wk.shape))
        return np.full((Wm.shape[0], b.shape[1]), 7.0)

    st = B.read_model_stages(None, rv, W, np.array([0]), ['2h'], glom, ('pn', 'bouton', 'claw'), claw_reader=reader)
    assert seen and seen[0][0] == (3, 2) and (st['claw_win'].to_numpy() == 7.0).all()
    plain = B.read_model_stages(None, rv, W, np.array([0]), ['2h'], glom, ('pn', 'bouton', 'claw'))
    assert not (plain['claw_win'].to_numpy() == 7.0).any()


def test_step_null_distance(monkeypatch):
    """The step's threshold per part is the 95th percentile of the resampled changes' RMS distances from the flies'
    change (hal's StepResample), the same as hal's step_distance threshold for a wiring with every entry."""
    change = {'strength/rate': pd.Series([0.0, 1.0, 2.0], index=list('abc')),
              'breadth/bins': pd.Series([np.nan, np.nan], index=list('xy'))}
    rng = np.random.default_rng(0)
    res = {'strength/rate': change['strength/rate'].to_numpy()[None, :] + rng.normal(0, 0.1, (500, 3)),
           'breadth/bins': np.full((500, 2), np.nan)}
    fake = P.StepResample(change=change, resampled=res, n_up=4, n_down=5)
    monkeypatch.setattr(P, 'step_resample', lambda *a, **k: fake)
    thr = B.step_null_distance([], [], [], None)
    assert set(thr) == {'strength/rate'}
    model = {'w': {'strength/rate': change['strength/rate'] + 0.05}}
    ref = P.step_distance(model, fake)
    assert thr['strength/rate'] == pytest.approx(float(ref.threshold.iloc[0]), rel=1e-12)
