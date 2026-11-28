"""Tests of E147's input flies, points and score helpers on planted values (no simulation)."""
import numpy as np
import pandas as pd
import pytest

from experiments import E147_megamat_input_flies as M
from kcq import harness as H


def planted_long(seed=0, missing=(('2023-04-22', 2, 'aa'),), shift5=0.0):
    """A long table like olfdb dataset 1: nine flies, megamat's 17 odors at -3 on the lab input's 36 glomeruli, three
    trials each (each fly a gain times the lab's megamat input plus trial noise), and one diagnostics row per fly that
    sets its largest |trial mean|."""
    rng = np.random.default_rng(seed)
    meg = H.abbrev_columns(M.lab_input('megamat'))
    keys = [('2023-04-22', 2), ('2023-04-22', 3), ('2023-04-26', 2), ('2023-04-26', 3), ('2023-05-08', 1),
            ('2023-05-08', 3), ('2023-05-09', 1), ('2023-05-10', 1), ('2023-06-22', 1)]
    rows = []
    for i, (d, f) in enumerate(keys):
        gain = 1.0 + 0.1 * i
        for o in meg.columns:
            if (d, f, o) in missing:
                continue
            for g in meg.index:
                base = gain * float(meg.loc[g, o]) / 100.0 + (shift5 if o not in H.MEGAMAT12 else 0.0)
                for t in range(3):
                    rows.append((d, f, str(g), 'megamat', f'{o} @ -3', t, base + 0.01 * rng.standard_normal()))
        rows.append((d, f, 'D', 'diagnostics', 'x @ -5', 0, 10.0))
    return pd.DataFrame(rows, columns=['date', 'fly_num', 'roi', 'panel', 'odor1', 'repeat', 'value'])


@pytest.fixture(scope='module')
def planted():
    return M.orn_flies(planted_long())


def test_mean_over_flies_is_the_lab_input(planted):
    for panel in M.PANELS:
        lab = M.lab_input(panel).to_numpy(dtype=float)
        ins = [M.input_fly(f, panel, F=planted).to_numpy(dtype=float) for f in planted['labels']]
        assert np.allclose(np.mean(ins, axis=0), lab, atol=1e-9)


def test_unmeasured_and_kept_entries_are_the_lab_input(planted):
    meg = M.input_fly('f0', 'megamat', F=planted)
    lab = M.lab_input('megamat')
    j = [H.odor_abbrev(c) for c in lab.columns].index('aa')
    assert (meg.iloc[:, j].to_numpy() == lab.iloc[:, j].to_numpy()).all()
    a16 = M.input_fly('f3', 'apl16', F=planted)
    lab16 = M.lab_input('apl16')
    for k, c in enumerate(lab16.columns):
        same = (a16.iloc[:, k].to_numpy() == lab16.iloc[:, k].to_numpy()).all()
        assert same == (H.odor_abbrev(c) in M.APL16_KEPT_AT_LAB)
    assert a16.index.equals(lab16.index) and a16.columns.equals(lab16.columns)


def test_deviation_is_reliability_times_trial_mean_minus_mean(planted):
    from experiments import E85_orn_fly_population as E85
    x = planted['x']
    m = np.nanmean(x, axis=3)
    j = planted['odors'].index('eb')
    k = 5
    M_ = np.nanmean(m[:, j, k])
    r = M.reliabilities(planted)[k]
    dev = M.megamat_deviations(F=planted)
    assert np.allclose(dev[:, j, k], r * (m[:, j, k] - M_))
    i12 = [planted['odors'].index(o) for o in H.MEGAMAT12]
    assert np.allclose(M.reliabilities(planted),
                       np.clip(E85.trial_stats(x[:, i12])['r_g'] * M.e85_part_a()['near_check']['c'], 0, 1))


def test_shared_odors_do_not_depend_on_the_other_five():
    a = M.orn_flies(planted_long(shift5=0.0))
    b = M.orn_flies(planted_long(shift5=0.3))
    da, db = M.megamat_deviations(F=a), M.megamat_deviations(F=b)
    i12 = [a['odors'].index(o) for o in H.MEGAMAT12]
    assert np.allclose(da[:, i12], db[:, i12])
    assert a['slope'] == pytest.approx(b['slope'])
    j = a['odors'].index('aa')
    assert not np.allclose(da[:, j], 0.0)


def test_same_deviation_on_both_panels(planted):
    meg = H.abbrev_columns(M.input_fly('f4', 'megamat', F=planted) - M.lab_input('megamat'))
    a16 = H.abbrev_columns(M.input_fly('f4', 'apl16', F=planted) - M.lab_input('apl16'))
    for o in H.MEGAMAT12:
        assert np.allclose(meg[o].to_numpy(), a16[o].to_numpy())


def test_zero_deviation_input_is_the_lab_input():
    for panel in M.PANELS:
        z = M.zero_deviation_input(panel)
        assert z.equals(M.lab_input(panel).astype(float))


def test_tags_and_points():
    assert M.parse_tag(M.point_tag('A1', 'apl16', 'f8', 'M3')) == ('A1', 'apl16', 'f8', 'M3')
    with pytest.raises(ValueError):
        M.point_tag('A1', 'validation2', 'f0', 'P3')
    p1, p2 = M.planned_points(1), M.planned_points(2)
    assert len(p1) == 72 and len(p2) == 108
    assert [p.wiring for p in p1[:36]] == [w for w in M.T_WIRINGS for _ in range(9)]
    assert all(p.wiring in M.T_WIRINGS for p in p2)
    assert len({M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in p1 + p2}) == 180
    assert all(p.setting.endswith('free') for p in M.planned_points(2, free=True))
    assert M.step_of('familyfree', 'megamat') == 1 and M.step_of('family', 'apl16') == 2


def test_terms_on_a_planted_table():
    w = np.array([0.0, 0.1, -0.1, 0.2])[:, None]
    f = np.linspace(-0.3, 0.3, 9)[None, :]
    E = pd.DataFrame(1.0 + w + f, index=list(M.T_WIRINGS), columns=list(M.FLY_LABELS))
    t = M.terms(E)
    assert t['input_fly'] == pytest.approx(np.std(f, ddof=1))
    assert t['wiring'] == pytest.approx(np.std(w, ddof=1))
    assert t['remainder'] == pytest.approx(0.0, abs=1e-12)
    assert M.term_class(0.4) == 'small' and M.term_class(0.7) == 'a part' and M.term_class(1.0).startswith('as large')
    E.iloc[0, 0] = np.nan
    assert not M.terms(E)['complete']


def test_classes_of_planted_values():
    assert M.stand_in_class(1.00, 1.04, 0.05) == 'a fair stand-in'
    assert M.stand_in_class(0.90, 1.00, 0.05) == 'moves it nearer'
    assert M.stand_in_class(1.10, 1.00, 0.05) == 'moves it further'
    assert M.flies_class(0.5) == 'vary too little' and M.flies_class(1.0) == 'as the flies'
    assert M.flies_class(1.3) == 'vary too much'
    row = dict(n_on=1900, target_on=1905, n_off=4100, target_off=4100, at_bound_thr=False, at_bound_apl='False')
    assert M.is_held(row)
    assert not M.is_held(dict(row, n_on=1880))
    assert not M.is_held(dict(row, at_bound_apl='True'))
    assert not M.is_held({})


def planted_side(n_w=3, n_f=3, seed=0):
    from kcq import population_fit as P
    rng = np.random.default_rng(seed)
    odors = ['a', 'b', 'c', 'd']
    ref = P.FlyReference(parts={}, mean={'strength/rate': pd.Series(1.0, index=odors)},
                         one_fly={'strength/rate': 0.5}, flies=['x', 'y'], odors=odors)
    tags, ws, fs, rows = [], [], [], []
    for i in range(n_w):
        for j in range(n_f):
            tags.append(f't{i}{j}'); ws.append(f'w{i}'); fs.append(f'f{j}')
            rows.append(1.0 + 0.1 * i + 0.05 * j + 0.01 * rng.standard_normal(len(odors)))
    return M.Side(('s', 'p', 'on'), ref, tags, ws, fs, {'strength/rate': np.array(rows)})


def test_against_others_uses_other_wirings_and_input_flies():
    sd = planted_side()
    om = M.others_means(sd)['strength/rate']
    A = sd.arrays['strength/rate']
    # t00 against the mean of the model flies with another wiring and another input fly (t11, t12, t21, t22)
    idx = [sd.tags.index(t) for t in ('t11', 't12', 't21', 't22')]
    assert np.allclose(om[0], A[idx].mean(axis=0))
    sc = M.against_others(sd)
    assert sc.loc['t00', 'part_strength_rate'] == pytest.approx(np.sqrt(np.mean((A[0] - om[0]) ** 2)) / 0.5)
    pop = M.population_score(sd)
    assert pop['part_strength_rate'] == pytest.approx(np.sqrt(np.mean((A.mean(axis=0) - 1.0) ** 2)) / 0.5)


def test_recompare_and_reclass_on_a_planted_record():
    c = dict(n=4, complete=True, changes=[-0.05, -0.06, -0.04, -0.05], change=-0.05, se=0.01, spread=0.02,
             improves=True, worse_than_spread=False)
    r = M.recompare(c, 0.06)
    assert r['se'] == pytest.approx(0.03) and not r['improves']
    assert M.recompare(c, None) == c
    ok = dict(n=4, complete=True, changes=[0.0, 0.01, -0.01, 0.0], change=0.0, se=0.01, spread=0.02, improves=False,
              worse_than_spread=False)
    bad = dict(ok, changes=[0.06, 0.07, 0.05, 0.06], change=0.06, worse_than_spread=True)
    rec = dict(arm='X', megamat=c, stages={s: ok for s in M.LOSS_STAGES}, panels={'apl16/apl16/on': bad},
               panels_beside={}, guard={}, control_fall=0.04, reasons=[], cls='megamat only', panel_missing=[])
    assert M.reclass_record(rec, 'E136_stage_population_fit')['cls'] == 'megamat only'
    rec2 = M.respread(rec, {'apl16/apl16/on': 0.08})
    assert M.reclass_record(rec2, 'E136_stage_population_fit')['cls'] == 'holds'
    rec3 = M.respread(rec, {'stage': 0.06})
    assert M.reclass_record(rec3, 'E136_stage_population_fit')['cls'] == 'does not improve'


def test_points_and_holds_from_planted_rows():
    rows = pd.DataFrame([
        dict(tag='family_megamat_f0_P3', n_on=1900, target_on=1900, n_off=4000, target_off=4000, at_bound_thr=False,
             at_bound_apl=False),
        dict(tag='family_megamat_f1_P3', n_on=1700, target_on=1900, n_off=4000, target_off=4000, at_bound_thr=False,
             at_bound_apl=False),
        dict(tag='familyfree_megamat_f1_P3', n_on=1700, target_on=1900, n_off=4000, target_off=4000),
        dict(tag='not_a_tag', n_on=1, target_on=1, n_off=1, target_off=1)])
    pts = M.points_from_rows(rows)
    assert list(pts.tag) == ['family_megamat_f0_P3', 'family_megamat_f1_P3', 'familyfree_megamat_f1_P3']
    assert list(pts.held) == [True, False, True]
    assert list(pts.step) == [1, 1, 1] and list(pts.group) == ['family/megamat'] * 2 + ['familyfree/megamat']
    assert not M.unheld_rule(8, 72) and M.unheld_rule(9, 72) and M.unheld_rule(5, 36)


def test_a1_fall_and_apl16_shares_on_planted_changes():
    idx = pd.MultiIndex.from_product([list(M.T_WIRINGS), list(M.FLY_LABELS)], names=['wiring', 'fly'])
    fam = pd.Series(1.0, index=idx)
    rng = np.random.default_rng(1)
    a1 = fam - 0.1 + 0.01 * rng.standard_normal(len(idx))
    r = M.a1_fall(M.paired(a1, fam))
    assert r['n'] == 36 and r['fall'] == pytest.approx(0.1, abs=0.01) and r['stands']
    assert r['se'] == pytest.approx(max(r['se_over_input_flies'], r['se_over_wirings']))
    assert not M.a1_fall(M.paired(fam - 0.03, fam))['stands']          # under the control's 0.040
    worse = pd.Series([0.2 if f in ('f0', 'f1', 'f2') else 0.0 for _, f in idx], index=idx)
    sh = M.a1_apl16_shares(M.paired(fam + worse, fam), 0.051)
    assert sh['share_worse'] == pytest.approx(3 / 9) and sh['share_better'] == 0.0


def test_without_odor_matches_e138():
    from experiments import E138_strength_by_odor_entries as E138
    from kcq import population_fit as P
    rng = np.random.default_rng(2)
    odors = ['aa', 'b', 'c', 'd', 'e']

    def pop(name, scale):
        v = pd.DataFrame(rng.gamma(2.0, scale, size=(40, len(odors))), columns=odors)
        return P.Population.from_criterion(v, 3.0, name=name)

    flies = [pop(f'fly{i}', 1.0 + 0.1 * i) for i in range(4)]
    ref = P.fly_reference(flies, odors, P.APL16_BREADTH_EDGES, parts=['strength/rate', 'strength/mean'])
    models = [pop(f'm{i}', 1.2) for i in range(3)]
    arrays = P.stack_views([P.reference_views(m, ref) for m in models], ref)
    sd = M.Side(('kc_soma', 'megamat', 'on'), ref, ['m0', 'm1', 'm2'], ['P3', 'P4', 'M2'], ['f0', 'f0', 'f1'], arrays)
    mine = M.without_odor_rows('x', sd).set_index('row')
    e = E138.odor_entries(models, ref, names=['m0', 'm1', 'm2'])
    e = e[e.odor == 'aa'].set_index('row')
    for row, theirs in (('m0', 'm0'), ('m2', 'm2'), ('population', 'common')):
        assert mine.loc[row, 'view_error'] == pytest.approx(float(e.loc[theirs, 'view_error']))
        assert mine.loc[row, 'view_error_without'] == pytest.approx(float(e.loc[theirs, 'view_error_without']))


def test_score_on_planted_reads(tmp_path, monkeypatch):
    """The score end to end on planted reads: every input fly's read is the lab input's saved run of its wiring
    (E134's held megamat runs, the wirings' own APL16 fits), and A1's reads are the family's, so the population is the
    lab's, the input-fly term is 0 and A1 does not fall."""
    import shutil
    from experiments import E134_population_baseline as B
    shutil.copy(M.LOCK_COPY, tmp_path / M.LOCK_NAME)
    monkeypatch.setitem(M.DIRS, 'results', tmp_path)
    monkeypatch.setitem(M.DIRS, 'local', tmp_path / 'local')
    monkeypatch.setitem(M.DIRS, 'archive', tmp_path / 'archive')
    cache = {}

    def load(tag, step=None):
        setting, panel, fly, w = M.parse_tag(tag)
        if (panel, w) not in cache:
            cache[(panel, w)] = (B.load_run(B.tag_of(B.Point('megamat', 'held', w))) if panel == 'megamat'
                                 else M.lab_apl16_runs([w])[f'lab_apl16_{w}'][2])
        return cache[(panel, w)]

    monkeypatch.setattr(M, 'load_point', load)
    monkeypatch.setattr(M, 'available', lambda pts: pts)
    rows = []
    for f in ('f0', 'f1'):
        for w in M.WIRINGS:
            rows.append(dict(tag=M.point_tag('family', 'megamat', f, w)))
        for w in M.T_WIRINGS:
            rows += [dict(tag=M.point_tag(s_, p_, f, w)) for s_, p_ in (('A1', 'megamat'), ('family', 'apl16'),
                                                                         ('A1', 'apl16'))]
    rows = pd.DataFrame([dict(r, n_on=100, target_on=100, n_off=200, target_off=200, at_bound_thr=False,
                              at_bound_apl=False) for r in rows])
    cl = M.part_score(rows)
    pop = pd.read_csv(tmp_path / 'score_populations.csv')
    fam = pop[(pop.population == 'family/megamat') & (pop.reading == 'all')].set_index('stage')
    assert np.allclose(fam.population_error, fam.lab_common, atol=1e-12)
    assert (fam.stand_in == 'a fair stand-in').all()
    assert np.allclose(fam.term_input_fly, 0.0, atol=1e-12) and (fam.input_fly_term == 'small').all()
    assert fam.loc['megamat_stage_error', 'lab_common'] == pytest.approx(M.E134_IDENTITY['common'])
    assert cl['a1_fall']['all']['fall'] == pytest.approx(0.0, abs=1e-12) and not cl['a1_fall']['all']['stands']
    g = pd.read_csv(tmp_path / 'score_apl16_groups.csv').set_index(['group', 'condition'])
    assert np.isclose(g.loc[('family/apl16', 'on'), 'strength_population'],
                      g.loc[('family/apl16', 'on'), 'strength_lab_common'])
    assert (tmp_path / 'classes.json').exists() and (tmp_path / 'tables.md').exists()


def test_check_inputs_passes_a_bitwise_differing_input_within_tolerance(tmp_path, monkeypatch):
    """A row whose hash differs from `input_fly`'s (the same input built on another host) passes when its saved input
    is within INPUT_ATOL of `input_fly`, and is bad beyond it or with no saved input."""
    from experiments.E77_input_screen import frame_sha256
    monkeypatch.setitem(M.DIRS, 'local', tmp_path / 'local')
    monkeypatch.setitem(M.DIRS, 'archive', tmp_path / 'archive')
    w = M.input_fly('f3', 'megamat')
    rows = []
    for tag, wiring, eps, save in (('family_megamat_f3_P3', 'P3', 0.0, False), ('family_megamat_f3_P4', 'P4', 1e-13, True),
                                   ('family_megamat_f3_M2', 'M2', 1e-6, True), ('family_megamat_f3_M3', 'M3', 1e-13, False)):
        x = w + eps
        sha = frame_sha256(w) if eps == 0.0 else frame_sha256(x) + 'x'
        if save:
            d = M.run_dirs(1, 'family', wiring, 'f3')[0]
            d.mkdir(parents=True, exist_ok=True)
            y = x.copy()
            y.columns = [c[-1] if isinstance(c, tuple) else c for c in y.columns]
            y.to_parquet(d / f'{tag}_orn.parquet')
        rows.append(dict(tag=tag, step=1, setting='family', panel='megamat', fly='f3', wiring=wiring, input_sha256=sha))
    bad, info = M.check_inputs_detail(pd.DataFrame(rows))
    assert sorted(M.check_inputs(pd.DataFrame(rows))) == sorted(bad)
    assert sorted(bad) == ['family_megamat_f3_M2', 'family_megamat_f3_M3']
    assert info['n_bitwise'] == 1 and info['n_within_atol'] == 1
    assert 0 < info['max_abs_diff'] <= M.INPUT_ATOL


def planted_flies_with_a_gap(seed=3):
    """`planted_long` with one glomerulus missing in one fly on one odor (2023-04-26 fly 2, eb, the lab input's first
    glomerulus) beside the aa missing in fly f0."""
    df = planted_long(seed)
    g0 = str(M.lab_input('megamat').index[0])
    gap = (df.date == '2023-04-26') & (df.fly_num == 2) & (df.odor1 == 'eb @ -3') & (df.roi == g0)
    return M.orn_flies(df[~gap].reset_index(drop=True)), g0


def test_measured_and_permuted_deviations_on_planted_values():
    from experiments import E85_orn_fly_population as E85
    F, g0 = planted_flies_with_a_gap()
    x = F['x']
    m = np.nanmean(x, axis=3)
    measured = np.isfinite(m)
    fa, ja, jb, k0 = 0, F['odors'].index('aa'), F['odors'].index('eb'), F['gloms'].index(g0)
    assert not measured[fa, ja].any() and not measured[2, jb, k0] and measured[2, jb].sum() == len(F['gloms']) - 1
    repro = M.megamat_deviations('repro', F)
    meas = M.megamat_deviations(M.VARIANT_MEASURED, F)
    perm = M.megamat_deviations(M.VARIANT_PERMUTED, F)
    for d in (repro, meas, perm):
        assert np.allclose(d.sum(axis=0), 0.0, atol=1e-9)               # the mean over the flies is the lab's input
        assert (d[~measured] == 0.0).all()
    # as measured: the trial mean minus the mean over the flies with the entry, no reliability factor
    mu = np.nanmean(m, axis=0)
    assert np.allclose(meas[measured], (m - mu[None])[measured])
    r = M.reliabilities(F)
    assert np.allclose(repro, r[None, None, :] * meas)
    # permuted: E85's rule (each fly's reproducible values permuted over its measured glomeruli within each odor, the
    # generator seeded by the fly's index), then re-centred per entry over the flies that have it
    Q = E85.permute_deviations(repro, measured)
    for f in range(x.shape[0]):
        for j in range(x.shape[1]):
            idx = measured[f, j]
            assert np.allclose(np.sort(Q[f, j, idx]), np.sort(repro[f, j, idx]))
            assert (Q[f, j, ~idx] == 0.0).all()
    shift = np.where(measured, Q - perm, np.nan)                          # the re-centring, one value per entry
    with np.errstate(invalid='ignore'):
        assert np.nanmax(np.nanmax(shift, axis=0) - np.nanmin(shift, axis=0)) < 1e-9
    assert np.allclose(perm, E85.recenter(Q, measured))
    assert not np.allclose(perm, repro)
    assert np.array_equal(perm, M.megamat_deviations(M.VARIANT_PERMUTED, F))   # seeded: the same draws every time
    # each input fly's input of either variant averages to the lab's input
    lab = M.lab_input('megamat').to_numpy(dtype=float)
    for v in (M.VARIANT_MEASURED, M.VARIANT_PERMUTED):
        ins = [M.input_fly(f, 'megamat', v, F=F).to_numpy(dtype=float) for f in F['labels']]
        assert np.allclose(np.mean(ins, axis=0), lab, atol=1e-9)
    with pytest.raises(ValueError):
        M.megamat_deviations('gainfree', F)


def test_step6a_points_tags_and_folders(tmp_path, monkeypatch):
    monkeypatch.setitem(M.DIRS, 'local', tmp_path / 'local')
    monkeypatch.setitem(M.DIRS, 'archive', tmp_path / 'archive')
    pts = M.planned_points('6a')
    assert len(pts) == 144 and {p.step for p in pts} == {M.STEP6A} and {p.panel for p in pts} == {'megamat'}
    assert {p.wiring for p in pts[:72]} == set(M.T_WIRINGS) and {p.wiring for p in pts[72:]} == set(M.S_WIRINGS)
    assert [(p.setting, p.fly, p.wiring) for p in pts[:18]] == \
        [(s, f, 'P3') for s in M.STEP6A_SETTINGS for f in M.FLY_LABELS]
    assert [p.wiring for p in pts[::18]] == list(M.T_WIRINGS + M.S_WIRINGS)
    for s in M.STEP6A_SETTINGS:
        assert sum(p.setting == s for p in pts[:72]) == 36 and sum(p.setting == s for p in pts[72:]) == 36
    tags = [M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in pts]
    old = {M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in M.planned_points(1) + M.planned_points(2)}
    assert len(set(tags)) == 144 and not set(tags) & old
    assert tags == M.step6a_lock_points()
    assert M.parse_tag('familyperm_megamat_f8_M1') == ('familyperm', 'megamat', 'f8', 'M1')
    assert M.step_of('familymeas', 'megamat') == '6a' and M.step_of('family', 'megamat') == 1
    assert M.step_key('6a') == '6a' and M.step_key('1') == 1 and M.step_key(2) == 2
    assert M.planned_points('1') == M.planned_points(1)
    with pytest.raises(ValueError):
        M.planned_points('6a', free=True)
    assert M.setting_parts('familymeas') == ('family', False, M.VARIANT_MEASURED)
    assert M.setting_parts('familyperm') == ('family', False, M.VARIANT_PERMUTED)
    assert M.setting_parts('family') == ('family', False, 'repro') and M.setting_parts('A1free') == ('A1', True, 'repro')
    assert set(M.SETTING_VARIANT) == set(M.SETTINGS) == set(M.SETTING_BASE)
    assert M.HELD_SETTINGS == ('family', 'A1')
    d = M.run_dirs('6a', 'familymeas', 'P3', 'f0')
    assert d[0] == tmp_path / 'local' / 'step6a' / 'familymeas' / 'P3' / 'f0'
    assert d[1] == tmp_path / 'archive' / 'step6a' / 'familymeas' / 'P3' / 'f0'
    assert M.run_dirs(1, 'family', 'P3', 'f0')[1] == tmp_path / 'archive' / 'step1' / 'family' / 'P3' / 'f0'
    rows = pd.DataFrame([dict(tag='familymeas_megamat_f0_P3', n_on=1, target_on=1, n_off=2, target_off=2),
                         dict(tag='family_megamat_f0_P3', n_on=1, target_on=1, n_off=2, target_off=2)])
    p = M.points_from_rows(rows)
    assert list(p.step) == ['6a', 1] and list(p.group) == ['familymeas/megamat', 'family/megamat']
    assert list(p.held) == [True, True]


def test_the_lock_names_step6a():
    """The committed lock carries the step-6a revision after the first one, the 144 tags and the step-6a block as this
    module has them, and step 1 and 2's points as before."""
    import json
    lock = json.loads(M.LOCK_COPY.read_text())
    assert [r['time'] for r in lock['revisions']] == [r['time'] for r in M.REVISIONS]
    assert lock['revisions'] == [dict(r) for r in M.REVISIONS]
    assert lock['points']['step6a'] == M.step6a_lock_points()
    assert lock['step6a'] == json.loads(json.dumps(M.STEP6A_LOCK))
    assert lock['points']['step1'] == [M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in M.planned_points(1)]
    assert lock['points']['step2'] == [M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in M.planned_points(2)]
    assert all(lock['settings'][s] == M.STEP6A_SETTING_TEXT[s] for s in M.STEP6A_SETTINGS)


def test_step6a_readings_on_planted_rows():
    def row(stand_in, against=0.5, model_flies=None):
        return dict(stand_in=stand_in, against_others_mean=against, population_error=0.9, lab_common=0.95,
                    spread_as_scored=0.02, model_flies=model_flies or M.flies_class(against))

    near, fair, far = 'moves it nearer', 'a fair stand-in', 'moves it further'
    R = row(near)
    assert M.step6a_reading(R, row(near), row(fair))['structure_or_variance'] == "the structure's"
    assert M.step6a_reading(R, row(near), row(near))['structure_or_variance'] == "the variance's"
    assert M.step6a_reading(R, row(near), row(far))['structure_or_variance'].startswith('neither')
    assert M.step6a_reading(row(fair), row(near), row(near))['structure_or_variance'] == 'not read'
    assert M.step6a_reading(row(far), row(near), row(fair))['structure_or_variance'] == 'not read'
    assert M.step6a_reading(R, row(near), row('no class'))['structure_or_variance'] == 'no reading'
    assert M.step6a_reading(R, row(near), None)['structure_or_variance'] == 'no reading'
    # (e) at the upper bound: under M a model fly against the others stays under 0.75 where R's vary too little
    x = M.step6a_reading(row(near, 0.45), row(near, 0.70), row(fair))
    assert x['e_applies'] and x['e_at_upper_bound'] is True
    x = M.step6a_reading(row(near, 0.45), row(near, 0.86), row(fair))
    assert x['e_applies'] and x['e_at_upper_bound'] is False
    x = M.step6a_reading(row(near, 0.95), row(near, 0.70), row(fair))
    assert not x['e_applies'] and x['e_at_upper_bound'] is None
    x = M.step6a_reading(row(near, 0.45), row('no class', 0.70, 'no class'), row(fair))
    assert x['e_applies'] and x['e_at_upper_bound'] is None
    assert x['repro_against_others_mean'] == 0.45 and x['measured_stand_in'] == 'no class'


def test_check_inputs_reads_the_settings_variant(tmp_path, monkeypatch):
    """A step-6a row passes on the hash of its setting's variant's input, is bad on another variant's, and is bad when
    its 'variant' column is not its setting's."""
    from experiments.E77_input_screen import frame_sha256
    monkeypatch.setitem(M.DIRS, 'local', tmp_path / 'local')
    monkeypatch.setitem(M.DIRS, 'archive', tmp_path / 'archive')
    sha = {v: frame_sha256(M.input_fly('f2', 'megamat', v)) for v in ('repro', 'measured', 'permuted')}
    assert len(set(sha.values())) == 3
    rows = [dict(tag='familymeas_megamat_f2_P3', setting='familymeas', input_sha256=sha['measured'], variant='measured'),
            dict(tag='familyperm_megamat_f2_P3', setting='familyperm', input_sha256=sha['permuted'], variant=np.nan),
            dict(tag='familymeas_megamat_f2_P4', setting='familymeas', input_sha256=sha['repro']),
            dict(tag='familyperm_megamat_f2_P4', setting='familyperm', input_sha256=sha['permuted'], variant='measured'),
            dict(tag='family_megamat_f2_P3', setting='family', input_sha256=sha['repro'])]
    rows = pd.DataFrame([dict(r, step=M.step_of(r['setting'], 'megamat'), panel='megamat', fly='f2',
                              wiring=r['tag'].split('_')[-1]) for r in rows])
    bad, info = M.check_inputs_detail(rows)
    assert sorted(bad) == ['familymeas_megamat_f2_P4', 'familyperm_megamat_f2_P4'] and info['n_bitwise'] == 3


def test_score_with_and_without_step6a_reads(tmp_path, monkeypatch):
    """The score on planted reads (every read the lab input's saved run of its wiring, as in
    `test_score_on_planted_reads`) without step-6a reads writes nothing of step 6a; with them (f0 and f1 of both
    variants on T, read as T's set) it writes step 6a's populations and readings and leaves every class of steps 1 and
    2 as it was."""
    import json
    import shutil
    from experiments import E134_population_baseline as B
    shutil.copy(M.LOCK_COPY, tmp_path / M.LOCK_NAME)
    monkeypatch.setitem(M.DIRS, 'results', tmp_path)
    monkeypatch.setitem(M.DIRS, 'local', tmp_path / 'local')
    monkeypatch.setitem(M.DIRS, 'archive', tmp_path / 'archive')
    cache = {}

    def load(tag, step=None):
        setting, panel, fly, w = M.parse_tag(tag)
        if (panel, w) not in cache:
            cache[(panel, w)] = (B.load_run(B.tag_of(B.Point('megamat', 'held', w))) if panel == 'megamat'
                                 else M.lab_apl16_runs([w])[f'lab_apl16_{w}'][2])
        return cache[(panel, w)]

    monkeypatch.setattr(M, 'load_point', load)
    monkeypatch.setattr(M, 'available', lambda pts: pts)
    monkeypatch.setattr(M, 'step6a_sets', lambda g, setting: [('T', M.T_WIRINGS)])
    held = dict(n_on=100, target_on=100, n_off=200, target_off=200, at_bound_thr=False, at_bound_apl=False)
    rows = []
    for f in ('f0', 'f1'):
        rows += [dict(tag=M.point_tag('family', 'megamat', f, w)) for w in M.WIRINGS]
        rows += [dict(tag=M.point_tag('A1', 'megamat', f, w)) for w in M.T_WIRINGS]
    rows = pd.DataFrame([dict(r, **held) for r in rows])
    cl0 = M.part_score(rows)
    before = json.loads((tmp_path / 'classes.json').read_text())
    assert 'step6a' not in cl0 and 'step6a' not in before
    assert not list(tmp_path.glob('score_step6a*')) and 'Step 6a' not in (tmp_path / 'tables.md').read_text()
    six = [dict(tag=M.point_tag(s, 'megamat', f, w), variant=M.SETTING_VARIANT[s], **held)
           for s in M.STEP6A_SETTINGS for f in ('f0', 'f1') for w in M.T_WIRINGS]
    six[0] = dict(six[0], n_on=50)                       # one unheld model fly of 'familymeas', left out and counted
    cl = M.part_score(pd.concat([rows, pd.DataFrame(six)], ignore_index=True))
    after = json.loads((tmp_path / 'classes.json').read_text())
    for k in ('megamat', 'apl16', 'a1_fall', 'a1_apl16', 'recomputed_changed', 'outcomes', 'spreads_model_flies_megamat'):
        assert after[k] == before[k], k
    assert after['step6a'] == json.loads(json.dumps(cl['step6a'], default=M._json_default))
    pop = pd.read_csv(tmp_path / 'score_step6a_populations.csv')
    a = pop[(pop.reading == 'all') & (pop.wirings_set == 'T')]
    assert set(a.population) == set(M.STEP6A_POPULATIONS.values())
    for name, g in a.groupby('population'):
        assert (g.wirings == ' '.join(M.T_WIRINGS)).all() and len(g) == 5
        if name != 'familymeas/megamat':                                   # M lacks one of P3's two model flies
            assert np.allclose(g.population_error, g.lab_common, atol=1e-12)    # every read its wiring's lab run
    m_, r_, p_ = (a[a.population == STEP] for STEP in ('familymeas/megamat', 'family/megamat', 'familyperm/megamat'))
    # one unheld of 8 is more than a ninth: M is read with no class, its unheld left out and counted
    assert (m_.n_unheld == 1).all() and (m_.n_model_flies == 7).all() and (m_.stand_in == 'no class').all()
    assert (r_.n_model_flies == 8).all() and (r_.stand_in == 'a fair stand-in').all()
    assert (p_.n_unheld == 0).all() and (p_.stand_in == 'a fair stand-in').all()
    rd = after['step6a']['readings']['T']
    assert set(rd) == {*M.LOSS_STAGES, 'megamat_stage_error'}
    for st, x in rd.items():
        assert x['structure_or_variance'] == 'not read'                     # R is a fair stand-in everywhere here
        assert x['e_applies'] == (x['repro_model_flies'] == 'vary too little')
        assert x['e_at_upper_bound'] is None and x['measured_stand_in'] == 'no class'
    assert after['step6a']['gates'] is False and after['step6a']['n_unheld'] == {'familymeas': 1, 'familyperm': 0}
    rt = pd.read_csv(tmp_path / 'score_step6a_readings.csv')
    assert set(rt.reading) == {'all', 'held'} and len(rt) == 2 * 5
    assert 'Step 6a' in (tmp_path / 'tables.md').read_text()
    mf = pd.read_csv(tmp_path / 'score_step6a_model_flies.csv')
    assert set(mf.population) == {'familymeas/megamat', 'familyperm/megamat'}
