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
