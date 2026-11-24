"""Tests of experiments/E134_archive_score.py on a planted archive and planted recordings (no model output is read)."""
import json
import shutil

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from kcq import harness as H
from kcq import population_fit as P
from kcq import stage_readings as SR
from experiments import E134_archive_score as E

ODORS = list(H.APL16_ODORS)
CRIT = 1.0
QS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


class Gen:
    """Planted recordings of one genotype (as in test_population_fit.py): rectified loadings on odor factors shared by
    every fly, each fly with its own odor strengths, units and noise."""

    def __init__(self, seed, gain=1.0):
        rng = np.random.default_rng(seed)
        self.F = rng.normal(size=(4, len(ODORS)))
        self.s = rng.uniform(0.6, 1.6, size=len(ODORS)) * gain
        self.rng = rng

    def values(self, n_units=300):
        r = self.rng
        s = self.s * np.exp(r.normal(0, 0.15, size=len(self.s)))
        W = r.normal(size=(n_units, 4))
        V = 2 * np.maximum(0, (W @ self.F) * s[None, :] * 0.7 - 0.3) + r.normal(0, 0.5, size=(n_units, len(ODORS)))
        return pd.DataFrame(V, columns=ODORS)

    def degraded(self, q, seed):
        """A fly-like table with a fraction q of its units' responses shuffled across odors."""
        v = self.values().to_numpy(copy=True)
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(v), int(round(q * len(v))), replace=False)
        for i in idx:
            v[i] = rng.permutation(v[i])
        return pd.DataFrame(v, columns=ODORS)


CALLS = []


def flylike(counts, condition='on', curve_counts=None):
    """A planted reading with stage_readings.SomaReading's call (the planted counts are already in the flies' units)."""
    CALLS.append((condition, curve_counts is not None))
    return P.Population.from_criterion(counts, CRIT)


READERS = {'flylike': flylike}
_ARCHIVED_VARIANTS = E.archived_variants   # the real one, which the end-to-end test calls through a patch


@pytest.fixture(scope='module')
def planted(tmp_path_factory):
    root = tmp_path_factory.mktemp('archive')
    g_on, g_off = Gen(1), Gen(2, gain=1.4)
    flies = {'on': [P.Population.from_criterion(g_on.values(), CRIT, name=f'c{i}') for i in range(11)],
             'off': [P.Population.from_criterion(g_off.values(), CRIT, name=f't{i}') for i in range(8)]}
    rows = []
    d = root / 'E901_planted'
    d.mkdir()
    for i, q in enumerate(QS):
        tag = f'q{i}'
        g_on.degraded(q, 10 + i).to_parquet(d / f'{tag}_on_spike_counts.parquet')
        g_off.degraded(q, 50 + i).to_parquet(d / f'{tag}_acute_off_spike_counts.parquet')
        rows.append(dict(exp='E901', folder='E901_planted', tag=tag, realism_flag='' if i != 9 else 'odors per KC, x',
            rate_on=0.07, rate_off=0.15, **{f'{E.AE.V}__apl_effect_rmse_centered': 0.05 + q / 10}))
    # an APL-on file without its APL-off pair, a variant of other odors, and a withheld subfolder
    g_on.values().to_parquet(d / 'lonely_on_spike_counts.parquet')
    g_on.values().iloc[:, :12].to_parquet(d / 'mm_on_spike_counts.parquet')
    g_on.values().iloc[:, :12].to_parquet(d / 'mm_acute_off_spike_counts.parquet')
    w = d / 'withheld_2026-11-15'
    w.mkdir()
    shutil.copy(d / 'q0_on_spike_counts.parquet', w / 'hidden_on_spike_counts.parquet')
    shutil.copy(d / 'q0_acute_off_spike_counts.parquet', w / 'hidden_acute_off_spike_counts.parquet')
    # the same fit saved by another experiment, and a held variant
    c = root / 'E902_copy'
    c.mkdir()
    for n in ('on', 'acute_off'):
        shutil.copy(d / f'q1_{n}_spike_counts.parquet', c / f'q1copy_{n}_spike_counts.parquet')
    h = root / 'E4_substrate'
    h.mkdir()
    for n in ('on', 'acute_off'):
        shutil.copy(d / f'q2_{n}_spike_counts.parquet', h / f'fafb_x_{n}_spike_counts.parquet')
    rows += [dict(exp='E901', folder='E901_planted', tag=t, realism_flag='') for t in ('lonely', 'mm', 'hidden')]
    rows += [dict(exp='E902', folder='E902_copy', tag='q1copy', realism_flag=''),
             dict(exp='E4', folder='E4_substrate', tag='fafb_x', realism_flag='')]
    bv = pd.DataFrame(rows)
    bv.loc[bv[f'{E.AE.V}__apl_effect_rmse_centered'].isna(), f'{E.AE.V}__apl_effect_rmse_centered'] = 0.2
    return dict(root=root, flies=flies, bv=bv)


@pytest.fixture(scope='module')
def read(planted):
    variants = E.archived_variants(planted['bv'], planted['root'])
    refs = {c: P.fly_reference(f, ODORS) for c, f in planted['flies'].items()}
    energy = E.energy_setup(planted['flies'], 'pooled', max_units=150)
    meta, data = E.variant_views(variants, 'flylike', energy=energy, readers=READERS)
    return dict(variants=variants, refs=refs, energy=energy, meta=meta, data=data)


def test_archived_variants(planted, read):
    v = read['variants']
    assert set(v.tag) == {f'q{i}' for i in range(10)} | {'mm', 'q1copy'}
    assert v.set_index('tag').loc['q9', 'inside_bands'] == False   # noqa: E712
    meta = read['meta']
    assert 'mm' not in set(meta.tag) and meta.attrs['skipped'][0]['tag'] == 'mm'
    assert meta.set_index('tag').loc['q1copy', 'duplicate'] and not meta.duplicate[meta.tag != 'q1copy'].any()


def test_scores_order_the_ladder(planted, read):
    sc = E.score_variants(read['meta'], read['data'], planted['flies'], read['refs'], read['energy'])
    assert set(sc.condition) == {'on', 'off'} and len(sc) == 2 * len(read['meta'])
    for cond in ('on', 'off'):
        s = sc[(sc.condition == cond) & (sc.folder == 'E901_planted')].set_index('tag').loc[[f'q{i}' for i in range(10)]]
        assert stats.spearmanr(QS, s['error']).statistic > 0.9
        assert 0.5 < s.loc['q0', 'error'] < 1.6
        assert s.loc['q9', 'error'] > 2 * s.loc['q0', 'error']
        assert stats.spearmanr(QS, s['dist_rho']).statistic < -0.8
        assert stats.spearmanr(QS, s['energy_ratio']).statistic > 0.5
        assert 0.5 < s.loc['q0', 'energy_ratio'] < 2.0


def test_scores_match_the_module(planted, read):
    sc = E.score_variants(read['meta'], read['data'], planted['flies'], read['refs'])
    p = flylike(E.read_counts(planted['root'] / 'E901_planted', 'q3', 'on'))
    direct = P.score(p, read['refs']['on'])
    row = sc[(sc.condition == 'on') & (sc.tag == 'q3')].iloc[0]
    assert row['error'] == pytest.approx(direct.error)
    assert row['part_sharing_who'] == pytest.approx(direct.parts['sharing/who'])
    assert 'energy_ratio' not in sc


def test_consistency(planted, read):
    cons = E.consistency(read['meta'], read['data'], planted['flies'], read['refs'], read['energy'], n_splits=12)
    g = cons.set_index(['condition', 'measure', 'variant_set'])
    for cond in ('on', 'off'):
        assert g.loc[(cond, 'error', 'inside_bands'), 'split_half_rho'] > 0.8
        # q1copy is a duplicate of q1, and q9 is outside the bands
        assert g.loc[(cond, 'error', 'inside_bands'), 'n_variants'] == 9
        assert g.loc[(cond, 'error', 'all'), 'n_variants'] == 10
        for m in ('view_strength', 'view_sharing', 'view_breadth', 'dist_rho', 'rate_rho', 'energy_ratio'):
            assert (cond, m, 'all') in g.index
    kept = E.consistency(read['meta'], read['data'], planted['flies'], read['refs'], None, n_splits=4, units_from='all')
    assert 'energy_ratio' not in set(kept.measure)


def test_tradeoff(planted, read):
    sc = E.score_variants(read['meta'], read['data'], planted['flies'], read['refs'])
    keep = read['meta'][read['meta'].inside_bands & ~read['meta'].duplicate]
    cal = pd.DataFrame(dict(folder=keep.folder, tag=keep.tag, cal_pattern_part=[float(t[1:]) for t in keep.tag]))
    per, best = E.tradeoff(sc, read['meta'], cal)
    assert set(per.tag) == {f'q{i}' for i in range(9)}
    b = best.set_index('measure')
    assert b.loc['cal_pattern_part', 'best_folder_tag'] == 'E901_planted/q0'
    assert b.loc['error_on', 'best_family'] == 'E901_planted'
    assert b.loc['error_on', 'rho_with_cal_pattern_part'] > 0.8
    assert per.loc[per.tag == 'q0', 'rank_cal_pattern_part'].item() == 1


def test_rowwise_spearman():
    rng = np.random.default_rng(0)
    M, d = rng.normal(size=(6, 20)), rng.normal(size=20)
    M[2, 3] = np.nan
    got = E._rowwise_spearman(M, d)
    for i in range(6):
        ok = np.isfinite(M[i])
        assert got[i] == pytest.approx(stats.spearmanr(M[i, ok], d[ok]).statistic)


def test_main_refuses_before_the_lock(monkeypatch):
    def no_lock():
        raise RuntimeError('the E134 lock is not committed')

    def never(*a, **k):
        raise AssertionError('read before the lock')

    monkeypatch.setattr(E, 'lock_committed', no_lock)
    monkeypatch.setattr(E, 'archived_variants', never)
    monkeypatch.setattr(E, 'fly_populations', never)
    with pytest.raises(RuntimeError, match='lock'):
        E.main(['--jobs', '1'])


def test_lock_committed_needs_a_tracked_lock(monkeypatch, tmp_path):
    from data_prep import stage_population_data as S

    def unwritten():
        raise NotImplementedError

    monkeypatch.setattr(S, 'require_lock', unwritten)
    monkeypatch.setattr(S, 'LOCK_COPY', tmp_path / 'E134_lock.json')
    with pytest.raises(RuntimeError, match='not committed'):
        E.lock_committed()
    (tmp_path / 'E134_lock.json').write_text('{}')
    with pytest.raises(RuntimeError, match='not committed'):
        E.lock_committed()


def test_main_end_to_end(planted, monkeypatch, tmp_path):
    bv = tmp_path / 'by_variant.csv'
    planted['bv'].to_csv(bv, index=False)
    monkeypatch.setattr(E, 'lock_committed', lambda: {'planted': True})
    monkeypatch.setattr(E, 'fly_populations', lambda cohort: planted['flies']['on' if cohort == 'control' else 'off'])
    monkeypatch.setattr(E, 'recorded_populations', lambda cohort: as_recorded(
        planted['flies']['on' if cohort == 'control' else 'off'], cohort))
    monkeypatch.setattr(E, 'archived_variants', lambda: _ARCHIVED_VARIANTS(None, planted['root']))
    monkeypatch.setattr(E, 'BY_VARIANT', bv)
    monkeypatch.setattr(E, 'RESULTS', tmp_path / 'results')
    monkeypatch.setattr(E, 'OUTPUTS', tmp_path / 'outputs')
    monkeypatch.setattr(E, 'readings', lambda: {'flylike': flylike})
    E.main(['--reading', 'flylike', '--jobs', '1', '--n-splits', '5', '--no-calibrated', '--energy-max-units', '100'])
    r = tmp_path / 'results'
    sc = pd.read_csv(r / 'archive_scores_flylike.csv')
    assert {'error', 'view_strength', 'part_strength_mean', 'energy_ratio', 'dist_rho', 'rate_rho'} <= set(sc.columns)
    cons = pd.read_csv(r / 'archive_consistency_flylike.csv')
    assert {'on', 'off'} == set(cons.condition) and set(cons.version) == set(E.VERSIONS)
    # the planted reading does not depend on the flies: a half read again is the shared reading
    g = cons.set_index(['version', 'condition', 'measure', 'variant_set'])['split_half_rho']
    for v, shared in zip(E.HALF_VERSIONS, E.SHARED_VERSIONS):
        pd.testing.assert_series_equal(g.loc[v], g.loc[shared])
    assert not np.allclose(g.loc['shared_curve'].loc['on'].to_numpy(),
        g.loc['shared_curve_all_units'].loc['on'].reindex(g.loc['shared_curve'].loc['on'].index).to_numpy(),
        equal_nan=True)
    splits = pd.read_csv(r / 'archive_consistency_splits_flylike.csv')
    assert set(splits.split) == set(range(5)) and set(splits.version) == set(E.VERSIONS)
    drawn = pd.read_csv(r / 'archive_consistency_draw_flylike.csv', dtype={'tag': str})
    assert len(drawn) == 10
    per = pd.read_csv(r / 'archive_tradeoff_flylike.csv')
    assert per['cal_pattern_part'].isna().all() and per['error_on'].notna().all()
    prov = json.loads((r / 'archive_provenance_flylike.json').read_text())
    assert prov['n_read'] == 11 and prov['n_distinct'] == 10 and prov['n_skipped'] == 1
    sh = prov['split_half']
    assert sh['n_splits'] == 5 and sh['halves_read_with_own_reading'] is False and sh['draw']['n_drawn'] == 10
    assert sh['reduced'] and not sh['meets_lock_minimum'] and sh['full_run_estimate_hours'] >= 0
    views = pd.read_parquet(tmp_path / 'outputs' / 'archive_views_flylike.parquet')
    assert set(views.part) == set(P.ALL_PARTS) - {'matrix/values'}
    E.main(['--reading', 'flylike', '--jobs', '1', '--n-splits', '3', '--no-calibrated', '--no-energy', '--parts',
        'probe'])
    probe = pd.read_csv(r / 'archive_scores_flylike_probe.csv')
    assert 'part_strength_mean' not in probe and 'energy_ratio' not in probe
    assert (probe['view_strength'] == probe['part_strength_rate']).all()
    # the split-half test alone, drawn, to another directory
    out = tmp_path / 'scratch'
    E.main(['--reading', 'flylike', '--jobs', '1', '--n-splits', '3', '--no-energy', '--consistency-only',
        '--draw-variants', '4', '--out', str(out)])
    assert sorted(x.name for x in out.iterdir()) == ['archive_consistency_draw_flylike.csv',
        'archive_consistency_flylike.csv', 'archive_consistency_splits_flylike.csv',
        'archive_provenance_consistency_flylike.json']
    d = pd.read_csv(out / 'archive_consistency_draw_flylike.csv')
    assert d.inside_bands.sum() == 4 and len(d) == 4   # the one variant outside: round(4 x 1 / 9) = 0


def test_calibrated_pattern_parts(planted, tmp_path):
    try:
        H.yang_apl_effect()
    except Exception as ex:   # the KC data are not on this machine
        pytest.skip(f'no KC data: {ex!r}')
    rng = np.random.default_rng(3)
    d = tmp_path / 'E903'
    d.mkdir()
    on = pd.DataFrame(rng.poisson(0.4, size=(400, 16)) * rng.integers(0, 2, size=(400, 1)), columns=ODORS).astype(float)
    off = on * 2 + rng.poisson(0.3, size=on.shape)
    on.to_parquet(d / 'v_on_spike_counts.parquet')
    off.to_parquet(d / 'v_acute_off_spike_counts.parquet')
    v = pd.DataFrame([dict(folder='E903', tag='v', archive_dir=str(d))])
    got = E.calibrated_pattern_parts(v, jobs=1, target=3.0).iloc[0]
    K = H.calibrate_readout(on, 'saturating', 3.0, stat='participation_ratio', responders_min=2.0)
    s = H.apl_effect_stats_readout(on, off, 'saturating', K)
    assert got['K'] == pytest.approx(K) and got['cal_pattern_part'] == pytest.approx(s['apl_effect_rmse_centered'])
    # the worker pool (its initializer must pickle for the forkserver start method; 2026-11-24 it did not)
    v2 = pd.concat([v, v.assign(tag='v')], ignore_index=True)
    got2 = E.calibrated_pattern_parts(v2, jobs=2, target=3.0)
    assert 'error' not in got2.columns or got2['error'].isna().all()
    assert got2['cal_pattern_part'].tolist() == pytest.approx([s['apl_effect_rmse_centered']] * 2)


def test_readings_are_stage_readings_soma_readings():
    r = E.readings()
    assert set(r) == {x.name for x in SR.readings_of('kc_soma')}
    assert all(isinstance(v, SR.SomaReading) for v in r.values())
    assert E.DEFAULT_READING == SR.DEFAULT_READINGS['kc_soma'].name


def test_off_is_read_with_the_on_counts(planted):
    v = E.archived_variants(planted['bv'], planted['root']).iloc[:2]
    CALLS.clear()
    E.variant_views(v, 'flylike', readers=READERS)
    assert CALLS == [('on', False), ('off', True)] * 2


def as_recorded(pops, cohort, trials=None):
    """The data module's recorded population of planted populations (P.Population), with planted single trials."""
    from data_prep import stage_population_data as D
    names = tuple(p.name for p in pops)
    return D.Population(stage='kc_soma', panel='apl16', cohort=cohort, odors=tuple(ODORS), flies=names,
        tables={p.name: p.values[ODORS] for p in pops}, criterion=CRIT, units='planted', trials=trials,
        has_trials=trials is not None)


def planted_recorded(gen, cohort, n, seed, n_trials=3):
    """Planted recordings with single trials (trial = the fly's mean + noise), as the data module records them."""
    rng = np.random.default_rng(seed)
    pops, trials = [], {}
    cols = pd.MultiIndex.from_tuples([(o, r) for o in ODORS for r in range(n_trials)], names=['odor', 'repeat'])
    for i in range(n):
        m = gen.values(200).to_numpy()
        T = m[:, :, None] + rng.normal(0, 0.6, size=m.shape + (n_trials,))
        name = f'{cohort[0]}{i}'
        trials[name] = pd.DataFrame(T.reshape(len(m), -1), columns=cols)
        pops.append(P.Population.from_criterion(pd.DataFrame(T.mean(axis=2), columns=ODORS), CRIT, name=name))
    return as_recorded(pops, cohort, trials), pops


@pytest.fixture(scope='module')
def soma(planted):
    """The archive read with stage_readings' default reading against planted recordings with trials."""
    rec_on, on = planted_recorded(Gen(5), 'control', 8, 1)
    rec_off, off = planted_recorded(Gen(6, gain=1.4), 'tnt', 6, 2)
    recorded, flies = {'on': rec_on, 'off': rec_off}, {'on': on, 'off': off}
    refs = {c: P.fly_reference(f, ODORS) for c, f in flies.items()}
    variants = E.archived_variants(planted['bv'], planted['root'])
    shared = SR.SomaReading('curve_noise', flies=recorded)
    meta, data = E.variant_views(variants, 'curve_noise', readers={'curve_noise': shared})
    sel = E.draw_variants(meta, 0)
    splits = E.paired_splits({c: len(r.flies) for c, r in refs.items()}, 3, seed=7)
    own = E.half_split_measures(sel, recorded, flies, refs, 'curve_noise', splits)
    return dict(recorded=recorded, flies=flies, refs=refs, meta=meta, data=data, sel=sel, splits=splits, own=own)


def test_paired_splits_are_the_first_runs():
    sp = E.paired_splits({'on': 11, 'off': 8}, 6, seed=3)
    on, off = P.split_halves(11, 6, 3), P.split_halves(8, 6, 4)
    for k, s in enumerate(sp):
        assert list(s) == ['on', 'off']
        for h in (0, 1):
            assert (s['on'][h] == on[k][h]).all() and (s['off'][h] == off[k][h]).all()


def test_draw_variants():
    meta = pd.DataFrame(dict(key=[f'v{i}' for i in range(40)], duplicate=[i % 10 == 9 for i in range(40)],
        inside_bands=[i % 4 != 0 for i in range(40)]))
    d = E.draw_variants(meta, 12, seed=1)
    distinct = meta[~meta.duplicate]
    n_in, n_out = int(distinct.inside_bands.sum()), int((~distinct.inside_bands).sum())
    assert d.inside_bands.sum() == 12 and (~d.inside_bands).sum() == round(12 * n_out / n_in)
    assert not d.duplicate.any() and list(d.index) == sorted(d.index)
    pd.testing.assert_frame_equal(d, E.draw_variants(meta, 12, seed=1))
    assert not d.index.equals(E.draw_variants(meta, 12, seed=2).index)
    pd.testing.assert_frame_equal(E.draw_variants(meta, 0), distinct)
    pd.testing.assert_frame_equal(E.draw_variants(meta, 1000), distinct)


def test_a_half_is_read_with_its_own_curve_noise_and_units(planted, soma):
    measures, columns = soma['own']
    sel, refs, rec = soma['sel'], soma['refs'], soma['recorded']
    assert set(measures) == set(E.HALF_VERSIONS)
    assert measures['half']['off'].shape[:3] == (len(sel), 3, 2)
    i = list(sel.tag).index('q3')
    on = E.read_counts(planted['root'] / 'E901_planted', 'q3', 'on')
    off = E.read_counts(planted['root'] / 'E901_planted', 'q3', 'acute_off')
    k, h = 1, 1
    idx = {c: soma['splits'][k][c][h] for c in refs}
    names = {c: [refs[c].flies[j] for j in idx[c]] for c in refs}
    halves = {c: SR.subset_flies(rec[c], names[c]) for c in refs}
    reader = SR.SomaReading('curve_noise', flies=halves)
    pops = {'on': reader(on), 'off': reader(off, condition='off', curve_counts=on)}
    for cond, pop in pops.items():
        cols = columns['half'][cond]
        own = P.score(pop, P.subset_reference(refs[cond], idx[cond]))
        kept = P.score(pop, P.subset_reference(refs[cond], idx[cond], keep_units=True))
        assert measures['half'][cond][i, k, h, cols.index('error')] == pytest.approx(own.error, rel=1e-10)
        assert measures['all_units'][cond][i, k, h, cols.index('error')] == pytest.approx(kept.error, rel=1e-10)
        assert measures['half'][cond][i, k, h, cols.index('part_sharing_who')] == pytest.approx(
            own.parts['sharing/who'], rel=1e-10)
        fr, fd = E.fly_summaries(soma['flies'][cond])
        rho = stats.spearmanr(pop.responds.mean(axis=0), np.nanmean(fr[idx[cond]], axis=0)).statistic
        assert measures['half'][cond][i, k, h, cols.index('rate_rho')] == pytest.approx(rho)
    # the halves are paired: APL off read with the other control half's curve is another population
    other = SR.subset_flies(rec['on'], [refs['on'].flies[j] for j in soma['splits'][k]['on'][1 - h]])
    swapped = SR.SomaReading('curve_noise', flies={'on': other, 'off': halves['off']})(off, condition='off',
        curve_counts=on)
    assert not np.allclose(swapped.values.to_numpy(), pops['off'].values.to_numpy())
    # and the half's own reading is not the one set on all recordings
    shared = SR.SomaReading('curve_noise', flies=rec)(on)
    assert not np.allclose(shared.values.to_numpy(), pops['on'].values.to_numpy())


def test_split_half_consistency_versions(soma):
    cons, per, own = E.split_half_consistency(soma['sel'], soma['data'], soma['flies'], soma['refs'],
        soma['recorded'], 'curve_noise', soma['splits'])
    assert own and list(dict.fromkeys(cons.version)) == list(E.VERSIONS)
    g = cons.set_index(['version', 'condition', 'measure', 'variant_set'])
    for v in E.VERSIONS:
        for cond in ('on', 'off'):
            assert g.loc[(v, cond, 'error', 'inside_bands'), 'n_splits'] == 3
            assert g.loc[(v, cond, 'error', 'inside_bands'), 'n_variants'] == 9
            assert g.loc[(v, cond, 'error', 'inside_bands'), 'split_half_rho'] > 0.5
    assert 'energy_ratio' not in set(cons.measure)
    m = per.groupby(['version', 'condition', 'measure', 'variant_set'])['rho'].mean()
    assert np.allclose(m.reindex(g.index).to_numpy(), g['split_half_rho'].to_numpy(), equal_nan=True)
    # the shared-curve versions are consistency()'s on the same splits
    for units, v in (('half', 'shared_curve'), ('all', 'shared_curve_all_units')):
        c = E.consistency(soma['sel'], soma['data'], soma['flies'], soma['refs'], splits=soma['splits'],
            units_from=units).set_index(['version', 'condition', 'measure', 'variant_set'])
        assert np.allclose(c['split_half_rho'], g.loc[c.index, 'split_half_rho'], equal_nan=True)


def test_raw_reads_a_half_as_all_recordings(soma):
    assert E.half_reading('raw', {}, {}) is None and E.half_reading('flylike', {}, {}) is None
    half = {'on': ['c0', 'c3'], 'off': ['t1', 't2']}
    r = E.half_reading('curve_noise', soma['recorded'], half)
    assert isinstance(r, SR.SomaReading)
    assert {c: tuple(r.flies(c).flies) for c in half} == {c: tuple(v) for c, v in half.items()}
    assert E.half_split_measures(soma['sel'], soma['recorded'], soma['flies'], soma['refs'], 'raw',
        soma['splits']) is None


def test_half_split_measures_in_worker_processes(soma):
    sel = soma['sel'].iloc[:3]
    one = E.half_split_measures(sel, soma['recorded'], soma['flies'], soma['refs'], 'curve_noise', soma['splits'][:2])
    two = E.half_split_measures(sel, soma['recorded'], soma['flies'], soma['refs'], 'curve_noise', soma['splits'][:2],
        jobs=2)
    for v in E.HALF_VERSIONS:
        for cond in ('on', 'off'):
            assert one[1][v][cond] == two[1][v][cond]
            assert np.allclose(one[0][v][cond], two[0][v][cond], equal_nan=True)
            assert np.allclose(one[0][v][cond], soma['own'][0][v][cond][:3, :2], equal_nan=True)
