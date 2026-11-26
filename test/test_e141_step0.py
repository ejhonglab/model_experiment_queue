#!/usr/bin/env python3
"""Tests for E141's step 0 (experiments/E141_step0.py): the KC-KC bins and A1's evaluations binned by them on made-up
histories, the strength split's rows on made-up populations (each wiring's row is E136's `strength_split`, the summary
rows are means of the wirings', the split adds up to the part's error), the lead table's pivot, what the holds reach on
a made-up run (parameters, rates by KC type, odors per KC, bounds, E90's parameters from the wiring's fit), the halves'
cache; with the shared directory, the somata of one variant reproduce E139's saved score."""

import dataclasses

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
S = pytest.importorskip('experiments.E141_step0')
F = pytest.importorskip('experiments.E136_stage_population_fit')
P = pytest.importorskip('kcq.population_fit')
R = pytest.importorskip('kcq.stage_readings')


def _shared() -> bool:
    try:
        H.data_targets_path().read_bytes()
        return (H.results_dir('E139_normalization_ladder') / 'score_panels.csv').exists()
    except (FileNotFoundError, OSError):
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='the shared directory is not available')


def test_kckc_bins():
    f = np.array([0.5, 1.0, 1.2, 2.0, 3.24, 3.25, 7.3])
    idx, edges = S.kckc_bins(f, 6)
    assert len(edges) == 8 and edges[0] == pytest.approx(0.5) and edges[1] == pytest.approx(1.0)
    assert edges[-1] == pytest.approx(7.3)
    assert np.allclose(np.diff(np.log(edges[1:])), np.log(7.3) / 6)
    assert idx[0] == 0 and idx[1] == 1 and idx[-1] == 6
    assert all(edges[i] <= x <= edges[i + 1] for i, x in zip(idx, f))
    idx2, edges2 = S.kckc_bins(np.array([1.0, 2.0, 4.0]), 2)
    assert len(edges2) == 3 and list(idx2) == [0, 1, 1]
    with pytest.raises(ValueError):
        S.kckc_bins(np.array([0.5, 1.0]))


def test_a1_by_kckc(tmp_path):
    zk = np.log(np.array([1.0, 3.2, 3.3, 7.3, 0.5]))
    hist = pd.DataFrame(dict(eval=range(5), z_aplpn=[0.0, -1.86, -1.8, -3.0, 0.0], z_kckc=zk,
                             z_inhadd=[0.0, 3.47, 3.4, 4.0, 0.0], z_inhsc=[0.0, 0.09, 0.1, -0.3, 0.0]))
    hist.to_csv(tmp_path / 'fit_history_A1.csv', index=False)
    rows = []
    for i in range(5):
        for half, add in (('fit', 0.0), ('held', 0.05)):
            v = 1.0 - 0.01 * i + add - (0.02 if i == 2 else 0.0)
            rows.append(dict(arm='A1', eval=i, measure='strength', half=half, stage_error=v,
                             stage_error_downstream=v + 0.1, pn_dendrite_error=v, pn_bouton_error=v, kc_claw_error=v,
                             kc_soma_error=v + 0.2))
            rows.append(dict(rows[-1], measure='all', stage_error=9.0))
    pd.DataFrame(rows).to_csv(tmp_path / 'rescore_A1.csv', index=False)
    t = S.a1_by_kckc(tmp_path, n_bins=3)
    assert set(t.columns) >= {'bin', 'bin_lo', 'bin_hi', 'bin_center', 'half', 'n_evals', 'statistic', 'eval',
                              'kckc_factor', 'aplpn_factor', 'inhsc_gain', 'pn_inhadd', 'stage_error',
                              'stage_error_downstream', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma'}
    best = t[t.statistic == 'best'].set_index(['bin', 'half'])
    # bins: [0.5, 1), [1, 1.94), [1.94, 3.76), [3.76, 7.3]
    assert sorted(set(t.bin)) == [0, 1, 2, 3]
    assert best.loc[(2, 'fit'), 'eval'] == 2 and best.loc[(2, 'held'), 'eval'] == 2
    assert best.loc[(2, 'fit'), 'n_evals'] == 2
    assert best.loc[(2, 'held'), 'stage_error'] == pytest.approx(1.0 - 0.02 + 0.05 - 0.02)
    assert best.loc[(2, 'fit'), 'kckc_factor'] == pytest.approx(3.3)
    assert best.loc[(3, 'fit'), 'kckc_factor'] == pytest.approx(7.3)
    assert best.loc[(0, 'fit'), 'bin_hi'] == pytest.approx(1.0)
    med = t[(t.statistic == 'median') & (t.bin == 2) & (t.half == 'fit')].iloc[0]
    assert med.stage_error == pytest.approx(np.median([0.99, 0.96]))
    assert med.kc_soma == pytest.approx(np.median([0.99, 0.96]) + 0.2)
    assert med.kckc_factor == pytest.approx(np.exp(np.median(zk[1:3])))


def _pop(rates, n_units=40, seed=0, name=''):
    """A made-up population: each odor's response rate `rates` (units respond at random), graded values."""
    rng = np.random.default_rng(seed)
    odors = [f'o{i}' for i in range(len(rates))]
    resp = pd.DataFrame((rng.random((n_units, len(rates))) < np.asarray(rates)).astype(float), columns=odors)
    vals = resp * (1.0 + rng.random(resp.shape)) + 0.01 * rng.random(resp.shape)
    return P.Population(vals, resp, name=name)


def _flies_and_wirings(n_odors=6):
    rng = np.random.default_rng(1)
    base = np.linspace(0.1, 0.5, n_odors)
    flies = [_pop(np.clip(base + 0.05 * rng.standard_normal(n_odors), 0.02, 0.9), seed=10 + i, name=f'f{i}')
             for i in range(4)]
    wirings = list(F.WIRINGS)
    pops = [_pop(np.clip(base[::-1] * (1 + 0.1 * i), 0.02, 0.9), seed=20 + i, name=f'C2/{w}')
            for i, w in enumerate(wirings)]
    odors = list(flies[0].values.columns)
    ref = P.fly_reference(flies, odors, P.APL16_BREADTH_EDGES, parts=['strength/rate', 'strength/mean'])
    return pops, ref, wirings


def test_split_rows():
    pops, ref, wirings = _flies_and_wirings()
    df = S.split_rows(pops, ref, wirings, F.S_WIRINGS, F.T_WIRINGS)
    assert list(df.columns) == ['wiring', 'statistic', 'odor', 'value']
    panel = df[df.odor.isna()].pivot(index='wiring', columns='statistic', values='value')
    assert set(panel.index) == set(wirings) | set(S.SUMMARY_ROWS)
    sp = F.strength_split(S._renamed(pops, wirings), ref)
    for w in wirings + ['common']:
        for c in sp.columns:
            assert panel.loc[w, c] == pytest.approx(sp.loc[w, c], nan_ok=True)
    for c in ('strength_rate_sd_ratio', 'strength_rate_spread', 'error', 'part_strength_rate'):
        assert panel.loc['mean_T', c] == pytest.approx(panel.loc[list(F.T_WIRINGS), c].mean())
        assert panel.loc['mean_S', c] == pytest.approx(panel.loc[list(F.S_WIRINGS), c].mean())
        assert panel.loc['mean', c] == pytest.approx(panel.loc[wirings, c].mean())
    # level^2 + spread^2 + order^2 is the part's error^2, on every wiring and the wirings' mean views
    for w in wirings + ['common']:
        for p in ('rate', 'mean'):
            s2 = sum(panel.loc[w, f'strength_{p}_{t}'] ** 2 for t in ('level', 'spread', 'order'))
            assert s2 == pytest.approx(panel.loc[w, f'part_strength_{p}'] ** 2, rel=1e-9)
    # the spread as a multiple of the flies': the SD over odors of the model's rate entries over the flies' mean's
    v = P.population_views(S._renamed(pops, wirings)[4], ref.odors, ref.breadth_edges, list(ref.parts))
    assert panel.loc[wirings[4], 'strength_rate_sd_ratio'] == pytest.approx(
        float(np.std(v['strength/rate'].reindex(ref.mean['strength/rate'].index)))
        / float(np.std(ref.mean['strength/rate'])))
    odor = df[df.odor.notna()]
    assert set(odor.statistic) == set(S.ODOR_STATS)
    assert set(odor.wiring) == set(wirings) | set(S.ODOR_SUMMARY_ROWS)
    assert set(odor.odor) == set(ref.odors)
    sh = odor[(odor.statistic == 'view_share') & (odor.wiring == 'common')].value
    assert sh.sum() == pytest.approx(1.0)
    e = odor[(odor.statistic == 'rate_e') & (odor.odor == ref.odors[0])].set_index('wiring').value
    assert e['mean_T'] == pytest.approx(e[list(F.T_WIRINGS)].mean())
    short = S.split_rows(pops, ref, wirings, F.S_WIRINGS, F.T_WIRINGS, odor_entries=False)
    assert len(short) == len(df[df.odor.isna()])


def test_spread_table():
    rows = []
    for v, x in (('family', 1.22), ('A1', 1.16), ('K_inhadd1000', 1.24)):
        for w, add in (('mean_T', 0.0), ('mean', 0.5)):
            rows.append(dict(variant=v, input='validation2', panel='validation2', condition='on', half='all', wiring=w,
                             statistic='strength_rate_sd_ratio', odor=None, value=x + add))
        rows.append(dict(rows[-1], wiring='mean_T', odor='paa', value=9.0))
        rows.append(dict(rows[-1], odor=None, half='fit', value=8.0))
    t = S.spread_table(pd.DataFrame(rows))
    assert list(t.columns) == ['family', 'K_inhadd1000', 'A1']
    assert t.loc[('validation2', 'validation2', 'on'), 'A1'] == pytest.approx(1.16)


def _run(source='E134:step3', info=None, n_kc=12, odors=('2h', 'IaA', 'pa', 'eb'), seed=0):
    rng = np.random.default_rng(seed)
    types = np.array(['g', "a'b'", 'ab'] * (n_kc // 3), dtype=object)
    idx = pd.MultiIndex.from_arrays([np.arange(n_kc), types], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    cols = [f'{o} @ -3' for o in odors]
    conds = {}
    for c, lam in (('on', 1.0), ('off', 2.5)):
        cnt = pd.DataFrame(rng.poisson(lam, (n_kc, len(odors))).astype(float), index=idx, columns=cols)
        conds[c] = dict(spike_counts=cnt, responses=(cnt >= 2).astype(float))
    return R.SavedRun(source=source, tag='v_megamat_P3', model='C2', wiring='P3', input_panel='megamat',
                      conditions=conds, info=dict(info or {}))


def test_panel_holds():
    run = _run()
    r, c = run.conditions['on']['responses'], run.conditions['on']['spike_counts']
    h = S.panel_holds(r, c, ['2h', 'IaA', 'eb'])
    sub = (c.iloc[:, [0, 1, 3]].to_numpy() >= 2)
    assert h['n_odors'] == 3 and h['n_kcs'] == 12 and h['rate'] == pytest.approx(sub.mean())
    types = r.index.get_level_values(H.KC_TYPE_LEVEL)
    assert h['rate_apbp'] == pytest.approx(sub[np.asarray(types == "a'b'")].mean()) and h['n_kcs_g'] == 4
    assert h['max_count'] == c.iloc[:, [0, 1, 3]].to_numpy().max()
    n = sub.sum(axis=1)
    assert [h[f'opk_{b}'] for b in ('0', '1', '2_3', '4_7', '8plus')] == pytest.approx(
        [np.mean(n == 0), np.mean(n == 1), np.mean((n >= 2) & (n <= 3)), 0.0, 0.0])
    assert sum(h[f'opk_{b}'] for b in ('0', '1', '2_3', '4_7', '8plus')) == pytest.approx(1.0)


def test_run_holds(monkeypatch):
    prm = dict(fixed_thr=50.0, wAPLKC_scale=1.0, wKCAPL_scale=1e-3, wAPLPN_scale=0.05, wPNAPL_scale=2e-3,
               wKCKC_scale=100.0)
    off = dict(prm, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    run = _run(info=dict(params=prm, params_off=off))
    monkeypatch.setattr(R, 'run_panels', lambda run, stage: [('megamat', ('2h', 'IaA', 'pa', 'eb'))])
    reads = {'E139': {'v_megamat_P3': dict(at_bound_thr='False', at_bound_apl='True')}}
    rows = S.run_holds(F, run, 'megamat', (0.4, 0.2), reads)
    assert [r['condition'] for r in rows] == ['on', 'off']
    on, offr = rows
    assert on['fixed_thr'] == 50.0 and on['wKCKC_scale'] == 100.0 and on['wAPLKC_scale'] == 1.0
    assert offr['wAPLKC_scale'] == 0.0 and offr['wKCKC_scale'] == 100.0
    assert on['frac_target'] == 0.2 and offr['frac_target'] == 0.4
    assert on['rate_held'] == pytest.approx(run.conditions['on']['responses'].to_numpy().mean())
    assert on['hold_miss'] == pytest.approx(abs(on['rate_held'] - 0.2) / 0.2)
    assert on['at_bound_thr'] == 0.0 and on['at_bound_apl'] == 1.0
    assert 'opk_tv' not in on
    # the APL16 panel: the odors-per-KC distance from the flies beside its band
    from kcq import apl_effect as A
    monkeypatch.setattr(R, 'run_panels', lambda run, stage: [('apl16', ('2h', 'IaA', 'pa', 'eb'))])
    monkeypatch.setattr(A, 'opk_data', lambda c, units='z': np.array([0.5, 0.2, 0.2, 0.1, 0.0]))
    monkeypatch.setattr(A, 'opk_band', lambda c, units='z': 0.15 if c == 'on' else 0.22)
    rows = S.run_holds(F, run, 'apl16', None, {})
    d = np.array([rows[0][f'opk_{b}'] for b in A.OPK_BINS])
    assert rows[0]['opk_tv'] == pytest.approx(0.5 * np.abs(d - [0.5, 0.2, 0.2, 0.1, 0.0]).sum())
    assert rows[0]['opk_tv_band'] == 0.15 and rows[1]['opk_tv_band'] == 0.22
    assert np.isnan(rows[0]['frac_target']) and np.isnan(rows[0]['at_bound_thr'])


def test_params_and_bounds(monkeypatch):
    from experiments import E83_cross_wiring_fit as E83
    p = dict(fixed_thr=77.9, wAPLKC_scale=0.2, wKCAPL_scale=1e-4, wAPLPN_scale=0.36, wPNAPL_scale=5e-4,
             wKCKC_scale=42.7)
    monkeypatch.setattr(E83, 'saved_fit', lambda w: dict(p=p))
    pr = S._params(F, _run(source='E90:B'))
    assert pr['on'] == pytest.approx(p)
    assert pr['off']['wAPLKC_scale'] == 0.0 and pr['off']['wAPLPN_scale'] == 0.0 and pr['off']['fixed_thr'] == 77.9
    assert pr['off']['wKCKC_scale'] == 42.7
    assert all(np.isnan(x) for x in S._params(F, _run(source='other'))['on'].values())
    assert S._bounds(_run(info=dict(at_bound_thr=True, at_bound_apl=False)), {}) == dict(at_bound_thr=1.0,
                                                                                         at_bound_apl=0.0)
    hold = dict(off=dict(at_end=False), on=dict(at_end=True))
    assert S._bounds(_run(info=dict(hold=hold)), {}) == dict(at_bound_thr=0.0, at_bound_apl=1.0)
    b = S._bounds(_run(), {})
    assert np.isnan(b['at_bound_thr']) and np.isnan(b['at_bound_apl'])


@dataclasses.dataclass
class _Flies:
    flies: tuple
    tables: dict
    trials: dict = None


class _Base:
    def __init__(self, pops):
        self.pops = pops

    def get(self, key, default=None):
        return self.pops.get(key, default)


def test_half_flies(monkeypatch):
    fl = _Flies(flies=('a', 'b', 'c', 'd'), tables={k: k for k in 'abcd'})
    monkeypatch.setattr(S, 'panel_half_names', lambda f: dict(fit=['b', 'd'], held=['a', 'c']))
    h = S.HalfFlies(_Base({('kc_soma', 'x', 'on'): fl}), 'fit')
    sub = h.get(('kc_soma', 'x', 'on'))
    assert sub.flies == ('b', 'd') and set(sub.tables) == {'b', 'd'}
    assert h.get(('kc_soma', 'x', 'on')) is sub
    assert h.get(('kc_soma', 'x', 'off')) is None and h.get(('kc_soma', 'x', 'off'), 0) == 0
    with pytest.raises(ValueError):
        S.HalfFlies(_Base({}), 'all')


@needs_shared
def test_soma_splits_reproduce_e139():
    # reads one variant's eight saved validation2 runs (about 5 s), so it runs under QUICK too
    sp = pd.read_csv(H.results_dir('E139_normalization_ladder') / 'score_panels.csv')
    sp = sp[(sp.arm == 'K_inhadd1000') & (sp.input == 'validation2')]
    df = S.soma_splits(variants=('K_inhadd1000',), halves=('all',), inputs=('validation2',), odor_entries=False)
    mine = df.set_index(['wiring', 'statistic']).value
    n = 0
    for r in sp.itertuples():
        for c in [c for c in sp.columns if c.startswith('split_')]:
            if np.isfinite(getattr(r, c)):
                assert mine[(r.wiring, c[len('split_'):])] == pytest.approx(getattr(r, c), abs=1e-12)
                n += 1
    assert n > 50
