"""E106 (the gap in ρ on odor-odor distances read by what could close it): the splits' algebra on made-up matrices (the
subsets of pairs, the rank residuals and the shares of the gap by odor, the split between the conditions, the odor part
and the pair part with both hybrids, leave-one-odor-out, the data's values put in), and the ρ columns kcq/apl_effect.py
adds to by_variant.csv (on a made-up results folder)."""
import numpy as np
import pandas as pd
import pytest

from kcq import harness as H
from kcq import apl_effect as AE
from experiments import E106_rho_gap_parts as E

ODORS = list(H.APL16_ODORS)


def pairs_series(values, odors=ODORS) -> pd.Series:
    idx = pd.MultiIndex.from_tuples(H.pair_list(odors), names=['odor_a', 'odor_b'])
    return pd.Series(np.asarray(values, dtype=float), index=idx)


def rand_pairs(seed, n=120) -> pd.Series:
    return pairs_series(np.random.default_rng(seed).uniform(0.1, 0.9, n))


def additive(s: pd.Series, r: pd.Series) -> pd.Series:
    return pd.Series([s[a] + s[b] for a, b in r.index], index=r.index) + r


def orthogonal_rest(seed) -> pd.Series:
    """Random pair values with zero projection onto the odor design."""
    x = rand_pairs(seed) - 0.5
    _, p, _ = E.odor_pair_parts(x)
    return p


def test_subset_sizes_and_public_selection():
    x = rand_pairs(0)
    n = {k: len(E.subset_pairs(x.index, k)) for k in AE.PAIR_SUBSETS}
    assert n == {'all120': 120, 'megamat66': 66, 'v2odors54': 54, 'no_paa_bmyr91': 91}
    assert list(E.select_pairs(x.index, odors=H.MEGAMAT12)) == list(E.subset_pairs(x.index, 'megamat66'))
    assert list(E.select_pairs(x.index, without=('paa', 'B-myr'))) == list(E.subset_pairs(x.index, 'no_paa_bmyr91'))
    # pairs given in either orientation, in any order
    some = [('IaA', '2h'), ('B-myr', 'paa'), ('eb', 'ep')]
    assert set(E.select_pairs(x.index, some)) == {('2h', 'IaA'), ('paa', 'B-myr'), ('eb', 'ep')}
    assert len(E.select_pairs(x.index)) == 120


def test_rho_pairs_and_two_rhos():
    m, d = rand_pairs(1), rand_pairs(2)
    for k in AE.PAIR_SUBSETS:
        idx = E.subset_pairs(m.index, k)
        want = H.spearman(m.loc[idx], d.loc[idx])
        assert E.rho(m, d, k) == pytest.approx(want, abs=1e-12)
    # a model with its keys reversed and shuffled reads the same
    rev = pd.Series(m.to_numpy()[::-1], index=pd.MultiIndex.from_tuples([(b, a) for a, b in m.index[::-1]]))
    assert E.rho_pairs(rev, d, odors=H.MEGAMAT12) == pytest.approx(E.rho(m, d, 'megamat66'), abs=1e-12)
    assert E.rho_pairs(m, d, without=('paa', 'B-myr')) == pytest.approx(E.rho(m, d, 'no_paa_bmyr91'), abs=1e-12)
    pl = list(E.subset_pairs(m.index, 'v2odors54'))
    assert E.rho_pairs(m, d, pl) == pytest.approx(E.rho(m, d, 'v2odors54'), abs=1e-12)
    t = E.two_rhos(m, d, d, m, odors=H.MEGAMAT12)
    assert t['on'] == pytest.approx(t['off'], abs=1e-12)
    # a missing value drops the pair
    m2 = m.copy(); m2.iloc[0] = np.nan
    assert E.rho_pairs(m2, d) == pytest.approx(H.spearman(m.iloc[1:], d.iloc[1:]), abs=1e-12)


def test_rank_residuals_sum_to_the_gap():
    m, d = rand_pairs(3), rand_pairs(4)
    r = E.rank_residuals(m, d)
    assert r['contribution'].sum() == pytest.approx(1 - H.spearman(m, d), abs=1e-12)
    sub = m.loc[E.subset_pairs(m.index, 'v2odors54')]
    r = E.rank_residuals(sub, d)
    assert len(r) == 54 and r['contribution'].sum() == pytest.approx(1 - H.spearman(sub, d), abs=1e-12)


def test_gap_shares_by_odor():
    m, d = rand_pairs(5), rand_pairs(6)
    g = E.gap_shares_by_odor(m, d)
    assert list(g.index) == ODORS and (g['n_pairs'] == 15).all()
    assert g['share_of_gap'].sum() == pytest.approx(2.0, abs=1e-12)
    loo = E.leave_one_odor_out(m, d)
    assert np.allclose(g['loo_rho'].to_numpy(), loo.reindex(g.index).to_numpy())
    assert np.allclose(g['loo_delta'], g['loo_rho'] - E.rho(m, d))
    # over a subset: the 12 megamat odors, 11 pairs each, ranks over the 66 pairs
    g66 = E.gap_shares_by_odor(m, d, odors=H.MEGAMAT12)
    assert list(g66.index) == H.MEGAMAT12 and (g66['n_pairs'] == 11).all()
    assert g66['rho'].iloc[0] == pytest.approx(E.rho(m, d, 'megamat66'), abs=1e-12)
    assert g66['share_of_gap'].sum() == pytest.approx(2.0, abs=1e-12)


def test_condition_split_algebra():
    d_on, d_off = rand_pairs(7), rand_pairs(8)
    # the model is the data: every reading is 1
    s = E.condition_split(d_on, d_off, d_on, d_off)
    for k in ('rho_on', 'rho_off', 'on_from_flies_off_model_effect', 'on_from_model_off_flies_effect',
              'off_from_flies_on_model_effect', 'off_from_model_on_flies_effect'):
        assert s[k] == pytest.approx(1.0), k
    # a model with the flies' exact APL effect on its own APL-off distances: APL on from the flies' APL-off distances and the
    # model's effect is the flies' APL on; from the model's APL-off distances and the flies' effect it is the model's APL on
    m_off = rand_pairs(9)
    m_on = m_off + (d_on - d_off)
    s = E.condition_split(m_on, m_off, d_on, d_off)
    assert s['on_from_flies_off_model_effect'] == pytest.approx(1.0)
    assert s['on_from_model_off_flies_effect'] == pytest.approx(s['rho_on'])
    assert s['off_from_flies_on_model_effect'] == pytest.approx(1.0)
    assert s['off_from_model_on_flies_effect'] == pytest.approx(s['rho_off'])
    assert s['flies_on_vs_flies_off'] == pytest.approx(H.spearman(d_on, d_off))
    # the residual agreement: ρ between the two rank-residual vectors
    r_on = E.rank_residuals(m_on, d_on)['rank_resid']
    r_off = E.rank_residuals(m_off, d_off)['rank_resid']
    assert s['resid_on_vs_off'] == pytest.approx(H.spearman(r_on, r_off))
    # over a subset, the same quantities on that subset's pairs alone
    s54 = E.condition_split(m_on, m_off, d_on, d_off, 'v2odors54')
    assert s54['rho_on'] == pytest.approx(E.rho(m_on, d_on, 'v2odors54'))


def test_odor_pair_parts_recover_an_additive_matrix():
    rng = np.random.default_rng(10)
    s = pd.Series(rng.uniform(0, 0.3, 16), index=ODORS)
    r = orthogonal_rest(11) * 0.2
    x = additive(s, r)
    o, p, s_fit = E.odor_pair_parts(x)
    assert np.allclose((o + p).to_numpy(), x.to_numpy())
    assert np.allclose(s_fit.reindex(ODORS).to_numpy(), s.to_numpy())
    assert np.allclose(p.to_numpy(), r.to_numpy())
    X, _ = H.pair_design(x.index, ODORS)
    assert np.allclose(X.T @ p.to_numpy(), 0.0, atol=1e-10)


def test_odor_pair_hybrids():
    rng = np.random.default_rng(12)
    r = orthogonal_rest(13) * 0.2
    t = pd.Series(rng.uniform(0, 0.3, 16), index=ODORS)
    u = pd.Series(rng.uniform(0, 0.3, 16), index=ODORS)
    data, model = additive(t, r), additive(u, r)
    h = E.odor_pair_hybrids(model, data)
    # the same pair part: the flies' odor part with the model's pair part is the data; the model's odor part with the
    # flies' pair part is the model
    assert h['flies_odor_model_pair'] == pytest.approx(1.0)
    assert h['model_odor_flies_pair'] == pytest.approx(h['model'])
    assert h['pair_parts_rho'] == pytest.approx(1.0)
    assert h['flies_odor_model_pair_best'] >= h['flies_odor_model_pair'] - 1e-12 or h['flies_odor_model_pair_best'] == pytest.approx(1.0)
    # the model is the data: every hybrid is 1
    h = E.odor_pair_hybrids(data, data)
    for k in ('model', 'flies_odor_model_pair', 'model_odor_flies_pair', 'flies_odor_model_pair_matched',
              'model_odor_flies_pair_matched', 'flies_odor_model_pair_quantile', 'model_odor_flies_pair_quantile',
              'odor_shares_rho', 'pair_parts_rho'):
        assert h[k] == pytest.approx(1.0), k
    # a model twice the data: the raw hybrids mix scales, the matched and quantile ones do not
    h = E.odor_pair_hybrids(2 * data, data)
    for k in ('flies_odor_model_pair_matched', 'model_odor_flies_pair_matched', 'flies_odor_model_pair_quantile',
              'model_odor_flies_pair_quantile'):
        assert h[k] == pytest.approx(1.0), k


def test_quantile_map_and_optimal_odor_part():
    x, y = rand_pairs(14), rand_pairs(15) ** 3
    q = E.quantile_map(x, y)
    assert H.spearman(q, x) == pytest.approx(1.0)
    assert np.allclose(np.sort(q.to_numpy()), np.sort(y.to_numpy()))
    # the Pearson-optimal odor part is in the odor design's span and beats other odor parts added to the same pair part
    p = orthogonal_rest(16) * 0.2
    d = rand_pairs(17)
    a = E.optimal_odor_part(p, d)
    _, rest, _ = E.odor_pair_parts(a)
    assert np.allclose(rest.to_numpy(), 0.0, atol=1e-10)
    yr = pd.Series(np.argsort(np.argsort(d.to_numpy())).astype(float), index=d.index)
    best = np.corrcoef((a + p).to_numpy(), yr.to_numpy())[0, 1]
    rng = np.random.default_rng(18)
    for _ in range(50):
        other = additive(pd.Series(rng.normal(0, 0.2, 16), index=ODORS), p * 0)
        assert np.corrcoef((other + p).to_numpy(), yr.to_numpy())[0, 1] <= best + 1e-9


def test_leave_one_odor_out_and_data_values():
    m, d = rand_pairs(19), rand_pairs(20)
    loo = E.leave_one_odor_out(m, d)
    keep = [p for p in m.index if 'paa' not in p]
    assert loo['paa'] == pytest.approx(H.spearman(m[keep], d[keep]))
    # every odor at the data's values: the data
    assert E.rho_pairs(E.at_data_values(m, d, ODORS, 'raw'), d) == pytest.approx(1.0)
    assert E.rho_pairs(E.at_data_values(m, d, ODORS, 'quantile'), d) == pytest.approx(1.0)
    # one odor's pairs only: the other pairs keep the model's values
    x = E.at_data_values(m, d, ['paa', 'B-myr'], 'raw')
    hit = np.array([('paa' in p) or ('B-myr' in p) for p in m.index])
    assert hit.sum() == 29
    assert np.allclose(x[hit], d[hit]) and np.allclose(x[~hit], m[~hit])
    xq = E.at_data_values(m, d, ['paa', 'B-myr'], 'quantile')
    assert np.allclose(xq[~hit], m[~hit])
    assert H.spearman(xq[hit], d[hit]) == pytest.approx(1.0)


def test_dist_rho_columns_and_held():
    m_on, m_off, d_on, d_off = rand_pairs(21), rand_pairs(22), rand_pairs(23), rand_pairs(24)
    c = AE.dist_rho_columns(m_on, m_off, d_on, d_off)
    assert len(c) == 8
    for cond, m, d in (('on', m_on, d_on), ('off', m_off, d_off)):
        for k in AE.PAIR_SUBSETS:
            assert c[f'{AE.V}__dist_rho_{cond}_{k}'] == pytest.approx(E.rho(m, d, k), abs=1e-12)
    assert AE.held('E97_glomerulus_claw_shares', 'B2_mean_real') and AE.held('E97_glomerulus_claw_shares', 'B1_flywire_L_P3')
    assert AE.held('E97_glomerulus_claw_shares', 'B3_perm_malecns_R_s1_M2') and AE.held('E97_glomerulus_claw_shares', 'B5_mean_x0.5_P1')
    # E97's re-plan of 2026-11-17 (the hemibrain and the male CNS): its points are read
    for t in ('B2_mean2_real', 'B1_mcnsL_P1', 'B3_permmean2_s1_M3', 'B3_invmean2_P4', 'B5_mean2_x2_real', 'B6_mean2_F_M1'):
        assert not AE.held('E97_glomerulus_claw_shares', t), t
    assert AE.held('E4_substrate', 'fafb-left-wd12_sp0.05') and not AE.held('E4_substrate', 'uniform_sp0.05_seed3')
    assert not AE.held('E72_wiring_population', 'popfit_real')


def _write_long(folder, tags_pairs):
    rows = []
    for tag, (on, off) in tags_pairs.items():
        for cond, s in (('on', on), ('acute_off', off)):
            for (a, b), v in s.items():
                rows.append(dict(tag=tag, condition=cond, odor_a=a, odor_b=b, distance=v))
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(folder / 'model_pair_distances_long.csv', index=False)


def test_score_all_appends_the_rho_columns(tmp_path, monkeypatch):
    monkeypatch.setattr(AE, 'RESULTS', tmp_path / 'results')
    monkeypatch.setattr(AE, 'ARCHIVE', tmp_path / 'archive')
    d_on, d_off = rand_pairs(30), rand_pairs(31)
    data = H.apl_effect(d_on, d_off)
    v = {'a': (rand_pairs(32), rand_pairs(33)), 'b': (rand_pairs(34), rand_pairs(35))}
    _write_long(tmp_path / 'results' / 'E999_made_up', v)
    _write_long(tmp_path / 'results' / 'E97_glomerulus_claw_shares', {'B2_mean_real': (rand_pairs(36), rand_pairs(37))})
    df = AE.score_all(data, d_on, d_off)
    rho_cols = [c for c in df.columns if '__dist_rho_' in c]
    assert len(rho_cols) == 8 and list(df.columns[-8:]) == rho_cols
    assert 'realism_zcrit' in df.columns and list(df.columns).index('realism_zcrit') < list(df.columns).index(rho_cols[0])
    for _, r in df.iterrows():
        if r.folder == 'E97_glomerulus_claw_shares':
            assert r[rho_cols].isna().all()
        else:
            on, off = v[r.tag]
            want = AE.dist_rho_columns(on, off, d_on, d_off)
            assert all(r[k] == pytest.approx(x, abs=1e-12) for k, x in want.items())
            assert r[f'{AE.V}__apl_effect_rmse'] == pytest.approx(H.apl_effect_stats(on, off, data)['apl_effect_rmse'])


def test_archive_dedupes_and_drops_held():
    cols = E.STAT_COLS + [E.rho_col(c, s) for c in E.CONDITIONS for s in E.SUBSETS]
    rng = np.random.default_rng(40)
    base = pd.DataFrame(rng.uniform(0, 1, (3, len(cols))), columns=cols)
    df = pd.concat([base, base.iloc[[0]], base.iloc[[1]]], ignore_index=True)
    df.insert(0, 'exp', 'E1'); df.insert(1, 'folder', ['E1_x'] * 5); df.insert(2, 'tag', list('abcde'))
    df['realism_flag'] = [np.nan, 'flagged', np.nan, np.nan, 'flagged']
    held = df.iloc[[2]].copy(); held['folder'] = 'E97_glomerulus_claw_shares'; held['tag'] = 'B2_mean_real'
    held[cols] = rng.uniform(0, 1, len(cols))
    a = E.archive(pd.concat([df, held], ignore_index=True))
    assert list(a.tag) == ['a', 'b', 'c']
    assert list(a.n_copies) == [2, 2, 1]
    assert a.copies.iloc[0] == 'E1_x:a ; E1_x:d'
    assert list(a.inside) == [True, False, True]


def test_registry_inputs_are_registry_keys():
    from kcq import data_targets as DT
    for keys in E.REGISTRY_INPUTS.values():
        assert set(keys) <= set(DT.INPUTS)
    for name in ('ctrl_ceiling_54', 'tnt_rec_rho_loo_91', 'family_rho_on_sd_wirings_120', 'family_rho_off_boot_sd_66',
                 'paired_rho_off_sd_e93_fuller_input', 'e85_population_agreement_on_66'):
        assert E.registry_inputs(name)


def test_probe_alternates_name_known_probes():
    assert set(E.PROBE_ALTS) <= set(E.PROBES)


def _made_up_results(d, seed=50):
    """The tables the figures read, with made-up values, in folder `d`."""
    rng = np.random.default_rng(seed)
    u = lambda: float(rng.uniform(0.5, 0.95))
    pd.DataFrame(dict(name=['ctrl_ceiling', 'tnt_ceiling'], value=[0.98, 0.97])).to_csv(d / 'yardsticks.csv', index=False)
    cc = ['rho_on', 'rho_off', 'on_from_flies_off_model_effect', 'on_from_model_off_flies_effect',
          'off_from_flies_on_model_effect', 'off_from_model_on_flies_effect', 'flies_on_vs_flies_off']
    rows = [dict(family=f, reading=r, subset=s, **{c: u() for c in cc}) for f in ('family', 'ceiling')
            for r in ('population', 'real', 'P1') for s in ('all120', 'megamat66')]
    pd.DataFrame(rows).to_csv(d / 'conditions.csv', index=False)
    pc = [f'{k}{r}' for k in ('flies_odor_model_pair', 'model_odor_flies_pair') for r in ('', '_matched', '_quantile', '_best')]
    pc += ['model', 'model_pair_optimal_odor', 'model_odor_r2', 'data_odor_r2']
    rows = [dict(family=f, reading=r, condition=c, **{k: u() for k in pc}) for f in ('family', 'ceiling')
            for r in ('population', 'real', 'P1') for c in E.CONDITIONS]
    pd.DataFrame(rows).to_csv(d / 'odor_pair_parts.csv', index=False)
    rows = [dict(family=f, condition=c, odors=o, rho_without=u(), at_data_raw=u(), at_data_quantile=u())
            for f in ('family', 'ceiling') for c in E.CONDITIONS for o in ('paa and B-myr', 'the four validation2 odors')]
    pd.DataFrame(rows).to_csv(d / 'odors.csv', index=False)
    rows = [dict(input=i, subset=s, rho_on=u(), rho_off=u()) for i in ('orn_pearson', 'input_cosine')
            for s in ('all120', 'megamat66')]
    pd.DataFrame(rows).to_csv(d / 'input.csv', index=False)
    rows = [dict(input=i, condition=c, model_odor_r2=u()) for i in ('orn_pearson', 'input_cosine') for c in E.CONDITIONS]
    pd.DataFrame(rows).to_csv(d / 'input_parts.csv', index=False)
    rows = [dict(family=f, condition=c, odor=o, share_of_gap=float(rng.uniform(0, 0.3))) for f in ('family', 'ceiling')
            for c in E.CONDITIONS for o in ODORS]
    pd.DataFrame(rows).to_csv(d / 'residuals_odors.csv', index=False)


def test_figures_from_the_saved_tables(tmp_path):
    _made_up_results(tmp_path)
    rows = E.gap_figure_rows(tmp_path)
    cond = pd.read_csv(tmp_path / 'conditions.csv')
    cp = cond[(cond.family == 'family') & (cond.reading == 'population') & (cond.subset == 'all120')].iloc[0]
    fam = [r for r in rows if r['label'].startswith('the family (mean matrix')][0]
    assert (fam['on'], fam['off']) == (cp.rho_on, cp.rho_off)
    w = cond[(cond.family == 'family') & (cond.reading != 'population') & (cond.subset == 'all120')]
    assert fam['on_range'] == (w.rho_on.min(), w.rho_on.max())
    hyb = [r for r in rows if r['label'] == "the flies' odor part + the family's pair part"][0]
    pp = pd.read_csv(tmp_path / 'odor_pair_parts.csv').query("family == 'family' and reading == 'population'").set_index('condition')
    vals = [pp.loc['on', f'flies_odor_model_pair{r}'] for r in ('', '_matched', '_quantile', '_best')]
    assert hyb['on'] == vals[0] and hyb['on_range'] == (min(vals), max(vals))
    assert rows[0]['on'] == 0.98 and rows[0]['off'] == 0.97
    for p in (E.figure_gap(tmp_path), E.figure_odor_parts(tmp_path)):
        assert p.exists() and p.stat().st_size > 10_000
