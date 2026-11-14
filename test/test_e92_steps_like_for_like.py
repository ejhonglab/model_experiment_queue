#!/usr/bin/env python3
"""Tests for E92 part A (experiments/E92_steps_like_for_like.py): the groups have their sizes; Pearson by centring on
synthetic data (the distances, the noise-removed estimate on centred trials, the expected distances through noise
against Monte Carlo and against `zcell_pair_distances`); the known-truth check on synthetic pseudo-recordings with a
share of the noise common to a repeat's odors, and on the real wiring's run (slow, skipped with QUICK=1); the verdict
functions and the line on made-up tables; and on the shared data: the flies' trial-mean rows equal E90's
stages_flies.csv, the cosine noise-removed values equal E81's noise_free_pairs.csv, the model's calibrated-readout rows
equal E90's stages_model.csv (and the linear map in baseline-noise units without noise is raw counts)."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E92_steps_like_for_like')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return (H.shared_dir() / 'results').exists()
    except Exception:
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='needs the shared directory')


def test_groups_have_their_sizes():
    m = X.group_masks(X.PAIRS66)
    sizes = {g: int(v.sum()) for g, v in m.items()}
    assert sizes == dict(ester=15, alcohol=6, ester_alcohol=24, other=21, near=6, all66=66, near_ester=4,
                         near_15ol_16ol=1, near_6al_t2h=1)
    assert not (m['near_ester'] & ~m['near']).any()
    m120 = X.group_masks(X.PAIRS120)
    assert int(m120['all66'].sum()) == 66 and int(m120['ester'].sum()) == 15


def _structured(rng, n_units=1500, n_odors=16, k=4):
    """Non-negative, sparse, correlated responses (units x odors)."""
    base = rng.standard_normal((n_units, k)) @ rng.standard_normal((k, n_odors))
    return np.maximum(base + 0.5 * rng.standard_normal((n_units, n_odors)), 0.0) + 0.2


def test_pearson_by_centring_synthetic():
    rng = np.random.default_rng(0)
    Xm = _structured(rng)
    np.testing.assert_allclose(X.distance_matrix(Xm, 'pearson'), 1 - np.corrcoef(Xm.T), atol=1e-12)
    n = np.linalg.norm(Xm, axis=0)
    np.testing.assert_allclose(X.distance_matrix(Xm, 'cosine'), 1 - (Xm.T @ Xm) / np.outer(n, n), atol=1e-12)
    # the noise-removed estimate on centred trials: single-trial Pearson r between different repeats
    T = Xm[None] + 0.7 * rng.standard_normal((3, *Xm.shape))
    est = X.noise_removed(T, 'pearson')
    o = list(X.ODORS16)
    r = lambda a, b: np.corrcoef(a, b)[0, 1]
    for a, b in [('2h', 'IaA'), ('eb', 'ms'), ('1-5ol', 'B-myr')]:
        ia, ib = o.index(a), o.index(b)
        cab = np.mean([r(T[i, :, ia], T[j, :, ib]) for i in range(3) for j in range(3) if i != j])
        caa = np.mean([r(T[i, :, ia], T[j, :, ia]) for i in range(3) for j in range(3) if i != j])
        cbb = np.mean([r(T[i, :, ib], T[j, :, ib]) for i in range(3) for j in range(3) if i != j])
        assert est[tuple(sorted((a, b)))] == pytest.approx(1 - cab / np.sqrt(caa * cbb), abs=1e-12)


@pytest.mark.parametrize('measure', ['pearson', 'cosine'])
def test_expected_distances_through_noise(measure):
    rng = np.random.default_rng(1)
    Xm = _structured(rng, n_units=400, n_odors=5)
    sd = rng.uniform(0.2, 1.0, Xm.shape[0])
    G, N2 = np.zeros((5, 5)), np.zeros(5)
    n_mc = 4000
    for _ in range(n_mc):
        Y = Xm + sd[:, None] * rng.standard_normal(Xm.shape)
        Yc = Y - Y.mean(axis=0) if measure == 'pearson' else Y
        g = Yc.T @ Yc
        G += g
        N2 += np.diag(g)
    G, N2 = G / n_mc, N2 / n_mc
    mc = 1 - G / np.sqrt(np.outer(N2, N2))
    ex = X.distance_matrix(Xm, measure, sd)
    off = ~np.eye(5, dtype=bool)
    np.testing.assert_allclose(ex[off], mc[off], atol=2e-3)


def test_expected_cosine_equals_zcell_pair_distances():
    rng = np.random.default_rng(2)
    counts = pd.DataFrame(np.round(_structured(rng, n_units=600) * 3), columns=list(X.ODORS16))
    z = H.session_zscore(counts, q=0.13, noise=0.5, amp='linear', noise_ratio=0.4, drift=1.0)
    sd = z['trial_sd'][0] / np.sqrt(z['n_trials'])
    ours = X._square_pairs(X.distance_matrix(z['center'][0], 'cosine', sd), X.ODORS16, X.PAIRS120)
    ref = X.E90._cs(H.zcell_pair_distances(z, 'graded')).reindex(X.PAIRS120)
    np.testing.assert_allclose(ours.to_numpy(), ref.to_numpy(), atol=1e-12)


@pytest.mark.parametrize('share', [0.0, 0.5])
def test_known_truth_check_synthetic(share):
    """The estimate on pseudo-recordings (three trials, noise at a reliability like the flies', a share of each unit's
    noise common to the odors of a repeat) returns the noise-free group distances within CHECK_TOL; the trial means do
    not."""
    rng = np.random.default_rng(3)
    truth = _structured(rng, n_units=2000)
    sig = 1.2 * truth.std()
    res = {m: [] for m in X.MEASURES}
    tms = {m: [] for m in X.MEASURES}
    for _ in range(10):
        xi = rng.standard_normal((3, *truth.shape))
        eta = rng.standard_normal((3, truth.shape[0], 1))
        T = truth[None] + sig * (np.sqrt(1 - share) * xi + np.sqrt(share) * eta)
        for m in X.MEASURES:
            res[m].append(X.noise_removed(T, m))
            tms[m].append(X._square_pairs(X.distance_matrix(T.mean(axis=0), m), X.ODORS16, X.PAIRS120))
    for m in X.MEASURES:
        t = X.group_means(X._square_pairs(X.distance_matrix(truth, m), X.ODORS16, X.PAIRS120))
        e = X.group_means(pd.concat(res[m], axis=1).mean(axis=1))
        tm = X.group_means(pd.concat(tms[m], axis=1).mean(axis=1))
        for g in X.OUTCOME_GROUPS + ('near', 'all66'):        # the groups of several pairs (single pairs are leads)
            assert abs(e[g] - t[g]) <= X.CHECK_TOL, (m, g, e[g], t[g])
        if share == 0.0:
            assert tm['all66'] - t['all66'] > 0.02            # the trial means are inflated by the noise


def test_check_table_and_passes():
    rows = []
    for w in ('a', 'b'):
        for c in X.CONDS:
            for lv in X.CHECK_LEVELS:
                for g in ('ester', 'alcohol'):
                    err = 0.05 if (g == 'alcohol' and lv == 'L3x' and c == 'on') else 0.01
                    rows.append(dict(wiring=w, cond=c, level=lv, units='z', measure='cosine', group=g, truth=0.5,
                                     estimate=0.5 + err, trial_mean=0.6, multiplier=1.0, share=0.0))
    df = X.check_table(pd.DataFrame(rows))
    p = X.check_passes(df)
    assert p[('ester', 'cosine', 'z')] and not p[('alcohol', 'cosine', 'z')]
    assert set(df.wiring) == {'a', 'b', 'population'}


@needs_shared
@pytest.mark.skipif(QUICK, reason='slow (E89 pseudo-recordings of one wiring)')
def test_known_truth_check_real_wiring():
    from experiments import E88_responder_sets as E88
    from experiments import E89_apl_breadth_weighting as E89
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    df = pd.DataFrame(X.check_wiring(('real', rel, gap)))
    d = df[df.group.isin(X.TABLE_GROUPS)]
    assert d.error.abs().max() <= X.CHECK_TOL


def _flies(mean):
    return pd.DataFrame({g: dict(mean=mean, lo=mean - 0.02, hi=mean + 0.02) for g in X.GROUPS}).T


def test_compare_and_verdicts_on_made_up_tables():
    fl = _flies(0.5)
    pop = {g: 0.5 for g in X.GROUPS}
    pop['ester'] = 0.42                  # too alike by 0.08, all wirings on that side
    pop['alcohol'] = 0.53                # outside the interval but below 0.05
    pop['other'] = 0.58                  # 0.08 more distinct but only 5 of 8 wirings on that side
    w = pd.DataFrame([{g: pop[g] for g in X.GROUPS} for _ in range(8)])
    w.loc[:2, 'other'] = 0.45
    r = {x['group']: x for x in X.compare(pop, w, fl)}
    assert r['ester']['differs'] and r['ester']['diff'] < 0
    assert not r['alcohol']['differs'] and r['alcohol']['outside'] and not r['alcohol']['big']
    assert not r['other']['differs'] and r['other']['n_same_side'] == 5
    # a population of four needs 3 on the side
    r4 = {x['group']: x for x in X.compare(pop, w.iloc[:4].assign(other=[0.45, 0.6, 0.6, 0.6]), fl, n_same=3)}
    assert r4['other']['differs'] and r4['ester']['differs']
    # conditions
    E = lambda d, n=4, s=(): dict(differs=d == n and n > 0, mixed=0 < d < n, n_differ=d, n_readings=n, signs=list(s))
    assert X.condition_verdict(E(4, s=[-1]), E(0)) == 'with APL on'
    assert X.condition_verdict(E(0), E(4, s=[1])) == 'without APL'
    assert X.condition_verdict(E(4), E(4)) == 'both'
    assert X.condition_verdict(E(0), E(0)) == 'neither'
    assert X.condition_verdict(E(4), E(2)) == 'not resolved'
    assert X.condition_verdict(E(0, 0), E(0)) == 'not read'
    # outcomes
    soma = {g: dict(on=E(0), off=E(0)) for g in X.GROUPS}
    soma['ester']['on'] = E(4, s=[-1])
    soma['alcohol']['on'] = E(4, s=[1])
    step = {g: dict(on=E(0), off=E(0)) for g in X.GROUPS}
    o = X.outcomes(soma, step)
    assert o['a'] and not o['b'] and not o['c']
    soma['other']['off'] = E(4, s=[-1])
    o = X.outcomes(soma, step)
    assert not o['a'] and o['b']
    soma['near_6al_t2h']['off'] = E(2, s=[1])       # a single pair: a lead, not read for the outcomes
    assert not X.outcomes(soma, step)['c']
    soma['ester_alcohol']['on'] = E(2, s=[-1])
    assert X.outcomes(soma, step)['c']


def test_verdict_rows_from_made_up_tables():
    """`soma_verdict_rows` reads the population and the eight wirings per reading and drops nothing that is read; a
    noise-removed comparison of a group that failed the check is kept with read False."""
    rows = []
    readings = [f'{p}_{ms}' for p, _, _, ms in X.like_for_like()] + [f'{r}_{ms}' for r in ('raw', 'cal') for ms in X.MEASURES]
    for rd in readings:
        for c in X.CONDS:
            v = 0.40 if c == 'on' else 0.5
            rows.append(dict(kind='population', part='', run='all8', family='C2', wiring='all8', cond=c, reading=rd,
                             **{g: v for g in X.GROUPS}))
            for w in X.ALL_WIRINGS:
                rows.append(dict(kind='run', part='B', run=w, family='C2', wiring=w, cond=c, reading=rd,
                                 **{g: v for g in X.GROUPS}))
    sm = pd.DataFrame(rows)
    passes = {(g, m, u): True for g in X.GROUPS for m in X.MEASURES for u in X.UNITS}
    passes[('alcohol', 'cosine', 'z')] = False
    vr = X.soma_verdict_rows(sm, passes, flies=lambda u, g, r, m: _flies(0.5))
    lf = vr[vr.set == 'like_for_like']
    assert len(lf) == 16 * 2 * len(X.GROUPS)
    assert lf[(lf.cond == 'on')].differs.all() and not lf[(lf.cond == 'off')].differs.any()
    x = lf[(lf.group == 'alcohol') & (lf.measure == 'cosine') & (lf.flies_units == 'z')]
    assert (~x[x.flies_reading == 'noise_removed'].read).all() and x[x.flies_reading == 'trial_mean'].read.all()
    v = X.verdicts(lf)
    assert v['ester']['verdict'] == 'with APL on' and v['alcohol']['on']['n_readings'] == 14


def test_line_and_place_on_line():
    a = np.linspace(-0.3, -0.1, 13)
    e = 0.5 * a - 0.05 + np.array([0.01, -0.01] * 5 + [0.02, -0.02, 0.0])
    ln = X.fit_line(e, a)
    assert ln['n'] == 13 and ln['slope'] == pytest.approx(0.5, abs=0.05) and ln['scatter'] > 0
    fp = {f'{rd}_{h}': dict(ester=0.1, alcohol=-0.25) for rd in X.FLY_READINGS for h in X.FLY_POINTS}
    ln = X._with_flies(ln, fp)
    assert all(v['resid'] > 0 for v in ln['flies'].values())
    on_line = X.place_on_line(ln['intercept'] + ln['slope'] * -0.2, -0.2, ln)
    assert abs(on_line['resid']) < 1e-12 and not on_line['off_toward_flies']
    up = X.place_on_line(ln['intercept'] + ln['slope'] * -0.2 + 3 * ln['scatter'], -0.2, ln)
    assert up['off_toward_flies'] and up['off_toward_flies_sd'] and up['toward_trial_mean_A'] > 0
    # beyond the SD but within the extent of the line's own residuals: off only by the SD reading
    mid = X.place_on_line(ln['intercept'] + ln['slope'] * -0.2 + 0.5 * (ln['scatter'] + ln['max_abs_resid']), -0.2, ln)
    if ln['max_abs_resid'] > ln['scatter']:
        assert mid['off_toward_flies_sd'] and not mid['off_toward_flies']
    # none of the line's own points is off it
    for (aa, ee) in zip(a, e):
        assert not X.place_on_line(ee, aa, ln)['off_toward_flies']
    down = X.place_on_line(ln['intercept'] + ln['slope'] * -0.2 - 3 * ln['scatter'], -0.2, ln)
    assert not down['off_toward_flies'] and down['toward_trial_mean_B'] < 0
    # relative to a reference residual
    assert not X.place_on_line(ln['intercept'] + ln['slope'] * -0.2 + 3 * ln['scatter'], -0.2, ln,
                               ref_resid=3 * ln['scatter'])['off_toward_flies']
    s = pd.Series(0.5, index=X.PAIRS120)
    m = X.group_masks(X.PAIRS120)
    s[m['ester']] = 0.3
    s[m['alcohol']] = 0.2
    est, alc = X.relative_distances(s)
    mu = s.mean()
    assert est == pytest.approx(0.3 - mu) and alc == pytest.approx(0.2 - mu)


@needs_shared
def test_flies_trial_means_equal_e90():
    S = X.flies_samples('z')
    f = pd.read_csv(H.results_dir('E90_model_stages') / 'stages_flies.csv')
    for key, (gk, ms) in {'soma_ctrl': ('ctrl', 'pearson'), 'soma_ctrl_cos': ('ctrl', 'cosine'),
                          'soma_tnt': ('tnt', 'pearson'), 'soma_tnt_cos': ('tnt', 'cosine')}.items():
        e = f[(f.key == key) & (f.kind == 'pair')].set_index(['odor_a', 'odor_b'])['mean']
        m = S[(gk, 'trial_mean', ms)].mean(axis=0).reindex(e.index)
        np.testing.assert_allclose(m.to_numpy(), e.to_numpy(), atol=1e-12)
        eg = f[(f.key == key) & (f.kind == 'group')].set_index('group')['mean']
        mg = X.group_means(S[(gk, 'trial_mean', ms)].mean(axis=0))
        for g in X.E90.GROUPS:
            assert mg[g] == pytest.approx(float(eg[g]), abs=1e-12)


@needs_shared
@pytest.mark.parametrize('units', ['z', 'noise'])
def test_cosine_noise_removed_equals_e81(units):
    nf = pd.read_csv(H.results_dir('E81_repeat_structure') / 'noise_free_pairs.csv')
    nf = nf[nf.units == units]
    nf.index = pd.MultiIndex.from_tuples([tuple(sorted(p)) for p in zip(nf.odor_a, nf.odor_b)])
    S = X.flies_samples(units)
    for gk, col in (('ctrl', 'ctrl_noise_free'), ('tnt', 'tnt_noise_free'), ('ctrl', 'ctrl_trial_mean'),
                    ('tnt', 'tnt_trial_mean')):
        rd = 'trial_mean' if 'trial_mean' in col else 'noise_removed'
        m = S[(gk, rd, 'cosine')].mean(axis=0)
        np.testing.assert_allclose(m.to_numpy(), nf[col].reindex(m.index).to_numpy(), atol=1e-12)
        a, b = X.group_means(m), X.group_means(nf[col])
        for g in X.GROUPS:
            assert a[g] == pytest.approx(b[g], abs=1e-12)


@needs_shared
def test_model_calibrated_readout_rows_equal_e90():
    o = X.e90_run('B', 'real')
    r = X.somata_readings(o['on'][1], o['off'][1])
    sm = pd.read_csv(H.results_dir('E90_model_stages') / 'stages_model.csv')
    e = sm[(sm.stage == 'soma') & (sm.wiring == 'real') & (sm.subset == 'all')]
    assert len(e) == 8
    for _, row in e.iterrows():
        mine = r[(r.cond == row.cond) & (r.reading == row.reading)].iloc[0]
        for g in X.E90.GROUPS:
            assert mine[g] == pytest.approx(float(row[g]), abs=1e-12)
    # the linear map in baseline-noise units without noise is raw counts (one noise level for all KCs)
    for c in X.CONDS:
        for ms in X.MEASURES:
            a = r[(r.cond == c) & (r.reading == f'linear_noise_none_{ms}')].iloc[0]
            b = r[(r.cond == c) & (r.reading == f'raw_{ms}')].iloc[0]
            for g in X.GROUPS:
                assert a[g] == pytest.approx(float(b[g]), abs=1e-9)


def test_effect_parts_keep_every_pair_in_e83_orientation():
    """Series in `_cs`'s sorted orientation (E83's targets, E90's pairs) lose no pair in the harness's alignment, and
    the pattern part equals E83's."""
    rng = np.random.default_rng(4)
    on = X.E90._cs(pd.Series(rng.uniform(0.3, 0.9, 120), index=X.PAIRS120))
    off = X.E90._cs(on - rng.uniform(0.0, 0.3, 120))
    data = X.E90._cs(pd.Series(rng.normal(0.18, 0.08, 120), index=X.PAIRS120))
    flipped = sum(1 for a, b in on.index if list(X.ODORS16).index(a) > list(X.ODORS16).index(b))
    assert flipped > 0
    ef = X.effect_parts(on, off, data)
    assert int(ef['apl_effect_n_pairs']) == 120
    from experiments import E83_cross_wiring_fit as E83
    sc = E83.effect_scores({'on': on, 'off': off}, {'t': data})
    assert ef['apl_effect_rmse_centered'] == pytest.approx(sc['t_pattern'], abs=1e-12)
    assert ef['apl_effect_amplitude'] ** 2 + ef['apl_effect_shape'] ** 2 == pytest.approx(sc['t_pattern'] ** 2, abs=1e-12)
    with pytest.raises((RuntimeError, KeyError)):
        X.effect_parts(on, off, data.iloc[:100])
    # a silent odor's pairs (NaN distances) are left out, and only those
    on2 = on.copy()
    on2[[p for p in on2.index if 'ms' in p]] = np.nan
    assert int(X.effect_parts(on2, off, data)['apl_effect_n_pairs']) == 105
