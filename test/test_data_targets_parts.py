"""The data-target registry's floors of the pattern part's amplitude and shape parts and its Pearson ceiling
(kcq/data_targets.py, `apl_effect_resamples(..., parts=)`; how the main target is read, 2026-11-13 16:22): filling
`parts` leaves the draws and the returned values unchanged, each resample's parts add in squares to its pattern part
and match the harness's functions, and on the KC data the floors and the ceiling are computed from the same resamples
and splits as apl_effect_rmse_floor and apl_effect_rho_ceiling (computed here, not read from the saved registry)."""
import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
DT = pytest.importorskip('kcq.data_targets')

ODORS = list(H.APL16_ODORS)
PAIRS = pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(ODORS) for b in ODORS[i + 1:]], names=['odor_a', 'odor_b'])


def _recs(seed: int = 0) -> dict:
    """Made-up per-recording pair distances (rows recordings, columns the 120 pairs), the APL-on genotype's larger."""
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.3, 0.9, len(PAIRS))
    effect = rng.normal(0.18, 0.05, len(PAIRS))
    return {'TNT3c': pd.DataFrame(np.clip(base + effect + rng.normal(0, 0.06, (11, len(PAIRS))), 0, 2), columns=PAIRS),
            'TNT': pd.DataFrame(np.clip(base + rng.normal(0, 0.06, (8, len(PAIRS))), 0, 2), columns=PAIRS)}


def _same(a: tuple, b: tuple) -> bool:
    return a[0].equals(b[0]) and a[1] == b[1] and a[2] == b[2] and np.array_equal(a[3], b[3], equal_nan=True)


def test_parts_leave_the_resamples_unchanged():
    recs = _recs()
    parts = {}
    assert _same(DT.apl_effect_resamples(recs, seed=3), DT.apl_effect_resamples(recs, seed=3, parts=parts))
    assert [len(parts[k]) for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson')] == [1000] * 5
    assert len(parts['pearson_split_half']) == 500
    p, a, s = (np.array(parts[k]) for k in ('pattern', 'amplitude', 'shape'))
    assert np.allclose(a ** 2 + s ** 2, p ** 2, rtol=1e-10, atol=1e-15)
    assert np.all(np.abs(parts['pearson']) <= 1) and np.all(np.abs(parts['pearson_split_half']) <= 1)


def test_parts_are_the_harness_measures_of_each_resample():
    """Replaying the generator: each resample's parts are what `apl_effect_stats` / `apl_effect_spread_shape` give a
    model whose effect is the resampled data's, and each split's Pearson r is that of the halves' effects."""
    recs = _recs(1)
    parts = {}
    effect, *_ = DT.apl_effect_resamples(recs, seed=5, parts=parts)
    rng = np.random.default_rng(5)
    for i in range(1000):
        m = H.resample_recording_means(recs, rng)
        if i in (0, 999):
            st, sp = H.apl_effect_stats(m['TNT3c'], m['TNT'], effect), H.apl_effect_spread_shape(m['TNT3c'], m['TNT'], effect)
            assert parts['pattern'][i] == st['apl_effect_rmse_centered'] and parts['pearson'][i] == st['apl_effect_pearson']
            assert all(parts[k][i] == sp[f'apl_effect_{k}'] for k in ('amplitude', 'shape', 'sd_ratio'))
    halves = [{}, {}]
    for g, df in recs.items():
        idx = rng.permutation(len(df))
        halves[0][g], halves[1][g] = df.iloc[idx[:len(df) // 2]].mean(), df.iloc[idx[len(df) // 2:]].mean()
    a = H.apl_effect(halves[0]['TNT3c'], halves[0]['TNT'])
    b = H.apl_effect(halves[1]['TNT3c'], halves[1]['TNT']).reindex(a.index)
    assert parts['pearson_split_half'][0] == H.pearson(a, b)


def test_floors_and_ceiling_on_the_kc_data():
    """On Yang's APL16 recordings: the parts' floors lie between 0 and the loss's floor, add in squares to about the
    pattern part's, and the Pearson ceiling is the Spearman-Brown-corrected split-half r, near the ρ ceiling (hal's
    read-only check, 2026-11-13, with its own 2000 resamples and 1000 splits: pattern 0.029, amplitude 0.009, shape
    0.026, Pearson ceiling 0.837)."""
    H.shared_dir()
    recs = {g: H.yang_recording_rdms(g, 'cosine') for g in ('TNT3c', 'TNT')}
    parts = {}
    effect, boot_mean, boot_rmse, sh_rho = DT.apl_effect_resamples(recs, parts=parts)
    floor = {k: float(np.mean(parts[k])) for k in ('pattern', 'amplitude', 'shape')}
    assert 0 < floor['amplitude'] < floor['shape'] < floor['pattern'] < float(np.mean(boot_rmse))
    assert abs(floor['pattern'] - 0.029) < 0.004 and abs(floor['amplitude'] - 0.009) < 0.003 and abs(floor['shape'] - 0.026) < 0.004
    r = float(np.nanmean(parts['pearson_split_half']))
    ceiling = np.sqrt(2 * r / (1 + r))
    assert abs(ceiling - 0.837) < 0.02
    r_rho = float(np.nanmean(sh_rho))
    assert abs(ceiling - np.sqrt(2 * r_rho / (1 + r_rho))) < 0.03
    sd = float(np.sqrt(np.mean((effect - effect.mean()) ** 2)))
    assert np.isclose(sd, H.apl_effect_spread_shape(recs['TNT3c'].mean(), recs['TNT'].mean(), effect)['apl_effect_sd_data'])
