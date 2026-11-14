"""The pattern part of the APL-effect loss split into its amplitude and shape parts
(`kcq.harness.apl_effect_spread_shape`), on made-up distances: the parts add in squares to
`apl_effect_stats`' pattern part, and the data's spread is the KC data's."""
import numpy as np
import pandas as pd

from kcq import harness as H


def _pairs(values: np.ndarray) -> pd.Series:
    idx = pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(H.APL16_ODORS) for b in H.APL16_ODORS[i + 1:]],
                                    names=['odor_a', 'odor_b'])
    return pd.Series(values, index=idx)


def test_parts_add_to_pattern_part():
    rng = np.random.default_rng(0)
    on = _pairs(rng.uniform(0.5, 1.0, 120))
    off = _pairs(on.to_numpy() - rng.uniform(0.0, 0.3, 120))
    st = H.apl_effect_stats(on, off)
    sp = H.apl_effect_spread_shape(on, off)
    assert np.isclose(sp['apl_effect_amplitude'] ** 2 + sp['apl_effect_shape'] ** 2,
                      st['apl_effect_rmse_centered'] ** 2, rtol=1e-10)
    d = H.yang_apl_effect()
    assert np.isclose(sp['apl_effect_sd_data'], float(np.sqrt(np.mean((d - d.mean()) ** 2))), rtol=1e-12)
    assert np.isclose(sp['apl_effect_sd_ratio'], sp['apl_effect_sd_model'] / sp['apl_effect_sd_data'])


def _counts(seed: int, p_resp: float) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    u = rng.random((800, len(H.APL16_ODORS))) * rng.gamma(3.0, 1 / 3.0, (800, 1))
    n = np.where(u < p_resp, 2 + rng.geometric(0.12, u.shape), np.where(u < p_resp + 0.05, 1, 0))
    return pd.DataFrame(n.astype(float), columns=list(H.APL16_ODORS))


def test_readout_stats_carry_the_parts():
    """`apl_effect_stats_readout` keeps every key and value it had (`apl_effect_stats` of the read-out distances, the
    mean distances) and adds the amplitude and shape parts of the same read-out distances, which add in squares to its
    pattern part."""
    on, off = _counts(1, 0.07), _counts(2, 0.14)
    for kind, param in (('raw', None), ('saturating', 4.0)):
        out = H.apl_effect_stats_readout(on, off, kind, param)
        d_on, d_off = H.readout_pair_distances(on, kind, param), H.readout_pair_distances(off, kind, param)
        old = dict(H.apl_effect_stats(d_on, d_off), mean_dist_on=float(d_on.mean()), mean_dist_off=float(d_off.mean()))
        new = H.apl_effect_spread_shape(d_on, d_off)
        assert list(out) == list(old) + list(new)
        assert all(out[k] == v or (np.isnan(out[k]) and np.isnan(v)) for k, v in old.items())
        assert all(out[k] == v for k, v in new.items())
        assert np.isclose(out['apl_effect_amplitude'] ** 2 + out['apl_effect_shape'] ** 2,
                          out['apl_effect_rmse_centered'] ** 2, rtol=1e-10)
        assert np.isclose(out['apl_effect_shape'] ** 2, 2 * out['apl_effect_sd_model'] * out['apl_effect_sd_data']
                          * (1 - out['apl_effect_pearson']), rtol=1e-8)
