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
