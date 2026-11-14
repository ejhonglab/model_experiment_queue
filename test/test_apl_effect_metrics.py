"""analyses/apl_effect_metrics.py's split-half extras (split_half_subsets.csv, best_pick_stability.csv,
odor_share_agreement.csv) on made-up recordings and variants: they use the halves of split_half_reliability.csv
(its 'all distinct' consistency is the existing table's), and the pick counts are within their bounds."""
import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
A = pytest.importorskip('analyses.apl_effect_metrics')

ODORS = list(H.APL16_ODORS)
PAIRS = pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(ODORS) for b in ODORS[i + 1:]], names=['odor_a', 'odor_b'])


def _setup(seed: int = 0, n: int = 60):
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.3, 0.9, len(PAIRS))
    effect = rng.normal(0.18, 0.05, len(PAIRS))
    ctrl = pd.DataFrame(base + effect + rng.normal(0, 0.06, (11, len(PAIRS))), columns=PAIRS)
    tnt = pd.DataFrame(base + rng.normal(0, 0.06, (8, len(PAIRS))), columns=PAIRS)
    d = (ctrl.mean() - tnt.mean()).to_numpy()
    M = 0.1 + rng.uniform(0.3, 1.5, (n, 1)) * (d - d.mean()) + rng.normal(0, 0.04, (n, len(PAIRS)))
    M[5] = M[4]
    X, _ = H.pair_design(PAIRS, H.APL16_ODORS)
    se = A.pair_se(ctrl, tnt)
    df = A.metrics(M, d, 1 / se ** 2, X, kendall=False)
    df.insert(0, 'folder', [f'E{i % 3}_x' for i in range(n)])
    df['duplicate'] = False
    df.loc[5, 'duplicate'] = True
    df['realistic'] = np.arange(n) % 4 != 0
    return df, M, ctrl, tnt, X


def test_extras_use_the_halves_of_split_half_reliability():
    df, M, ctrl, tnt, X = _setup()
    uniq = ~df.duplicate.to_numpy()
    rel = A.split_half_reliability(M[uniq], ctrl, tnt, X, 12, np.random.default_rng(A.SEED)).set_index('metric')
    sub, picks, odors = A.split_half_extras(df, M, ctrl, tnt, X, 12)
    alld = sub[sub.subset == 'all distinct'].set_index('metric')
    assert set(alld.index) == set(rel.index)
    assert np.allclose(alld.split_half_rho, rel.split_half_rho.reindex(alld.index), atol=1e-12)
    assert set(sub.subset) == {'all distinct', 'inside the bands', 'within experiments'}
    n_in = int((uniq & df.realistic.to_numpy()).sum())
    assert (sub[sub.subset == 'inside the bands'].n == n_in).all()
    assert list(picks.metric) == A.PICK and (picks.n_experiments == 3).all()
    assert (picks.differs_between_halves <= 3).all() and (picks.differs_from_pattern_part_all_recordings <= 3).all()
    assert picks.set_index('metric').loc['pattern_part', 'differs_from_pattern_part_all_recordings'] == 0
    assert picks.set_index('metric').loc['pattern_part', 'differs_from_pattern_part_one_half'] == 0
    assert ((picks.best_rank_by_other_half_median >= 1) & (picks.best_rank_by_other_half_median <= n_in)).all()
    assert list(odors.n_variants) == [n_in, A.TOP_K]
    assert ((odors.median_r_between_halves >= -1) & (odors.median_r_between_halves <= 1)).all()
