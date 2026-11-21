#!/usr/bin/env python3
"""Tests for E122 (one feedback neuron that survives the block in place of KC-KC inhibition;
experiments/E122_feedback_neuron.py) on made-up matrices: the arms' matrices (total weight, structure, rank), the
within-type shuffle, step 0's predictors and scores, the class rule on T, and the refits' tags and commands. No model
simulation."""

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E122_feedback_neuron')


def _toy(seed=0):
    rng = np.random.default_rng(seed)
    types = np.array(['g'] * 7 + ["a'b'"] * 4 + ['ab'] * 9, dtype=object)
    n = len(types)
    ref = rng.poisson(0.6, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.4)
    np.fill_diagonal(ref, 0.0)
    ref *= n / ref.sum()
    return ref, types


@pytest.mark.parametrize('arm', ['U', 'TY', 'G', 'SH', 'GAB', 'HB'])
def test_total_weight_matched(arm):
    ref, types = _toy()
    M = X.arm_array(arm, ref, types)
    assert M.shape == ref.shape
    assert np.isclose(M.sum(), ref.sum())
    assert (M >= 0).all()


def test_structure_of_pooled_arms():
    ref, types = _toy()
    n = len(types)
    U = X.arm_array('U', ref, types)
    assert np.allclose(U, U[0, 0]) and np.linalg.matrix_rank(U) == 1
    TY = X.arm_array('TY', ref, types)
    for t in X.TYPES:
        m = types == t
        blk = TY[np.ix_(m, m)]
        assert np.allclose(blk, blk[0, 0]) and np.linalg.matrix_rank(blk) == 1
        # each type's total incoming weight is the hemibrain's
        assert np.isclose(TY[m].sum(), ref[m].sum())
        assert np.allclose(TY[np.ix_(m, ~m)], 0.0)
    for arm, out in (('G', 'g'), ('GAB', 'ab')):
        G = X.arm_array(arm, ref, types)
        m = types == out
        assert np.allclose(G[~m], 0.0) and np.allclose(G[m], G[m][0, 0])
        assert np.linalg.matrix_rank(G) == 1
        assert np.isclose(G[m].sum(axis=1)[0], n / m.sum())


def test_shuffle_keeps_blocks():
    ref, types = _toy()
    S = X.arm_array('SH', ref, types, seed=3)
    assert np.allclose(np.diag(S), 0.0)
    assert not np.array_equal(S, ref)
    for tp in X.TYPES:
        for tq in X.TYPES:
            a = ref[np.ix_(types == tp, types == tq)]
            b = S[np.ix_(types == tp, types == tq)]
            assert np.isclose(a.sum(), b.sum())
            assert sorted(a.ravel()) == sorted(b.ravel())
    assert np.array_equal(S, X.arm_array('SH', ref, types, seed=3))


def test_arm_matrix_frame_and_unknown():
    ref, types = _toy()
    ids = pd.Index(np.arange(100, 100 + len(types)))
    df = pd.DataFrame(ref, index=ids, columns=ids)
    M = X.arm_matrix('U', df, types)
    assert M.index.equals(ids) and M.columns.equals(ids)
    with pytest.raises(ValueError):
        X.arm_array('T', ref, types)


def test_block_table_and_predictors():
    ref, types = _toy()
    bt = X.block_table(ref, types)
    assert np.isclose(bt['weight_share'].sum(), 1.0)
    assert np.isclose(bt['pair_share'].sum(), 1.0)
    P = X.predictors(ref, types)
    assert set(P) == {'U', 'TY', 'G', 'GAB', 'blocks', 'own_amounts'}
    # own amounts keep every KC's row sum from each pre type; blocks keep each block's total
    for tq in X.TYPES:
        m = types == tq
        assert np.allclose(P['own_amounts'][:, m].sum(axis=1), ref[:, m].sum(axis=1))
        for tp in X.TYPES:
            mp = types == tp
            assert np.isclose(P['blocks'][np.ix_(mp, m)].sum(), ref[np.ix_(mp, m)].sum())
    # a KC with the same activity everywhere: the own amounts predict the input exactly
    S = np.ones((len(types), 3))
    sc = X.prediction_scores(ref @ S, P['own_amounts'] @ S)
    assert np.isclose(sc['r2'], 1.0)


def test_prediction_scores_perfect_and_scaled():
    rng = np.random.default_rng(1)
    a = rng.random((20, 5))
    s = X.prediction_scores(a, a)
    assert np.isclose(s['r2'], 1.0) and np.isclose(s['pearson_r'], 1.0) and np.isclose(s['mean_r_within_odor'], 1.0)
    s2 = X.prediction_scores(a, 2 * a)
    assert s2['r2'] < 1.0 and np.isclose(s2['pearson_r'], 1.0)


def _rd(kind_vals):
    """A readings table: kind -> {unit: pattern} on the readout, every target the same; mean effect and raw loss 0."""
    rows = []
    for kind, units in kind_vals.items():
        for u, v in units.items():
            for ro in ('cal', 'raw'):
                rows.append(dict(kind=kind, unit=u, readout=ro, mean_effect=0.0, all_loss=0.0,
                                 **{f'{t}_pattern': v for t in X.TARGETS}))
    return pd.DataFrame(rows)


def _runs(kind, in_bands=True):
    return pd.DataFrame([dict(kind=kind, wiring=w, in_bands=in_bands) for w in X.T_WIRINGS])


@pytest.mark.parametrize('delta,cls', [(0.002, 'replaces KC-KC'), (-0.004, 'below the family'), (0.005, 'neither')])
def test_classify(monkeypatch, delta, cls):
    monkeypatch.setattr(X, 'carried', lambda arm: {})
    fam = {u: 0.050 for u in ('T4',) + X.T_WIRINGS}
    arm = {u: 0.050 + delta for u in ('T4',) + X.T_WIRINGS}
    c = X.classify(_rd({'family': fam, 'T_U': arm}), _runs('T_U'), 'U')
    assert c['cls'] == cls


def test_classify_three_of_four_and_bands(monkeypatch):
    monkeypatch.setattr(X, 'carried', lambda arm: {})
    fam = {u: 0.050 for u in ('T4',) + X.T_WIRINGS}
    arm = {'T4': 0.052, 'P3': 0.052, 'P4': 0.052, 'M2': 0.052, 'M3': 0.060}
    assert X.classify(_rd({'family': fam, 'T_U': arm}), _runs('T_U'), 'U')['cls'] == 'replaces KC-KC'
    arm['M2'] = 0.060
    assert X.classify(_rd({'family': fam, 'T_U': arm}), _runs('T_U'), 'U')['cls'] == 'neither'
    arm = {u: 0.052 for u in ('T4',) + X.T_WIRINGS}
    assert X.classify(_rd({'family': fam, 'T_U': arm}), _runs('T_U', False), 'U')['cls'] == 'neither'
    assert X.classify(_rd({'family': fam}), _runs('T_G').iloc[:0], 'G') == dict(arm='G', complete=False)


def test_refit_order():
    g = dict(passing=['SH', 'G', 'U'], arms={'SH': dict(pattern_cal_S4=0.050), 'G': dict(pattern_cal_S4=0.054),
                                              'U': dict(pattern_cal_S4=0.052)})
    assert X.refit_order(g) == ['U', 'G', 'SH']


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


@pytest.mark.skipif(not _shared(), reason='needs the shared directory (the data-target registry)')
def test_fit_command_and_tag():
    spec = X.FitSpec('TY_P1', 'TY', 'P1')
    cmd = X.fit_command(spec, python='py')
    i = cmd.index('--kckc-matrix')
    assert cmd[i + 1].endswith('matrices/kckc_TY.parquet') and cmd[cmd.index('--kckc-matrix-tag') + 1] == 'E122TY'
    assert cmd[cmd.index('--init-params') + 1] == 'results/E122_feedback_neuron/starts.json:S1_TY_P1'
    assert '--kckc-fixed' not in cmd and cmd[cmd.index('--experiment') + 1] == 'E122_feedback_neuron/fits/TY_P1'
    assert '--wiring' in cmd
    tag = X.fit_tag(spec)
    assert tag.startswith('wiring-claw_permutation_within_type-seed') and tag.endswith('_aplboth_kcm-E122TY_k2_cap3_R')
