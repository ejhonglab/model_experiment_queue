"""Tests of E128's write-up checks' tables (analyses/E128_writeup_*.py); synthetic data only, no shared-directory reads."""
import json

import numpy as np
import pandas as pd
import pytest

from analyses import E128_writeup_equal_weight as EW
from analyses import E128_writeup_leads as EL
from analyses import E128_writeup_thresholds as ET
from experiments import E128_bouton_groups as E


def _flies(seed=0):
    rng = np.random.default_rng(seed)
    t = np.array([[20.0] + [1.0] * 11, [4.0, 5, 6, 5, 4, 5, 6, 5, 4, 5, 6, 5], [6.0, 5, 4, 5, 6, 5, 4, 5, 6, 5, 4, 5],
                  [1.0] * 5 + [15.0] + [1.0] * 6])
    M = np.vstack([t[0] + rng.normal(0, 0.05, (30, 12)), t[1] + rng.normal(0, 0.05, (10, 12)),
                   t[2] + rng.normal(0, 0.05, (10, 12)), t[3] + rng.normal(0, 0.05, (2, 12)),
                   rng.normal(0, 0.3, (5, 12))])
    return {'fly_a': np.stack([M, M, M])}


def test_group_sizes_table():
    G = EW.group_sizes_table(_flies())
    r = G.iloc[0]
    assert (r.n_grouped, r.n_groups, r.n_groups_small, r.remainder) == (52, 4, 1, True)
    assert r.pct_grouped_top3_by_size == pytest.approx(100 * 50 / 52)
    assert r.pct_grouped_small == pytest.approx(100 * 2 / 52)
    assert r.pct_weight_small == pytest.approx(100 * 1 / 5)        # one small group among four groups and the remainder
    p = np.array([30, 10, 10, 2]) / 52
    assert r.n_effective == pytest.approx(1 / (p ** 2).sum())
    # the three of largest variance share hold the dominant group, so at least its 30 ROIs
    assert 100 * 30 / 52 <= r.pct_grouped_top3_by_share <= 100


def test_equal_weight_and_below():
    T = EW.equal_weight_table(_flies())
    assert list(T.columns[:3]) == ['acq', 'measured', 'equal'] and T.acq.iloc[0] == 'fly_a'
    lab = E.fly_groups(_flies()['fly_a'].mean(axis=0))
    assert T.equal.iloc[0] == pytest.approx(E.readings(_flies()['fly_a'].mean(axis=0), lab)['equal'])
    B = EW.below_table(T)
    assert set(B.reading) == set(T.columns) - {'acq'}
    assert (B.n_flies_below == [int(T[c].iloc[0] < E.HALF_WAY) for c in B.reading]).all()


def test_ranges():
    T = pd.DataFrame(dict(acq=['a', 'b', 'c'], x=[3.0, 1.0, 2.0], ok=[True, False, True]))
    R = EW.ranges(T, 'acq')
    assert list(R.statistic) == ['x']
    assert R.iloc[0][['lo', 'hi', 'acq_lo', 'acq_hi']].tolist() == [1.0, 3.0, 'b', 'a']


def test_amplitude_tables():
    sp = dict(totals_cv=(0.3, 0.4), resp_cv=(0.4, 0.6), top1=(0.05, 0.07))
    D = pd.DataFrame(dict(variant=['all'] * 4 + ['big'] * 4, wiring=['w1', 'w1', 'w2', 'w2'] * 2,
                          totals_cv=[0.35, 0.2, 0.35, 0.35] * 2, resp_cv=[0.7, 0.5, 0.5, 0.5] * 2,
                          top1=[0.06, 0.06, 0.08, 0.06] * 2, n_gloms=[10, 12, 14, 16] * 2, n_rois=[100] * 8))
    D['tin'] = D.totals_cv.between(*sp['totals_cv'])
    D['resp_above'] = D.resp_cv > sp['resp_cv'][1]
    D['top_above'] = D.top1 > sp['top1'][1]
    D['passes'] = D.tin & ~D.resp_above & ~D.top_above
    A, S = EW.amplitude_tables(D, sp)
    a = A[A.variant == 'all'].set_index('wiring')
    assert a.loc['w1', 'frac_totals_in'] == 0.5 and a.loc['w2', 'frac_totals_in'] == 1.0
    assert a.loc['w1', 'p_resp'] == 0.0 and a.loc['w1', 'p_top'] == 0.5 and a.loc['w2', 'passes'] == 0.5
    s = S.set_index('variant').loc['all']
    assert s.n_deals_totals_in == 3 and s.frac_resp_above == pytest.approx(1 / 3) and s.frac_either_above == pytest.approx(2 / 3)
    assert s.n_wirings_totals_in_ge_5pct == 2 and s.n_gloms_median == 13


def test_thresholds_summary():
    T = pd.DataFrame({'grouped': [0.02, 0.03], 'all_roi_rate': [0.03, 0.04], 'e119_fixed': [0.0418, 0.0546],
                      'rcv_grouped': [0.56, 0.57], 'rcv_all_roi_rate': [0.47, 0.48], 'rcv_e119_fixed': [0.43, 0.44]},
                     index=pd.Index(['M2', 'P3'], name='wiring'))
    S = ET.summary_table(T).set_index('rule')
    assert S.loc['e119_fixed', 'n_wirings_ge_5pct'] == 1 and S.loc['grouped', 'n_wirings_ge_5pct'] == 0
    assert S.loc['e119_fixed', 'pct_pass_min'] == pytest.approx(4.18) and S.loc['e119_fixed', 'wiring_max'] == 'P3'
    assert S.loc['all_roi_rate', 'resp_cv_median_hi'] == 0.48


def test_model_sample():
    nb = pd.Series([16, 14, 14, 10, 4], index=list('abcde'))
    s = EL.sample_table(nb).iloc[0]
    assert (s.n_glomeruli, s.n_boutons, s.n_boutons_top3) == (5, 58, 44)
    assert s.pct_boutons_top3 == pytest.approx(100 * 44 / 58)


def test_save_writeup(tmp_path):
    p = E.save_writeup(dict(x=pd.DataFrame(dict(a=[1, 2]))), 'analyses/E128_writeup_test.py', tmp_path, note='n')
    assert p.name == 'provenance_writeup_test.json'
    assert pd.read_csv(tmp_path / 'writeup_x.csv').a.tolist() == [1, 2]
    d = json.loads(p.read_text())
    assert d['tables'] == ['writeup_x.csv'] and d['note'] == 'n' and 'queue_commit' in d
