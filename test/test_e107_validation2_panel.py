#!/usr/bin/env python3
"""Tests for E107 (the four validation2 odors in the panel their ORN input was recorded in;
experiments/E107_validation2_panel.py), on made-up inputs: the distance bands and their edges; the outcome classifier
on made-up ρ's (each of (a) to (d) and 'none', the middle band read on order alone, the cosine and correlation bands
disagreeing); the odor map through olfdb's abbreviations and aliases; the stimulus table's pairs and flies by
concentration; the quality-control identity; the per-fly reader (excluded stimuli, the drop list, pooled names by
concentration, one odor twice refused); the mean over flies with its least number of flies; the ceiling from halves of
the flies on the pairs both halves hold; the complete odor sets of the lock's rule for the 136 and the 153; and, with the
shared directory, the committed lock and the per-fly reading against the manuscript's file."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
M = pytest.importorskip('experiments.E107_validation2_panel')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return M.remy_dir().exists() and M.LOCK_PATH.exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# the bands and the classifier

def test_distance_band_edges():
    assert M.distance_band(0.70) == 'close'
    assert M.distance_band(0.7000001) == 'middle'
    assert M.distance_band(0.849) == 'middle'
    assert M.distance_band(0.85) == 'apart'
    assert M.distance_band(0.2) == 'close' and M.distance_band(1.1) == 'apart'


def _q(pb, all_, mm):
    return dict(pb=pb, all=all_, megamat=mm)


def test_classify_a_when_close():
    c = M.classify(0.6, 0.65, _q(0.9, 0.9, 0.9), _q(0.3, 0.3, 0.9))
    assert c['class'] == 'a' and c['band'] == 'close'
    assert c['order_class'] == 'b'                      # on order alone it would be (b)


def test_classify_b_c_d_none_when_apart():
    # (d): the family orders the pairs of paa and B-myr and the set within 0.05 of its megamat level
    assert M.classify(0.9, 0.95, _q(0.5, 0.7, 0.8), _q(0.76, 0.76, 0.8))['class'] == 'd'
    # (b): the ORN input orders them, the family does not
    c = M.classify(0.9, 0.95, _q(0.76, 0.8, 0.8), _q(0.6, 0.8, 0.8))
    assert c['class'] == 'b' and c['orn_orders'] and not c['fam_orders_pb']
    # (c): neither orders them
    assert M.classify(0.9, 0.95, _q(0.6, 0.8, 0.8), _q(0.6, 0.8, 0.8))['class'] == 'c'
    # 'none': the family orders the pairs of paa and B-myr but not the set
    assert M.classify(0.9, 0.95, _q(0.6, 0.8, 0.8), _q(0.8, 0.6, 0.8))['class'] == 'none'


def test_classify_order_tolerance_edge():
    # exactly the megamat level minus ORDER_TOL counts as ordering
    mm = 0.8
    c = M.classify(0.9, 0.9, _q(mm - M.ORDER_TOL, 0.8, mm), _q(0.5, 0.8, mm))
    assert c['orn_orders'] and c['class'] == 'b'
    assert c['margins']['orn_pb'] == pytest.approx(0.0)


def test_classify_middle_and_disagreeing_bands_read_on_order():
    c = M.classify(0.8, 0.82, _q(0.6, 0.8, 0.8), _q(0.6, 0.8, 0.8))
    assert c['band'] == 'middle' and c['on_order_alone'] and c['class'] == 'c'
    # cosine close, correlation apart: the band is 'middle', so not (a)
    c = M.classify(0.65, 0.9, _q(0.6, 0.8, 0.8), _q(0.8, 0.8, 0.8))
    assert c['band'] == 'middle' and c['class'] == 'd'
    assert (c['band_cos'], c['band_corr']) == ('close', 'apart')


# ---------------------------------------------------------------------------
# part 0's pieces on made-up tables

def _olfdb():
    return pd.DataFrame([dict(id=1, name='phenylacetaldehyde', abbrev='paa', aliases=['PAA'], cid=998),
                         dict(id=2, name='1-propanol', abbrev='1-prop', aliases=['1-3ol'], cid=1031),
                         dict(id=3, name='beta-myrcene', abbrev='B-myr', aliases=['b-myrcene'], cid=31253),
                         dict(id=4, name='eucalyptol', abbrev='euc', aliases=[], cid=2758)])


def test_odor_map_by_abbreviation_and_alias():
    mp = M.odor_map(['PAA', '1-prop', 'B-myr', 'euc'], _olfdb(), ['paa', '1-3ol', 'B-myr', 'euc'])
    assert {k: v['input'] for k, v in mp.items()} == {'PAA': 'paa', '1-prop': '1-3ol', 'B-myr': 'B-myr', 'euc': 'euc'}
    assert mp['PAA']['olfdb_id'] == 1 and mp['1-prop']['cid'] == 1031
    # case-insensitive only when the exact match fails
    assert M.odor_map(['B-MYR'], _olfdb(), ['B-myr'])['B-MYR']['input'] == 'B-myr'
    # no input odor for it: None
    assert M.odor_map(['euc'], _olfdb(), ['paa'])['euc']['input'] is None


def test_odor_map_refuses_unknown_and_ambiguous():
    with pytest.raises(ValueError):
        M.odor_map(['geos'], _olfdb(), ['paa'])
    o = pd.concat([_olfdb(), pd.DataFrame([dict(id=9, name='x', abbrev='y', aliases=['PAA'], cid=None)])], ignore_index=True)
    with pytest.raises(ValueError):
        M.odor_map(['PAA'], o, ['paa'])


def _stims():
    rows = []
    for fly, held in (('f1', {'PAA': -4.0, 'B-myr': -2.0, 'euc': -2.0}), ('f2', {'PAA': -4.0, 'B-myr': -2.0, 'euc': -1.75}),
                      ('f3', {'PAA': -4.0, 'euc': -2.0})):
        for name in ('PAA', 'B-myr', 'euc'):
            c = held.get(name, {'PAA': -4.0, 'B-myr': -2.0, 'euc': -2.0}[name])
            pres = name in held
            rows.append(dict(fly=fly, stim=f'{name} @ {c}', name=name, conc=c, presented=pres, in_rdm=pres,
                             qc_dropped=False, fly_in_manuscript=True))
    return pd.DataFrame(rows)


def _names():
    return {'PAA': dict(input='paa'), 'B-myr': dict(input='B-myr'), 'euc': dict(input='euc')}


def test_flies_by_pair_and_concentrations():
    fp = M.flies_by_pair(_stims(), _names())
    a = fp[fp.conc_a == 'all'].set_index(['odor_a', 'odor_b'])
    assert a.loc[('paa', 'B-myr'), 'n_flies'] == 2
    assert a.loc[('paa', 'euc'), 'n_flies'] == 3
    # euc's pairs by concentration: paa / euc at -2 (f1, f3) and at -1.75 (f2)
    by = fp[(fp.odor_a == 'paa') & (fp.odor_b == 'euc') & (fp.conc_a != 'all')]
    assert sorted(zip(by.conc_b, by.n_flies)) == [('-1.75', 1), ('-2', 2)]
    oc = M.odor_concentrations(_stims(), _names())
    e = oc[oc.odor == 'euc']
    assert list(zip(e.conc, e.n_in_rdm)) == [(-2.0, 2), (-1.75, 1)]


def test_qc_identity():
    s = _stims()
    s.loc[0, 'qc_dropped'] = True
    s.loc[0, 'in_rdm'] = False
    drops = pd.DataFrame(dict(fly=['f1', 'f3'], stim=['PAA @ -4.0', 'B-myr @ -2.0']))   # f3 was not shown B-myr
    q = M.qc_identity(s, drops)
    assert q['equal'] and q['in_rdm_is_presented_not_dropped'] and q['n_dropped'] == 1
    q = M.qc_identity(s, drops.iloc[1:])
    assert not q['equal'] and q['table_only'] == [('f1', 'PAA @ -4.0')]


# ---------------------------------------------------------------------------
# the readers' algebra

def test_pair_series_orientation():
    s = M.pair_series({('B-myr', 'paa'): 0.9, ('sab', 'euc'): 0.5})
    assert list(s.index) == [('paa', 'B-myr'), ('euc', 'sab')]          # APL16 order first, then by name
    assert M.pkey('zz', 'aa') == ('aa', 'zz')


def test_fly_rows_excluded_dropped_and_pooled():
    stims = ['PAA @ -4.0', 'B-myr @ -2.0', 'euc @ -2.0', 'euc @ -1.75', 'pfo @ 0.0']
    names = dict(_names(), pfo=dict(input='pfo'))
    X = np.full((2, 5, 5), np.nan)
    for k, held in enumerate(([0, 1, 2, 4], [0, 1, 3])):
        for i in held:
            for j in held:
                X[k, i, j] = 0.0 if i == j else 0.1 * (i + j) + k
    rows = M._fly_rows(X, ['f1', 'f2'], stims, names, drops={('f2', 'B-myr @ -2.0')})
    assert ('paa', 'B-myr') in rows.columns and ('paa', 'pfo') not in rows.columns     # the solvent is excluded
    assert rows.loc['f1', ('paa', 'B-myr')] == pytest.approx(0.1)
    assert np.isnan(rows.loc['f2', ('paa', 'B-myr')])                                  # dropped for f2
    assert rows.loc['f2', ('paa', 'euc')] == pytest.approx(0.3 + 1)                    # euc at -1.75, pooled
    bc = M._fly_rows(X, ['f1', 'f2'], stims, names, drops=set(), by_conc=True)
    assert ('paa', 'euc @ -2') in bc.columns and ('paa', 'euc @ -1.75') in bc.columns
    X2 = X.copy()
    X2[0, 3, 3] = 0.0                                                                  # f1 shown euc twice
    with pytest.raises(RuntimeError):
        M._fly_rows(X2, ['f1', 'f2'], stims, names, drops=set())


def test_fly_mean_least_flies():
    rows = pd.DataFrame({('a', 'b'): [1.0, 2.0, 3.0, np.nan], ('a', 'c'): [1.0, np.nan, np.nan, np.nan]})
    m = M.fly_mean(rows, min_flies=2)
    assert list(m.index) == [('a', 'b')] and m.iloc[0] == pytest.approx(2.0)


def test_split_half_ceiling():
    rng = np.random.default_rng(1)
    truth = rng.uniform(0, 1, 30)
    pairs = [('o%02d' % i, 'p') for i in range(30)]
    rows = pd.DataFrame([truth for _ in range(8)], columns=pairs)
    r = M.split_half_ceiling(rows, n_splits=20)
    assert r['split_half'] == pytest.approx(1.0) and r['ceiling'] == pytest.approx(1.0) and r['n_flies'] == 8
    noisy = pd.DataFrame([truth + rng.normal(0, 0.3, 30) for _ in range(10)], columns=pairs)
    a, b = M.split_half_ceiling(noisy, n_splits=50, seed=3), M.split_half_ceiling(noisy, n_splits=50, seed=3)
    assert a == b and 0 < a['split_half'] < 1 and a['spearman_brown_full'] > a['split_half']
    # a pair held by fewer than min_flies flies is left out; the halves read the pairs both hold
    sparse = noisy.copy()
    sparse.iloc[1:, 0] = np.nan
    assert M.split_half_ceiling(sparse, n_splits=10)['n_pairs'] == 29
    assert M.split_half_ceiling(noisy, pairs=pairs[:10], n_splits=10)['n_pairs'] == 10


def test_complete_odors_rule():
    odors = ['a', 'b', 'c', 'd']
    rows = []
    for x in odors:
        for y in odors:
            v = np.nan if {x, y} in ({'a', 'd'}, {'b', 'd'}) else 1.0
            rows.append(dict(a=x, b=y, col=v))
    t = pd.DataFrame(rows)
    assert M.complete_odors(t, 'col', odors) == ['a', 'b', 'c']           # d misses two pairs, removed first
    t2 = t.copy()
    t2.loc[(t2.a == 'a') & (t2.b == 'b'), 'col'] = np.nan                 # one order of a pair missing is enough
    assert M.complete_odors(t2, 'col', odors) == ['a', 'b', 'c']


def test_subsets():
    pairs = H.pair_list(['paa', 'B-myr', 'euc', 'sab', 'long'])
    assert len(M.subset_of(pairs, 'pb')) == 7
    assert len(M.subset_of(pairs, 'four')) == 9
    assert len(M.subset_of(pairs, 'nopb')) == 3
    assert M.subset_of(pairs, 'all') == pairs


# ---------------------------------------------------------------------------
# with the shared directory

@pytest.mark.skipif(not _shared(), reason="Remy's folder or the lock not on hand")
def test_lock_committed_and_matches_module():
    lock = M.require_lock()
    assert lock['definitions'] == M.lock_definitions()
    assert lock['validation2']['n_flies'] == 20 and len(lock['validation2']['odors']) == 22
    assert lock['validation2']['qc_identity']['equal']
    assert lock['odor_map']['validation2']['PAA']['input'] == 'paa'


@pytest.mark.skipif(QUICK or not _shared(), reason='QUICK, or no shared directory')
def test_per_fly_correlation_is_the_manuscripts():
    a, b = M.per_fly('correlation'), M.manuscript_per_fly()
    a = a.reindex(index=b.index, columns=b.columns)
    assert np.nanmax(np.abs(a.to_numpy() - b.to_numpy())) == 0.0
    assert (a.isna() == b.isna()).all().all()
