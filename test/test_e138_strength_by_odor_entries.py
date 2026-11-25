"""Tests of experiments/E138_strength_by_odor_entries.py (E138, the odors that carry the strength error and the glomeruli
that feed them): the odor entries and their shares against hal's scores on planted populations, the strength split
against E136's, the halves' agreement, the class rules ("concentrated", (d)), the feeders' rule, the counts of a planted
wiring, the IR class from Task et al. 2022's table, the classes lock's guard, and the held design item 6 (refused
without Tom's or Yang's word; its function on planted data only)."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments import E138_strength_by_odor_entries as E
from kcq import population_fit as P

ODORS = [f'o{i}' for i in range(8)]


def planted_pop(name, seed, n_units=300, scale=None, odors=ODORS):
    rng = np.random.default_rng(seed)
    sc = np.ones(len(odors)) if scale is None else np.asarray(scale, float)
    v = pd.DataFrame(rng.gamma(1.0, 1.0, size=(n_units, len(odors))) * sc[None, :],
                     index=[str(i) for i in range(n_units)], columns=odors)
    return P.Population.from_criterion(v, 1.5, name=name)


def planted_case(n_wirings=4, n_flies=4):
    flies = [planted_pop(f'f{i}', 100 + i, scale=np.linspace(0.6, 1.6, len(ODORS))) for i in range(n_flies)]
    sc = np.linspace(0.6, 1.6, len(ODORS))
    sc[2] = 4.0                                     # one odor far too strong in the model
    pops = [planted_pop(f'w{i}', i, scale=sc) for i in range(n_wirings)]
    ref = E.strength_reference(flies, ODORS)
    return pops, flies, ref


def test_imports_and_constants():
    assert E.EXPERIMENT == 'E138_strength_by_odor_entries'
    assert E.WIRINGS == E.S_WIRINGS + E.T_WIRINGS and len(E.WIRINGS) == 8
    assert E.STRENGTH_PARTS == ('strength/rate', 'strength/mean')
    assert E.CLASSES == ('ir', 'not_dendrite28', 'few_hemibrain')
    assert set(E.MAIN_READING) == {s.panel for s in E.PANEL_SPECS}


def test_entries_match_hal_scores_and_shares():
    pops, flies, ref = planted_case()
    names = [p.name for p in pops]
    ent = E.odor_entries(pops, ref, names)
    tab = P.wirings_as_flies(pops, ref)
    for row in names + ['common']:
        x = ent[ent.row == row]
        assert len(x) == len(ODORS)
        assert np.isclose(x.view_error.iloc[0], tab.loc[row, 'view_strength'])
        assert np.isclose(x.rate_error.iloc[0], tab.loc[row, 'part_strength_rate'])
        assert np.isclose(x.mean_error.iloc[0], tab.loc[row, 'part_strength_mean'])
        assert np.isclose(x.view_share.sum(), 1.0) and np.isclose(x.rate_share.sum(), 1.0)
        # the error without an odor's entry: its term removed, the number of odors kept
        k = x.view_share.idxmax()
        assert np.isclose(x.loc[k, 'view_error_without'] ** 2, x.loc[k, 'view_error'] ** 2 * (1 - x.loc[k, 'view_share']))
        assert np.allclose(x.rate_resid, x.rate_model - x.rate_flies)
        assert (x.rate_flies_min <= x.rate_flies + 1e-12).all() and (x.rate_flies <= x.rate_flies_max + 1e-12).all()
    c = ent[ent.row == 'common']
    assert E.top_odors(ent, 'common', 1) == ['o2']
    assert c.set_index('odor').view_share['o2'] > 0.5


def test_split_matches_e136():
    from experiments import E136_stage_population_fit as F
    pops, flies, ref = planted_case()
    mine = E.strength_split_of(pops, ref).set_index('row')
    theirs = F.strength_split(pops, ref)
    assert list(mine.index) == list(theirs.index)
    pd.testing.assert_frame_equal(mine[theirs.columns].astype(float), theirs.astype(float), check_names=False)
    # level^2 + spread^2 + order^2 is the part's error^2
    tab = P.wirings_as_flies(pops, ref)
    for row in ('w0', 'common'):
        r = mine.loc[row]
        assert np.isclose(r.strength_rate_level ** 2 + r.strength_rate_spread ** 2 + r.strength_rate_order ** 2,
                          tab.loc[row, 'part_strength_rate'] ** 2)


def test_partial_reference_entries():
    """A partial-panel reference (validation2's 'fly_scale' reading beside E134's): the entries use hal's
    `reference_views`, and their view error is hal's score."""
    rng = np.random.default_rng(5)
    flies = []
    for i in range(5):
        f = planted_pop(f'f{i}', 200 + i, scale=np.linspace(0.6, 1.6, len(ODORS)))
        drop = list(rng.choice(ODORS[3:], size=2, replace=False))
        f = P.Population(f.values.assign(**{o: np.nan for o in drop}), f.responds, name=f.name)
        flies.append(f)
    ref = P.partial_reference(flies, rule='fly_scale', odors=ODORS, breadth_edges=P.APL16_BREADTH_EDGES,
                              parts=list(E.STRENGTH_PARTS))
    pops = [planted_pop(f'w{i}', i) for i in range(3)]
    ent = E.odor_entries(pops, ref)
    sc = P.score_many(pops, ref)
    for p in pops:
        assert np.isclose(ent[ent.row == p.name].view_error.iloc[0], sc.loc[p.name, 'view_strength'])


def test_concentrated_rule():
    def ent(shares):
        return pd.DataFrame(dict(row='common', odor=list(shares), view_share=list(shares.values()),
                                 rate_share=list(shares.values()), mean_share=list(shares.values())))
    a = ent({'aa': 0.6, 'va': 0.2, 'eb': 0.1, 'ms': 0.1})
    yes = E.concentrated(a, {'fit': ent({'aa': 0.5, 'va': 0.1, 'eb': 0.3, 'ms': 0.1}),
                             'held': ent({'aa': 0.4, 'va': 0.3, 'eb': 0.2, 'ms': 0.1})})
    assert yes['odors'] == ['aa', 'va'] and yes['cls'] == 'concentrated'
    assert np.isclose(yes['share_halves']['fit'], 0.6) and yes['halves_own']['fit']['odors'] == ['aa', 'eb']
    # the same two odors, below half on one half: not concentrated
    no = E.concentrated(a, {'fit': ent({'aa': 0.5, 'va': 0.1, 'eb': 0.3, 'ms': 0.1}),
                            'held': ent({'aa': 0.2, 'va': 0.2, 'eb': 0.4, 'ms': 0.2})})
    assert no['cls'] == 'not concentrated'
    assert E.concentrated(a, {'fit': a})['cls'] == 'not readable'


def test_halves_agreement():
    pops, flies, ref = planted_case()
    a = E.odor_entries(pops, E.strength_reference(flies[:2], ODORS))
    b = E.odor_entries(pops, E.strength_reference(flies[2:], ODORS))
    ag = E.halves_agreement(a, b)
    assert set(ag) == {'rate_resid', 'mean_resid', 'rate_share', 'mean_share', 'view_share'}
    assert ag['rate_resid']['n'] == len(ODORS) and ag['rate_resid']['pearson'] > 0.8
    same = E.halves_agreement(a, a)
    assert np.isclose(same['view_share']['pearson'], 1.0) and np.isclose(same['view_share']['spearman'], 1.0)


def test_feeders_rule():
    x = pd.Series({'DC4': 90.0, 'DM1': 40.0, 'DP1m': 30.0, 'DM4': 20.0, 'VA2': -10.0, 'DM2': 20.0})
    f = E.feeders(x, share=2 / 3)
    # positive total 200; DC4 0.45, + DM1 0.65, + DP1m 0.80 >= 2/3
    assert list(f.glomerulus) == ['DC4', 'DM1', 'DP1m'] and np.isclose(f.cum_share.iloc[-1], 0.8)
    assert list(E.feeders(x, share=0.45).glomerulus) == ['DC4']
    assert E.feeders(pd.Series({'A': -1.0, 'B': 0.0})).empty


def planted_wiring(boutons):
    """A claw x bouton matrix from {glomerulus: [claws per bouton]}."""
    cols, blocks = [], []
    for g, claws in boutons.items():
        for j, n in enumerate(claws):
            cols.append((g, 0, j))
            blocks.append(n)
    n_claws = sum(blocks)
    W = np.zeros((n_claws, len(cols)))
    k = 0
    for j, n in enumerate(blocks):
        W[k:k + n, j] = 1
        k += n
    return pd.DataFrame(W, columns=pd.MultiIndex.from_tuples(cols, names=['glomerulus', 'pn_id', 'bouton_id']))


def test_wiring_counts_and_few():
    W = planted_wiring({'A': [3, 2], 'B': [5], 'C': [1, 1, 1], 'D': [4]})
    c = E.wiring_counts(W)
    assert c.loc['A', 'n_boutons'] == 2 and c.loc['A', 'n_claws'] == 5 and c.loc['C', 'n_boutons'] == 3
    assert np.isclose(c.share_claws.sum(), 1) and np.isclose(c.loc["B", "share_claws"], 5 / 17)
    few = E.few_hemibrain(c, ['A', 'B', 'C', 'D', 'E'])
    # boutons A2 B1 C3 D1 E0 (median 1): E; claws A5 B5 C3 D4 E0 (median 4): C, E
    assert few['few_boutons'] == ['E'] and few['few_claws'] == ['C', 'E'] and few['glomeruli'] == ['C', 'E']
    assert few['few_both'] == ['E']


def test_d_reading():
    hemi = E.wiring_counts(planted_wiring({'A': [2], 'B': [2], 'C': [2, 2]}))
    male_more = E.wiring_counts(planted_wiring({'A': [2, 2], 'B': [2], 'C': [2]}))
    male_less = E.wiring_counts(planted_wiring({'A': [1], 'B': [3], 'C': [3, 3]}))
    counts = {E.HEMIBRAIN: hemi, 'malecns-left': male_more, 'malecns-right': male_more}
    d = E.d_reading(['A'], counts)
    assert d['cls'] == 'more' and d['more']['malecns-left']['boutons']
    counts['malecns-right'] = male_less
    assert E.d_reading(['A'], counts)['cls'] == 'not more'
    assert E.d_reading([], counts)['cls'] == 'not readable'


def test_ir_class_from_task_table():
    ir = E.ir_glomeruli()
    model36 = {'DC4', 'DL2d', 'DL2v', 'DP1l', 'DP1m', 'VL2a', 'VL2p'}
    assert model36 <= set(ir['glomeruli'])
    assert not {'DM1', 'DA2', 'VA2', 'V', 'VM6v'} & set(ir['glomeruli'])
    assert ir['receptors']['DC4'] == ['Ir64a'] and len(ir['sha256']) == 64 and 'eLife' in ir['citation']


def test_classes_lock_guard(tmp_path, monkeypatch):
    missing = tmp_path / 'E138_classes_lock.json'
    with pytest.raises(PermissionError):
        E.require_classes_lock(copy=missing, shared=missing)
    missing.write_text(json.dumps({'x': 1}))
    with pytest.raises(PermissionError):                 # not tracked by git
        E.require_classes_lock(copy=missing, shared=missing)
    monkeypatch.setattr(E, 'CLASSES_LOCK_COPY', tmp_path / 'none.json')
    with pytest.raises(PermissionError):
        E.require_all()
    with pytest.raises(PermissionError):
        E.part_entries()


def test_regression_refused_without_word():
    with pytest.raises(SystemExit):
        E.main(['--part', 'regression'])
    with pytest.raises(SystemExit):
        E.part_regression('')


def test_class_regression_planted():
    """The held regression's function on planted numbers only (never on the real residuals)."""
    od = [f'x{i}' for i in range(12)]
    rng = np.random.default_rng(0)
    share = pd.Series(rng.uniform(0, 1, 12), index=od)
    res = 0.2 + 2.0 * share
    rd_share = pd.Series(rng.uniform(0, 1, 10), index=[f'y{i}' for i in range(10)])
    rd = 0.1 + 2.0 * rd_share + rng.normal(0, 0.01, 10)
    rand = [(rng.uniform(0, 1, 12), rng.uniform(0, 1, 10)) for _ in range(200)]
    out = E.class_regression(res, share, rd, rd_share, rand)
    assert np.isclose(out['slope'], 2.0) and np.isclose(out['intercept'], 0.2) and out['r'] > 0.99
    assert out['passes'] and out['quantile'] > 0.95
