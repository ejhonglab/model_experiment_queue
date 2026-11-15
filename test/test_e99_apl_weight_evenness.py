#!/usr/bin/env python3
"""Tests for E99 (the evenness of APL's weights among KCs and over a KC's claws; experiments/E99_apl_weight_evenness.py):
the variant names and the blocks; the claw-count groups and their pooling; the eight settings' defining properties on a
made-up table (every mean 1, the two products of the identities, K1's group means, K2's and K4's constancy, K3 at none
and all of the inequality, K4's zeros, W2 for a KC without a kept synapse, the sums per KC, T1 and T2 as E10's new modes)
and on the cached real tables (part 0's identities, the entry's probe values); the variance split; the two descriptions;
part C's shares, classes, rule by wiring, bound and the KCs with none; the choice, the equal rule and part D's rule on
made-up rows; A2's file kept once part B has a row; A3's reading by fifths and of the claws; E10_joint_fit.py's new
--aplkc-flat modes; part D's commands and tags."""

import json
import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
G = pytest.importorskip('experiments.E99_apl_weight_evenness')
J = pytest.importorskip('experiments.E10_joint_fit')

QUICK = os.environ.get('QUICK') == '1'
KC = H.KC_ID_LEVEL


def _shared() -> bool:
    try:
        return H.shared_dir().exists() and (H.shared_dir() / 'results' / 'E67_apl_effect_fits' / 'fits.csv').exists()
    except Exception:
        return False


def _cache() -> bool:
    try:
        return (H.outputs_dir('E96_offclaw_apl') / 'W' / 'real_kc.parquet').exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# names

def test_the_blocks_and_their_sizes():
    n = len(G.S_WIRINGS)
    assert len(G.BLOCKS['I']) * n == 24 and len(G.BLOCKS['II']) * n == 40          # the entry's 64 points
    assert set(G.BLOCKS['I']) == {'K1_offclaw_d1', 'K2_offclaw_d1', 'W1_d1', 'W2_d1', 'T1_d1', 'T2_d1'}
    assert set(G.BLOCKS['II']) == {'K1_offclaw_d2', 'K2_offclaw_d2', 'K1_offclawped_d1', 'K2_offclawped_d1', 'K3_offclaw_d1',
                                   'K4_offclaw_d1', 'T1_d0.5', 'T1_d2', 'T2_d0.5', 'T2_d2'}
    assert not set(G.BLOCKS['I']) & set(G.BLOCKS['II'])
    assert G.APL_BRACKET == (0.01, 10.0) and G.D_DOSE == 0.1 and G.GROUP_MIN == 20 and G.K3_FRAC == 0.5


def test_variant_fields():
    p = G.parse_variant('K1_offclaw_d1')
    assert (p['shape'], p['counts'], p['size'], p['dose'], p['part'], p['screen']) == ('K1', 'offclaw', 1.0, 0.1, 'B', 'K1')
    assert p['construction'] and not p['equal']
    p = G.parse_variant('K2_offclawped_d1')
    assert p['counts'] == 'offclaw_ped' and p['construction']
    p = G.parse_variant('W2_d1')
    assert p['construction'] and not p['equal'] and p['counts'] is None
    for v in ('W1_d1', 'T1_d0.5', 'T2_d2'):
        p = G.parse_variant(v)
        assert p['equal'] and not p['construction'] and p['part'] == 'B'
    assert G.parse_variant('T1_d0.5')['size'] == 0.5
    p = G.parse_variant('D1_uniform_d2')
    assert (p['part'], p['screen'], p['mode'], p['size']) == ('e95', 'e95', 'uniform', 2.0)
    p = G.parse_variant('B6_allsyn_d1')
    assert (p['part'], p['screen'], p['mode'], p['equal']) == ('e96', 'e96', 'allsyn', False)
    p = G.parse_variant('B4_offclawped_d1_permall1')
    assert p['permute'] == 'all' and p['part'] == 'e96'
    assert G.parse_variant('ref')['screen'] == 'ref'
    assert G.parse_variant('D0')['control']
    for bad in ('K5_offclaw_d1', 'K1_d1', 'W1_offclaw_d1', 'T3_d1', 'T1_d1.0', 'T1_d0', 'K1_offclaw_d01', 'B7_offclaw_d1',
                'U_uniform_d1_pool1', 'B3_pool1x'):
        with pytest.raises(ValueError):
            G.parse_variant(bad)
    for shape in G.SHAPES:
        counts = 'offclaw' if shape in G.K_SHAPES else None
        for d in G.SIZES:
            v = G.vname(shape, d, counts)
            assert G.parse_variant(v)['shape'] == shape and G.parse_variant(v)['size'] == d
    with pytest.raises(ValueError):
        G.vname('K1', 1)
    with pytest.raises(ValueError):
        G.vname('W1', 1, 'offclaw')


def test_the_corners_and_the_controls():
    assert G.CORNERS['E95'] == ('D1_uniform_d0.5', 'D1_uniform_d1', 'D1_uniform_d2', 'D2_connectome_d1')
    assert set(G.CORNERS['E96']) == {'B1_offclaw_d1', 'B1_offclaw_d2', 'B2_offclawped_d1', 'B6_allsyn_d0.5', 'B6_allsyn_d1',
                                     'B6_allsyn_d2', 'B4_offclawped_d1_permclaw1', 'B4_offclawped_d1_permall1'}
    for vs in G.CORNERS.values():
        for v in vs:
            assert G.parse_variant(v)['part'] in ('e95', 'e96')
    assert G.control_of99('T1_d2') == 'D1_uniform_d2' and G.control_of99('W1_d1') == 'D1_uniform_d1'
    for v in ('K1_offclaw_d1', 'W2_d1', 'B6_allsyn_d1', 'ref'):
        with pytest.raises(ValueError):
            G.control_of99(v)


# ---------------------------------------------------------------------------
# the groups

def test_claw_groups_pool_small_groups_with_the_next_lower_count():
    n = np.repeat([1, 2, 3, 4, 5, 6], [25, 30, 5, 3, 40, 10])
    g = G.claw_groups(n)
    lab = dict(zip(n, g))
    assert lab[1] == 1 and lab[5] == 5
    assert lab[6] == 5                                    # 10 KCs with 6 claws join those with 5
    assert lab[4] == lab[3] == lab[2] == 2                # 3 + 5 still short of 20: on to the 30 with 2
    n2 = np.repeat([1, 2], [5, 30])
    assert len(set(G.claw_groups(n2))) == 1               # a lowest group still short joins the next higher one
    n3 = np.repeat([3, 4], [20, 19])
    assert len(set(G.claw_groups(n3))) == 1 and len(set(G.claw_groups(np.repeat([3, 4], [20, 20])))) == 2


# ---------------------------------------------------------------------------
# the settings on a made-up table

def _toy(seed: int = 0):
    """A made-up wiring: 25 KCs with 2 claws, 25 with 3, 10 with 4 (pooled with 3); off-claw counts with some zeros, the
    peduncle's, kept counts per claw with some zeros and one KC without any."""
    rng = np.random.default_rng(seed)
    n = np.repeat([2, 3, 4], [25, 25, 10])
    ids = np.arange(1000, 1000 + n.size)
    no_claw = rng.poisson(4, n.size).astype(float)
    no_claw[[0, 7, 30]] = 0
    ped = rng.poisson(3, n.size).astype(float)
    kc_of = np.repeat(ids, n)
    claw_id = np.concatenate([np.arange(k) for k in n])
    kept = rng.poisson(1.5, kc_of.size).astype(float)
    kept[kc_of == ids[5]] = 0                              # a KC without a synapse on a kept claw
    claws = pd.DataFrame({KC: kc_of, 'claw_id': claw_id, 'kept': kept})
    kc = pd.DataFrame({KC: ids, 'n_claws': n, 'no_claw': no_claw, 'peduncle': ped,
                       'kept': pd.Series(kept).groupby(kc_of).sum().reindex(ids).to_numpy()})
    claws['wnorm'] = kept / kept.sum() * 3.0               # the connectome's normalized weights (proportional to kept)
    return kc, claws


def test_every_setting_has_mean_one_and_is_not_negative():
    kc, claws = _toy()
    for s in G.SHAPES:
        for c in G.COUNTS:
            x = G.shape99(s, kc, claws, c)
            assert x.shape == (len(claws),) and np.isclose(x.mean(), 1.0, rtol=0, atol=1e-12) and (x >= 0).all()


def test_the_products_of_the_identities():
    kc, claws = _toy()
    kc_of = claws[KC].to_numpy()
    k = kc.set_index(KC)
    for c in G.COUNTS:
        off = J.aplkc_flat_shape(c, kc_of, J.aplkc_flat_offclaw(c, k))         # E96's off-claw component
        assert np.allclose(G.shape99('K1', kc, claws, c) * G.shape99('K2', kc, claws, c), off, rtol=1e-12, atol=0)
        assert np.allclose(G.shape99('K3', kc, claws, c, frac=1.0), off, rtol=1e-12, atol=0)
        assert np.allclose(G.shape99('K3', kc, claws, c, frac=0.0), G.shape99('K2', kc, claws, c), rtol=1e-12, atol=0)
    kept = claws['kept'].to_numpy()
    assert np.allclose(G.shape99('W1', kc, claws) * G.shape99('W2', kc, claws), kept / kept.mean(), rtol=1e-12, atol=0)


def test_k1_group_means_k2_k4_constancy_and_k4_zeros():
    kc, claws = _toy()
    kc_of = claws[KC].to_numpy()
    grp = G.group_of_kc(kc).reindex(kc_of).to_numpy()
    assert set(grp) == {2, 3}                                                  # the 10 KCs with 4 claws pooled with 3
    k1 = G.shape99('K1', kc, claws)
    assert np.allclose(pd.Series(k1).groupby(grp).mean(), 1.0, rtol=0, atol=1e-12)
    k2 = G.shape99('K2', kc, claws)
    assert (pd.Series(k2).groupby(grp).agg(lambda s: s.max() - s.min()) < 1e-12).all()
    # K2 is mu_G (the group's off-claw synapses per claw) over the mean: the group with 3 and 4 claws pooled
    a, n = kc['no_claw'].to_numpy(), kc['n_claws'].to_numpy()
    mu = {g: a[np.isin(n, ns)].sum() / n[np.isin(n, ns)].sum() for g, ns in ((2, [2]), (3, [3, 4]))}
    mean = a.sum() / n.sum()
    assert np.isclose(k2[grp == 2][0], mu[2] / mean) and np.isclose(k2[grp == 3][0], mu[3] / mean)
    k4 = G.shape99('K4', kc, claws)
    zero_kcs = set(kc[KC][kc.no_claw == 0])
    assert set(kc_of[k4 == 0]) == zero_kcs and np.ptp(k4[k4 > 0]) < 1e-12
    assert np.isclose(k4[k4 > 0][0], len(claws) / (k4 > 0).sum())
    for x in (k2, k4):
        assert (pd.Series(x).groupby(kc_of).agg(lambda s: s.max() - s.min()) < 1e-12).all()
    # K3: halfway between K2 and the off-claw component, claw by claw (before the common normalization, which is 1)
    k3 = G.shape99('K3', kc, claws)
    off = G.shape99('K3', kc, claws, frac=1.0)
    assert np.allclose(k3, (k2 + off) / 2, rtol=1e-12, atol=1e-15)


def test_w2_and_the_sums_per_kc():
    kc, claws = _toy()
    kc_of = claws[KC].to_numpy()
    ids = kc[KC].to_numpy()
    w2 = G.shape99('W2', kc, claws)
    none = kc[KC].iloc[5]
    assert np.allclose(w2[kc_of == none], 1.0)                                 # no synapse on a kept claw: 1 on each claw
    assert np.allclose(G.kc_totals(w2, kc_of, ids), kc['n_claws'].to_numpy())  # W2's sum per KC is n_k
    c_k = kc['kept'].to_numpy()
    for s, t in (('W1', c_k), ('T1', c_k + kc.no_claw.to_numpy() + kc.peduncle.to_numpy()), ('T2', c_k + kc.no_claw.to_numpy())):
        tot = G.kc_totals(G.shape99(s, kc, claws), kc_of, ids)
        assert np.allclose(tot * t.sum() / len(claws), t, rtol=1e-12)
        assert (pd.Series(G.shape99(s, kc, claws)).groupby(kc_of).agg(lambda q: q.max() - q.min()) < 1e-12).all()   # even


def test_shape_identities_pass_on_the_made_up_table():
    kc, claws = _toy()
    claws = claws.copy()
    k = kc.set_index(KC)
    kc_of = claws[KC].to_numpy()
    for c in G.COUNTS:
        claws[f'comp_{c}'] = J.aplkc_flat_shape(c, kc_of, J.aplkc_flat_offclaw(c, k))
    r = G.shape_identities(kc, claws)
    assert r['passes'], r
    assert r['w2_kcs_without_kept'] == int((kc.kept == 0).sum()) >= 1
    assert np.isclose(r['k4_zero_kcs_share']['offclaw'], (kc.no_claw == 0).mean())
    bad = claws.copy()
    bad.loc[int(np.argmax(bad.comp_offclaw > 0)), "comp_offclaw"] *= 1.01     # a cache that is not the counts fails
    assert not G.shape_identities(kc, bad)['passes']


def test_t1_t2_are_e10s_new_modes():
    kc, claws = _toy()
    idx = pd.MultiIndex.from_frame(claws[[KC, 'claw_id']])
    per_claw = pd.Series(claws['kept'].to_numpy(), index=idx)
    wn = pd.Series(claws['wnorm'].to_numpy(), index=idx)
    k = kc.set_index(KC)
    for s, mode in G.T_MODE.items():
        comp = J.aplkc_flat_component(wn, mode, (k, per_claw))
        assert np.array_equal(comp, float(wn.to_numpy().mean()) * G.shape99(s, kc, claws))
        assert np.isclose(comp.mean(), wn.mean())
    assert np.allclose(J.aplkc_flat_per_kc('allsyn_even', k), k.kept + k.no_claw + k.peduncle)
    assert np.allclose(J.aplkc_flat_per_kc('allsyn_even_calyx', k), k.kept + k.no_claw)
    assert np.allclose(J.aplkc_flat_per_kc('offclaw', k), k.no_claw)          # the old modes unchanged
    assert np.allclose(J.aplkc_flat_per_kc('allsyn', k), k.no_claw + k.peduncle)


def test_the_tables_are_checked():
    kc, claws = _toy()
    with pytest.raises(ValueError):
        G.shape99('K1', kc.assign(n_claws=kc.n_claws + 1), claws)
    with pytest.raises(ValueError):
        G.shape99('K1', kc.iloc[1:], claws)
    with pytest.raises(ValueError):
        G.shape99('K9', kc, claws)


# ---------------------------------------------------------------------------
# the variance split

def test_the_variance_split_adds_up():
    kc, claws = _toy()
    kc_of = claws[KC].to_numpy()
    grp = G.group_of_kc(kc).reindex(kc_of).to_numpy()
    for s in G.SHAPES:
        x = G.shape99(s, kc, claws)
        v = G.variance_split(x, kc_of, grp)
        assert np.isclose(v['parts_sum'], v['total'], rtol=1e-10) and np.isclose(v['total'], x.var())
    v = G.variance_split(G.shape99('K2', kc, claws), kc_of, grp)
    assert v['among'] < 1e-20 and v['within'] < 1e-20 and v['between'] > 0
    v = G.variance_split(G.shape99('K1', kc, claws), kc_of, grp)
    assert v['between'] < 1e-20 and v['within'] < 1e-20
    v = G.variance_split(G.shape99('W2', kc, claws), kc_of, grp)
    assert v['between'] < 1e-20 and v['among'] < 1e-20 and v['within'] > 0
    v = G.variance_split(np.ones(len(claws)), kc_of, grp)
    assert v['total'] == 0
    v2 = G.variance_split(3 * G.shape99('K1', kc, claws), kc_of, grp)          # relative to the squared mean
    assert np.isclose(v2['among'], G.variance_split(G.shape99('K1', kc, claws), kc_of, grp)['among'])


@pytest.mark.skipif(not _cache(), reason="needs E96's cached per-KC and per-claw tables (viki's outputs)")
def test_the_real_tables_identities_and_the_entrys_probe():
    kc, claws = G.cached_tables('real')
    g = pd.Series(G.claw_groups(kc.n_claws.to_numpy())).value_counts()
    assert g[11] == 45 and 12 not in g.index                                   # the 10 KCs with 12 claws join those of 11
    r = G.shape_identities(kc, claws)
    assert r['passes'] and np.isclose(r['k4_zero_kcs_share']['offclaw'], 0.069, atol=0.0005)
    rows = {}
    for name in G.SAVED_SHAPES:
        rows[name] = G.shape_row(G.saved_shape(name, kc, claws, 'real'), kc, claws)
    for key, (sh, c) in G.NEW_SHAPES.items():
        rows[key] = G.shape_row(G.shape99(sh, kc, claws, c or 'offclaw'), kc, claws)
    probe = {'offclaw_ped': (0.04, 0.23, 0, 0.67), 'perm_claw': (0.04, 0.23, 0, 0.67), 'perm_all': (0.26, 0.37, 0, 0.0),
             'allsyn': (0.01, 0.12, 0.18, 0.87), 'offclaw': (0.06, 0.39, 0, 0.47), 'connectome': (0.02, 0.22, 0.77, 1.10),
             'K1': (0, 0.31, 0, 0.98), 'K2': (0.06, 0, 0, 0.48), 'K3': (0.06, 0.10, 0, 0.51), 'K4': (0.00, 0.06, 0, 1.00),
             'W1': (0.02, 0.22, 0, 1.10), 'W2': (0, 0, 0.91, 1.00), 'T1': (0.01, 0.12, 0, 0.86), 'T2': (0.01, 0.14, 0, 0.91)}
    for k, (b, a, w, s) in probe.items():
        r = rows[k]
        assert abs(r['between'] - b) <= 0.006 and abs(r['among'] - a) <= 0.006 and abs(r['within'] - w) <= 0.006, (k, r)
        assert abs(r['loglog_slope_n'] - s) <= 0.051, (k, r['loglog_slope_n'])


# ---------------------------------------------------------------------------
# the two descriptions

def test_the_descriptions_recover_exact_coefficients():
    rng = np.random.default_rng(3)
    var = pd.DataFrame(rng.uniform(0, 0.8, (7, 3)), columns=list(G.PARTS3), index=list(G.SAVED7))
    var['total'] = var.sum(axis=1)
    means = 0.18 + 0.05 * var.between - 0.09 * var.among - 0.05 * var.within
    d = G.fit_descriptions(var, means)
    assert d['ii']['rms'] < 1e-12 and np.isclose(d['ii']['coef']['among'], -0.09) and d['ii']['n_coef'] == 4
    assert d['i']['rms'] > 0 and d['i']['n_points'] == 7
    p1, p2 = G.predict(d, dict(between=0.0, among=0.31, within=0.0, total=0.31))
    assert np.isclose(p2, 0.18 - 0.09 * 0.31)
    assert np.isclose(p1, d['i']['coef']['intercept'] + d['i']['coef']['total'] * 0.31)


# ---------------------------------------------------------------------------
# part C: the shares, the classes, the rule by wiring, the bound, the KCs with none

def test_the_classes():
    assert G.class_between(0.1, 0.7) == 'inequality' and G.class_between(0.7, 0.1) == 'slope'
    assert G.class_between(0.1, 0.2) == 'neither' and G.class_between(0.3, 0.7) == 'both'
    assert G.class_between(0.6, 0.6) == 'both' and G.class_between(0.24, 0.6) == 'inequality'
    assert G.class_between(0.25, 0.6) == 'both' and G.class_between(np.nan, 0.6) == 'not read'
    assert G.class_within(0.7, 0.1) == 'placement' and G.class_within(0.1, 0.7) == 'amounts'
    assert G.class_within(0.1, 0.1) == 'neither' and G.class_within(0.5, 0.5) == 'both'
    s = G.shares_between(U=0.181, O=0.147, K1=0.150, K2=0.182)
    assert np.isclose(s['gap'], 0.034) and np.isclose(s['slope_share'], 0.003 / 0.034) and s['cls'] == 'inequality'
    assert G.shares_between(0.15, 0.16, 0.15, 0.15)['cls'] == 'no gap'
    assert G.shares_between(np.nan, 0.16, 0.15, 0.15)['cls'] == 'not read'
    w = G.shares_within(U=0.181, D=0.124, W1=0.164, W2=0.132)
    assert np.isclose(w['placement_share'], 0.040 / 0.057) and w['cls'] == 'placement'


def test_the_rule_by_wiring():
    b = G.by_wiring({'real': 'inequality', 'P1': 'inequality', 'P2': 'both', 'M1': 'inequality'}, 'inequality')
    assert b['holds_by_wiring'] and b['n_same'] == 3 and b['reported_as'] == 'by wiring'
    b = G.by_wiring({'real': 'inequality', 'P1': 'both', 'P2': 'both', 'M1': 'inequality'}, 'inequality')
    assert not b['holds_by_wiring'] and b['reported_as'] == "S's alone"
    assert not G.by_wiring({w: 'not read' for w in ('real', 'P1', 'P2', 'M1')}, 'not read')['holds_by_wiring']


def test_the_bound_and_the_kcs_with_none():
    U, O, K2 = 0.181, 0.147, 0.180
    LO = K2 - O
    sq = G.evenness_bound(U, O, K2, K2 - LO / 4, cv_O=0.6, cv_K3=0.3)        # the loss with the square of the inequality
    assert sq['applies'] and np.isclose(sq['r3'], 0.75) and np.isclose(sq['exponent'], 2.0)
    assert np.isclose(sq['cv_at_loss'], 0.6 * np.sqrt(0.01 / LO)) and not sq['extrapolated']
    lin = G.evenness_bound(U, O, K2, K2 - LO / 2, cv_O=0.6, cv_K3=0.3)       # in proportion to it
    assert np.isclose(lin['r3'], 0.5) and np.isclose(lin['exponent'], 1.0) and np.isclose(lin['cv_at_loss'], 0.6 * 0.01 / LO)
    small = G.evenness_bound(U, O, O + 0.005, O + 0.004, cv_O=0.6, cv_K3=0.3)
    assert not small['applies']                                               # K2 recovers less than a quarter
    flat = G.evenness_bound(U, O, K2, K2 + 0.001, cv_O=0.6, cv_K3=0.3)        # K3 above K2: no power law
    assert not np.isfinite(flat['cv_at_loss']) and 'note' in flat
    assert G.none_verdict(0.181, 0.147, 0.160)['verdict'] == 'the KCs without an off-claw synapse carry the loss'
    assert G.none_verdict(0.181, 0.147, 0.176)['verdict'] == 'presence alone holds the mean'
    assert G.none_verdict(0.181, 0.147, 0.168)['verdict'] == 'between'


def _fake_readings(patterns, flagged=(), at_bound=(), T=None, raw=None, mean=None, per_target=None, units=('S',),
                   unit_means=None):
    """Readings of made-up settings ({variant: pattern part on S, the same against every target unless `per_target`
    gives {variant: {target: pattern}}}), with T's ({variant: pattern}); raw-count losses and mean effects default to the
    reference's (as E96's test); `unit_means` {(variant, unit): mean effect} for single wirings."""
    rows = []
    for v, p in patterns.items():
        pv = G.parse_variant(v)
        for unit, pat in tuple((u, p) for u in units) + ((('T', T[v]),) if T and v in T else ()):
            me = (unit_means or {}).get((v, unit), (mean or {}).get(v, 0.177))
            base = dict(variant=v, arm=pv['arm'], part=pv['part'], screen=pv['screen'], shape=pv['shape'], counts=pv['counts'],
                        mode=pv['mode'], size=pv['size'], pool=pv['pool'], dose=pv['dose'], control=pv['control'],
                        construction=pv['construction'], equal=pv['equal'], unit=unit, n_wirings=4,
                        n_flagged=int(v in flagged), n_at_bound=int(v in at_bound), mean_effect=me)
            for t in ('halfA', 'halfB', 'halfC', 'halfD', 'z8', 'z6'):
                pt = (per_target or {}).get(v, {}).get(t, pat) if unit == 'S' else pat
                base[f'{t}_pattern'] = pt
                base[f'{t}_loss'] = pt + 0.01
                base.update({f'{t}_amplitude': 0.6 * pt, f'{t}_shape': 0.8 * pt, f'{t}_sd_ratio': 0.9,
                             f'{t}_pearson': 0.5, f'{t}_rho': 0.45})
            cal, rw = dict(base, readout='cal'), dict(base, readout='raw', z8_loss=(raw or {}).get(v, 0.078))
            rows += [cal, rw]
            if unit == 'T':
                for w in G.T_WIRINGS:
                    rows += [dict(cal, unit=w), dict(rw, unit=w)]
    return pd.DataFrame(rows)


def test_part_c_verdicts_on_made_up_readings(tmp_path):
    means = {'D1_uniform_d1': 0.181, 'B1_offclaw_d1': 0.147, 'K1_offclaw_d1': 0.150, 'K2_offclaw_d1': 0.180,
             'D1_uniform_d2': 0.180, 'B1_offclaw_d2': 0.136, 'K1_offclaw_d2': 0.140, 'K2_offclaw_d2': 0.176,
             'B2_offclawped_d1': 0.161, 'K1_offclawped_d1': 0.170, 'K2_offclawped_d1': 0.175,
             'D2_connectome_d1': 0.124, 'W1_d1': 0.170, 'W2_d1': 0.130, 'K3_offclaw_d1': 0.172, 'K4_offclaw_d1': 0.178,
             'B6_allsyn_d0.5': 0.159, 'B6_allsyn_d1': 0.156, 'B6_allsyn_d2': 0.151, 'T1_d0.5': 0.170, 'T1_d1': 0.169,
             'T1_d2': 0.165, 'T2_d0.5': 0.171, 'T2_d1': 0.168, 'T2_d2': 0.166}
    um = {}
    for w in G.S_WIRINGS:
        for v, m in means.items():
            um[(v, w)] = m
    um[('K1_offclaw_d1', 'P2')] = 0.170                                       # P2 gives 'both'
    rd = _fake_readings({v: 0.054 for v in means}, mean=means, units=('S',) + G.S_WIRINGS, unit_means=um,
                        at_bound=('K3_offclaw_d1',))
    a1 = pd.DataFrame([dict(wiring='real', setting='offclaw', within_cv_median=0.58),
                       dict(wiring='real', setting='K3', within_cv_median=0.29)])
    cv, cs = G.c_verdicts(rd, a1=a1, out=tmp_path)
    b = cv['between_calyx_d1']
    assert b['S']['cls'] == 'inequality' and b['per_wiring']['P2'] == 'both' and b['holds_by_wiring']
    assert cv['between_calyx_d2']['agrees_with_class'] and not cv['between_peduncle_d1']['agrees_with_class']
    assert cv['within_d1']['S']['cls'] == 'placement'
    assert cv['none']['verdict'] == 'presence alone holds the mean'
    assert not np.isfinite(cv['bound']['K3'])                                 # K3 at a bound: not part of a share
    t = {r['setting']: r for r in cv['t_vs_b6']}
    assert np.isclose(t['T1_d1']['d_mean_cal'], 0.013)
    assert set(cs.kind) == {'between', 'within'} and set(cs.unit) == {'S'} | set(G.S_WIRINGS)


def test_the_choice_and_the_equal_rule():
    ref = G.REF
    c = G.choose99(_fake_readings({ref: 0.050, 'T1_d1': 0.0505}), 'A')
    assert c['chosen'] == ref and not c['confirm']
    rd = _fake_readings({ref: 0.050, 'K2_offclaw_d1': 0.030, 'W2_d1': 0.031, 'D1_uniform_d1': 0.030, 'B6_allsyn_d1': 0.030,
                         'T1_d2': 0.040, 'T2_d1': 0.047, 'W1_d1': 0.0475}, flagged=('T1_d2',))
    c = G.choose99(rd, 'B')
    assert c['chosen'] == 'T2_d1' and c['confirm'] and np.isclose(c['gain'], 0.003) and c['control'] == 'D1_uniform_d1'
    assert 'T1_d2' in json.loads(c['excluded'])                               # constructions and corners never chosen
    eq = {e['variant']: e for e in G.equal_rule99(_fake_readings(
        {ref: 0.050, 'T1_d1': 0.052, 'T2_d1': 0.0545, 'W1_d1': 0.0495, 'K1_offclaw_d1': 0.050, 'W2_d1': 0.050},
        mean={'W1_d1': 0.160}))}
    assert set(eq) == {'T1_d1', 'T2_d1', 'W1_d1'}
    assert eq['T1_d1']['equal'] and not eq['T2_d1']['equal'] and not eq['W1_d1']['equal']    # W1 loses the mean


def test_part_d_rule():
    pats = {'D1_uniform_d0.5': 0.0535, 'D1_uniform_d1': 0.0545, 'D1_uniform_d2': 0.0539,
            'T1_d0.5': 0.0540, 'T1_d1': 0.0550, 'T1_d2': 0.0545, 'T2_d0.5': 0.0544, 'T2_d1': 0.0561, 'T2_d2': 0.0570}
    mean = {'T1_d0.5': 0.170, 'T1_d1': 0.169, 'T1_d2': 0.160, 'T2_d0.5': 0.171, 'T2_d1': 0.172, 'T2_d2': 0.168}
    r = G.d_rule99(_fake_readings(pats, mean=mean))
    meets = {x['variant']: x['meets'] for x in r['rows']}
    assert meets['T1_d0.5'] and meets['T1_d1'] and not meets['T1_d2']          # T1_d2 loses the mean
    assert meets['T2_d0.5'] and not meets['T2_d1'] and not meets['T2_d2']       # T2_d1, d2 above the uniform's + 0.0015
    assert r['runs'] and r['best'] == {'T1': 'T1_d0.5', 'T2': 'T2_d0.5'}
    assert r['component'] == 'T2' and r['mode'] == 'allsyn_even_calyx'        # 0.0544 within 0.0005 of 0.0540
    r = G.d_rule99(_fake_readings(dict(pats, **{'T2_d0.5': 0.0546}), mean=mean))
    assert r['component'] == 'T1' and r['mode'] == 'allsyn_even'              # 0.0006 above: the lower
    r = G.d_rule99(_fake_readings(pats, mean=mean, flagged=('T1_d0.5', 'T1_d1'), at_bound=('T2_d0.5',)))
    assert not r['runs'] and r['component'] is None
    r = G.d_rule99(_fake_readings(pats, mean=mean, flagged=('T1_d0.5', 'T1_d1')))
    assert r['component'] == 'T2'                                               # only T2 meets it


# ---------------------------------------------------------------------------
# A2's file is written once, before part B

def test_a2_is_kept_once_part_b_has_rows(tmp_path, monkeypatch):
    kc, claws = _toy()
    rows = []
    for name in ('uniform', 'connectome', 'offclaw', 'offclaw_ped', 'allsyn', 'perm_claw', 'perm_all'):
        x = (np.ones(len(claws)) if name == 'uniform' else G.shape99('W1', kc, claws) * G.shape99('W2', kc, claws)
             if name == 'connectome' else G.shape99('K3', kc, claws, 'offclaw_ped' if name != 'offclaw' else 'offclaw', 1.0)
             if name in ('offclaw', 'offclaw_ped', 'perm_claw') else G.shape99('T1', kc, claws) if name == 'allsyn'
             else G.shape99('K1', kc, claws))
        rows.append(dict(wiring='real', setting=name, **G.shape_row(x, kc, claws)))
    for key, (sh, c) in G.NEW_SHAPES.items():
        rows.append(dict(wiring='real', setting=key, **G.shape_row(G.shape99(sh, kc, claws, c or 'offclaw'), kc, claws)))
    a1 = pd.DataFrame(rows)
    means = dict(zip(G.SAVED7, (0.181, 0.161, 0.159, 0.158, 0.156, 0.147, 0.124)))
    rd96 = pd.DataFrame([dict(variant=G.SAVED_SHAPES[s], unit='S', readout='cal', mean_effect=m) for s, m in means.items()])
    p = tmp_path / 'readings96.csv'
    rd96.to_csv(p, index=False)
    monkeypatch.setattr(G, '_e96_readings', lambda: rd96)
    monkeypatch.setattr(G.H, 'results_dir', lambda e: tmp_path if e != G.EXPERIMENT else tmp_path)
    (tmp_path / 'readings.csv').write_text(rd96.to_csv(index=False))
    a = G.a2_predictions(tmp_path, a1)
    assert a['written_at'] and set(a['predictions']) == set(G.PRED_KEYS)
    assert {'between', 'replicate_peduncle', 'within', 'none', 'bound', 'd_rule_mean'} <= set(a['classes']['ii'])
    d = pd.read_csv(tmp_path / 'a2_predictions.csv')
    assert set(d.kind) == {'saved', 'new'} and np.isclose(d[d.kind == 'saved'].observed, list(means.values())).all()
    (tmp_path / 'runs_B.csv').write_text('tag\nx\n')
    first = (tmp_path / 'a2_predictions.json').read_text()
    a2 = G.a2_predictions(tmp_path, a1.assign(total=a1.total * 2))
    assert a2['written_at'] == a['written_at'] and (tmp_path / 'a2_predictions.json').read_text() == first


# ---------------------------------------------------------------------------
# A3: who responds

def test_the_reading_by_fifths_and_of_the_claws():
    kc, claws = _toy()
    kc = kc.assign(kc_type=np.where(np.arange(len(kc)) % 2, 'ab', 'g'))
    ids = kc[KC].to_numpy()
    odors = [f'o{i}' for i in range(8)]
    rng = np.random.default_rng(5)
    on = pd.DataFrame(rng.poisson(1.0, (len(ids), 8)), index=pd.MultiIndex.from_arrays([ids, kc.kc_type], names=[KC, 'kc_type']),
                      columns=odors)
    off = on * 2
    f = G.kc_response_frame(on, off, kc, claws, G.shape99('K3', kc, claws, frac=1.0))
    assert set(f.fifth_offclaw) == {1, 2, 3, 4, 5} and (f.fifth_offclaw.value_counts() == 12).all()
    assert np.allclose(f.rel_offclaw, (kc.no_claw / kc.n_claws).to_numpy()
                       / G.group_means(G.offclaw_amounts(kc), kc.set_index(KC).n_claws, G.group_of_kc(kc)).to_numpy())
    assert (f.fifth_own == f.fifth_offclaw).all()                              # the off-claw component ranks the same
    rows = pd.DataFrame(G.summarize_responses(f, dict(variant='x', unit='real')))
    for by in ('fifth_offclaw', 'claw_class', 'kc_type'):
        q = rows[rows.by == by]
        assert np.isclose(q.share_on.sum(), 1.0) and q.n_kcs.sum() == len(kc)
    assert np.isclose(rows[rows.by == 'all'].odors_on.iloc[0], (on.to_numpy() > G.RAW_CRIT).sum(axis=1).mean())
    fu = G.kc_response_frame(on, off, kc, claws, np.ones(len(claws)))
    assert (fu.fifth_own == 0).all() and 'fifth_own' not in set(pd.DataFrame(G.summarize_responses(fu, {})).by)
    # the claws: the most and least weighted claw of each KC with an uneven component
    x = G.shape99('W2', kc, claws)
    ix = pd.MultiIndex.from_frame(claws[[KC, 'claw_id']])
    st_on = pd.DataFrame(np.outer(x, np.ones(3)), index=ix, columns=['a', 'b', 'c'])
    st_off = st_on * 2 + 1
    c = G.claw_extremes(st_on, st_off, x, claws)
    assert set(c.which) == {'most', 'least'} and (c.groupby('kc_id').size() == 2).all()
    most = c[c.which == 'most'].set_index('kc_id').weight
    least = c[c.which == 'least'].set_index('kc_id').weight
    assert (most > least).all() and np.allclose(c.mean_on, c.weight)
    s = pd.DataFrame(G.summarize_claws(c, {}))
    assert set(s.which) == {'most', 'least'} and 'all' in set(s.claw_class)


# ---------------------------------------------------------------------------
# E10_joint_fit.py's new --aplkc-flat modes

def test_e10_new_modes():
    assert J.aplkc_flat_tag(J.aplkc_flat_spec('allsyn_even', 3.8)) == 'flatt1_'
    assert J.aplkc_flat_tag(J.aplkc_flat_spec('allsyn_even_calyx', 3.8)) == 'flatt2_'
    assert len(set(J.APLKC_FLAT_TAGS.values())) == len(J.APLKC_FLAT_KINDS)
    assert set(G.T_MODE.values()) <= set(J.APLKC_FLAT_COUNT_KINDS)
    assert J.aplkc_flat_spec(None) is None                                      # off by default


def test_e10_new_modes_set_the_simulators_weights():
    kc, claws = _toy()
    idx = pd.MultiIndex.from_frame(claws[[KC, 'claw_id']])
    per_claw = pd.Series(claws['kept'].to_numpy(), index=idx)
    wn = pd.Series(claws['wnorm'].to_numpy(), index=idx)
    wk = pd.Series(np.random.default_rng(2).uniform(0, 2, len(idx)), index=idx)
    J._FLAT_COUNTS.clear()
    key = tuple(map(tuple, idx.to_frame(index=False)[[KC, 'claw_id']].to_numpy().tolist()))
    J._FLAT_COUNTS[key] = (kc.set_index(KC), per_claw)
    try:
        for s, mode in G.T_MODE.items():
            held = wk.to_numpy().reshape(-1, 1).copy()

            def make_sim(w):
                wAPLKC_unscaled = w

                def simulate(params):
                    return wAPLKC_unscaled
                return simulate
            sim = make_sim(held)
            p = J.aplkc_flat_param(wn, wk, J.aplkc_flat_spec(mode, 2.0, 0.0))
            p['set'](sim, 0.0)
            assert np.array_equal(held.ravel(), wk.to_numpy())                  # size 0: the family's weights exactly
            p['set'](sim, 1.5)
            assert np.allclose(held.ravel(), wk.to_numpy() + 1.5 * 2.0 * wn.mean() * G.shape99(s, kc, claws))
    finally:
        J._FLAT_COUNTS.clear()


# ---------------------------------------------------------------------------
# part D's commands and tags

@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the fits of C2's family)")
def test_part_d_commands_and_tags():
    c = G.d_command('real', 'T2', 'own', 3.84)
    assert c[c.index('--aplkc-flat') + 1] == 'allsyn_even_calyx' and c[c.index('-k') + 1] == '2'
    assert np.isclose(float(c[c.index('--aplpn-hold') + 1]), 0.1 * G.E94.aplpn_fit_multiple('real'))
    assert c[c.index('--init-params') + 1] == f'results/{G.EXPERIMENT}/d_start_params.json:real_T2_own'
    for a, v in (('--rate-target', 'ctrl_frac'), ('--loss', 'apl'), ('--apl-readout', 'both'), ('--tag', '_g1'),
                 ('--maxfev', '200'), ('--experiment', G.EXPERIMENT), ('--aplkc-flat-max', '8')):
        assert c[c.index(a) + 1] == v
    assert '--claw-inh-pool' not in c
    assert G.d_command('P1', 'T1', 'e96', 2.5)[G.d_command('P1', 'T1', 'e96', 2.5).index('--tag') + 1] == '_g1e'
    assert G.d_tag('real', 'T1', 'own') == 'rate0.068_aplboth_aplpnhold_flatt1_k2_g1'
    assert G.d_tag('P1', 'T2', 'e96') == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_aplpnhold_flatt2_k2_g1e'


@pytest.mark.skipif(not _shared() or not (H.results_dir('E96_offclaw_apl') / 'd_scores.csv').exists(),
                    reason="needs E96's refits in the shared directory")
def test_the_e96_start():
    k = G.e96_kept_refit('real')
    r = k['row']
    s = G.e96_start('real', float(r['aplkc_flat_unit']))
    assert s['aplkc_flat'] == float(r['aplkc_flat_size']) and s['wKCKC_scale'] == float(r['kc_kc_inh']) and s['start'] == 'e96'
    with pytest.raises(RuntimeError):
        G.e96_start('real', float(r['aplkc_flat_unit']) * 1.01)               # another unit is refused


# ---------------------------------------------------------------------------
# the runner's rows on saved data (no simulation; E96's part 0 cache and the shared directory)

def _saved_model(monkeypatch):
    """A stand-in for a StaticModel on the real wiring at its fit (E96's test's), with E96's `reading96` replaced by a
    reading of the saved APL-on run with the fit's stage matrices, the counts from E96's cache."""
    from types import SimpleNamespace
    E83 = G.E83
    f = E83.saved_fit('real')
    W, _ = G.E94.wiring_W('real')
    kc, claws = G.cached_tables('real')
    idx = pd.MultiIndex.from_frame(claws[[KC, 'claw_id']])
    wn = pd.Series(claws.wnorm.to_numpy(), index=idx)
    per_claw = pd.Series(claws.kept.to_numpy(), index=idx)
    per_kc = kc.set_index(KC)[['n_claws', 'kept', 'kept_unassigned', 'dropped_claw', 'no_claw', 'no_claw_nan',
                               'peduncle']].sort_index()
    m = object.__new__(G.E95.StaticModel)
    m.label, m.folder, m.tag, m.p, m.saved, m.K = 'real', f['folder'], f['tag'], dict(f['p']), f['saved'], f['K']
    m.target_on, m.target_off, m.W = f['target_on'], f['target_off'], W
    m.wAPLKC = pd.Series(claws.w_fam.to_numpy(), index=idx)
    m.rv = SimpleNamespace(kc=SimpleNamespace(thr=np.ones(len(kc))))
    seen = {}

    def fake_reading(model, pv, wk, label=''):
        seen[label] = (dict(pv), None if wk is None else wk.copy())
        params = dict(model.p, wAPLPN_scale=model.p['wAPLPN_scale'] * pv['dose'])
        return dict(on=f['saved']['on'], off=f['saved']['off'], factor=1.2, at_bound=False, n_sims=3, params=params,
                    seconds=0.0, history=[], stages=dict(G.E94.fit_stage('real', 'on')),
                    gain=dict(L=np.zeros(f['saved']['on'][1].shape), row={}), bracket=G.APL_BRACKET, ext=None, whole=False)
    key = tuple(map(tuple, idx.to_frame(index=False)[[KC, 'claw_id']].to_numpy().tolist()))
    J._FLAT_COUNTS[key] = (per_kc, per_claw)
    monkeypatch.setattr(G.E96, 'reading96', fake_reading)
    monkeypatch.setattr(G.E96, 'wiring_wnorm', lambda w: wn)
    monkeypatch.setattr(G.E96, 'wiring_counts', lambda w='real': (per_kc, per_claw))
    monkeypatch.setattr(G.E95, 'd_unit', lambda model, L=None: 3.0)
    return m, seen, (wn, per_kc, per_claw, kc, claws)


@pytest.mark.skipif(not _shared() or not _cache(), reason="needs the shared directory and E96's per-KC tables")
def test_run_point_rows_on_saved_data(monkeypatch):
    m, seen, (wn, per_kc, per_claw, kc, claws) = _saved_model(monkeypatch)
    saved = {}
    monkeypatch.setattr(G, 'save_run99', lambda tag, o, sub='runs', **k: saved.update({tag: k}))
    try:
        vs = [v for b in G.BLOCKS.values() for v in b] + ['D1_uniform_d1', 'D2_connectome_d1', 'B1_offclaw_d1',
                                                          'B6_allsyn_d1']
        for v in vs:
            row, _ = G.run_point99(m, v, part='test', cc={'queue_commit': 'x'})
            pv = G.parse_variant(v)
            tag = f'{v}_real'
            assert row['variant'] == v and row['tag'] == tag and row['component_unit'] == 3.0 and row['source'] == 'E99'
            assert row['dose'] == 0.1 and row['size'] == pv['size'] and row['shape'] == pv['shape']
            assert all(k in row for k in ('rate_on', 'realism_flag', 'max_count_on', 'K_run', 'bouton_share', 'comp_share'))
            rpv, wk = seen[tag]
            assert rpv['pool'] == 0.0 and rpv['dose'] == 0.1
            add = wk.to_numpy() - m.wAPLKC.to_numpy()
            assert np.isclose(add.mean(), pv['size'] * 3.0 * wn.mean())          # d x unit x the connectome's mean
            sh = saved[tag]['component_claw']['shape'].to_numpy()
            assert np.allclose(add, pv['size'] * 3.0 * wn.mean() * sh, rtol=1e-10)
            comp = saved[tag]['component']
            assert np.isclose(comp.component_sum.sum(), add.sum())
            if pv['part'] == 'B':
                assert np.allclose(sh, G.setting_shape(v, kc, claws), rtol=1e-12)    # the live counts are the cache's
        # the saved corners through this path are E96's components
        for v, mode in (('B1_offclaw_d1', 'offclaw'), ('B6_allsyn_d1', 'allsyn'), ('D2_connectome_d1', 'connectome')):
            assert np.allclose(saved[f'{v}_real']['component_claw']['shape'].to_numpy(), claws[f'comp_{mode}'].to_numpy(),
                               rtol=1e-12)
    finally:
        J._FLAT_COUNTS.clear()


# ---------------------------------------------------------------------------
# the scoring path on saved runs (no simulation)

@pytest.mark.skipif(QUICK or not _shared() or not (H.outputs_dir('E95_static_apl') / 'runs').exists(),
                    reason="needs E95's and E96's saved runs (viki's outputs) and the shared directory; not with QUICK=1")
def test_the_readings_on_saved_runs(tmp_path):
    """The scoring path on E95's and E96's saved runs standing in for this experiment's points on S."""
    stand = {'K1_offclaw_d1': ('E95', 'D1_uniform_d1'), 'T1_d1': ('E96', 'B6_allsyn_d1')}
    rows = []
    for v, (src, e) in stand.items():
        for w in G.S_WIRINGS:
            rows.append(dict(tag=f'{e}_{w}', variant=v, wiring=w, source=src, at_bound=False, realism_flag='', factor=1.0))
    pd.DataFrame(rows).to_csv(tmp_path / 'runs_B.csv', index=False)
    rd = G.readings99(tmp_path, readouts=('cal', 'raw'))
    s = rd[(rd.unit == 'S') & (rd.readout == 'cal')].set_index('variant')
    assert set(s.index) >= {G.REF, G.D0, 'K1_offclaw_d1', 'T1_d1'} | set(G.CORNERS['E95']) | set(G.CORNERS['E96'])
    assert np.isclose(s.loc[G.REF, 'z8_pattern'], G.UNCHANGED_S['all'], atol=5e-5)
    e96 = pd.read_csv(H.results_dir('E96_offclaw_apl') / 'readings.csv')
    e96 = e96[(e96.unit == 'S') & (e96.readout == 'cal')].drop_duplicates('variant').set_index('variant')
    for v in G.CORNERS['E95'] + G.CORNERS['E96']:
        assert np.isclose(s.loc[v, 'mean_effect'], e96.loc[v, 'mean_effect'], rtol=1e-9), v
    assert np.isclose(s.loc['T1_d1', 'z8_pattern'], e96.loc['B6_allsyn_d1', 'z8_pattern'], rtol=1e-9)
    cv, _ = G.c_verdicts(rd, out=tmp_path)
    b = cv['between_calyx_d1']['S']
    assert np.isclose(b['U'], e96.loc['D1_uniform_d1', 'mean_effect']) and np.isclose(b['gap'], 0.0348, atol=0.001)
    assert np.isclose(b['slope_share'], 1.0)                                  # K1 stands in with the uniform's runs
    v = G.verdict99(tmp_path, rd)
    assert set(v['choice']) == set(G.HALVES) and 'D_rule' in v and not v['D_rule']['runs']
    assert v['classes']['between'] == 'not read' and not v['classes']['between_by_wiring']      # K2 has no stand-in
    text = G.comparison_table(tmp_path)
    assert '| K1_offclaw_d1 |' in text and 'Verdict' in text and "Part C" in text


def test_the_figure_on_made_up_readings(tmp_path):
    vs = {G.REF: 0.050, 'D1_uniform_d1': 0.0545, 'B1_offclaw_d1': 0.0541, 'K1_offclaw_d1': 0.054, 'T1_d1': 0.053, 'D0': 0.055}
    _fake_readings(vs, mean={'D1_uniform_d1': 0.181, 'B1_offclaw_d1': 0.147}, at_bound=('T1_d1',)).to_csv(
        tmp_path / 'readings.csv', index=False)
    p = G.figure(tmp_path)
    assert p is not None and p.exists() and p.stat().st_size > 1000
