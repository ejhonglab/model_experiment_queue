#!/usr/bin/env python3
"""Tests for E133's transplant builders (experiments/E133_builders.py; model_experiment_queue) on made-up calyces: the
identity (donor = host returns the tables bit for bit, every arm and option), each arm changing only its own feature,
C's claw counts and kept KC sums, S's shares, A's and B's kept totals and transplanted profiles, K's block shares and
mean row sum, the step-3 options touching only their subset, ALL as the chain, the rank rule, the arm strings, step 2′'s
permutation with APL's weights and E59's shuffle, and the constants shared with E59 and E97. No simulations. Needs the
shared directory (E97 imports E83, which reads the data-target registry at import); skipped without it."""

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
B = pytest.importorskip('experiments.E133_builders')

GLOMS = ('DM1', 'DM2', 'DP1m', 'VA2', 'VM7d', 'DL1', 'VC3')
WN = ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')


def make_calyx(seed, n_kcs=None, claws=None, gloms=GLOMS, boutons=2, calyx='test-host', kc_id0=1000, skew=0.0):
    """A made-up calyx's tables (E133_calyx_transplants.calyx_tables' keys): three KC types, one-hot claw x bouton
    wiring, positive APL weights with some zeros, a KC x KC matrix with mean row sum 1."""
    rng = np.random.default_rng(seed)
    n_kcs = n_kcs or {'g': 14, "a'b'": 6, 'ab': 12}
    claws = claws or {'g': (5, 10), "a'b'": (2, 5), 'ab': (3, 7)}
    cols = pd.MultiIndex.from_tuples([(g, 100 + i, b) for i, g in enumerate(gloms) for b in range(boutons)],
                                     names=['glomerulus', 'pn_id', 'bouton_id'])
    p = np.exp(skew * np.arange(len(cols)) / len(cols))
    p = p / p.sum()
    rows, cs = [], []
    kid = kc_id0
    for t, n in n_kcs.items():
        for _ in range(n):
            lo, hi = claws[t]
            for c in range(int(rng.integers(lo, hi + 1))):
                rows.append((kid, c, t))
                cs.append(int(rng.choice(len(cols), p=p)))
            kid += 1
    idx = pd.MultiIndex.from_tuples(rows, names=[H.KC_ID_LEVEL, 'claw_id', H.KC_TYPE_LEVEL])
    v = np.zeros((len(rows), len(cols)))
    v[np.arange(len(rows)), cs] = 1.0
    W = pd.DataFrame(v, index=idx, columns=cols)
    akc = rng.gamma(2.0, 1.0, len(rows)) * (rng.random(len(rows)) > 0.15)
    kca = rng.gamma(3.0, 1.0, len(rows)) * (rng.random(len(rows)) > 0.1)
    wn = {'wAPLKC': pd.Series(akc / akc.mean(), index=idx), 'wKCAPL': pd.Series(kca / kca.mean(), index=idx),
          'wAPLPN': pd.Series(rng.gamma(2.0, 1.0, len(cols)), index=cols),
          'wPNAPL': pd.Series(rng.gamma(2.0, 1.0, len(cols)), index=cols)}
    # one KC without APL->KC synapses on any claw
    first = W.index.get_level_values(H.KC_ID_LEVEL) == kc_id0 + 1
    wn['wAPLKC'][first] = 0.0
    ids = pd.unique(W.index.get_level_values(H.KC_ID_LEVEL))
    A = rng.poisson(0.6, (len(ids), len(ids))).astype(float)
    np.fill_diagonal(A, 0.0)
    M = pd.DataFrame(A / (A.sum() / len(ids)), index=pd.Index(ids, name='kc_id_post'), columns=pd.Index(ids, name='kc_id_pre'))
    breadth = pd.Series(rng.integers(0, 30, len(ids)).astype(float), index=pd.Index(ids, name=H.KC_ID_LEVEL))
    return dict(calyx=calyx, seed=0, label=calyx, orn=None, base={'fixed_thr': 1.0}, wnorm=wn, W=W, breadth=breadth, kckc=M,
                meta=dict(wiring='own', wiring_seed=0))


@pytest.fixture(scope='module')
def pair():
    host = make_calyx(1, calyx='test-host')
    donor = make_calyx(2, n_kcs={'g': 10, "a'b'": 8, 'ab': 15}, claws={'g': (3, 7), "a'b'": (3, 6), 'ab': (4, 8)},
                       calyx='test-donor', kc_id0=5000, skew=1.5)
    return host, donor


def _changed(T, host):
    out = set()
    if not B._same(T['W'], host['W']):
        out.add('W')
    for k in WN:
        if not B._same(T['wnorm'][k], host['wnorm'][k]):
            out.add(k)
    if not B._same(T['kckc'], host['kckc']):
        out.add('kckc')
    return out


def _kc_sums(T, k):
    return T['wnorm'][k].groupby(level=H.KC_ID_LEVEL, sort=False).sum()


def _small1(monkeypatch):
    """E97's small-glomerulus rule at 1 claw (the made-up calyces are small)."""
    monkeypatch.setattr(B, 'SMALL_CLAWS', 1)
    orig = B.glomerulus_shares
    monkeypatch.setattr(B, 'glomerulus_shares', lambda h, d, seed=0, glomeruli=None, small=1: orig(h, d, seed, glomeruli, small))
    monkeypatch.setitem(B.BUILDERS, 'S', B.glomerulus_shares)


# ---------------------------------------------------------------------------
# the identity

ALL_ARMS = B.ARMS + ('C-g', 'C-ab', 'S-three', 'A-three', 'B-three', 'A-aplkc', 'A-kcapl', 'B-aplpn', 'B-pnapl')


@pytest.mark.parametrize('arm', ALL_ARMS)
@pytest.mark.parametrize('seed', (0, 1, 3))
def test_identity(arm, seed, pair, monkeypatch):
    """Donor equal to host: every arm and option returns the tables bit for bit (W, the four weight vectors, kckc) and
    keeps the breadth."""
    _small1(monkeypatch)
    monkeypatch.setattr(B, 'three_glomeruli', lambda t: ['DM1', 'VA2'])
    host, _ = pair
    assert B.identity_holds(arm, host, seed) == {}


def test_identity_without_a_matrix(pair, monkeypatch):
    """A hemibrain host has no matrix in its tables (olfsysm loads it): K with itself as donor leaves it None."""
    host, _ = pair
    hb = dict(host, calyx='hemibrain-argmax', kckc=None)
    other = make_calyx(2, calyx='hemibrain-x', kc_id0=5000)
    mats = {id(host['W']): host['kckc'], id(other['W']): other['kckc']}
    monkeypatch.setattr(B, 'hemibrain_kckc', lambda W: mats[id(W)])
    assert B.identity_holds('K', hb, 0) == {}
    T = B.build('K', hb, dict(other, kckc=None), 0)
    assert T['kckc'] is not None and T['meta']['kckc_sha256'] and T['kckc'].index.equals(host['kckc'].index)


# ---------------------------------------------------------------------------
# each arm's own feature

@pytest.mark.parametrize('arm,expected', [('C', {'W', 'wAPLKC', 'wKCAPL'}), ('S', {'W'}), ('A', {'wAPLKC', 'wKCAPL'}),
                                          ('B', {'wAPLPN', 'wPNAPL'}), ('K', {'kckc'})])
def test_arm_changes_only_its_feature(arm, expected, pair, monkeypatch):
    _small1(monkeypatch)
    host, donor = pair
    T = B.build(arm, host, donor, seed=0)
    assert _changed(T, host) == expected
    assert T['orn'] is host['orn'] and T['base'] is host['base']
    assert (T['breadth'] is None) == ('W' in expected)
    assert T['meta']['transplant']['arm'] == arm and arm in T['meta']['transplant']['parts']


def test_builders_do_not_change_their_inputs(pair, monkeypatch):
    _small1(monkeypatch)
    host, donor = pair
    snap = {k: host['wnorm'][k].copy() for k in WN}
    W0, M0 = host['W'].copy(), host['kckc'].copy()
    B.build('ALL', host, donor, seed=2)
    assert host['W'].equals(W0) and host['kckc'].equals(M0)
    assert all(host['wnorm'][k].equals(snap[k]) for k in WN)


# ---------------------------------------------------------------------------
# C

def test_C_counts_and_kept_sums(pair):
    host, donor = pair
    T = B.build('C', host, donor, seed=0)
    n0, t0, ids0 = B.claw_counts(host['W'])
    n1, t1, ids1 = B.claw_counts(T['W'])
    nd, td, _ = B.claw_counts(donor['W'])
    assert np.array_equal(ids0, ids1) and np.array_equal(t0, t1)
    for t in set(t0):
        m, md = t0 == t, td == t
        # each type's counts are the donor's at the same rank
        want = np.sort(np.sort(nd[md])[((2 * np.arange(m.sum()) + 1) * md.sum()) // (2 * m.sum())])
        assert np.array_equal(np.sort(n1[m]), want)
        # the rank order is kept: a KC with more claws before has at least as many after (ties in random order)
        a0, a1 = n0[m], n1[m]
        for c in np.unique(a0)[:-1]:
            assert a1[a0 == c].max() <= a1[a0 > c].min()
    # each KC keeps its summed APL weights; the rows keep their boutons
    for k in ('wAPLKC', 'wKCAPL'):
        assert np.allclose(_kc_sums(T, k).to_numpy(), _kc_sums(host, k).to_numpy(), rtol=1e-12, atol=0)
    v0, v1 = host['W'].to_numpy(), T['W'].to_numpy()
    assert np.array_equal(v0 > 0, v1 > 0) and ((v1 > 0).sum(axis=1) == 1).all()
    off = np.isclose(v1.max(axis=1), B.SWITCHED_OFF_WEIGHT, rtol=0, atol=0)
    assert off.any() and (T['wnorm']['wAPLKC'].to_numpy()[off] == 0).all()
    assert set(np.unique(v1[v1 > 0])) <= {B.SWITCHED_OFF_WEIGHT} | set(range(1, 20))
    size = T['meta']['transplant']['parts']['C']
    assert size['claws_off'] == int(off.sum()) and size['claws_copied'] == int((v1[v1 >= 1] - 1).sum())
    assert size['off_share'] == pytest.approx(off.sum() / n0.sum())


def test_C_split_is_the_full_arm_on_one_type(pair):
    host, donor = pair
    full = B.build('C', host, donor, seed=1)
    g = B.build('C-g', host, donor, seed=1)
    t = np.asarray(host['W'].index.get_level_values(H.KC_TYPE_LEVEL))
    m = t == 'g'
    assert np.array_equal(g['W'].to_numpy()[m], full['W'].to_numpy()[m])
    assert np.array_equal(g['W'].to_numpy()[~m], host['W'].to_numpy()[~m])
    for k in ('wAPLKC', 'wKCAPL'):
        assert np.array_equal(g['wnorm'][k].to_numpy()[m], full['wnorm'][k].to_numpy()[m])
        assert np.array_equal(g['wnorm'][k].to_numpy()[~m], host['wnorm'][k].to_numpy()[~m])


def test_C_more_copies_than_claws():
    host = make_calyx(3, n_kcs={'g': 4, "a'b'": 2, 'ab': 2}, claws={'g': (2, 2), "a'b'": (1, 1), 'ab': (1, 1)})
    donor = make_calyx(4, n_kcs={'g': 4, "a'b'": 2, 'ab': 2}, claws={'g': (7, 7), "a'b'": (1, 1), 'ab': (1, 1)}, kc_id0=9000)
    T = B.build('C', host, donor, seed=0)
    n1, t1, _ = B.claw_counts(T['W'])
    assert (n1[t1 == 'g'] == 7).all()
    # 5 copies on 2 claws: each copied twice, one three times
    v = T['W'].to_numpy()
    g = np.asarray(host['W'].index.get_level_values(H.KC_TYPE_LEVEL)) == 'g'
    assert sorted(set(v[g].max(axis=1))) == [3.0, 4.0]


# ---------------------------------------------------------------------------
# S

def test_S_moves_to_the_donor_shares(pair, monkeypatch):
    from experiments import E97_glomerulus_claw_shares as E97
    _small1(monkeypatch)
    host, donor = pair
    T = B.build('S', host, donor, seed=0)
    r, sm, tab = B.shares_ratios(host['W'], donor['W'], small=1)
    assert not sm.any()
    # within each KC type the targets are the host's counts times r_g, rescaled to the type's total (E97)
    c0, c1 = E97.type_glomerulus_counts(host['W']), E97.type_glomerulus_counts(T['W'])
    assert c0.sum(axis=1).equals(c1.sum(axis=1))
    assert c1.equals(E97.type_targets(c0, r.to_dict()).reindex_like(c1))
    # each KC keeps its claws and the APL weights stay at their claws
    assert np.array_equal(B.claw_counts(T['W'])[0], B.claw_counts(host['W'])[0])
    size = T['meta']['transplant']['parts']['S']
    assert size['claws_moved'] > 0 and size['share_gap_to_donor_max'] < size['share_gap_to_donor_before_max']


def test_S_small_rule_keeps_everything(pair):
    """With every glomerulus below the small-glomerulus count, r_g = 1 and nothing moves."""
    host, donor = pair
    T = B.glomerulus_shares(host, donor, 0, None, small=10 ** 6)
    assert T['W'] is host['W'] and T['meta']['transplant']['parts']['S']['n_glomeruli_small'] == len(GLOMS)


def test_S_three_keeps_the_others_in_proportion(pair, monkeypatch):
    _small1(monkeypatch)
    sub = ['DM1', 'VA2']
    monkeypatch.setattr(B, 'three_glomeruli', lambda t: sub)
    host, donor = pair
    T = B.build('S-three', host, donor, seed=0)
    g0 = np.asarray(host['W'].columns.get_level_values('glomerulus'))[host['W'].to_numpy().argmax(axis=1)]
    g1 = np.asarray(T['W'].columns.get_level_values('glomerulus'))[T['W'].to_numpy().argmax(axis=1)]
    mv = g0 != g1
    assert mv.any()
    # every moved claw comes from or goes to one of the subset's glomeruli
    assert (np.isin(g0[mv], sub) | np.isin(g1[mv], sub)).all()
    assert set(T['meta']['transplant']['parts']['S']['glomeruli'].split()) == set(sub)


# ---------------------------------------------------------------------------
# A

def test_A_totals_rank_and_proportions(pair):
    host, donor = pair
    T = B.build('A', host, donor, seed=0)
    types = B._kc_types(host['W'])
    dtypes = B._kc_types(donor['W'])
    for k in ('wAPLKC', 'wKCAPL'):
        # the calyx's total is kept
        assert T['wnorm'][k].sum() == pytest.approx(host['wnorm'][k].sum(), rel=1e-12)
        h0, h1, d = _kc_sums(host, k), _kc_sums(T, k), _kc_sums(donor, k)
        matched = {}
        for t in ('g', "a'b'", 'ab'):
            m, md = (types == t).to_numpy(), (dtypes == t).to_numpy()
            ds = np.sort(d.to_numpy()[md])
            matched[t] = ds[((2 * np.arange(m.sum()) + 1) * md.sum()) // (2 * m.sum())]
        allm = np.concatenate(list(matched.values()))
        c = h0.sum() / allm.sum()
        for t in matched:
            m = (types == t).to_numpy()
            assert np.allclose(np.sort(h1.to_numpy()[m]), np.sort(matched[t] * c), rtol=1e-12, atol=1e-15)
        lit = T['meta']['transplant']['parts']['A'][('aplkc' if k == 'wAPLKC' else 'kcapl') + '_literal_total_ratio']
        assert lit == pytest.approx(allm.mean() / d.mean(), rel=1e-12)
        # within a KC the weights keep their proportions (KCs with weight before)
        f = (T['wnorm'][k] / host['wnorm'][k]).replace([np.inf], np.nan)
        kc = T['wnorm'][k].index.get_level_values(H.KC_ID_LEVEL)
        spread = f.groupby(kc).agg(lambda x: np.nanmax(x) - np.nanmin(x) if np.isfinite(x).any() else 0.0)
        assert (spread.fillna(0) < 1e-12 * f.abs().max()).all()
    # the KC with no APL->KC weight takes the donor's amount, evenly over its claws
    z = host['wnorm']['wAPLKC'].index.get_level_values(H.KC_ID_LEVEL) == 1001
    after = T['wnorm']['wAPLKC'].to_numpy()[z]
    assert np.allclose(after, after[0])


def test_A_and_B_splits(pair, monkeypatch):
    sub = ['DM1', 'VA2']
    monkeypatch.setattr(B, 'three_glomeruli', lambda t: sub)
    host, donor = pair
    fa, fb = B.build('A', host, donor, seed=2), B.build('B', host, donor, seed=2)
    a3, b3 = B.build('A-three', host, donor, seed=2), B.build('B-three', host, donor, seed=2)
    W = host['W']
    rg = np.asarray(W.columns.get_level_values('glomerulus'))[W.to_numpy().argmax(axis=1)]
    m = np.isin(rg, sub)
    for k in ('wAPLKC', 'wKCAPL'):
        assert np.array_equal(a3['wnorm'][k].to_numpy()[m], fa['wnorm'][k].to_numpy()[m])
        assert np.array_equal(a3['wnorm'][k].to_numpy()[~m], host['wnorm'][k].to_numpy()[~m])
    bm = np.isin(np.asarray(W.columns.get_level_values('glomerulus')), sub)
    for k in ('wAPLPN', 'wPNAPL'):
        assert np.array_equal(b3['wnorm'][k].to_numpy()[bm], fb['wnorm'][k].to_numpy()[bm])
        assert np.array_equal(b3['wnorm'][k].to_numpy()[~bm], host['wnorm'][k].to_numpy()[~bm])
    for arm, k, other in (('A-aplkc', 'wAPLKC', 'wKCAPL'), ('A-kcapl', 'wKCAPL', 'wAPLKC'), ('B-aplpn', 'wAPLPN', 'wPNAPL'),
                          ('B-pnapl', 'wPNAPL', 'wAPLPN')):
        T = B.build(arm, host, donor, seed=2)
        full = fa if arm[0] == 'A' else fb
        assert B._same(T['wnorm'][k], full['wnorm'][k]) and T['wnorm'][other] is host['wnorm'][other]
        assert _changed(T, host) == {k}


# ---------------------------------------------------------------------------
# B

def test_B_profile_total_and_missing_glomerulus(pair):
    host, donor = pair
    # the donor lacks VC3: its boutons keep their weights; the others take the donor's profile
    keep = np.asarray(donor['W'].columns.get_level_values('glomerulus')) != 'VC3'
    d2 = dict(donor, W=donor['W'].loc[:, keep],
              wnorm=dict(donor['wnorm'], wAPLPN=donor['wnorm']['wAPLPN'][keep], wPNAPL=donor['wnorm']['wPNAPL'][keep]))
    d2['W'] = d2['W'][d2['W'].to_numpy().sum(axis=1) > 0]
    d2['wnorm'] = dict(d2['wnorm'], wAPLKC=donor['wnorm']['wAPLKC'].reindex(d2['W'].index),
                       wKCAPL=donor['wnorm']['wKCAPL'].reindex(d2['W'].index))
    T = B.build('B', host, d2, seed=0)
    gh = np.asarray(host['W'].columns.get_level_values('glomerulus'))
    gd = np.asarray(d2['W'].columns.get_level_values('glomerulus'))
    for k in ('wAPLPN', 'wPNAPL'):
        w0, w1, wd = host['wnorm'][k].to_numpy(), T['wnorm'][k].to_numpy(), d2['wnorm'][k].to_numpy()
        assert w1.sum() == pytest.approx(w0.sum(), rel=1e-12)
        assert np.array_equal(w1[gh == 'VC3'], w0[gh == 'VC3'])
        rel1 = {g: w1[gh == g].mean() / w1.mean() for g in set(gd)}
        reld = {g: wd[gd == g].mean() / wd.mean() for g in set(gd)}
        # one common factor between the host's new profile and the donor's
        q = np.array([rel1[g] / reld[g] for g in sorted(reld)])
        assert np.allclose(q, q[0], rtol=1e-12)
        # within a glomerulus one factor
        for g in set(gd):
            f = w1[gh == g] / w0[gh == g]
            assert np.allclose(f, f[0], rtol=1e-13)
        assert 'VC3' in T['meta']['transplant']['parts']['B'][('aplpn' if k == 'wAPLPN' else 'pnapl') + '_kept_missing']


# ---------------------------------------------------------------------------
# K

def test_K_shares_and_row_sum(pair):
    host, donor = pair
    T = B.build('K', host, donor, seed=0)
    th, td = B._kc_types(host['W']), B._kc_types(donor['W'])
    b1, bd = B.block_sums(T['kckc'], th), B.block_sums(donor['kckc'], td)
    assert np.allclose((b1 / b1.to_numpy().sum()).to_numpy(), (bd / bd.to_numpy().sum()).reindex_like(b1).to_numpy(), rtol=1e-12)
    assert T['kckc'].to_numpy().sum() / len(T['kckc']) == pytest.approx(1.0, rel=1e-12)
    assert T['kckc'].index.equals(host['kckc'].index) and T['kckc'].columns.equals(host['kckc'].columns)
    # zeros stay zero; within a block one factor
    A0, A1 = host['kckc'].to_numpy(), T['kckc'].to_numpy()
    assert np.array_equal(A0 == 0, A1 == 0)
    tp = th.reindex(host['kckc'].index).to_numpy()
    for a in set(tp):
        for b in set(tp):
            blk = np.ix_(tp == a, tp == b)
            nz = A0[blk] > 0
            f = A1[blk][nz] / A0[blk][nz]
            assert np.allclose(f, f[0], rtol=1e-12)


def test_K_empty_block_stays_empty(pair):
    host, donor = pair
    th = B._kc_types(host['W'])
    M = host['kckc'].copy()
    g = (th.reindex(M.index) == 'g').to_numpy()
    ab = (th.reindex(M.columns) == 'ab').to_numpy()
    M.iloc[np.flatnonzero(g), np.flatnonzero(ab)] = 0.0
    M = M / (M.to_numpy().sum() / len(M))
    T = B.build('K', dict(host, kckc=M), donor, seed=0)
    assert (T['kckc'].to_numpy()[np.ix_(g, ab)] == 0).all()
    assert T['meta']['transplant']['parts']['K']['blocks_empty_in_host'] == 'g<ab'


# ---------------------------------------------------------------------------
# ALL, the arm strings, the rank rule, the size columns

def test_ALL_is_the_chain(pair, monkeypatch):
    _small1(monkeypatch)
    host, donor = pair
    T = B.build('ALL', host, donor, seed=3)
    U = host
    for a in B.ALL_ORDER:
        U = B.BUILDERS[a](U, donor, seed=3)
    assert B._same(T['W'], U['W']) and B._same(T['kckc'], U['kckc'])
    assert all(B._same(T['wnorm'][k], U['wnorm'][k]) for k in WN)
    assert list(T['meta']['transplant']['parts']) == list(B.ALL_ORDER)
    row = B.size_row(T)
    assert row['size_arm'] == 'ALL' and row['size_donor'] == 'test-donor' and 'size_K_rescale' in row
    assert 'size_C_off_share' in row and 'size_S_moved_share' in row and 'size_A_aplkc_literal_total_ratio' in row


def test_rank_match():
    rng = np.random.default_rng(0)
    h = rng.integers(0, 5, 17).astype(float)
    assert np.array_equal(B.rank_match(h, h, np.random.default_rng(1)), h)
    d = np.arange(34.0)
    out = B.rank_match(np.arange(17.0), d)
    assert np.array_equal(out, d[((2 * np.arange(17) + 1) * 34) // 34])
    assert np.array_equal(B.rank_match(np.arange(4.0)[::-1], np.array([10.0, 20.0])), np.array([20.0, 20.0, 10.0, 10.0]))


def test_parse_arm_and_donors():
    p = B.parse_arm('B-pnapl-dmr')
    assert p['base'] == 'B' and p['directions'] == ('wPNAPL',) and p['donor'] == 'malecns-right'
    assert B.parse_arm('C-g')['types'] == ('g',) and B.parse_arm('S-three')['subset'] == 'three'
    for bad in ('C-three', 'K-aplkc', 'A-aplpn', 'PA-dhb', 'X', 'C_g', 'S-g'):
        with pytest.raises(ValueError):
            B.parse_arm(bad)
    assert B.donor_of('hemibrain-argmax', 'C') == 'malecns-left' and B.donor_of('malecns-left', 'ALL') == 'hemibrain-argmax'
    assert B.donor_of('malecns-right', 'A') == 'hemibrain-argmax' and B.donor_of('malecns-left', 'S-dmr') == 'malecns-right'
    assert B.donor_of('malecns-left', 'PA') is None
    with pytest.raises(ValueError):
        B.donor_of('malecns-left', 'C-dml')
    # the runner's tags take the option strings
    from experiments import E133_calyx_transplants as E
    pt = E.Point(2, 'malecns-left', 2, 'B-pnapl-dmr')
    assert E.parse_tag(E.tag_of(pt)) == pt


# ---------------------------------------------------------------------------
# step 2′

def test_numpy_permutation_of_values_is_of_positions():
    """`within_type_permutation` replicates shuffle_wPNKC's `rng.permutation(values)` by `rng.permutation(n)`."""
    x = np.random.default_rng(5).integers(0, 50, 200)
    assert np.array_equal(np.random.default_rng(7).permutation(x), x[np.random.default_rng(7).permutation(len(x))])


def test_PA_is_the_plain_permutation_with_apl_moving(pair):
    from pyolfsysm.wpnkc_shuffles import shuffle_wPNKC
    host, _ = pair
    for s in (1, 2, 3):
        T = B.build('PA', host, None, seed=s)
        ref = shuffle_wPNKC(host['W'], 'claw_permutation', within_kc_type=True, seed=s)
        assert B._same(T['W'], ref)
        src = B.within_type_permutation(host['W'], s)
        v0 = host['W'].to_numpy()
        assert np.array_equal(T['W'].to_numpy(), v0[src])
        t = np.asarray(host['W'].index.get_level_values(H.KC_TYPE_LEVEL))
        assert np.array_equal(t[src], t)
        for k in ('wAPLKC', 'wKCAPL'):
            assert np.array_equal(T['wnorm'][k].to_numpy(), host['wnorm'][k].to_numpy()[src])
            assert T['wnorm'][k].sum() == pytest.approx(host['wnorm'][k].sum(), rel=1e-12)
        assert T['wnorm']['wAPLPN'] is host['wnorm']['wAPLPN'] and T['kckc'] is host['kckc'] and T['breadth'] is None
    with pytest.raises(ValueError):
        B.build('PA', dict(host, seed=2), None, seed=1)


def test_CI_on_a_made_up_calyx(pair):
    host, _ = pair
    T = B.coinput_shuffle(host, 1, match_kws=dict(n_steps=2000, t_start=5.0, t_end=0.05, type_layer_weight=4.0))
    assert np.array_equal(B.claw_counts(T['W'])[0], B.claw_counts(host['W'])[0])
    assert T['wnorm'] is host['wnorm'] and T['kckc'] is host['kckc']
    size = T['meta']['transplant']['parts']['CI']
    assert size['kc_n_claws_preserved'] and size['match_seed'] == 101 and np.isfinite(size['objective_end'])


def test_constants_shared_with_E59_and_E97():
    from experiments import E97_glomerulus_claw_shares as E97
    assert B.SMALL_CLAWS == E97.SMALL_CLAWS
    src = (Path(H.HERE) / 'experiments' / 'E59_coinput_shuffle.py').read_text()
    tree = ast.parse(src)
    kws = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'MATCH_KWS' for t in node.targets):
            kws = {kw.arg: ast.literal_eval(kw.value) for kw in node.value.keywords}
    assert kws == B.COINPUT_MATCH_KWS
    assert 'seed=100 + seed' in src and "within_kc_type=within_type, match_type_layer=within_type" in src
