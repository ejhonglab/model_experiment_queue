#!/usr/bin/env python3
"""Tests for E97 (the glomeruli's shares of the claws in other connectomes; experiments/E97_glomerulus_claw_shares.py)
on a made-up wiring and made-up tables: the builder's invariants of part 0 (each KC's claws and type, the per-type
totals, the targets, every moved claw on a bouton of its new glomerulus, every other row unchanged); ratios of 1 move
nothing; the seeds reproduce; the largest-remainder targets; the rule for small glomeruli; the mean proportions; the
name mapping; the permutation keeps its strata and the small glomeruli; A1's and A3's rules; the variants and the plan's
116 points."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
G = pytest.importorskip('experiments.E97_glomerulus_claw_shares')

QUICK = os.environ.get('QUICK') == '1'


# ---------------------------------------------------------------------------
# a made-up wiring: claw x bouton, one row per claw, (kc_id, claw_id, kc_type) rows, (glomerulus, bouton_id) columns

GLOMS = ('DA1', 'DC1', 'DL1', 'VM7d', 'VA4', 'V')


def made_up_wiring(seed: int = 0, n_kcs: int = 60) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    cols = [(g, 100 * i + b) for i, g in enumerate(GLOMS) for b in range(1 + (i % 3))]
    columns = pd.MultiIndex.from_tuples(cols, names=['glomerulus', 'bouton_id'])
    gl_w = np.array([3.0, 2.0, 1.5, 1.0, 0.6, 0.4])
    p_col = np.array([gl_w[GLOMS.index(g)] / (1 + (GLOMS.index(g) % 3)) for g, _ in cols])
    p_col = p_col / p_col.sum()
    rows, vals = [], []
    types = ['g'] * (n_kcs // 2) + ['ab'] * (n_kcs // 3) + ["a'b'"] * (n_kcs - n_kcs // 2 - n_kcs // 3)
    for k, t in enumerate(types):
        for c in range(int(rng.integers(3, 8))):
            rows.append((1000 + k, c, t))
            v = np.zeros(len(cols), dtype=np.int64)
            v[rng.choice(len(cols), p=p_col)] = 1
            vals.append(v)
    # one claw row without a bouton (left as it is)
    rows.append((1000, 99, 'g'))
    vals.append(np.zeros(len(cols), dtype=np.int64))
    index = pd.MultiIndex.from_tuples(rows, names=[H.KC_ID_LEVEL, 'claw_id', H.KC_TYPE_LEVEL])
    return pd.DataFrame(np.vstack(vals), index=index, columns=columns)


RATIOS = {'DA1': 0.7, 'DC1': 1.4, 'DL1': 0.9, 'VM7d': 1.8, 'VA4': 1.0, 'V': 0.6}


# ---------------------------------------------------------------------------
# the builder

def test_ratios_of_one_move_nothing():
    W = made_up_wiring()
    for r in ({}, {g: 1.0 for g in GLOMS}):
        W1, mv, tg = G.move_claws(W, r, [1, 0])
        assert W1.equals(W)
        assert len(mv) == 0
        assert (tg.base == tg.target).all()
        assert G.check_built(W, W1, mv, tg)['passes']


def test_builder_invariants():
    W = made_up_wiring()
    W1, mv, tg = G.move_claws(W, RATIOS, [1, 0])
    chk = G.check_built(W, W1, mv, tg)
    assert chk['passes'], chk
    assert len(mv) > 0
    # each KC's claws and type, the per-type totals
    t0, t1 = G.claw_table(W), G.claw_table(W1)
    assert (t0[H.KC_ID_LEVEL] == t1[H.KC_ID_LEVEL]).all() and (t0[H.KC_TYPE_LEVEL] == t1[H.KC_TYPE_LEVEL]).all()
    c0, c1 = G.type_glomerulus_counts(W), G.type_glomerulus_counts(W1)
    assert c0.sum(axis=1).equals(c1.sum(axis=1))
    # the glomeruli move the way of their ratios (VM7d gains, V and DA1 lose)
    assert c1['VM7d'].sum() > c0['VM7d'].sum() and c1['V'].sum() < c0['V'].sum() and c1['DA1'].sum() < c0['DA1'].sum()
    # a moved claw sits on a bouton of its new glomerulus; the row without a bouton is left as it is
    gl = np.asarray(W.columns.get_level_values('glomerulus'))
    assert (gl[W1.to_numpy()[mv.row.to_numpy()].argmax(axis=1)] == mv.to_glomerulus.to_numpy()).all()
    assert W1.loc[(1000, 99, 'g')].sum() == 0
    # each moved claw went to a bouton that has claws of its KC type in the base (drawn in proportion to them)
    for r in mv.itertuples():
        same_type = W.to_numpy()[(t0[H.KC_TYPE_LEVEL] == r.kc_type).to_numpy()]
        assert same_type[:, r.to_col].sum() > 0


def test_a_broken_wiring_fails_the_check():
    W = made_up_wiring()
    W1, mv, tg = G.move_claws(W, RATIOS, [1, 0])
    bad = W1.copy()
    i = int(mv.row.iloc[0])
    v = bad.to_numpy().copy()
    v[i] = 0
    v[i, int(mv.from_col.iloc[0])] = 1                  # put a moved claw back
    bad = pd.DataFrame(v, index=W1.index, columns=W1.columns)
    assert not G.check_built(W, bad, mv, tg)['passes']


def test_seeds_reproduce():
    W = made_up_wiring()
    a = G.move_claws(W, RATIOS, [1, 3])
    b = G.move_claws(W, RATIOS, [1, 3])
    c = G.move_claws(W, RATIOS, [2, 3])
    assert a[0].equals(b[0]) and a[1].equals(b[1])
    assert not a[0].equals(c[0])
    # the targets do not depend on the draw
    assert a[2].equals(c[2])


def test_type_targets_largest_remainder():
    counts = pd.DataFrame([[5, 3, 2]], index=['g'], columns=['A', 'B', 'C'])
    t = G.type_targets(counts, {'A': 2.0})
    # x = 10, 3, 2 (sum 15) scaled to 10: 6.67, 2, 1.33 -> 6, 2, 1 and the one left to A (the largest remainder)
    assert t.loc['g'].tolist() == [7, 2, 1]
    # a glomerulus without claws in a type gets none
    counts = pd.DataFrame([[5, 0, 5]], index=['g'], columns=['A', 'B', 'C'])
    t = G.type_targets(counts, {'B': 5.0, 'A': 1.3})
    assert t.loc['g', 'B'] == 0 and t.loc['g'].sum() == 10
    with pytest.raises(ValueError):
        G.type_targets(counts, {'A': 0.0})


# ---------------------------------------------------------------------------
# the derived tables' names, shares and ratios

def test_name_mapping():
    assert G.glomerulus_name('VC3l_adPN', 'hemibrain') == 'VC3'
    assert G.glomerulus_name('VC3m_lvPN', 'hemibrain') == 'VC3'
    assert G.glomerulus_name('VC5_lvPN', 'hemibrain') == 'VC5'
    assert G.glomerulus_name('VC3_adPN', 'flywire') == 'VC3'
    assert G.glomerulus_name('VC5_lvPN', 'flywire') == 'VC3'
    assert G.glomerulus_name('VM6_adPN', 'malecns') == 'VC5'
    assert G.glomerulus_name('VM7d_adPN', 'banc') == 'VM7d'
    for t in ('VP1d+VP4_l2PN1', 'VP3+_vPN', 'M_vPNml50', None, float('nan'), ''):
        assert G.glomerulus_name(t, 'flywire') is None
    assert G.kc_class('KCg-m') == 'g' and G.kc_class('KCab-c') == 'ab' and G.kc_class("KCa'b'-ap1") == "a'b'"
    assert G.kc_class('KCapbp-m') == "a'b'" and G.kc_class('KC') == 'unknown' and G.kc_class("a'b'") == "a'b'"


def test_ratios_small_rule_and_mean_proportions():
    gl = ['A', 'B', 'C', 'D']
    counts = pd.DataFrame({'hemibrain': [100, 200, 20, 80], 'flywire_L': [200, 200, 200, 25],
                           'flywire_R': [200, 200, 30, 100], 'malecns_L': [50, 200, 100, 100],
                           'malecns_R': [50, 200, 100, 100]}, index=gl)
    r, sm = G.ratios(counts, gl)
    # C has 20 in the hemibrain: 1 everywhere; D has 25 in flywire_L: 1 there only
    assert (r.loc['C'] == 1).all() and sm.loc['C'].all()
    assert r.loc['D', 'flywire_L'] == 1 and sm.loc['D', 'flywire_L'] and not sm.loc['D', 'flywire_R']
    s = counts / counts.sum()
    assert r.loc['A', 'flywire_R'] == pytest.approx(s.loc['A', 'flywire_R'] / s.loc['A', 'hemibrain'])
    m = G.mean_proportions(r)
    expect = 2 ** ((0 + np.log2(r.loc['A', ['flywire_L', 'flywire_R']]).mean()
                    + np.log2(r.loc['A', ['malecns_L', 'malecns_R']]).mean()) / 3)
    assert m['A'] == pytest.approx(expect)
    assert m['C'] == 1.0
    sh = G.shares(pd.Series({'A': 1, 'B': 3, 'X': 10}), ['A', 'B'])
    assert sh.tolist() == [0.25, 0.75]


def test_permutation_keeps_its_strata_and_the_small_ones():
    gl = [f'g{i}' for i in range(12)]
    totals = pd.Series(np.arange(12) * 10.0 + 5, index=gl)
    st = G.strata(totals, with_input=gl[:8])
    assert st.value_counts().to_dict() == {'input_1': 3, 'input_2': 3, 'input_3': 2, 'noinput_1': 2, 'noinput_2': 1,
                                           'noinput_3': 1}
    assert st['g0'] == 'input_1' and st['g7'] == 'input_3' and st['g11'] == 'noinput_3'
    r = pd.Series(np.linspace(0.5, 2.0, 12), index=gl)
    fixed = ['g0']
    p1 = G.permuted_ratios(r, st, 1, fixed)
    assert p1['g0'] == r['g0']
    for s in st.unique():
        gs = st.index[st == s]
        assert sorted(p1[gs]) == sorted(r[gs])
    assert p1.equals(G.permuted_ratios(r, st, 1, fixed))
    assert not p1.equals(G.permuted_ratios(r, st, 2, fixed))
    assert G.scaled_ratios(r, -1.0).mul(r).round(12).eq(1).all()
    assert G.scaled_ratios(r, 0.5).pow(2).sub(r).abs().max() < 1e-12


# ---------------------------------------------------------------------------
# the rules

def test_a1_rule():
    gl = [f'g{i}' for i in range(20)]
    model = pd.Series(np.linspace(50, 500, 20), index=gl)
    ok = G.a1_rule(model, model * 1.05)
    assert ok['sd'] == pytest.approx(0.0, abs=1e-12) and not ok['stops']
    rng = np.random.default_rng(0)
    off = model * 2 ** rng.normal(0, 0.4, 20)
    bad = G.a1_rule(model, off)
    assert bad['sd'] > G.A1_RULE_SD and bad['stops']
    assert bad['n_glomeruli'] == int((model >= G.A1_RULE_MIN_CLAWS).sum())
    gone = G.a1_rule(model, model.where(model.index != 'g19', 0))
    assert gone['stops'] and gone['n_infinite'] == 1


def test_a3_rule():
    null = np.linspace(0.3, 0.7, 101)
    assert G.a3_rule(0.9, null) == 'with'
    assert G.a3_rule(0.1, null) == 'against'
    assert G.a3_rule(0.5, null) == 'neither'
    assert G.sign_agreement([1, -1, 1, 0], [2, -3, -1, 1]) == 0.5


def test_counting_noise():
    n = np.full(40, 400.0)
    cn = G.counting_noise_sd(n, n, n_draws=200)
    assert cn['analytic'] == pytest.approx(np.sqrt(2 / 400) / np.log(2))
    assert cn['simulated'] == pytest.approx(cn['analytic'], rel=0.15)


# ---------------------------------------------------------------------------
# the variants and the plan

def test_variants_and_the_plan():
    p = G.plan(G.ARMS, G.ALL_WIRINGS)
    assert len(p) == 116
    sizes = {}
    for v, w in p:
        a = G.parse_variant(v)['arm']
        sizes[a] = sizes.get(a, 0) + 1
    assert sizes == {'B1': 32, 'B2': 8, 'B3': 56, 'B4': 4, 'B5': 8, 'B6': 8}
    assert len(G.plan(('C2',), G.ALL_WIRINGS)) == 4 and {w for _, w in G.plan(('C2',))} == set(G.T_WIRINGS)
    assert len(G.plan(G.ARMS, G.S_WIRINGS)) == 4 * 1 + 4 * 4 + 4 * 7 + 4 + 8 + 4
    assert G.parse_variant('B3_perm_malecns_L_s1') == dict(variant='B3_perm_malecns_L_s1', arm='B3',
                                                           proportions='malecns_L', control='perm', seed=1, scale=1.0,
                                                           draw=1, input='lab')
    assert G.parse_variant('B4_mean_d2')['draw'] == 2 and G.parse_variant('C2_mean_d2')['draw'] == 2
    assert G.parse_variant('B5_mean_x0.5')['scale'] == 0.5 and G.parse_variant('B3_invmean')['scale'] == -1.0
    assert G.parse_variant('B6_mean_F')['input'] == 'F'
    for bad in ('B1_banc_L', 'B2_x', 'B7_mean', 'B5_mean_x3'):
        with pytest.raises((ValueError, AttributeError)):
            G.parse_variant(bad)
    assert p[0] == ('B2_mean', 'real')
    assert G.draw_seed('B4_mean_d2', 'P3') == [G.DRAW_SEEDS[1], G.ALL_WIRINGS.index('P3')]


def test_variant_ratios_on_a_made_up_table():
    gl = list(GLOMS)
    rt = pd.DataFrame(dict(glomerulus=gl, stratum=['input_1', 'input_1', 'input_2', 'input_2', 'noinput_1', 'noinput_1'],
                           r_mean=[0.7, 1.4, 0.9, 1.8, 1.0, 0.6], fixed_mean=[False, False, False, False, True, False]))
    for m in G.OTHER_MBS:
        rt[f'r_{m}'] = [1.1, 0.9, 1.2, 0.8, 1.0, 1.3]
        rt[f'fixed_{m}'] = [False] * 6
    rt = rt.set_index('glomerulus')
    assert G.variant_ratios('B2_mean', rt).tolist() == rt.r_mean.tolist()
    inv = G.variant_ratios('B3_invmean', rt)
    assert np.allclose(inv * rt.r_mean, 1.0)
    p = G.variant_ratios('B3_permmean_s1', rt)
    assert p['VA4'] == 1.0
    assert sorted(p[['DA1', 'DC1']]) == [0.7, 1.4] and sorted(p[['DL1', 'VM7d']]) == [0.9, 1.8]
    assert np.allclose(G.variant_ratios('B5_mean_x2', rt), rt.r_mean ** 2)
    assert G.variant_ratios('B1_flywire_R', rt).tolist() == rt['r_flywire_R'].tolist()
    W = made_up_wiring()
    W1, mv, tg = G.move_claws(W, G.variant_ratios('B2_mean', rt).to_dict(), [1, 0])
    assert G.check_built(W, W1, mv, tg)['passes']


# ---------------------------------------------------------------------------
# part C's verdict on made-up rows

def _rows(b2=0.0, b2_s=None, b2_t=None, perm=(0.0, 0.0), inv=0.0, r2=0.0, r2p=0.0, six=None, flagged=0, raw=0.0,
          mean=0.178, mbs=(0.0, 0.0, 0.0, 0.0), b6=0.0):
    """readings rows: each variant's pattern part as the reference's (0.0500 on S, 0.0486 on T, 0.0482 on the eight)
    minus the given gain."""
    ref = {'S': G.REF_S, 'T': G.REF_T, 'all8': G.REF_8}
    gains = {G.REF: (0, 0, 0), 'B2_mean': (b2_s if b2_s is not None else b2, b2_t if b2_t is not None else b2, b2),
             'B3_permmean_s1': (perm[0],) * 3, 'B3_permmean_s2': (perm[1],) * 3, 'B3_invmean': (inv,) * 3,
             G.R2: (r2,) * 3, G.R2_PERM: (r2p,) * 3, f'B6_mean_{G.B6_INPUT}': (b2 + b6,) * 3}
    for m, g in zip(G.OTHER_MBS, mbs):
        gains[f'B1_{m}'] = (g,) * 3
    rows = []
    for v, gs in gains.items():
        for unit, g in zip(('S', 'T', 'all8'), gs):
            sx = (six if six is not None else g) if v == 'B2_mean' else g
            rows.append(dict(variant=v, unit=unit, readout='cal', z8_pattern=ref[unit] - g, z6_pattern=0.05 - sx,
                             mean_effect=mean if v == 'B2_mean' else 0.177, n_flagged=flagged if v == 'B2_mean' else 0))
            rows.append(dict(variant=v, unit=unit, readout='raw', z8_loss=0.073 + (raw if v == 'B2_mean' else 0.0),
                             n_flagged=0))
    return pd.DataFrame(rows)


def _boot(hi):
    return pd.DataFrame([dict(variant='B2_mean', unit='all8', ref_variant=G.REF, diff_hi=hi)])


def test_verdicts(tmp_path, monkeypatch):
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    v = G.verdict97(_rows(b2=0.001, r2=-0.0005), _boot(0.001))
    assert v['outcome'] == 'no move' and v['final'] == 'no move' and not v['needs_draw2']
    v = G.verdict97(_rows(b2=0.004, perm=(0.001, 0.0005), inv=-0.001, r2=0.002, r2p=0.0005), _boot(-0.001))
    assert v['outcome'] == 'moves' and v['needs_draw2'] and v['final'] == 'pending: C2'
    assert v['conditions']['interval_below_zero'] and v['conditions']['controls_perm']
    # the same gain without the interval below zero: not 'moves', but 'a part' needs a gain under 0.003
    v = G.verdict97(_rows(b2=0.004, perm=(0.001, 0.0005), inv=-0.001, r2=0.002), _boot(0.0002))
    assert not v['moves'] and v['outcome'] == "none of the entry's"
    v = G.verdict97(_rows(b2=0.002, perm=(0.0005, 0.0), inv=0.0), _boot(-0.0001))
    assert v['outcome'] == 'a part'
    # a control that gains as much fails both
    v = G.verdict97(_rows(b2=0.004, perm=(0.003, 0.0), inv=0.0, r2=0.002), _boot(-0.001))
    assert not v['moves'] and not v['a_part']
    v = G.verdict97(_rows(b2=-0.002), _boot(0.003))
    assert v['outcome'] == 'worse'
    v = G.verdict97(_rows(b2=0.0005, r2=0.002, r2p=0.0015), _boot(0.001))
    assert v['outcome'] == 'averaging' and not v['no_move']
    # S or T alone below 0.0015: not 'moves'
    v = G.verdict97(_rows(b2=0.004, b2_t=0.001, perm=(0.0, 0.0), inv=0.0, r2=0.002), _boot(-0.001))
    assert not v['moves']
    # the bands, the raw-count loss and the mean APL effect
    for kw in (dict(flagged=1), dict(raw=0.006), dict(mean=0.165), dict(six=0.001)):
        v = G.verdict97(_rows(b2=0.004, perm=(0.0, 0.0), inv=0.0, r2=0.002, **kw), _boot(-0.001))
        assert not v['moves'], kw
    # (f) and (h)
    v = G.verdict97(_rows(mbs=(0.004, -0.001, 0.0, 0.0), b6=0.002), _boot(0.0))
    assert v['disagree'] and v['B6_differs']
    assert (tmp_path / 'verdict.json').exists()


def test_an_incomplete_part_b_gives_no_verdict(tmp_path, monkeypatch):
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    rd = _rows(b2=0.004)
    rd = rd[~((rd.variant == G.R2) & (rd.unit == 'all8'))]
    v = G.verdict97(rd, _boot(-0.001))
    assert v['outcome'] == 'incomplete' and 'R2' in v['missing'] and not v['needs_draw2']
