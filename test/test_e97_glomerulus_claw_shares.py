#!/usr/bin/env python3
"""Tests for E97 (the glomeruli's shares of the claws in other connectomes; experiments/E97_glomerulus_claw_shares.py)
on a made-up wiring and made-up tables in the form of Pratyush's two per-synapse tables (the hemibrain's and the male
CNS's): the reader's identity with the loader's steps, the three calyx rules and the confidence reading; the names; the
rule for small glomeruli with a label missing on one side; the two-brain mean; no held table named or opened; the
builder's invariants of part 0 (each KC's claws and type, the per-type totals, the targets, every moved claw on a bouton
of its new glomerulus, every other row unchanged); ratios of 1 move nothing; the seeds reproduce; the largest-remainder
targets; the permutation keeps its strata and the small glomeruli; A1's and A3's rules; the variants and the plan's 84
points; part C's verdict and its "at a price" label; part (R)'s classes, its second draw and its
readings; B6's input check within a last-bit tolerance."""

import json
import os
from typing import Dict, Tuple

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
# shares and ratios

def test_kc_classes():
    assert G.kc_class('KCg-m') == 'g' and G.kc_class('KCab-c') == 'ab' and G.kc_class("KCa'b'-ap1") == "a'b'"
    assert G.kc_class('KCg') == 'g' and G.kc_class('KC') == 'unknown' and G.kc_class("a'b'") == "a'b'"
    assert G.kc_class_of(None) is None and G.kc_class_of('KC part due to gap') is None


def test_ratios_small_rule_and_the_two_brain_mean():
    gl = ['A', 'B', 'C', 'D', 'E']
    # E has no PN of its label on the right (missing there): it counts 0, so it is small there only
    counts = pd.DataFrame({'hemibrain': [100, 200, 20, 80, 60], 'malecns_L': [200, 200, 200, 25, 50],
                           'malecns_R': [50, 200, 100, 100, np.nan]}, index=gl)
    r, sm = G.ratios(counts, gl)
    assert (r.loc['C'] == 1).all() and sm.loc['C'].all()
    assert r.loc['D', 'malecns_L'] == 1 and sm.loc['D', 'malecns_L'] and not sm.loc['D', 'malecns_R']
    assert r.loc['E', 'malecns_R'] == 1 and sm.loc['E', 'malecns_R'] and not sm.loc['E', 'malecns_L']
    s = counts.fillna(0) / counts.fillna(0).sum()
    assert r.loc['A', 'malecns_R'] == pytest.approx(s.loc['A', 'malecns_R'] / s.loc['A', 'hemibrain'])
    m = G.mean_proportions(r)
    # the hemibrain's 1 and the male CNS's two sides averaged on log2 first: 2 ** (mean log2 / 2)
    expect = 2 ** ((0 + np.log2(r.loc['A', ['malecns_L', 'malecns_R']]).mean()) / 2)
    assert m['A'] == pytest.approx(expect)
    assert m['C'] == 1.0
    assert m['E'] == pytest.approx(2 ** (np.log2(r.loc['E', 'malecns_L']) / 4))
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
    counts = pd.DataFrame({'malecns_L': model * 1.02, 'malecns_R': model * 1.01})
    pns = pd.DataFrame({'malecns_L': 2.0, 'malecns_R': 1.0}, index=gl)
    ok = G.a1_rule(model, counts, pns)
    assert not ok['stops'] and ok['mbs']['malecns_L']['ratio_to_model'] == pytest.approx(1.02)
    assert ok['n_big'] == int((model >= G.A1_RULE_MIN_CLAWS).sum())
    for k, f in (('malecns_L', 0.79), ('malecns_R', 1.26)):
        bad = G.a1_rule(model, counts.assign(**{k: model * f}), pns)
        assert bad['stops'] and not bad['mbs'][k]['inside_band']
    gone = G.a1_rule(model, counts, pns.assign(malecns_R=pns.malecns_R.where(pns.index != 'g19', 0)))
    assert gone['stops'] and gone['mbs']['malecns_R']['big_without_pn'] == ['g19']
    # a small glomerulus without a PN does not stop it
    assert not G.a1_rule(model, counts, pns.assign(malecns_R=pns.malecns_R.where(pns.index != 'g0', 0)))['stops']


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
    assert len(p) == 84
    sizes = {}
    for v, w in p:
        a = G.parse_variant(v)['arm']
        sizes[a] = sizes.get(a, 0) + 1
    assert sizes == {'B1': 16, 'B2': 8, 'B3': 40, 'B4': 4, 'B5': 8, 'B6': 8}
    assert len(G.plan(('C2',), G.ALL_WIRINGS)) == 4 and {w for _, w in G.plan(('C2',))} == set(G.T_WIRINGS)
    assert len(G.plan(G.ARMS, G.S_WIRINGS)) == 4 * 1 + 4 * 2 + 4 * 5 + 4 + 8 + 4
    assert G.parse_variant('B3_perm_mcnsL_s1') == dict(variant='B3_perm_mcnsL_s1', arm='B3', proportions='malecns_L',
                                                       control='perm', seed=1, scale=1.0, draw=1, input='lab')
    assert G.parse_variant('B1_mcnsR')['proportions'] == 'malecns_R'
    assert G.parse_variant('B4_mean2_d2')['draw'] == 2 and G.parse_variant('C2_mean2_d2')['draw'] == 2
    assert G.parse_variant('B5_mean2_x0.5')['scale'] == 0.5 and G.parse_variant('B3_invmean2')['scale'] == -1.0
    assert G.parse_variant('B6_mean2_F')['input'] == 'F'
    for bad in ('B1_malecns_L', 'B2_mean', 'B2_x', 'B7_mean2', 'B5_mean2_x3', 'B3_invmean'):
        with pytest.raises((ValueError, AttributeError, KeyError)):
            G.parse_variant(bad)
    assert p[0] == ('B2_mean2', 'real')
    assert G.draw_seed('B4_mean2_d2', 'P3') == [G.DRAW_SEEDS[1], G.ALL_WIRINGS.index('P3')]


def test_no_new_tag_is_a_withheld_one():
    """No point of the re-plan takes the name of one of the withheld points of 2026-11-15 (their arms on the eight
    wirings, the five mushroom bodies of that design and its mean), and the run table is a new one."""
    old_mbs = ['flywire_L', 'flywire_R', 'malecns_L', 'malecns_R']
    old = (['B2_mean', 'B3_permmean_s1', 'B3_permmean_s2', 'B3_invmean', 'B4_mean_d2', 'B5_mean_x0.5', 'B5_mean_x2',
            'B6_mean_F', 'C2_mean_d2'] + [f'B1_{m}' for m in old_mbs] + [f'B3_perm_{m}_s1' for m in old_mbs])
    old_tags = {f'{v}_{w}' for v in old for w in G.ALL_WIRINGS}
    new_tags = {G.run_tag(v, w) for v, w in G.plan(G.ARMS + ('C2',), G.ALL_WIRINGS)}
    assert 'B2_mean_real' in old_tags and not (new_tags & old_tags)
    assert G.RUNS_TABLE != 'runs_B.csv' and G.RUNS_SUB != 'runs'


def test_variant_ratios_on_a_made_up_table():
    gl = list(GLOMS)
    rt = pd.DataFrame(dict(glomerulus=gl, stratum=['input_1', 'input_1', 'input_2', 'input_2', 'noinput_1', 'noinput_1'],
                           r_mean=[0.7, 1.4, 0.9, 1.8, 1.0, 0.6], fixed_mean=[False, False, False, False, True, False]))
    for m in G.OTHER_MBS:
        rt[f'r_{m}'] = [1.1, 0.9, 1.2, 0.8, 1.0, 1.3]
        rt[f'fixed_{m}'] = [False] * 6
    rt = rt.set_index('glomerulus')
    assert G.variant_ratios('B2_mean2', rt).tolist() == rt.r_mean.tolist()
    inv = G.variant_ratios('B3_invmean2', rt)
    assert np.allclose(inv * rt.r_mean, 1.0)
    p = G.variant_ratios('B3_permmean2_s1', rt)
    assert p['VA4'] == 1.0
    assert sorted(p[['DA1', 'DC1']]) == [0.7, 1.4] and sorted(p[['DL1', 'VM7d']]) == [0.9, 1.8]
    assert np.allclose(G.variant_ratios('B5_mean2_x2', rt), rt.r_mean ** 2)
    assert G.variant_ratios('B1_mcnsR', rt).tolist() == rt['r_malecns_R'].tolist()
    W = made_up_wiring()
    W1, mv, tg = G.move_claws(W, G.variant_ratios('B2_mean2', rt).to_dict(), [1, 0])
    assert G.check_built(W, W1, mv, tg)['passes']


# ---------------------------------------------------------------------------
# part C's verdict on made-up rows

def _rows(b2=0.0, b2_s=None, b2_t=None, perm=(0.0, 0.0), inv=0.0, r2=0.0, r2p=0.0, six=None, flagged=0, raw=0.0,
          mean=0.178, mbs=(0.0, 0.0), b6=0.0):
    """readings rows: each variant's pattern part as the reference's (0.0500 on S, 0.0486 on T, 0.0482 on the eight)
    minus the given gain."""
    ref = {'S': G.REF_S, 'T': G.REF_T, 'all8': G.REF_8}
    gains = {G.REF: (0, 0, 0), 'B2_mean2': (b2_s if b2_s is not None else b2, b2_t if b2_t is not None else b2, b2),
             'B3_permmean2_s1': (perm[0],) * 3, 'B3_permmean2_s2': (perm[1],) * 3, 'B3_invmean2': (inv,) * 3,
             G.R2: (r2,) * 3, G.R2_PERM: (r2p,) * 3, f'B6_mean2_{G.B6_INPUT}': (b2 + b6,) * 3}
    for m, g in zip(G.OTHER_MBS, mbs):
        gains[f'B1_{G.MB_SHORT[m]}'] = (g,) * 3
    rows = []
    for v, gs in gains.items():
        for unit, g in zip(('S', 'T', 'all8'), gs):
            sx = (six if six is not None else g) if v == 'B2_mean2' else g
            rows.append(dict(variant=v, unit=unit, readout='cal', z8_pattern=ref[unit] - g, z6_pattern=0.05 - sx,
                             mean_effect=mean if v == 'B2_mean2' else 0.177, n_flagged=flagged if v == 'B2_mean2' else 0))
            rows.append(dict(variant=v, unit=unit, readout='raw', z8_loss=0.073 + (raw if v == 'B2_mean2' else 0.0),
                             n_flagged=0))
    return pd.DataFrame(rows)


def _boot(hi):
    return pd.DataFrame([dict(variant='B2_mean2', unit='all8', ref_variant=G.REF, diff_hi=hi)])


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
    v = G.verdict97(_rows(mbs=(0.004, -0.001), b6=0.002), _boot(0.0))
    assert v['disagree'] and v['B6_differs']
    assert (tmp_path / 'verdict.json').exists()


def test_part_c_at_a_price(tmp_path, monkeypatch):
    """Part C's verdict stands; "at a price" when R1 lowers a ρ on odor-odor distances by more than 0.005 or response
    rate ρ by more than 0.05 (added 2026-11-18 with part R), and with the pending columns applied."""
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    rows = _rows(b2=0.004, perm=(0.001, 0.0005), inv=-0.001, r2=0.002, r2p=0.0005)
    extra = {(G.REF, 'raw'): dict(rho_on_all=0.690, rho_off_all=0.756, run_rate_rho_on=0.47, run_rate_rho_off=0.40,
                                  block=0.05, vcf_on=0.3, vcf_off=0.2, vcf_diff=0.1, struct_wt_on=0.4),
             ('B2_mean2', 'raw'): dict(rho_on_all=0.680, rho_off_all=0.754, run_rate_rho_on=0.46, run_rate_rho_off=0.40,
                                       block=0.05, vcf_on=0.3, vcf_off=0.2, vcf_diff=0.1, struct_wt_on=0.3)}

    def with_extra(ex):
        rd = rows.copy()
        for (v, ro), cols in ex.items():
            for c, x in cols.items():
                rd.loc[(rd.variant == v) & (rd.unit == 'all8') & (rd.readout == ro), c] = x
        return rd
    v = G.verdict97(with_extra(extra), _boot(-0.001))
    assert v['outcome'] == 'moves' and v['at_a_price'] and v['label'] == 'moves at a price'
    assert list(v['price']) == ['ρ on odor-odor distances, APL on'] and v['price'][
        'ρ on odor-odor distances, APL on'] == pytest.approx(-0.010)
    assert v['price_pending'] == {'structure vs wild-type somata, APL on': pytest.approx(-0.1)}
    ex2 = {k: dict(c) for k, c in extra.items()}
    ex2[('B2_mean2', 'raw')].update(rho_on_all=0.686, struct_wt_on=0.4)
    v = G.verdict97(with_extra(ex2), _boot(-0.001))
    assert not v['at_a_price'] and v['label'] == 'moves' and v['label_with_pending'] == 'moves'
    ex2[('B2_mean2', 'raw')].update(run_rate_rho_off=0.34)
    v = G.verdict97(with_extra(ex2), _boot(-0.001))
    assert v['at_a_price'] and list(v['price']) == ['response rate ρ, APL off']
    # no price on a verdict that is not 'a part' or 'moves'
    v = G.verdict97(_rows(b2=0.001, r2=-0.0005), _boot(0.001))
    assert v['label'] == 'no move' and not v['at_a_price']


# ---------------------------------------------------------------------------
# part (R), the second verdict on ρ on odor-odor distances, on made-up rows

RHO_BASE = {'on': 0.690, 'off': 0.756}


def _rho_rows(on=0.0, off=0.0, per_w=None, halves=None, perm=((0.0, 0.0), (0.0, 0.0)), inv=(0.0, 0.0),
              r2=(0.0, 0.0), r2p=(0.0, 0.0), pattern=0.0, half_pattern=None, raw=0.0, mean=0.0, flagged=0,
              rate=(0.0, 0.0), block=0.0, vcf=(0.0, 0.0), struct=0.0, d2=None, drop=()):
    """readings rows for part R: every reading of the reference at a base value, B2's (and the controls') the base plus
    the given change; per_w {wiring: (on, off)} B2's change on each wiring (default the population's), halves (on,
    off) on each of halves A to D (default the population's); d2 the second draw's changes (on, off), or None."""
    ch = {'on': on, 'off': off}
    hv = {'on': on, 'off': off} if halves is None else dict(zip(('on', 'off'), halves))
    rows = []

    def raw_row(v, unit, d_on, d_off, d_h=None, extra=None):
        d_h = d_h or {'on': d_on, 'off': d_off}
        r = dict(variant=v, unit=unit, readout='raw', n_flagged=0, z8_loss=0.073, run_rate_rho_on=0.47,
                 run_rate_rho_off=0.40, block=0.05, vcf_on=0.30, vcf_off=0.20, vcf_diff=0.10, struct_wt_on=0.40)
        for c, dc in (('on', d_on), ('off', d_off)):
            r[f'rho_{c}_all'] = RHO_BASE[c] + dc
            for h in G.HALVES:
                r[f'rho_{c}_{h}'] = RHO_BASE[c] - 0.02 + d_h[c]
        r.update(extra or {})
        rows.append(r)

    def cal_row(v, unit, extra=None):
        r = dict(variant=v, unit=unit, readout='cal', z8_pattern=0.0482, mean_effect=0.177,
                 **{f'half{h}_pattern': 0.055 for h in G.HALVES})
        r.update(extra or {})
        rows.append(r)
    for w in G.ALL_WIRINGS:
        raw_row(G.REF, w, 0.0, 0.0)
        pw = (per_w or {}).get(w, (on, off))
        raw_row('B2_mean2', w, *pw)
    raw_row(G.REF, 'all8', 0.0, 0.0)
    cal_row(G.REF, 'all8')
    hp = pattern if half_pattern is None else half_pattern
    raw_row('B2_mean2', 'all8', on, off, hv,
            dict(n_flagged=flagged, z8_loss=0.073 + raw, run_rate_rho_on=0.47 + rate[0], run_rate_rho_off=0.40 + rate[1],
                 block=0.05 + block, vcf_on=0.30 + vcf[0], vcf_off=0.20 + vcf[1], vcf_diff=0.10 + vcf[0] - vcf[1],
                 struct_wt_on=0.40 + struct))
    cal_row('B2_mean2', 'all8', dict(z8_pattern=0.0482 + pattern, mean_effect=0.177 + mean,
                                     **{f'half{h}_pattern': 0.055 + hp for h in G.HALVES}))
    for s, p in zip(G.PERM_SEEDS, perm):
        raw_row(f'B3_permmean2_s{s}', 'all8', *p)
    raw_row('B3_invmean2', 'all8', *inv)
    raw_row(G.R2, 'all8', *r2)
    raw_row(G.R2_PERM, 'all8', *r2p)
    if d2 is not None:
        for w in G.ALL_WIRINGS:
            raw_row(f'B4_mean2_d2' if w in G.S_WIRINGS else 'C2_mean2_d2', w, *d2)
        raw_row(G.D2, 'all8', *d2)
        cal_row(G.D2, 'all8')
    rd = pd.DataFrame(rows)
    for v, u in drop:
        rd = rd[~((rd.variant == v) & (rd.unit == u))]
    return rd


def test_rho_classes():
    c = G.rho_class(_rho_rows(on=0.002, off=-0.001))
    assert c['cls'] == 'no change' and c['no_change'] and not c['raises']
    c = G.rho_class(_rho_rows(on=0.015, off=0.0))
    assert c['cls'] == 'raises' and c['gains'] == {'on': True, 'off': False} and c['n_wirings_up_on'] == 8
    # the other ρ lower by more than 0.005: no gain
    c = G.rho_class(_rho_rows(on=0.015, off=-0.006))
    assert not c['gains']['on'] and c['cls'] == 'below the gate'
    # only 5 of the 8 wirings gain 0.010
    pw = {w: ((0.015, 0.0) if i < 5 else (0.005, 0.0)) for i, w in enumerate(G.ALL_WIRINGS)}
    c = G.rho_class(_rho_rows(on=0.015, per_w=pw))
    assert c['n_wirings_up_on'] == 5 and c['cls'] == 'below the gate'
    # one half below the gain
    rd = _rho_rows(on=0.015)
    rd.loc[(rd.variant == 'B2_mean2') & (rd.unit == 'all8'), 'rho_on_C'] -= 0.01
    assert G.rho_class(rd)['cls'] == 'below the gate'
    # the controls: a permuted seed gaining half as much, or the inverted gaining 0.005
    assert G.rho_class(_rho_rows(on=0.015, perm=((0.008, 0.0), (0.0, 0.0))))['cls'] == 'below the gate'
    assert G.rho_class(_rho_rows(on=0.015, perm=((0.007, 0.0), (0.0, 0.0))))['cls'] == 'raises'
    assert G.rho_class(_rho_rows(on=0.015, inv=(0.0055, 0.0)))['cls'] == 'below the gate'
    # R2 below the eight as they are by more than 0.005 on either ρ
    assert G.rho_class(_rho_rows(on=0.015, r2=(0.0, -0.006)))['cls'] == 'below the gate'
    # the not-worse margins: a price, each named
    for kw, name in ((dict(pattern=0.0031), 'pattern part (all)'), (dict(half_pattern=0.0031), 'pattern part (A)'),
                     (dict(raw=0.0051), 'raw-count loss'), (dict(mean=-0.011), 'mean APL effect'),
                     (dict(flagged=1), 'realism (runs outside the bands)'), (dict(rate=(-0.051, 0.0)), 'response rate ρ, APL on')):
        c = G.rho_class(_rho_rows(on=0.015, **kw))
        assert c['cls'] == 'raises at a price' and name in c['fails'], kw
    c = G.rho_class(_rho_rows(on=0.015, pattern=0.0029, raw=0.0049, mean=-0.009, rate=(-0.049, -0.049)))
    assert c['cls'] == 'raises' and not c['fails']
    # the pending columns enter no class, but give the class as it would read with them applied
    for kw, name in ((dict(block=-0.011), 'the block'), (dict(vcf=(-0.051, 0.0)), 'natural-source alignment, APL on'),
                     (dict(vcf=(0.0, 0.06)), 'natural-source alignment, APL on minus off'),
                     (dict(struct=-0.051), 'structure vs wild-type somata, APL on')):
        c = G.rho_class(_rho_rows(on=0.015, **kw))
        assert c['cls'] == 'raises' and c['cls_with_pending'] == 'raises at a price' and name in c['pending_fails'], kw
    # lowers, no change, averaging
    c = G.rho_class(_rho_rows(off=-0.012))
    assert c['cls'] == 'lowers' and c['lowers_off'] and not c['lowers_on']
    c = G.rho_class(_rho_rows(on=0.003, r2=(0.012, 0.0), r2p=(0.007, 0.0)))
    assert c['cls'] == 'averaging'
    c = G.rho_class(_rho_rows(on=0.003, r2=(0.012, 0.0), r2p=(0.005, 0.0)))
    assert c['cls'] == 'below the gate'
    c = G.rho_class(_rho_rows(on=0.003, r2=(0.006, 0.0)))
    assert c['cls'] == 'below the gate' and not c['no_change']


def test_rho_verdict_second_draw_and_incomplete(tmp_path, monkeypatch):
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    v = G.verdict_rho97(_rho_rows(on=0.015))
    assert v['cls'] == 'raises' and v['final'] == 'pending: C2' and v['holds_on_draw2'] is None
    v = G.verdict_rho97(_rho_rows(on=0.015, d2=(0.014, 0.0)))
    assert v['holds_on_draw2'] and v['final'] == 'raises' and v['draw2']['n_wirings_up_on'] == 8
    v = G.verdict_rho97(_rho_rows(on=0.015, d2=(0.004, 0.0)))
    assert v['holds_on_draw2'] is False and v['final'] == 'raises on the first draw only'
    v = G.verdict_rho97(_rho_rows(on=0.002))
    assert v['final'] == 'no change' and v['draw2'] is None
    v = G.verdict_rho97(_rho_rows(on=0.015, drop=[(G.R2, 'all8')]))
    assert v['cls'] == 'incomplete' and 'R2' in v['missing']
    v = G.verdict_rho97(_rho_rows(on=0.015, drop=[('B2_mean2', 'P3')]))
    assert v['cls'] == 'incomplete' and 'B2_mean2 on P3' in v['missing']
    assert (tmp_path / 'verdict_rho.json').exists()


def test_rho_scores_on_made_up_distances(monkeypatch):
    """`rho_scores` reads E106's ρ against E108's distance targets, over the 120 pairs and the subsets."""
    E108 = pytest.importorskip('experiments.E108_paa_extrapolation')
    odors = list(H.APL16_ODORS)
    pairs = pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(odors) for b in odors[i + 1:]],
                                      names=['odor_a', 'odor_b'])
    rng = np.random.default_rng(1)
    tg = {t: {c: pd.Series(rng.random(len(pairs)), index=pairs) for c in ('on', 'off')} for t in G.DIST_TARGETS}
    monkeypatch.setattr(E108, 'dist_targets', lambda: tg)
    d = {c: pd.Series(rng.random(len(pairs)), index=pairs) for c in ('on', 'off')}
    r = G.rho_scores(d)
    assert r['rho_on_all'] == pytest.approx(H.spearman(d['on'], tg['all']['on']))
    assert r['rho_off_B'] == pytest.approx(H.spearman(d['off'], tg['B']['off']))
    mm = [p for p in pairs if p[0] in H.MEGAMAT12 and p[1] in H.MEGAMAT12]
    assert len(mm) == 66
    assert r['rho_on_megamat66'] == pytest.approx(H.spearman(d['on'].loc[mm], tg['all']['on'].loc[mm]))
    assert 'block' not in r


def test_fuller_input_accepts_a_last_bit_difference(tmp_path, monkeypatch):
    """B6's input check (2026-11-18): E93's fuller input rebuilt with a last-bit difference from its saved file is
    accepted and says so; a real difference is refused."""
    E93 = pytest.importorskip('experiments.E93_missing_glomeruli_input')
    E108 = pytest.importorskip('experiments.E108_paa_extrapolation')
    from experiments.E77_input_screen import frame_sha256
    orn = pd.DataFrame(np.random.default_rng(0).random((4, 16)) * 200, index=['DA1', 'DC1', 'DL1', 'VM7d'],
                       columns=list(H.APL16_ODORS))
    res = tmp_path / 'results' / E93.EXPERIMENT / 'inputs'
    res.mkdir(parents=True)
    lay = H.abbrev_columns(orn)
    lay.index.name = 'glomerulus'
    lay.to_csv(res / 'F.csv')
    (res / 'sha256.json').write_text(json.dumps({'F': dict(file=f'results/{E93.EXPERIMENT}/inputs/F.csv',
                                                           frame_sha256=frame_sha256(orn))}))
    monkeypatch.setattr(H, 'results_dir', lambda e: tmp_path / 'results' / e)
    monkeypatch.setattr(H, 'shared_dir', lambda: tmp_path)
    built = {'x': orn}
    monkeypatch.setattr(E93, 'input_of', lambda name: built['x'].copy())
    _, _, rec = G.fuller_input()
    assert 'recorded_frame_sha256' not in rec
    built['x'] = orn + 1e-12
    _, _, rec = G.fuller_input()
    assert rec['recorded_frame_sha256'] == frame_sha256(orn) and rec['frame_sha256'] == frame_sha256(built['x'])
    assert 0 < rec['frame_max_abs_diff'] <= G.INPUT_ATOL
    built['x'] = orn + 1e-6
    with pytest.raises(RuntimeError, match='differs from its saved record'):
        G.fuller_input()
    assert E108.INPUT_ATOL == G.INPUT_ATOL


def test_an_incomplete_part_b_gives_no_verdict(tmp_path, monkeypatch):
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    rd = _rows(b2=0.004)
    rd = rd[~((rd.variant == G.R2) & (rd.unit == 'all8'))]
    v = G.verdict97(rd, _boot(-0.001))
    assert v['outcome'] == 'incomplete' and 'R2' in v['missing'] and not v['needs_draw2']


def test_the_lock(tmp_path):
    p = tmp_path / 'runs_B.csv'
    with G._Lock(p):
        assert (tmp_path / 'runs_B.csv.lock').is_dir()
        with pytest.raises(TimeoutError):
            with G._Lock(p, timeout=0.6):
                pass
    assert not (tmp_path / 'runs_B.csv.lock').exists()
    G._append_locked(p, dict(tag='a', x=1))
    G._append_locked(p, dict(tag='b', x=2))
    assert pd.read_csv(p).tag.tolist() == ['a', 'b']


def test_brain_decomposition_on_made_up_effects(monkeypatch, tmp_path):
    """The decomposition's terms add up to the total variance of the per-pair effect over the grid, and the SDs of the
    pattern part are taken within a wiring and over the wirings as they are."""
    rng = np.random.default_rng(0)
    pairs = pd.MultiIndex.from_tuples([(f'a{i}', f'b{i}') for i in range(30)], names=['odor_a', 'odor_b'])
    tg = pd.Series(rng.normal(0.18, 0.05, 30), index=pairs)
    vs = [G.REF] + [f'B1_{G.MB_SHORT[m]}' for m in G.OTHER_MBS]
    eff = {(v, w): pd.Series(tg.to_numpy() + rng.normal(0, 0.03, 30), index=pairs) for v in vs for w in G.ALL_WIRINGS}
    runs = pd.DataFrame([dict(variant=v, wiring=w) for v, w in eff])
    monkeypatch.setattr(G, 'runs_table', lambda: runs)
    monkeypatch.setattr(G, 'targets', lambda: {'z8': tg})
    monkeypatch.setattr(G, '_dist', lambda v, w, ro: {'on': eff[(v, w)], 'off': pd.Series(0.0, index=pairs)})
    monkeypatch.setattr(G, '_effect', lambda d: d['on'] - d['off'])
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    df = G.brain_decomposition()
    e = np.stack([[eff[(v, w)].to_numpy() for v in vs] for w in G.ALL_WIRINGS])
    tot = float(((e - e.mean(axis=(0, 1))) ** 2).mean())
    r = df.iloc[0]
    assert r.var_wiring + r.var_brain + r.var_interaction == pytest.approx(tot)
    assert (tmp_path / 'brain_decomposition.csv').exists()
    monkeypatch.setattr(G, 'runs_table', lambda: runs.iloc[1:])
    assert G.brain_decomposition() is None


# ---------------------------------------------------------------------------
# the re-plan of 2026-11-17: the reader of Pratyush's per-synapse tables (the hemibrain and the male CNS), on made-up
# tables in the form of his two files

def _syn(kc, claw, pn, ptype, ktype, n, *, bouton=0, roi='CA(R)', roi_pre=None, side='Right', region='CA', conf=0.9,
         dist=1000.0, kinst=None):
    return [dict(bodyId_pre=pn, bodyId_post=kc, roi_pre=roi if roi_pre is None else roi_pre, roi_post=roi,
                 x_pre=100 + i, y_pre=200, z_pre=300, x_post=101 + i, y_post=201, z_post=301, confidence_pre=0.9,
                 confidence_post=conf, instance_pre=f'{ptype}_R', type_pre=ptype,
                 instance_post=kinst or f'{ktype}_R', type_post=ktype, claw=claw, dist=dist,
                 bouton=bouton, side=side, region=region) for i in range(n)]


def made_up_hemibrain() -> Tuple[pd.DataFrame, set]:
    """A made-up hemibrain table (the loader's file's columns) and its counted claws under the loader's rule."""
    r = []
    r += _syn(1, 0, 10, 'DM1_lPN', 'KCg-m', 5)                                   # counted
    r += _syn(1, 1, 11, 'VC3l_adPN', 'KCg-m', 3, bouton=0) + _syn(1, 1, 11, 'VC3l_adPN', 'KCg-m', 3, bouton=1)  # VC3
    r += _syn(1, 2, 12, 'VC3m_lvPN', 'KCg-m', 4)                                 # below 5
    r += _syn(2, -1, 10, 'DM1_lPN', 'KCab-s', 5)                                 # -1 is a claw
    r += _syn(2, 0, 10, 'DM1_lPN', 'KCab-s', 5) + _syn(2, 0, 13, 'DA1_lPN', 'KCab-s', 6)   # two PNs: dropped
    r += _syn(2, 1, 14, 'VP1d+VP4_l2PN1', 'KCab-s', 8)                            # '+': dropped
    r += _syn(3, 0, 15, 'DL1_adPN', "KCa'b'-m", 4) + _syn(3, 0, 15, 'DL1_adPN', "KCa'b'-m", 1, roi_pre='SLP(R)')
    r += _syn(3, 1, 15, 'DL1_adPN', "KCa'b'-m", 4) + _syn(3, 1, 15, 'DL1_adPN', "KCa'b'-m", 2, conf=0.3)
    r += _syn(4, 0, 10, 'DM1_lPN', None, 6, kinst='KCy(half)')                    # no KC type: dropped
    r += _syn(5, np.nan, 10, 'DM1_lPN', 'KCg-m', 6)                              # no claw label: dropped
    r += _syn(5, 0, 16, 'M_adPNm4', 'KCg-m', 6)                                  # not a task glomerulus
    r += _syn(5, 1, 10, 'DM1_lPN', 'KCg-m', 4) + _syn(5, 1, 10, 'DM1_lPN', 'KCg-m', 2, bouton=-1)
    r += _syn(5, 2, 10, 'DM1_lPN', 'KCg-m', 4) + _syn(5, 2, 10, 'DM1_lPN', 'KCg-m', 1, dist=np.inf)
    df = pd.DataFrame(r)
    df = df.rename(columns=dict(claw='anatomical_claw_corrected', dist='dist_to_root', bouton='anatomical_bouton'))
    df['anatomical_bouton'] = df['anatomical_bouton'].astype(float)
    df = df.drop(columns=['side', 'region'])
    return df, {(1, 0, 'DM1'), (1, 1, 'VC3'), (2, -1, 'DM1'), (3, 1, 'DL1')}


def made_up_malecns() -> Tuple[pd.DataFrame, Dict[str, set]]:
    """A made-up male-CNS table (his file's columns) and its counted claws under each calyx rule."""
    r = []
    r += _syn(101, 0, 20, 'VC5_lvPN', 'KCg-m', 5, roi='CA(L)', side='Left')                 # the model's VC3
    r += _syn(101, 1, 21, 'VM6_adPN', 'KCg-m', 5, roi='CentralBrain-unspecified', side='Left')   # the model's VC5
    r += _syn(102, 0, 22, 'VC3_adPN', 'KCab-c', 5)                                           # VC3 stays VC3
    r += _syn(102, 1, 22, 'VC3_adPN', 'KCab-c', 5, region='ACA')
    r += _syn(103, 0, 22, 'VC3_adPN', 'KC', 6)                                               # untyped KC: dropped
    r += _syn(104, 0, 24, 'WEDPN12', "KCa'b'-ap1", 5, side='Left', roi='CA(L)')              # dropped
    df = pd.DataFrame(r)
    df = df.rename(columns=dict(claw='anatomical_claw_corrected_post', dist='dist_to_root_post',
                                bouton='anatomical_bouton_pre', side='side_post', region='region_post'))
    df['anatomical_claw_post'] = df['anatomical_claw_corrected_post']
    df['anatomical_bouton_pre'] = df['anatomical_bouton_pre'].astype(np.int64)
    want = {'none': {(101, 0, 'VC3'), (101, 1, 'VC5'), (102, 0, 'VC3'), (102, 1, 'VC3')},
            'official': {(101, 0, 'VC3'), (102, 0, 'VC3'), (102, 1, 'VC3')},
            'region': {(101, 0, 'VC3'), (101, 1, 'VC5'), (102, 0, 'VC3')}}
    return df, want


@pytest.fixture
def made_up_tables(tmp_path, monkeypatch):
    hb, mc = tmp_path / 'hemibrain', tmp_path / 'malecns'
    hb.mkdir()
    mc.mkdir()
    h, hw = made_up_hemibrain()
    m, mw = made_up_malecns()
    h.to_parquet(hb / G.HEMIBRAIN_PN_CLAW_FILE)
    m.to_parquet(mc / f'PNbouton-to-KCclaw_Connectivity_Male-CNS_{G.MALECNS_STAMP}.parquet')
    monkeypatch.setattr(G, 'hemibrain_dir', lambda: hb)
    monkeypatch.setattr(G, 'malecns_dir', lambda: mc)
    G._CACHE.clear()
    yield dict(hemibrain=hw, malecns=mw, root=tmp_path)
    G._CACHE.clear()


def _keys(c: pd.DataFrame) -> set:
    return set(zip(c.kc_id.astype(int), c.claw_id.astype(int), c.glomerulus))


def test_reader_follows_the_loader_on_made_up_tables(made_up_tables):
    c = G.counted_claws('hemibrain')
    assert _keys(c) == made_up_tables['hemibrain']
    assert set(c.mb) == {'hemibrain'} and set(c.side) == {'R'}
    r = c.set_index(['kc_id', 'claw_id']).loc[(1, 1)]
    assert r.n_synapses == 6 and r.n_boutons == 2 and r.bouton_ids == (0, 1) and r.main_bouton_id == 0
    assert set(c.kc_class) == {'g', 'ab', "a'b'"}
    # the ROI rule and the confidence reading beside
    assert (3, 0, 'DL1') in _keys(G.counted_claws('hemibrain', roi='none'))
    assert (3, 1, 'DL1') not in _keys(G.counted_claws('hemibrain', min_confidence=0.5))
    for roi, want in made_up_tables['malecns'].items():
        assert _keys(G.counted_claws('malecns', roi=roi)) == want, roi
    m = G.counted_claws('malecns')
    assert m.groupby('mb').size().to_dict() == {'malecns_L': 2, 'malecns_R': 2}
    n = G.claw_counts(m)
    assert n.set_index(['mb', 'glomerulus']).n_claws.to_dict() == {('malecns_L', 'VC3'): 1, ('malecns_L', 'VC5'): 1,
                                                                    ('malecns_R', 'VC3'): 2}
    with pytest.raises(ValueError):
        G.counted_claws('hemibrain', roi='region')


def test_the_names():
    assert G.glomerulus_name('VC5_lvPN', 'malecns') == 'VC3' and G.glomerulus_name('VM6_adPN', 'malecns') == 'VC5'
    assert G.glomerulus_name('VC3_adPN', 'malecns') == 'VC3' and G.glomerulus_name('VC3l_adPN', 'hemibrain') == 'VC3'
    assert G.glomerulus_name('VC3m_lvPN', 'hemibrain') == 'VC3' and G.glomerulus_name('VC5_lvPN', 'hemibrain') == 'VC5'
    for t in ('VP1d+VP4_l2PN1', 'VP3+_vPN', 'M_vPNml50', 'WEDPN12', None, float('nan'), ''):
        assert G.glomerulus_name(t, 'malecns') is None
    # one lookup on the type's own name: VM6 -> VC5 does not chain on to VC3
    types = pd.Series(['VM6_adPN', 'VC5_lvPN', 'VC3_adPN', 'M_adPNm4', 'VP2+_adPN'])
    got = G._glomeruli_of_types(types, 'malecns')
    assert got[:3].tolist() == ['VC5', 'VC3', 'VC3'] and got[3:].isna().all()
    with pytest.raises(ValueError):
        G.counted_claws('nobrain')
    with pytest.raises(RuntimeError):
        G._glomeruli_of_types(pd.Series(['DA1_lPN_x']), 'hemibrain')


HELD_WORDS = ('claude_APL_modeling_YZ4TC', 'Claude_APL_analysis', 'apl_synapses', 'flywire', 'fafb', 'banc',
              '2024-09-13', '_derived', 'withheld_')


def test_the_module_names_no_held_table():
    """The module (and so every code path of the experiment) names neither the folder of the tables that Yang's own
    Claude's project derived, nor FlyWire, FAFB or BANC, nor the withheld outputs (CLAUDE.md, "Connectome data")."""
    import inspect
    src = inspect.getsource(G).lower()
    found = [w for w in HELD_WORDS if w.lower() in src]
    assert not found, found


def test_the_reader_opens_no_held_file(made_up_tables, monkeypatch):
    """Every file the reader and part A's counting open is under the two allowed folders."""
    opened = []
    real_parquet, real_csv = pd.read_parquet, pd.read_csv

    def rp(path, *a, **k):
        opened.append(str(path))
        return real_parquet(path, *a, **k)

    def rc(path, *a, **k):
        opened.append(str(path))
        return real_csv(path, *a, **k)
    monkeypatch.setattr(pd, 'read_parquet', rp)
    monkeypatch.setattr(pd, 'read_csv', rc)
    for brain in G.BRAINS2:
        for roi in (None, 'official', 'none'):
            G.counted_claws(brain, roi=roi)
    G.mb_counts(5, True)
    G.pn_labels()
    root = str(made_up_tables['root'])
    assert opened and all(p.startswith(root) for p in opened), opened
    assert not any(w.lower() in p.lower() for p in opened for w in HELD_WORDS)


def test_the_real_folders_are_the_allowed_ones():
    assert str(G.malecns_dir()).endswith('data/Male-CNS_Tom')
    assert str(G.hemibrain_dir()).endswith('from_pratyush/2025-12-05')


def test_apl_effect_reads_no_withheld_point(tmp_path):
    """kcq/apl_effect.py: E97's points of 2026-11-15 are held (no new reading) and a withheld_* folder of the archive is
    never read; the re-plan's points are read."""
    AE = pytest.importorskip('kcq.apl_effect')
    assert AE.held(G.EXPERIMENT, 'B2_mean_real') and not AE.held(G.EXPERIMENT, 'B2_mean2_real')
    for v, w in G.plan(G.ARMS + ('C2',), G.ALL_WIRINGS):
        assert not AE.held(G.EXPERIMENT, G.run_tag(v, w))
    d = tmp_path / 'E97_x'
    (d / 'withheld_2026-11-15' / 'runs').mkdir(parents=True)
    (d / 'runs_mcns').mkdir()
    pd.DataFrame({'a': [1.0]}).to_parquet(d / 'withheld_2026-11-15' / 'runs' / 'B2_mean_real_on_responses.parquet')
    pd.DataFrame({'a': [1.0]}).to_parquet(d / 'runs_mcns' / 'B2_mean2_real_on_responses.parquet')
    idx = AE._archive_index(str(d))
    assert set(idx) == {'B2_mean2_real_on_responses.parquet'}


def test_the_arm_order_and_the_claims(tmp_path, monkeypatch):
    pa = G.plan(G.ARMS, G.ALL_WIRINGS, order='arm')
    assert sorted(pa) == sorted(G.plan(G.ARMS, G.ALL_WIRINGS)) and len(pa) == 84
    assert pa[:8] == [('B2_mean2', w) for w in G.ALL_WIRINGS]
    with pytest.raises(ValueError):
        G.plan(order='odor')
    monkeypatch.setattr(G, 'out_dir', lambda: tmp_path)
    assert G.claim_point('B2_mean2_real') and not G.claim_point('B2_mean2_real')
    # a claim of this host whose process is gone is taken over; another host's is not
    (tmp_path / 'claims' / 'B2_mean2_P1').mkdir()
    (tmp_path / 'claims' / 'B2_mean2_P1' / 'owner').write_text(f'{G._host()} 999999999')
    assert G.claim_point('B2_mean2_P1')
    (tmp_path / 'claims' / 'B2_mean2_P2').mkdir()
    (tmp_path / 'claims' / 'B2_mean2_P2' / 'owner').write_text('otherhost 999999999')
    assert not G.claim_point('B2_mean2_P2')
    G.release_point('B2_mean2_real')
    assert G.claim_point('B2_mean2_real')
