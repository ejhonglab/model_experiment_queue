#!/usr/bin/env python3
"""Tests for E92 parts B to F (APL's inhibition of the KCs as a gain on the claw's drive; experiments/E92_gain.py) and
E10_joint_fit.py's --apl-divisive: the variant names; the divisive flags set and reset on a simulator's olfsysm
parameters, also around a reading that fails; the claws' loss of drive and the divisor of a made-up run, and the
control's matrix from it (the loss claw state x i, summed by KC, one value per KC); the amplitude and shape parts beside
the pattern part; the choice, the gate and the equal rule on made-up tables (the entry's thresholds); D's plan and
commands with --apl-divisive; the E10 option; and, slow (skipped with QUICK=1), the divisive form on the real wiring's
simulator: its APL-off run is the saved fit's, the flag set in place equals a simulator built with it."""

import inspect
import json
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
G = pytest.importorskip('experiments.E92_gain')

ODORS = list(H.APL16_ODORS)
QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# names

def test_variant_names_round_trip():
    cases = [((('F', 2.0), 'sub', False), 'B_F_k2', 'ref'), ((('F', 1.0), 'div', False), 'B_F_k1_div', 'B1'),
             ((('F', 0.0), 'div', False), 'B_none_div', 'B1'), ((('F', 0.5), 'div', False), 'B_F_k0.5_div', 'B1'),
             ((('F', 2.0), 'bdiv', False), 'B_F_k2_bdiv', 'B2'), ((('F', 1.0), 'both', False), 'B_F_k1_both', 'B2'),
             ((('F', 1.0), 'div', True), 'B_F_k1_div_ctrl', 'control')]
    for (w, form, ctl), name, screen in cases:
        assert G.variant_name(w, form, ctl) == name
        p = G.parse_variant(name)
        assert p['exponent'] == w[1] and p['form_of'] == form and p['control'] == ctl and p['screen'] == screen
        assert p['form'] == ('sub' if ctl else form)                     # a control runs in the subtractive form
        assert (p['apl_inh_divisive'], p['bouton_inh_divisive']) == G.FORMS[p['form']]
        assert G.variant_name(p['weighting'], p['form_of'], p['control']) == name
    assert G.parse_variant('B_none_div')['measure'] == 'none'
    assert G.ref_variant() == 'B_F_k2' and G.control_of('B_F_k0.5_div') == 'B_F_k0.5_div_ctrl'
    assert G.b1_variants() == ['B_F_k2_div', 'B_F_k1_div', 'B_F_k0.5_div', 'B_none_div']
    assert G.b2_variants() == ['B_F_k2_bdiv', 'B_F_k2_both', 'B_F_k1_both']
    for bad in ('B_F_k2_bogus', 'B_F_k2_ctrl', 'B_F_k2_div_div', 'B1_F_k2', 'B_F'):
        with pytest.raises(ValueError):
            G.parse_variant(bad)
    with pytest.raises(ValueError):
        G.variant_name(('F', 2.0), 'sub', True)


def test_the_screens_sizes():
    assert len(G.B1_SETTINGS) * len(G.E83.S_WIRINGS) == 16 and len(G.B2_SETTINGS) * len(G.E83.S_WIRINGS) == 12
    lo, hi = G.APL_BRACKET
    assert lo <= (1 / 15) / 6 and hi >= G.E83.APL_BRACKET[1] and lo < G.E83.APL_BRACKET[0]   # 15x, linear, margin 6


# ---------------------------------------------------------------------------
# the flags on a simulator

def _mp():
    olfsysm = pytest.importorskip('olfsysm')
    mp = olfsysm.ModelParams()
    mp.pn.preset_wAPLPN = True
    return mp


def test_divisive_flags_set_and_reset():
    mp = _mp()
    assert G.get_form(mp) == 'sub'
    for form, flags in G.FORMS.items():
        G.set_form(mp, form)
        assert (mp.kc.apl_inh_divisive, mp.pn.bouton_inh_divisive) == flags and G.get_form(mp) == form
        G.set_form(mp, 'sub')
        assert not mp.kc.apl_inh_divisive and not mp.pn.bouton_inh_divisive
    assert G.form_kws('div') == {'apl_inh_divisive': True} and G.form_kws('sub') == {}
    assert G.form_kws('both') == {'apl_inh_divisive': True, 'bouton_inh_divisive': True}
    with pytest.raises(ValueError):
        G.set_form(mp, 'divisive')
    mp.kc.claw_inh_pool = 1.0                                           # olfsysm refuses the two together
    with pytest.raises(RuntimeError):
        G.set_form(mp, 'div')
    G.set_form(mp, 'bdiv')                                              # the bouton form does not touch the claws
    mp.kc.claw_inh_pool = 0.0
    mp.pn.preset_wAPLPN = False
    with pytest.raises(RuntimeError):
        G.set_form(mp, 'bdiv')


def _fake_model(mp, fail: bool = False):
    seen = []
    idx = pd.MultiIndex.from_arrays([np.arange(3), ['g'] * 3], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    counts = pd.DataFrame(np.array([[3.0, 0.0], [0.0, 2.0], [0.0, 0.0]]), index=idx, columns=ODORS[:2])
    mp.kc.N = 3

    class M:
        label, target_on, size, p = 'fake', 2, 6, dict(fixed_thr=1.0, wAPLKC_scale=1.0, wKCAPL_scale=1.0, wPNAPL_scale=1.0,
                                                       wAPLPN_scale=1.0, wKCKC_scale=0.0)
        order = np.array([0, 1])
        saved = {'off': ((counts > 1.5).astype(float), counts)}
        _kc0 = np.ones((3, 1))
        _kc = np.ones((3, 1))

        def run_on(self, x=1.0):
            seen.append((G.get_form(mp), np.asarray(self.rv.kc.ext_inh).shape))
            if fail:
                raise RuntimeError('a C++ check')
            return (counts > 1.5).astype(float), counts
    m = M()
    m.mp = mp
    m.mp.kc.ext_inh_scale = 1.0
    m.rv = SimpleNamespace(kc=SimpleNamespace(ext_inh=np.zeros((0, 0))))
    return m, seen


def test_b_reading_sets_and_resets_the_form(monkeypatch):
    E69 = pytest.importorskip('experiments.E69_kc_firing_ceiling')
    monkeypatch.setattr(E69, 'bisect_factor', lambda run, n, size, br, step, label: (0.5, run(0.5), [{}], False))
    monkeypatch.setattr(G, 'gain_stats', lambda rv, mp, odors, resp, form: dict(L=np.zeros((3, 2)), D=np.zeros((3, 2)),
                                                                              row={}, form_seen=G.get_form(mp)))
    mp = _mp()
    m, seen = _fake_model(mp)
    res = G.b_reading92(m, None, 'both', ext=np.array([1.0, 2.0, 0.0]), label='t')
    assert res['apl_factor'] == 0.5 and res['stats']['form_seen'] == 'both'
    assert all(f == 'both' and sh == (3, 2) for f, sh in seen)          # set (with the static inhibition) during the runs
    assert G.get_form(mp) == 'sub' and m.rv.kc.ext_inh.size == 0 and not mp.kc.save_inh_sims
    assert res['params']['wAPLKC_scale'] == 0.5 and res['params']['wAPLPN_scale'] == 1.0   # the tied scales only
    # reset also when the reading fails
    m2, _ = _fake_model(mp, fail=True)
    with pytest.raises(RuntimeError):
        G.b_reading92(m2, None, 'div', label='t')
    assert G.get_form(mp) == 'sub' and m2.rv.kc.ext_inh.size == 0


# ---------------------------------------------------------------------------
# the loss of drive, the divisor and the control's matrix

def test_claw_loss_of_the_two_forms():
    x = np.array([[10.0, 4.0, 0.0, 6.0]])
    i = np.array([[4.0, 4.0, 3.0, 0.0]])
    s, loss = G.claw_loss(x, i, True)
    assert np.allclose(s, [[2.0, 0.8, 0.0, 6.0]]) and np.allclose(loss, x - x / (1 + i)) and np.allclose(loss, s * i)
    s0, loss0 = G.claw_loss(x, i, False)
    assert np.allclose(s0, [[6.0, 0.0, 0.0, 6.0]]) and np.allclose(loss0, [[4.0, 4.0, 0.0, 0.0]])
    # with no inhibition the two forms are the same claw (APL off)
    assert np.allclose(G.claw_loss(x, 0 * i, True)[0], G.claw_loss(x, 0 * i, False)[0])


def _fake_run(monkeypatch):
    """A made-up run: 2 KCs (3 and 2 claws), 2 boutons, 2 odors, 5 time steps at t = -0.001 .. 0.003, the pulse patched
    to cover t >= 0; APL's activity constant per odor."""
    import pyolfsysm.mb_model as mbm
    t = np.array([-0.001, 0.0, 0.001, 0.002, 0.003])
    monkeypatch.setattr(mbm, 'get_time_index', lambda mp, rv: pd.Index(t))
    monkeypatch.setattr(G, 'PULSE', (0.0, 0.0035))
    Wp = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.0, 2.0], [1.0, 1.0]])   # claws x boutons
    c2k = np.array([0, 0, 0, 1, 1])
    wk = np.array([1.0, 1.0, 1.0, 0.5, 0.5])                                       # APL>claw weights (scaled)
    wpn = np.array([0.2, 0.4])                                                     # APL>bouton weights (scaled)
    bouton = [np.array([[4.0] * 5, [1.0] * 5]), np.array([[0.0] * 5, [3.0] * 5])]  # odors: boutons x time
    inh = [np.array([[2.0] * 5]), np.array([[1.0] * 5])]                           # odors: APL x time
    rv = SimpleNamespace(kc=SimpleNamespace(wPNKC=Wp, wAPLKC=wk[:, None], claw_to_kc=c2k, inh_sims=inh),
                         pn=SimpleNamespace(bouton_sims=bouton, wAPLPN=wpn[:, None]))
    mp = SimpleNamespace(kc=SimpleNamespace(N=2))
    return rv, mp, Wp, c2k, wk, wpn, bouton, inh


def test_gain_stats_and_the_control_of_a_made_up_run(monkeypatch):
    rv, mp, Wp, c2k, wk, wpn, bouton, inh = _fake_run(monkeypatch)
    resp = np.array([[True, False], [True, True]])
    st = G.gain_stats(rv, mp, [0, 1], resp, 'div')
    for j in range(2):
        x = Wp @ bouton[j][:, 0]                         # constant over time
        i = wk * inh[j][0, 0]
        state = np.maximum(x / (1 + i), 0)
        L = np.array([(state * i)[c2k == k].sum() for k in range(2)])
        D = np.array([x[c2k == k].clip(min=0).sum() for k in range(2)])
        assert np.allclose(st['L'][:, j], L) and np.allclose(st['D'][:, j], D)
        assert np.allclose(st['L'][:, j], D - np.array([state[c2k == k].sum() for k in range(2)]))   # = x - x / (1 + i)
    # odor 0: KC 0's claws 4, 1, 2 under i = 2 -> states 4/3, 1/3, 2/3, losses 8/3, 2/3, 4/3; KC 1's 2, 5 under 1
    assert np.allclose(st['L'][:, 0], [14 / 3, 3.5])
    # the divisor at the responders' claws: odor 0 KCs 0 and 1 (1 + 2 three times, 1 + 1 twice), odor 1 KC 1 (1 + 0.5)
    pool = np.r_[[3.0] * 3 * 4, [2.0] * 2 * 4, [1.5] * 2 * 4]
    row = st['row']
    assert np.isclose(row['divisor_kc_max'], 3.0) and np.isclose(row['divisor_kc_p50'], np.quantile(pool, 0.5))
    assert np.isclose(row['divisor_kc_p10'], np.quantile(pool, 0.1)) and row['divisor_kc_share_above10'] == 0.0
    assert np.isclose(row['inh_kc_max'], 2.0) and np.isnan(row['divisor_bouton_max'])
    # the bouton form: the divisor at the boutons feeding the responders' claws (both boutons), APL's activity 2 and 1
    stb = G.gain_stats(rv, mp, [0, 1], resp, 'both')
    assert np.isclose(stb['row']['divisor_bouton_max'], 1 + 0.4 * 2) and np.isclose(stb['row']['divisor_bouton_p10'], 1 + 0.2)
    # the subtractive form: no divisor; its loss min(x, i) per claw
    sts = G.gain_stats(rv, mp, [0, 1], resp, 'sub')
    assert np.isnan(sts['row']['divisor_kc_max']) and np.allclose(sts['L'][:, 0], [2 + 1 + 2, 1 + 1])
    # the control: per KC the mean over the odors of the loss, the same for every odor
    idx = pd.MultiIndex.from_arrays([np.arange(2), ['g'] * 2], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    loss = pd.DataFrame(st['L'], index=idx, columns=ODORS[:2])
    cm = G.control_matrix92(loss)
    assert cm.shape == loss.shape and cm.index.equals(loss.index) and np.allclose(cm.to_numpy(), cm.to_numpy()[:, [0]])
    assert np.allclose(cm.iloc[:, 0], st['L'].mean(axis=1))
    # without APL's activity saved there is nothing to read
    rv.kc.inh_sims = []
    with pytest.raises(RuntimeError):
        G.gain_stats(rv, mp, [0, 1], resp, 'div')


# ---------------------------------------------------------------------------
# the parts of the pattern part

def test_spread_row_splits_the_pattern_part():
    E83 = pytest.importorskip('experiments.E83_cross_wiring_fit')
    rng = np.random.default_rng(3)
    pairs = pd.MultiIndex.from_tuples([tuple(sorted((a, b))) for i, a in enumerate(ODORS) for b in ODORS[i + 1:]],
                                      names=['odor_a', 'odor_b'])
    d = {'on': E83._cs(pd.Series(rng.uniform(0.5, 1.0, len(pairs)), index=pairs)),
         'off': E83._cs(pd.Series(rng.uniform(0.3, 0.8, len(pairs)), index=pairs))}
    tg = {'x': E83._cs(pd.Series(rng.normal(0.18, 0.06, len(pairs)), index=pairs))}
    row = G.spread_row(d, tg)
    sc = E83.effect_scores(d, tg)
    assert np.isclose(np.hypot(row['x_amplitude'], row['x_shape']), sc['x_pattern'])
    me = (d['on'] - d['off']).reindex(tg['x'].index)
    assert np.isclose(row['x_sd_ratio'], me.std(ddof=0) / tg['x'].std(ddof=0))
    assert np.isclose(row['x_pearson'], np.corrcoef(me, tg['x'])[0, 1])
    assert np.isclose(row['x_pattern_h'], sc['x_pattern'])
    # the pitfall: E83's `_cs` orientation given to the harness as it is loses the pairs oriented the other way
    lost = H.apl_effect_stats(d['on'], d['off'], tg['x'])['apl_effect_n_pairs']
    assert lost < len(pairs) and G.spread_row(d, tg)['x_pattern_h'] == row['x_pattern_h']
    # every pair of the target must be aligned (H.canonical_pairs refuses a missing pair, spread_row a misaligned one)
    with pytest.raises((KeyError, RuntimeError)):
        G.spread_row(d, {'x': tg['x'].iloc[:-1]})
    # the data's effect stretched about its mean: shape 0, the stretch as the spread ratio
    de = me - me.mean()
    st = G.spread_row(d, {'y': E83._cs(0.5 * de + 0.18)})
    assert np.isclose(st['y_shape'], 0.0, atol=1e-9) and np.isclose(st['y_sd_ratio'], 2.0)


def test_score_row_checks_the_parts(monkeypatch):
    E83 = pytest.importorskip('experiments.E83_cross_wiring_fit')
    E89 = pytest.importorskip('experiments.E89_apl_breadth_weighting')
    rng = np.random.default_rng(4)
    pairs = pd.MultiIndex.from_tuples([tuple(sorted((a, b))) for i, a in enumerate(ODORS) for b in ODORS[i + 1:]],
                                      names=['odor_a', 'odor_b'])
    d = {c: E83._cs(pd.Series(rng.uniform(0.3, 1.0, len(pairs)), index=pairs)) for c in ('on', 'off')}
    tg = {t: E83._cs(pd.Series(rng.normal(0.18, 0.06, len(pairs)), index=pairs)) for t in ('z8', 'halfA')}
    monkeypatch.setattr(E89, 'targets89', lambda: tg)
    row = G.score_row92(d)
    for t in tg:
        assert np.isclose(row[f'{t}_pattern'] ** 2, row[f'{t}_amplitude'] ** 2 + row[f'{t}_shape'] ** 2)
        assert f'{t}_pattern_h' not in row
    # a misaligned part fails loudly
    real = G.spread_row
    monkeypatch.setattr(G, 'spread_row', lambda d, tg: {k: (v * 1.1 if k.endswith('_shape') else v) for k, v in real(d, tg).items()})
    with pytest.raises(RuntimeError):
        G.score_row92(d)


# ---------------------------------------------------------------------------
# the choice, the gate, the equal rule

def _fake_readings(patterns, flagged=(), at_bound=(), T=None, raw=None, mean=None):
    """Readings of made-up settings ({variant: pattern part on S, the same against every half and all recordings}),
    with T's ({variant: pattern}) when given; raw-count losses and mean effects default to the reference's."""
    rows = []
    for v, p in patterns.items():
        pv = G.parse_variant(v)
        for unit, pat in (('S', p),) + ((('T', T[v]),) if T and v in T else ()):
            base = dict(arm='B', variant=v, screen=pv['screen'], weighting=pv['wtag'], exponent=pv['exponent'],
                        form=pv['form'], form_of=pv['form_of'], control=pv['control'], permute_seed=np.nan, unit=unit,
                        n_wirings=4, n_flagged=int(v in flagged), n_at_bound=int(v in at_bound),
                        mean_effect=(mean or {}).get(v, 0.177))
            for t in ('halfA', 'halfB', 'halfC', 'halfD', 'z8', 'z6'):
                base[f'{t}_pattern'] = pat
                base[f'{t}_loss'] = pat + 0.01
                base.update({f'{t}_amplitude': 0.6 * pat, f'{t}_shape': 0.8 * pat, f'{t}_sd_ratio': 0.9,
                             f'{t}_pearson': 0.5, f'{t}_rho': 0.45})
            cal, rw = dict(base, readout='cal'), dict(base, readout='raw', z8_loss=(raw or {}).get(v, 0.078))
            rows += [cal, rw]
            if unit == 'T':
                for w in ('P3', 'P4', 'M2', 'M3'):
                    rows += [dict(cal, unit=w), dict(rw, unit=w)]
    return pd.DataFrame(rows)


def test_choice_needs_a_divisive_form_and_the_gain():
    rd = _fake_readings({'B_F_k2': 0.050, 'B_F_k1': 0.048, 'B_F_k2_div': 0.0495})
    c = G.choose92(rd, 'A')
    assert c['chosen'] == 'B_F_k1' and not c['confirm']                  # a clipped weighting wins: nothing to confirm
    rd = _fake_readings({'B_F_k2': 0.050, 'B_F_k1': 0.052, 'B_F_k2_div': 0.0490})
    c = G.choose92(rd, 'B')
    assert c['chosen'] == 'B_F_k2_div' and np.isclose(c['gain'], 0.001) and not c['confirm']   # below 0.0015
    rd = _fake_readings({'B_F_k2': 0.050, 'B_F_k1_div': 0.0484, 'B_none_div': 0.040, 'B_F_k0.5_div': 0.041,
                         'B_F_k2_bdiv': 0.030, 'B_F_k1_div_ctrl': 0.030}, flagged=('B_none_div',),
                        at_bound=('B_F_k0.5_div',))
    c = G.choose92(rd, 'C')
    assert c['chosen'] == 'B_F_k1_div' and c['confirm'] and np.isclose(c['gain'], 0.0016)
    assert {'B_none_div', 'B_F_k0.5_div'} <= set(json.loads(c['excluded']))   # outside the bands, at a bound
    assert np.isclose(c['chosen_shape'], 0.8 * 0.0484) and np.isclose(c['ref_amplitude'], 0.6 * 0.050)
    assert G.choose92(_fake_readings({'B_F_k1_div': 0.04}), 'A')['chosen'] is None   # no reference: runs missing


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the flies' mean APL effect)")
def test_gate_on_four_directions():
    v, ctl = 'B_F_k1_div', 'B_F_k1_div_ctrl'
    runs = pd.DataFrame(dict(variant=[v] * 4, wiring=['P3', 'P4', 'M2', 'M3'], flagged=[False] * 4))
    rd = _fake_readings({'B_F_k2': 0.050, v: 0.046, ctl: 0.049}, T={'B_F_k2': 0.0538, v: 0.0490, ctl: 0.0530})
    rows = [G.gate_direction92(rd, runs, h, v, 'B_F_k2', [ctl]) for h in G.HALVES]
    assert all(np.isclose(r['g'], 0.0048) and r['n_singles_lower'] == 4 and r['v'] is True for r in rows)
    assert np.isclose(rows[0]['chosen_other_shape'], 0.8 * 0.049) and np.isclose(rows[0]['ref_z8_amplitude'], 0.6 * 0.0538)
    assert rows[0]['chosen_z6_rho'] == 0.45 and np.isclose(rows[0]['control0_other_pattern'], 0.053)
    g = G.gate92(rows)
    assert g['i'] and g['ii'] and g['n_singles'] == 16 and g['iii'] and g['iv'] and g['control'] and g['passes']
    # (ii): 12 of 16 cells is enough, (i) the mean over the four directions
    rows2 = rows[:3] + [G.gate_direction92(rd, runs, 'D', None, 'B_F_k2', [])]
    g2 = G.gate92(rows2)
    assert g2['n_singles_lower'] == 12 and g2['ii'] and np.isclose(g2['g_mean'], 0.75 * 0.0048) and not g2['passes']
    # (iii): a run on T outside the bands; the raw-count loss 0.012 above the reference's
    rows3 = [G.gate_direction92(rd, runs.assign(flagged=[True, False, False, False]), h, v, 'B_F_k2', [ctl]) for h in G.HALVES]
    assert not G.gate92(rows3)['iii']
    rd4 = _fake_readings({'B_F_k2': 0.050, v: 0.046}, T={'B_F_k2': 0.0538, v: 0.0490}, raw={v: 0.090})
    assert not G.gate92([G.gate_direction92(rd4, runs, h, v, 'B_F_k2', []) for h in G.HALVES])['iii']
    # the control gaining half as much or more fails the control rule only
    rd5 = _fake_readings({'B_F_k2': 0.050, v: 0.046, ctl: 0.049}, T={'B_F_k2': 0.0538, v: 0.0490, ctl: 0.0510})
    g5 = G.gate92([G.gate_direction92(rd5, runs, h, v, 'B_F_k2', [ctl]) for h in G.HALVES])
    assert not g5['control'] and g5['passes_i_to_iv'] and not g5['passes']


def test_equal_rule():
    rd = _fake_readings({'B_F_k2': 0.050, 'B_F_k1_div': 0.052, 'B_F_k0.5_div': 0.0545, 'B_none_div': 0.049,
                         'B_F_k2_div': 0.0501, 'B_F_k1': 0.050, 'B_F_k1_both': 0.0505, 'B_F_k1_div_ctrl': 0.05},
                        raw={'B_F_k1_both': 0.090}, mean={'B_none_div': 0.150})
    eq = {e['variant']: e for e in G.equal_rule92(rd)}
    assert set(eq) == {'B_F_k1_div', 'B_F_k0.5_div', 'B_none_div', 'B_F_k1_both'}   # divisive, exponent <= 1, no controls
    assert eq['B_F_k1_div']['equal'] and not eq['B_F_k1_div']['better'] and eq['B_F_k1_div']['screen'] == 'B1'
    assert not eq['B_F_k0.5_div']['equal']                               # 0.0045 worse
    assert not eq['B_none_div']['equal'] and eq['B_none_div']['better']  # the mean APL effect 0.027 lower
    assert not eq['B_F_k1_both']['equal'] and eq['B_F_k1_both']['screen'] == 'B2'   # the raw-count loss 0.012 above
    assert set(e['variant'] for e in G.equal_rule92(rd, screens=('B1',))) == {'B_F_k1_div', 'B_F_k0.5_div', 'B_none_div'}


def test_d_plan(tmp_path):
    rd = _fake_readings({'B_F_k2': 0.050, 'B_F_k1_div': 0.0484, 'B_F_k0.5_div': 0.0540, 'B_F_k2_bdiv': 0.040})
    p = G.d_plan(tmp_path, rd)
    assert p['allowed'] and p['exponents'] == [1.0] and 'choice on half A' in p['why']['1'][0]   # B2 is not a D setting
    assert any('equal' in x for x in p['why']['1'])
    rd0 = _fake_readings({'B_F_k2': 0.050, 'B_F_k1_div': 0.0545, 'B_F_k2_bdiv': 0.040})
    assert not G.d_plan(tmp_path, rd0)['allowed']


# ---------------------------------------------------------------------------
# E10_joint_fit.py's option and D's commands

def test_e10_apl_divisive_option():
    J = pytest.importorskip('experiments.E10_joint_fit')
    assert J.apl_divisive_kws(False) == {} and J.apl_divisive_tag(False) == ''
    assert J.apl_divisive_kws(True) == {'apl_inh_divisive': True} and J.apl_divisive_tag(True) == 'div_'
    assert inspect.signature(J.fit_one).parameters['apl_divisive'].default is False
    src = inspect.getsource(J.main)
    assert "'--apl-divisive'" in src and 'apl_divisive_tag(args.apl_divisive)' in src
    assert "ap.error('--apl-divisive with --claw-inh-pool" in src and "ap.error('--apl-divisive with --block-proxy" in src
    fsrc = inspect.getsource(J.fit_one)
    assert '**pool_kws, **div_kws' in fsrc and "NotImplementedError('--block-proxy with --apl-divisive')" in fsrc


@pytest.mark.skipif(QUICK or not _shared(), reason='runs E10_joint_fit.py up to its option checks (data loading, ~1 min)')
def test_e10_refuses_the_pooled_share_with_the_divisor():
    import subprocess
    import sys
    r = subprocess.run([sys.executable, 'experiments/E10_joint_fit.py', '-k', '2', '--apl-divisive', '--claw-inh-pool', '1',
                        '--experiment', G.EXPERIMENT], cwd=H.HERE, capture_output=True, text=True, timeout=600)
    assert r.returncode == 2 and '--apl-divisive with --claw-inh-pool' in r.stderr


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the fits of C2's family)")
def test_d_commands():
    J = pytest.importorskip('experiments.E10_joint_fit')
    c = G.d_command('real', 2.0)
    assert c[c.index('-k') + 1] == '2' and '--apl-divisive' in c and '--claw-inh-pool' not in c
    assert c[c.index('--experiment') + 1] == G.EXPERIMENT and c[c.index('--maxfev') + 1] == '150'
    assert c[c.index('--init-params') + 1] == f'results/{G.EXPERIMENT}/d_start_params.json:real_k2'
    assert '--init-from' not in c and '--wiring' not in c and '--wiring-file' not in c
    ce = G.d_command('real', 2.0, start='exact')
    assert '--init-from-exact' in ce and '--init-params' not in ce
    assert '--init-from-exact' not in G.d_command('P1', 1.0, start='exact')
    for a in ('--aplpn-max', '--rate-target', '--loss', '--apl-readout', '--tag'):
        assert c[c.index(a) + 1] == {'--aplpn-max': '3', '--rate-target': 'ctrl_frac', '--loss': 'apl',
                                     '--apl-readout': 'both', '--tag': '_cap3'}[a]
    assert G.d_tag('real', 2.0) == 'rate0.068_aplboth_div_k2_cap3'
    c1 = G.d_command('P1', 1.0)
    assert c1[c1.index('--wiring-seed') + 1] == '1'
    assert G.d_tag('P1', 1.0) == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_div_k1_cap3'
    cm = G.d_command('M1', 2.0)
    assert cm[cm.index('--wiring-tag') + 1] == 'curveball_within_type_matched-seed1'
    assert J.apl_divisive_tag(True) in G.d_tag('M1', 2.0)
    with pytest.raises(ValueError):
        G.d_command('real', 2.0, start='bogus')


# ---------------------------------------------------------------------------
# the divisive form on a real simulator (slow)

SLOW = pytest.mark.skipif(QUICK or not _shared(), reason="simulations of C2's family on the real wiring (minutes)")


@SLOW
def test_divisive_form_on_the_real_wiring():
    E83 = pytest.importorskip('experiments.E83_cross_wiring_fit')
    E91 = pytest.importorskip('experiments.E91_summed_drive_inhibition')
    m = E91.WiringModel91('real', check=True)
    for form in ('div', 'both'):
        G.set_form(m.mp, form)
        try:
            off = m.sim(m.p_off())
            on = m.run_on(0.3)
        finally:
            G.set_form(m.mp, 'sub')
        assert np.array_equal(off[1].to_numpy(), m.saved['off'][1].to_numpy())         # APL off: x / (1 + 0) = x - 0
        assert (on[1].to_numpy() != m.saved['on'][1].to_numpy()).sum() > 0
    # back in the subtractive form the simulator is the fit again
    ns, nr = E83._same({'on': m.run_on(1.0), 'off': m.sim(m.p_off())}, m.saved)
    assert ns == 0 and nr == 0
    # the flag set in place equals a simulator built with it
    G.set_form(m.mp, 'div')
    try:
        on_in = m.run_on(0.3)
    finally:
        G.set_form(m.mp, 'sub')
    mk = E91.WiringModel91('real', check=False, apl_inh_divisive=True)
    assert mk.mp.kc.apl_inh_divisive
    assert np.array_equal(mk.run_on(0.3)[1].to_numpy(), on_in[1].to_numpy())


# ---------------------------------------------------------------------------
# part F: the summed form (E91) on part A's readings

def test_f_families_and_the_refits():
    if not _shared():
        pytest.skip('needs the shared directory')
    r = G.f_runs()
    assert r.groupby('family').size().to_dict() == {'C2': 8, 'E91_B1_k2_s1': 4, 'E91_B2_k1_s1': 8, 'E91_B2_k2': 8,
                                                    'E91_refit': 4}
    assert set(G.F_FAMILIES) == set(r.family) and set(G.F_STAGE_FAMILIES) <= set(G.F_FAMILIES)
    tags = G.f_refit_tags()
    fits = pd.read_csv(H.results_dir(G.F_SOURCE) / 'fits.csv')
    assert set(tags) == set(fits.tag) and sorted(tags.values()) == sorted(G.E83.S_WIRINGS)
    assert (fits.claw_inh_pool == G.F_POOL).all() and (fits.k == G.F_K).all()
    assert G.F_B1 == 'B1_F_k2_s1' and G.F_B2_REF == 'B2_F_k2' and G.F_B2_CHOICE == 'B2_F_k1_s1'
    # the base's own family on S is E89's C2 run (E91's reference on S), on T E91's
    x = r[r.family == 'E91_B2_k2'].set_index('wiring')
    assert x.loc['real', 'source'] == 'E89' and x.loc['real', 'tag'] == 'C2_F_k2_real' and x.loc['P3', 'tag'] == 'B2_F_k2_P3'


def _fake_stage_sim(n_gl=3, n_bt=4, n_odors=3, seed=0):
    """A made-up olfsysm run: 3 glomeruli, 4 boutons, 5 claws on 2 KCs, 3 odors, 200 steps of 0.01 s from -1 s."""
    rng = np.random.default_rng(seed)
    n_steps = 200
    mp = SimpleNamespace(time_dt=0.01, time_pre_start=-1.0, time_start=-0.5, time_end=1.0)
    gl = ['DA1', 'DA1', 'VA2', 'DM1'][:n_bt]
    cols = pd.MultiIndex.from_arrays([gl, [f'b{i}' for i in range(n_bt)]], names=['glomerulus', 'bouton'])
    idx = pd.MultiIndex.from_arrays([[10, 10, 10, 11, 11], range(5)], names=[H.KC_ID_LEVEL, 'claw'])
    W = pd.DataFrame(rng.uniform(0, 1, (5, n_bt)) * (rng.uniform(0, 1, (5, n_bt)) > 0.4), index=idx, columns=cols)
    rv = SimpleNamespace(orn=SimpleNamespace(sims=[np.zeros((2, n_steps))]),
                         pn=SimpleNamespace(pn_sims=[rng.uniform(0, 50, (n_gl, n_steps)) for _ in range(n_odors)],
                                            bouton_sims=[rng.uniform(0, 80, (n_bt, n_steps)) for _ in range(n_odors)]),
                         kc=SimpleNamespace(inh_sims=[rng.uniform(0, 30, (1, n_steps)) for _ in range(n_odors)],
                                            wAPLKC=rng.uniform(0, 1, (5, 1))))
    return mp, rv, W


def test_read_stages_odor_is_e90s_reading():
    E90 = pytest.importorskip('experiments.E90_model_stages')
    mp, rv, W = _fake_stage_sim()
    odors = ['2h', 'IaA', 'pa']
    full = E90.read_stages(mp, rv, W, np.arange(3), odors, E90.MODEL_WINDOW)
    recs = [G.read_stages_odor(mp, rv, W, i, E90.MODEL_WINDOW) for i in range(3)]
    st = G.stage_frames(recs, odors, W)
    for k in G.STAGE_KEYS:
        assert st[k].index.equals(full[k].index) and list(st[k].columns) == odors
        assert np.array_equal(st[k].to_numpy(), full[k].to_numpy())
    assert np.allclose(st['apl'].loc[odors, ['win', 'pre']].to_numpy(), full['apl'].loc[odors, ['win', 'pre']].to_numpy())
    # an odor run alone: the others' saved variables may be empty, and E90's reading of all odors fails where this holds
    rv.kc.inh_sims[1] = np.zeros((0,))
    rv.pn.bouton_sims[2] = np.zeros((0, 0))
    one = G.read_stages_odor(mp, rv, W, 0, E90.MODEL_WINDOW)
    assert all(np.array_equal(one[k], recs[0][k]) for k in G.STAGE_KEYS)
    with pytest.raises(Exception):
        E90.read_stages(mp, rv, W, np.arange(3), odors, E90.MODEL_WINDOW)
    with pytest.raises(RuntimeError):
        G.read_stages_odor(mp, rv, W.iloc[:, :3], 0, E90.MODEL_WINDOW)          # the wiring's boutons are not the run's


def test_stage_save_and_load(tmp_path, monkeypatch):
    E90 = pytest.importorskip('experiments.E90_model_stages')
    mp, rv, W = _fake_stage_sim()
    st = E90.read_stages(mp, rv, W, np.arange(3), ['2h', 'IaA', 'pa'], E90.MODEL_WINDOW)
    monkeypatch.setattr(G, '_f_dir', lambda: tmp_path)
    G.save_stage_f('F_x_real', 'off', st)
    assert not list(tmp_path.glob('*spike_counts*')) and not list(tmp_path.glob('*responses*'))   # the counts stay the source's
    got = G.load_stage_f('F_x_real', 'off')
    for k in G.STAGE_KEYS:
        assert np.array_equal(got[k].to_numpy(), st[k].to_numpy()) and len(got[k].index) == len(st[k].index)
    assert np.array_equal(np.asarray(got['claw_win'].index.get_level_values(H.KC_ID_LEVEL)), [10, 10, 10, 11, 11])
    assert G.load_stage_f('F_x_real', 'on') is None
    meta = G.units_meta(W)
    assert meta['kc_of_claw'] == [10, 10, 10, 11, 11] and meta['bouton_glomerulus'] == ['DA1', 'DA1', 'VA2', 'DM1']


def test_b2_odor_params_in_the_runs_order():
    row = pd.Series(dict(tag='B2_F_k1_s1_real', thr_per_odor=json.dumps({'b @ -3': 1.2, 'a @ -3': 0.8}),
                         apl_per_odor=json.dumps({'a @ -3': 0.5, 'b @ -3': 2.5})))
    assert G.b2_odor_params(row, ['a @ -3', 'b @ -3']) == ([0.8, 1.2], [0.5, 2.5])
    with pytest.raises(RuntimeError):
        G.b2_odor_params(row, ['a @ -3', 'c @ -3'])
    with pytest.raises(RuntimeError):
        G.b2_odor_params(row, ['a @ -3'])


def test_identity_counts_every_difference():
    idx = pd.Index([1, 2], name=H.KC_ID_LEVEL)
    r = pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=idx, columns=['x', 'y'])
    s = pd.DataFrame([[3.0, 1.0], [0.0, 5.0]], index=idx, columns=['x', 'y'])
    assert G._identity((r, s), (r, s))['exact']
    d = G._identity((r, s + np.array([[0, 1], [0, 0]])), (r, s))
    assert d['counts_differ'] == 1 and d['responses_differ'] == 0 and not d['exact']
    assert not G._identity((r, s.set_axis(['y', 'x'], axis=1)), (r, s))['exact']


def _fake_somata(values: dict, readings, wirings_of: dict):
    """somata_F_means-like rows: per family, reading, condition, the runs (`values[family][wiring]`, one value for every
    group) and the populations (mean over the wirings) of S and, with eight wirings, all8."""
    A = G._A()
    rows = []
    for fam, ws in wirings_of.items():
        for rd in readings:
            for c in A.CONDS:
                vals = {w: values[fam][w] for w in ws}
                rows += [dict(kind='run', family=fam, wiring=w, cond=c, reading=rd, **{g: v for g in A.GROUPS})
                         for w, v in vals.items()]
                for name, need in G.F_POPULATIONS:
                    if set(need) <= set(ws) and name != 'T':
                        m = float(np.mean([vals[w] for w in need]))
                        rows.append(dict(kind='population', family=fam, wiring=name, cond=c, reading=rd,
                                         **{g: m for g in A.GROUPS}))
    return pd.DataFrame(rows)


def _flies_groups(mean):
    A = G._A()
    return pd.DataFrame({g: dict(mean=mean, lo=mean - 0.02, hi=mean + 0.02) for g in A.GROUPS}).T


def test_f_soma_rows_population_rules():
    A = G._A()
    readings = [f'{p}_{ms}' for p, _, _, ms in A.like_for_like()] + [f'{r}_{ms}' for r in ('raw', 'cal') for ms in A.MEASURES]
    S, T = G.E83.S_WIRINGS, G.E83.T_WIRINGS
    # eight wirings: 6 of them 0.1 below the flies, 2 above -> the population differs by 6 of 8, and on S by 3 of 4
    eight = {w: 0.4 for w in S[:3] + T[:3]}
    eight.update({S[3]: 0.55, T[3]: 0.55})
    # four wirings: 2 of 4 below -> its S population is below by more than 0.05 but does not differ
    four = {S[0]: 0.2, S[1]: 0.2, S[2]: 0.55, S[3]: 0.55}
    sm = _fake_somata({'X': eight, 'Y': four}, readings, {'X': S + T, 'Y': S})
    passes = {(g, m, u): True for g in A.GROUPS for m in A.MEASURES for u in A.UNITS}
    rows = G.f_soma_rows(sm, passes, flies=lambda u, g, r, m: _flies_groups(0.5))
    assert set(map(tuple, rows[['family', 'population']].drop_duplicates().to_numpy())) == {('X', 'S'), ('X', 'all8'), ('Y', 'S')}
    lf = rows[rows.set == 'like_for_like']
    assert len(lf) == 3 * 16 * 2 * len(A.GROUPS) and len(rows[rows.set == 'beside']) == 3 * 8 * 2 * len(A.GROUPS)
    assert lf[lf.family == 'X'].differs.all()                        # 6 of 8; 3 of 4
    y = lf[lf.family == 'Y']
    assert (y['diff'] < -0.05).all() and not y.differs.any() and (y.n_same_side == 2).all()
    assert np.allclose(lf[(lf.family == 'X') & (lf.population == 'all8')].model, np.mean(list(eight.values())))
    with pytest.raises(ValueError):
        G._n_same(5)


def test_f_step_rows_are_somata_minus_claws():
    A = G._A()
    readings = [f'{p}_{ms}' for p, _, _, ms in A.like_for_like()]
    S = G.E83.S_WIRINGS
    sm = _fake_somata({'Y': {w: 0.6 for w in S}}, readings, {'Y': S})
    claws = pd.DataFrame([dict(family='Y', wiring=w, reading=rd, **{g: c for g in A.GROUPS})
                          for w, c in zip(S, (0.3, 0.3, 0.3, 0.5)) for rd in ('session_z', 'noisy_measured')])
    passes = {(g, m, u): True for g in A.GROUPS for m in A.MEASURES for u in A.UNITS}
    got = []

    def flies_step(u, g, fr, ms, fcr):
        got.append((g, fr, fcr))
        return _flies_groups(0.15)
    rows = G.f_step_rows(sm, claws, passes, claw_ok=False, flies_step=flies_step)
    assert len(rows) == 16 * 2 * len(A.GROUPS)
    assert np.allclose(rows.model, 0.6 - 0.35) and np.allclose(rows.claw_model, 0.35)
    # three wirings' steps are 0.3, one 0.1: 3 of 4 on the side, 0.1 above the flies' 0.15 -> differs
    assert rows[rows.read].differs.all() and (rows.n_same_side == 3).all()
    # the claws as part A reads them: the session z without noise, as measured through it; noise removed is read only
    # with E90's check of the claws' noise removal
    assert set(rows.claws_model) == {'session_z', 'noisy_measured'}
    assert (~rows[rows.flies_reading == 'noise_removed'].read).all() and rows[rows.flies_reading == 'trial_mean'].read.all()
    assert {fcr for _, fr, fcr in got if fr == 'noise_removed'} == {'noise_removed'}
    with pytest.raises(RuntimeError):
        G.f_step_rows(sm, claws[claws.wiring != 'M1'], passes, claw_ok=True, flies_step=flies_step)


def test_f_line_rows_name_the_refits(tmp_path, monkeypatch):
    A = G._A()
    folder = tmp_path / G.F_SOURCE
    folder.mkdir()
    refit = 'rate0.068_aplboth_pool1_k2_cap3'
    for tag in ('B2_F_k1_s1_real', refit):
        for c in ('on', 'acute_off'):
            for k in ('spike_counts', 'responses'):
                (folder / f'{tag}_{c}_{k}.parquet').write_bytes(b'')
    (folder / 'base_P3_acute_off_spike_counts.parquet').write_bytes(b'')
    monkeypatch.setattr(A, 'ARCHIVE', tmp_path)
    files = A._find_files(G.F_SOURCE)
    r = G._f_line_row(G.F_SOURCE, refit, 'E91_refit', 'real', files)
    assert r['included'] and r['path'].endswith(f'{refit}_on_spike_counts.parquet') and r['wiring'] == 'real'
    assert G._f_line_row(G.F_SOURCE, 'base_P3', 'base', 'P3', files)['reason']
    assert A._variant_key(G.F_SOURCE, 'B2_F_k1_s1_ctrl_M3', None) == ('B2_F_k1_s1_ctrl', 'M3')


@SLOW
def test_refit_stage_run_reproduces_e91():
    st, ident = G.f_stage_refit('real')
    assert ident['on']['exact'] and ident['off']['exact']
    assert (st['on']['claw_win'].to_numpy() >= 0).all()
    assert list(st['on']['claw_win'].columns) == [H.odor_abbrev(c) for c in st['on']['spike_counts'].columns]
