#!/usr/bin/env python3
"""Tests for E96 (APL's synapses on the KCs outside the claws; experiments/E96_offclaw_apl.py): the variant names and the
arms' sizes; the counts on a made-up table (a kept claw, a dropped claw, no claw, the peduncle, a KC with none); the
component (mean 1 over the claws, a KC's sum in proportion to a_k, the proportion in the mode with all synapses); the
permutations keep the values within a claw count; size 0 gives the weights of the path at a tenth; the uniform mode is
E95's; E10_joint_fit.py's new --aplkc-flat modes; the choice, the gate, the equal rule, part D's plan (its condition, the
choice of the component, the pooled share, B6's rule) and the rule for the two starts on made-up rows; the commands and
tags of part D; the scoring path on E95's saved runs."""

import json
import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
G = pytest.importorskip('experiments.E96_offclaw_apl')
J = pytest.importorskip('experiments.E10_joint_fit')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists() and (H.shared_dir() / 'results' / 'E67_apl_effect_fits' / 'fits.csv').exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# names

def test_the_arms_and_their_sizes():
    n = len(G.S_WIRINGS)
    assert sum(len(v) for v in G.ARMS.values()) * n == 64                  # the entry's 64 points
    assert {a: len(v) * n for a, v in G.ARMS.items()} == {'B1': 12, 'B2': 12, 'B3': 16, 'B4': 8, 'B5': 4, 'B6': 12}
    for arm, vs in G.ARMS.items():
        for v in vs:
            p = G.parse_variant(v)
            assert p['arm'] == arm and p['part'] == 'B' and not p['control']
            assert p['mode'] == (None if v == G.B3_REF else G.ARM_MODE[arm])
    assert G.APL_BRACKET == (0.01, 10.0) and G.D_DOSE == 0.1 and G.B_SIZES == (0.5, 1.0, 2.0)


def test_variant_fields():
    p = G.parse_variant('B1_offclaw_d0.5')
    assert (p['mode'], p['size'], p['pool'], p['dose'], p['equal'], p['screen']) == ('offclaw', 0.5, 0.0, 0.1, True, 'B1')
    p = G.parse_variant('B3_offclawped_d2_pool1')
    assert (p['mode'], p['size'], p['pool'], p['dose'], p['equal']) == ('offclaw_ped', 2.0, 1.0, 0.1, True)
    p = G.parse_variant('B3_pool1')
    assert (p['size'], p['pool'], p['screen'], p['equal']) == (0.0, 1.0, 'B3ref', False)
    p = G.parse_variant('B4_offclawped_d1_permclaw1')
    assert (p['permute'], p['seed'], p['equal'], p['screen']) == ('claw', 1, False, 'B4')
    p = G.parse_variant('B5_offclawped_d1_dose0')
    assert (p['dose'], p['equal']) == (0.0, True)
    p = G.parse_variant('B6_allsyn_d1')
    assert (p['mode'], p['equal']) == ('allsyn', True)
    p = G.parse_variant('B1_offclaw_d1_permclaw1')                       # part C's control of B1 at d = 1
    assert p['control'] and p['control_of'] == 'B1_offclaw_d1' and p['screen'] == 'control' and p['part'] == 'control'
    p = G.parse_variant('U_uniform_d1_pool1')
    assert p['control'] and (p['mode'], p['pool'], p['dose']) == ('uniform', 1.0, 0.1)
    p = G.parse_variant('D1_uniform_d2')
    assert (p['part'], p['screen'], p['mode'], p['size'], p['dose']) == ('e95', 'e95', 'uniform', 2.0, 0.1)
    p = G.parse_variant('D3_uniform_d1_dose0')
    assert p['dose'] == 0.0 and p['part'] == 'e95'
    p = G.parse_variant('D0')
    assert p['control'] and p['size'] == 0.0 and p['dose'] == 0.1
    assert G.parse_variant('ref')['screen'] == 'ref'
    for bad in ('B1_offclawped_d1', 'B2_offclaw_d1', 'B3_offclawped_d1', 'B2_offclawped_d1_pool1', 'B4_offclawped_d2_permclaw1',
                'B4_offclawped_d1', 'B5_offclawped_d2_dose0', 'B1_offclaw_d1.0', 'B1_offclaw_d0', 'B6_offclaw_d1',
                'B5_offclawped_d1_dose0_permclaw1', 'U_uniform_d1', 'C1_bck', 'B7_offclaw_d1'):
        with pytest.raises(ValueError):
            G.parse_variant(bad)


def test_the_controls():
    assert G.control_of96('B1_offclaw_d2') == 'D1_uniform_d2'
    assert G.control_of96('B6_allsyn_d0.5') == 'D1_uniform_d0.5'
    assert G.control_of96('B5_offclawped_d1_dose0') == 'D3_uniform_d1_dose0'
    assert G.control_of96('B3_offclawped_d1_pool1') == 'U_uniform_d1_pool1'
    assert G.beside_controls('B2_offclawped_d1') == ['B2_offclawped_d1_permclaw1', 'D0']
    for c in ('D1_uniform_d2', 'D3_uniform_d1_dose0', 'U_uniform_d1_pool1') + tuple(G.beside_controls('B1_offclaw_d1')):
        G.parse_variant(c)
    for bad in ('B4_offclawped_d1_permall1', 'B3_pool1', 'ref', 'D0'):
        with pytest.raises(ValueError):
            G.control_of96(bad)


# ---------------------------------------------------------------------------
# the counts and the component on made-up data

def _syn(kc, claw, roi='CA(R)'):
    return dict(bodyId_pre=0, bodyId_post=kc, roi_pre=roi, roi_post=roi, instance_pre='APL_R', instance_post='KC',
                anatomical_claw_corrected=claw)


def _toy():
    """A made-up wiring (KC 1: claws 0, 1; KC 2: claws 0, 1, 2; KC 3: claw 0, no APL synapse) and table: kept claws, a
    claw the wiring drops, no claw (-1 and none), the peduncle."""
    from pyolfsysm.mb_model import connectome_APL_KC_synapse_counts
    idx = pd.MultiIndex.from_tuples([(1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (3, 0)], names=[H.KC_ID_LEVEL, 'claw_id'])
    rows = ([_syn(1, 0)] * 2 + [_syn(1, 1)] + [_syn(1, 7)] + [_syn(1, -1)] * 3 + [_syn(1, np.nan)]
            + [_syn(1, np.nan, 'PED(R)')] * 2
            + [_syn(2, 2)] * 4 + [_syn(2, -1)] * 3 + [_syn(2, np.nan, 'PED(R)')] * 6)
    per_kc, per_claw = connectome_APL_KC_synapse_counts(idx, synapses=pd.DataFrame(rows))
    return idx, per_kc, per_claw


def test_the_counts_on_a_made_up_table():
    idx, per_kc, per_claw = _toy()
    assert per_kc.loc[1].to_dict() == dict(n_claws=2, kept=3, kept_unassigned=0, dropped_claw=1, no_claw=4, no_claw_nan=1,
                                           peduncle=2)
    assert per_kc.loc[2].to_dict() == dict(n_claws=3, kept=4, kept_unassigned=0, dropped_claw=0, no_claw=3, no_claw_nan=0,
                                           peduncle=6)
    assert (per_kc.loc[3].drop('n_claws') == 0).all()                    # a KC with none
    assert per_claw.tolist() == [2.0, 1.0, 0.0, 0.0, 4.0, 0.0]
    assert G.offclaw_counts('offclaw', per_kc).tolist() == [4.0, 3.0, 0.0]
    assert G.offclaw_counts('offclaw_ped', per_kc).tolist() == [6.0, 9.0, 0.0]
    assert G.offclaw_counts('allsyn', per_kc).tolist() == [6.0, 9.0, 0.0]


def test_the_component_mean_and_proportion():
    idx, per_kc, per_claw = _toy()
    wn = pd.Series([2.0, 1.0, 0.0, 0.0, 4.0, 0.0], index=idx)          # proportional to the kept counts
    kc_of = np.asarray(idx.get_level_values(H.KC_ID_LEVEL))
    for mode in G.MODES:
        sh = G.component_shape(mode, wn, (per_kc, per_claw))
        assert np.isclose(sh.mean(), 1.0)
        a = G.offclaw_counts(mode, per_kc)
        s = pd.Series(sh).groupby(kc_of).sum()
        if mode == 'allsyn':                                             # the proportion in which they are counted
            v = per_claw.to_numpy() + a.reindex(kc_of).to_numpy() / per_kc.n_claws.reindex(kc_of).to_numpy()
            assert np.allclose(sh, v / v.mean())
            tot = a + per_kc.kept
            assert np.allclose((s / tot)[tot > 0], (s / tot)[tot > 0].iloc[0])
        else:
            assert np.allclose(s[a > 0] / a[a > 0], (s / a)[a > 0].iloc[0])  # a KC's sum in proportion to a_k
            assert s[3] == 0
            # the same on each of a KC's claws
            for k in (1, 2):
                assert np.ptp(sh[kc_of == k]) == 0
        comp = J.aplkc_flat_component(wn, mode, (per_kc, per_claw))
        assert np.allclose(comp, wn.mean() * sh) and np.isclose(comp.mean(), wn.mean())
    # E10's shape refuses what it cannot use
    with pytest.raises(ValueError):
        J.aplkc_flat_shape('uniform', kc_of, G.offclaw_counts('offclaw', per_kc))
    with pytest.raises(ValueError):
        J.aplkc_flat_shape('allsyn', kc_of, G.offclaw_counts('offclaw', per_kc))
    with pytest.raises(ValueError):
        J.aplkc_flat_shape('offclaw', kc_of, pd.Series([1.0, 2.0], index=[1, 2]))     # KC 3 missing
    with pytest.raises(ValueError):
        J.aplkc_flat_shape('offclaw', kc_of, pd.Series([0.0, 0.0, 0.0], index=[1, 2, 3]))


def test_the_permutations_keep_the_values_within_a_claw_count():
    rng = np.random.default_rng(0)
    n = rng.integers(1, 6, size=200)
    per_kc = pd.DataFrame(dict(n_claws=n, no_claw=rng.integers(0, 20, size=200), peduncle=rng.integers(0, 9, size=200),
                               kept=0), index=pd.Index(np.arange(200) + 10, name=H.KC_ID_LEVEL))
    a = G.offclaw_counts('offclaw_ped', per_kc)
    pc = G.offclaw_counts('offclaw_ped', per_kc, 'claw', 1)
    pa = G.offclaw_counts('offclaw_ped', per_kc, 'all', 1)
    assert pc.index.equals(a.index) and not np.array_equal(pc, a)
    for k in np.unique(n):
        assert sorted(pc[n == k]) == sorted(a[n == k])
    assert sorted(pa) == sorted(a) and not np.array_equal(pa, a)
    assert np.array_equal(G.offclaw_counts('offclaw_ped', per_kc, 'claw', 1), pc)   # one seed, the same draw
    with pytest.raises(ValueError):
        G.offclaw_counts('offclaw_ped', per_kc, 'type', 1)


def test_size_zero_and_the_uniform_mode_are_e95s():
    idx, per_kc, per_claw = _toy()
    wn = pd.Series([2.0, 1.0, 0.0, 0.0, 4.0, 0.0], index=idx)
    w_fam = pd.Series([0.5, 0.2, 3.0, 0.0, 1.0, 9.0], index=idx, name='weight')
    for mode in G.MODES + ('uniform',):
        comp = G.component_values(mode, wn, (per_kc, per_claw))
        w0 = G.e96_weights(w_fam, wn, comp, 0.0, 3.8)
        assert w0.equals(w_fam) and w0 is not w_fam                     # size 0: the family's weights exactly
    cu = G.component_values('uniform', wn)
    for d in (0.5, 1.0, 2.0):
        assert G.e96_weights(w_fam, wn, cu, d, 3.8).equals(G.E95.d_weights(w_fam, wn, 'uniform', d, 3.8))
    assert G.e96_weights(w_fam, wn, G.component_values('connectome', wn), 1.0, 3.8).equals(
        G.E95.d_weights(w_fam, wn, 'connectome', 1.0, 3.8))
    # the permuted component keeps the mean
    cp = G.component_values('offclaw_ped', wn, (per_kc, per_claw), 'all', 1)
    assert np.isclose(cp.mean(), wn.mean())
    for bad in (dict(d=-1.0), dict(d=np.nan)):
        with pytest.raises(ValueError):
            G.e96_weights(w_fam, wn, cu, bad['d'], 3.8)
    with pytest.raises(ValueError):
        G.e96_weights(w_fam, wn, -cu, 1.0, 3.8)
    # the parsed settings: the bouton path at a tenth unless dose0
    assert all(G.parse_variant(v)['dose'] == 0.1 for a in ('B1', 'B2', 'B3', 'B4', 'B6') for v in G.ARMS[a])
    assert G.parse_variant(G.ARMS['B5'][0])['dose'] == 0.0


def test_synapse_share():
    _, per_kc, _ = _toy()
    assert np.isclose(G.synapse_share('offclaw', per_kc), 7 / (7 + 7))
    assert np.isclose(G.synapse_share('offclaw_ped', per_kc), 15 / (15 + 7))


# ---------------------------------------------------------------------------
# E10_joint_fit.py's new --aplkc-flat modes

def test_e10_new_modes():
    import inspect
    for kind, tag in (('uniform', 'flatu_'), ('connectome', 'flatc_'), ('offclaw', 'flato_'), ('offclaw_ped', 'flatp_'),
                      ('allsyn', 'flata_')):
        assert J.aplkc_flat_tag(J.aplkc_flat_spec(kind, 3.8)) == tag
    assert len(set(J.APLKC_FLAT_TAGS.values())) == len(J.APLKC_FLAT_KINDS)
    assert J.aplkc_flat_spec(None) is None                                  # off by default
    src = inspect.getsource(J.main)
    assert "choices=list(APLKC_FLAT_KINDS)" in src


def test_e10_new_modes_set_the_simulators_weights():
    idx, per_kc, per_claw = _toy()
    wn = pd.Series([2.0, 1.0, 0.0, 0.0, 4.0, 0.0], index=idx)
    wk = pd.Series([0.5, 0.2, 3.0, 0.0, 1.0, 9.0], index=idx)
    J._FLAT_COUNTS.clear()
    key = tuple(map(tuple, idx.to_frame(index=False)[[H.KC_ID_LEVEL, 'claw_id']].to_numpy().tolist()))
    J._FLAT_COUNTS[key] = (per_kc, per_claw)                           # the counts as if read for this wiring
    try:
        for mode in G.MODES:
            held = wk.to_numpy().reshape(-1, 1).copy()

            def make_sim(w):
                wAPLKC_unscaled = w

                def simulate(params):
                    return wAPLKC_unscaled
                return simulate
            sim = make_sim(held)
            p = J.aplkc_flat_param(wn, wk, J.aplkc_flat_spec(mode, 2.0, 0.0))
            p['set'](sim, 0.0)
            assert np.array_equal(held.ravel(), wk.to_numpy())
            p['set'](sim, 1.5)
            assert np.allclose(held.ravel(), wk.to_numpy() + 1.5 * 2.0 * wn.mean() * G.component_shape(mode, wn, (per_kc, per_claw)))
    finally:
        J._FLAT_COUNTS.clear()


@pytest.mark.skipif(not _shared(), reason='needs the shared directory and the model data')
def test_the_counts_on_the_real_wiring_are_the_loaders():
    from pyolfsysm.mb_model import connectome_APL_weights, connectome_wPNKC
    w0 = connectome_wPNKC(prat_claws=True, prat_boutons=True)
    per_kc, per_claw = J.aplkc_flat_counts(w0.index)
    dropped = {}
    wa = connectome_APL_weights(wPNKC=w0, prat_claws=True, prat_boutons=True, dropped_counts=dropped)[0]
    assert np.array_equal(per_claw.to_numpy(), wa.to_numpy(dtype=float))
    d = dropped['wAPLKC']
    assert (d.nan_claw + d.unassigned_claw == per_kc.no_claw).all() and (d.other_claw == per_kc.dropped_claw).all()
    assert len(per_kc) == 1732 and len(per_claw) == 9867
    for mode in G.MODES:
        comp = J.aplkc_flat_component(wa / wa.mean(), mode)
        assert np.isclose(comp.mean(), 1.0) and (comp >= 0).all()


# ---------------------------------------------------------------------------
# the choice, the gate, the equal rule, part D's plan (made-up rows)

TARGETS = ('halfA', 'halfB', 'halfC', 'halfD', 'z8', 'z6')


def _fake_readings(patterns, flagged=(), at_bound=(), T=None, raw=None, mean=None, per_target=None):
    """Readings of made-up settings ({variant: pattern part on S, the same against every target unless `per_target`
    gives {variant: {target: pattern}}}), with T's ({variant: pattern}); raw-count losses and mean effects default to the
    reference's (as E95's test)."""
    rows = []
    for v, p in patterns.items():
        pv = G.parse_variant(v)
        for unit, pat in (('S', p),) + ((('T', T[v]),) if T and v in T else ()):
            base = dict(variant=v, arm=pv['arm'], part=pv['part'], screen=pv['screen'], mode=pv['mode'], size=pv['size'],
                        pool=pv['pool'], dose=pv['dose'], control=pv['control'], equal=pv['equal'], unit=unit, n_wirings=4,
                        n_flagged=int(v in flagged), n_at_bound=int(v in at_bound), mean_effect=(mean or {}).get(v, 0.177))
            for t in TARGETS:
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


def test_the_choice():
    ref = G.REF
    c = G.choose96(_fake_readings({ref: 0.050, 'B1_offclaw_d1': 0.0505}), 'A')
    assert c['chosen'] == ref and not c['confirm'] and c['control'] == ''
    c = G.choose96(_fake_readings({ref: 0.050, 'B2_offclawped_d1': 0.049}), 'B')
    assert c['chosen'] == 'B2_offclawped_d1' and np.isclose(c['gain'], 0.001) and not c['confirm']     # below 0.0015
    rd = _fake_readings({ref: 0.050, 'B4_offclawped_d1_permclaw1': 0.030, 'B3_pool1': 0.031, 'D0': 0.030,
                         'D1_uniform_d1': 0.032, 'B1_offclaw_d1_permclaw1': 0.029, 'U_uniform_d1_pool1': 0.03,
                         'B6_allsyn_d2': 0.040, 'B5_offclawped_d1_dose0': 0.041, 'B3_offclawped_d0.5_pool1': 0.047},
                        flagged=('B6_allsyn_d2',), at_bound=('B5_offclawped_d1_dose0',))
    c = G.choose96(rd, 'C')
    assert c['chosen'] == 'B3_offclawped_d0.5_pool1' and c['confirm'] and np.isclose(c['gain'], 0.003)
    assert c['control'] == 'U_uniform_d0.5_pool1'
    assert {'B6_allsyn_d2', 'B5_offclawped_d1_dose0'} <= set(json.loads(c['excluded']))
    assert np.isclose(c['chosen_shape'], 0.8 * 0.047)
    assert G.choose96(_fake_readings({'B1_offclaw_d1': 0.04}), 'A')['chosen'] is None                 # no reference


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the flies' mean APL effect)")
def test_the_gate_with_the_uniform_control_and_the_beside_controls():
    v, ref = 'B2_offclawped_d1', G.REF
    u = G.control_of96(v)
    pc, d0 = G.beside_controls(v)
    runs = pd.DataFrame(dict(variant=[v] * 4, wiring=list(G.T_WIRINGS), flagged=[False] * 4))
    rd = _fake_readings({ref: 0.050, v: 0.046, u: 0.049, pc: 0.047, d0: 0.055},
                        T={ref: 0.0538, v: 0.0490, u: 0.0530, pc: 0.0500, d0: 0.0560})
    plan = G.plan_c(rd)
    assert all(plan[h]['chosen'] == v and plan[h]['confirm'] for h in G.HALVES)
    dirs = G.gate_rows(rd, runs, plan)
    g = G.E92G.gate92(dirs)
    assert g['passes']
    bg = json.loads(dirs[0]['beside_gains'])
    assert set(bg) == {pc, d0} and np.isclose(bg[pc], 0.0038) and np.isclose(bg[d0], -0.0022)
    assert dirs[0]['rule_v_control'] == u
    # the uniform component gains half as much: rule (v) fails; the beside controls do not enter it
    rd5 = _fake_readings({ref: 0.050, v: 0.046, u: 0.049, pc: 0.046}, T={ref: 0.0538, v: 0.0490, u: 0.0510, pc: 0.0490})
    g5 = G.E92G.gate92(G.gate_rows(rd5, runs, G.plan_c(rd5)))
    assert not g5['control'] and g5['passes_i_to_iv'] and not g5['passes']


def test_the_equal_rule():
    ref = G.REF
    rd = _fake_readings({ref: 0.050, 'B1_offclaw_d1': 0.052, 'B2_offclawped_d1': 0.0545, 'B5_offclawped_d1_dose0': 0.0495,
                         'B3_offclawped_d1_pool1': 0.0505, 'B6_allsyn_d1': 0.0501, 'B6_allsyn_d2': 0.0502,
                         'B3_pool1': 0.050, 'B4_offclawped_d1_permclaw1': 0.050, 'D1_uniform_d1': 0.050, 'D0': 0.050},
                        raw={'B3_offclawped_d1_pool1': 0.090}, mean={'B6_allsyn_d1': 0.160},
                        per_target={'B1_offclaw_d1': {'halfA': 0.0525}}, at_bound=('B6_allsyn_d2',))
    eq = {e['variant']: e for e in G.equal_rule96(rd)}
    assert set(eq) == {'B1_offclaw_d1', 'B2_offclawped_d1', 'B5_offclawped_d1_dose0', 'B3_offclawped_d1_pool1',
                       'B6_allsyn_d1', 'B6_allsyn_d2'}                 # no reference, permutation, E95 point or control
    assert eq['B1_offclaw_d1']['equal']                                 # 0.0025 off on half A
    assert not eq['B2_offclawped_d1']['equal']                          # 0.0045 worse
    assert eq['B5_offclawped_d1_dose0']['equal'] and eq['B5_offclawped_d1_dose0']['better']
    assert not eq['B3_offclawped_d1_pool1']['equal']                    # raw-count loss 0.012 above
    assert not eq['B6_allsyn_d1']['equal']                              # the mean APL effect lost
    assert not eq['B6_allsyn_d2']['equal']                              # at a bound


def test_part_d_plan():
    ref = G.REF
    base = {ref: 0.050, 'B1_offclaw_d0.5': 0.0560, 'B1_offclaw_d1': 0.0550, 'B1_offclaw_d2': 0.0570,
            'B2_offclawped_d0.5': 0.0565, 'B2_offclawped_d1': 0.0547, 'B2_offclawped_d2': 0.0560,
            'B3_offclawped_d0.5_pool1': 0.0565, 'B3_offclawped_d1_pool1': 0.0540, 'B3_offclawped_d2_pool1': 0.0530,
            'B6_allsyn_d1': 0.0600}
    p = G.d_plan(_fake_readings(base))
    assert p['runs'] and np.isclose(p['above'], 0.0047)
    assert p['mode'] == 'offclaw' and p['mode_tie'] and p['mode_size'] == 1.0     # within 0.0005: B1, which assumes less
    assert p['pool'] is None and np.isclose(p['pool_diff'], -0.0007)              # B3 only 0.0007 below B2 at d = 1
    assert not p['allsyn']
    p = G.d_plan(_fake_readings(dict(base, **{'B2_offclawped_d1': 0.0540, 'B3_offclawped_d1_pool1': 0.0520})))
    assert p['mode'] == 'offclaw_ped' and not p['mode_tie'] and p['pool'] == 1.0   # B2 lower by 0.001; B3 0.002 below
    assert G.d_arm(p['mode'], p['pool']) == 'B3'
    p = G.d_plan(_fake_readings(dict(base, **{'B6_allsyn_d1': 0.0501})))
    assert p['allsyn'] and p['b6_equals'] == ['B6_allsyn_d1']
    bad = {k: (v + 0.004 if k.startswith(('B1', 'B2')) else v) for k, v in base.items()}
    p = G.d_plan(_fake_readings(bad))
    assert not p['runs'] and np.isclose(p['above'], 0.0087)                       # 0.0087 above: no refit
    assert not G.d_plan(_fake_readings({ref: 0.050}))['runs']                     # B1 / B2 not run
    assert G.d_arm('offclaw', None) == 'B1' and G.d_arm('offclaw_ped', None) == 'B2' and G.d_arm('allsyn', 1.0) == 'B6'


def test_the_rule_for_the_two_starts():
    a = dict(start='own', fit_loss=4.60, tag='a')
    b = dict(start='e95', fit_loss=4.55, tag='b')
    assert G.pick_start([a, b]) is b and G.pick_start([b, a]) is b
    assert G.pick_start([a, dict(b, fit_loss=np.nan)]) is a                       # a fit without a final loss
    assert G.pick_start([a, None]) is a and G.pick_start([]) is None


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the fits of C2's family)")
def test_part_d_commands_and_tags():
    c = G.d_command('real', 'offclaw', None, 'own', 3.84)
    assert c[c.index('-k') + 1] == '2' and c[c.index('--aplkc-flat') + 1] == 'offclaw'
    assert np.isclose(float(c[c.index('--aplpn-hold') + 1]), 0.1 * G.E94.aplpn_fit_multiple('real'))
    assert float(c[c.index('--aplkc-flat-unit') + 1]) == 3.84 and c[c.index('--aplkc-flat-max') + 1] == '8'
    assert c[c.index('--init-params') + 1] == f'results/{G.EXPERIMENT}/d_start_params.json:real_offclaw_own'
    for a, v in (('--rate-target', 'ctrl_frac'), ('--loss', 'apl'), ('--apl-readout', 'both'), ('--tag', '_g1'),
                 ('--maxfev', '200'), ('--experiment', G.EXPERIMENT)):
        assert c[c.index(a) + 1] == v
    assert '--claw-inh-pool' not in c
    cp = G.d_command('P1', 'offclaw_ped', 1.0, 'e95', 2.5)
    assert cp[cp.index('--claw-inh-pool') + 1] == '1' and cp[cp.index('--tag') + 1] == '_g1e'
    assert cp[cp.index('--init-params') + 1].endswith(':P1_offclawped_pool1_e95')
    assert G.d_tag('real', 'offclaw', None, 'own') == 'rate0.068_aplboth_aplpnhold_flato_k2_g1'
    assert G.d_tag('P1', 'offclaw_ped', 1.0, 'e95') == ('wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_pool1_'
                                                        'aplpnhold_flatp_k2_g1e')
    assert G.d_tag('M1', 'allsyn', None, 'own').endswith('rate0.068_aplboth_aplpnhold_flata_k2_g1')
    assert '--wiring-file' in G.d_command('M1', 'allsyn', None, 'own', 3.2)


@pytest.mark.skipif(not _shared(), reason="needs E95's fits in the shared directory")
def test_the_e95_start():
    f = pd.read_csv(G.H.results_dir('E95_static_apl') / 'fits.csv').drop_duplicates('tag', keep='last').set_index('tag')
    r = f.loc[G.E95.f_tag('real', 'uniform')]
    s = G.e95_start('real', float(r['aplkc_flat_unit']))
    assert s['aplkc_flat'] == float(r['aplkc_flat_size']) and s['wKCKC_scale'] == float(r['kc_kc_inh'])
    assert s['fixed_thr'] == float(r['fixed_thr']) and s['start'] == 'e95'
    with pytest.raises(RuntimeError):
        G.e95_start('real', float(r['aplkc_flat_unit']) * 1.01)               # another unit is refused


# ---------------------------------------------------------------------------
# the scoring path on saved runs (no simulation)

@pytest.mark.skipif(QUICK or not _shared() or not (G.H.outputs_dir('E95_static_apl') / 'runs').exists(),
                    reason="needs E95's saved runs (viki's outputs) and the shared directory; not with QUICK=1")
def test_the_readings_on_saved_runs(tmp_path):
    """The scoring path on E95's saved part D runs standing in for this experiment's points on S."""
    stand = {'B1_offclaw_d1': 'D1_uniform_d1', 'B2_offclawped_d1': 'D2_connectome_d1'}
    rows = []
    for v, e in stand.items():
        for w in G.S_WIRINGS:
            rows.append(dict(tag=f'{e}_{w}', variant=v, wiring=w, source='E95', at_bound=False, realism_flag='', factor=1.0))
    pd.DataFrame(rows).to_csv(tmp_path / 'runs_B.csv', index=False)
    rd = G.readings96(tmp_path, readouts=('cal', 'raw'))
    s = rd[(rd.unit == 'S') & (rd.readout == 'cal')].set_index('variant')
    assert set(s.index) >= {G.REF, 'B1_offclaw_d1', 'B2_offclawped_d1', 'D1_uniform_d1', G.D0}
    assert np.isclose(s.loc[G.REF, 'z8_pattern'], G.UNCHANGED_S['all'], atol=5e-5)
    e95 = pd.read_csv(G.H.results_dir('E95_static_apl') / 'readings.csv')
    e95 = e95[(e95.unit == 'S') & (e95.readout == 'cal')].set_index('variant')
    assert np.isclose(s.loc['B1_offclaw_d1', 'z8_pattern'], e95.loc['D1_uniform_d1', 'z8_pattern'], rtol=1e-9)
    v = G.verdict96(tmp_path, rd)
    assert set(v['choice']) == set(G.HALVES) and 'D_plan' in v
    text = G.comparison_table(tmp_path)
    assert '| B1_offclaw_d1 |' in text and 'Verdict' in text


# ---------------------------------------------------------------------------
# run_point's rows on saved data (no simulation; part 0's cache and the shared directory)

def _saved_model(monkeypatch):
    """A stand-in for a StaticModel on the real wiring at its fit, from the saved fit, E94's wiring cache and part 0's
    per-KC tables (no simulator); reading96 replaced by a reading of the saved APL-on run with the fit's stage
    matrices."""
    from types import SimpleNamespace
    E83 = G.E83
    f = E83.saved_fit('real')
    W, _ = G.E94.wiring_W('real')
    kc, claws = G.load_kc_tables96('real')
    idx = pd.MultiIndex.from_frame(claws[[H.KC_ID_LEVEL, 'claw_id']])
    wn = pd.Series(claws.wnorm.to_numpy(), index=idx)
    per_claw = pd.Series(claws.kept.to_numpy(), index=idx)
    per_kc = kc.set_index(H.KC_ID_LEVEL)[['n_claws', 'kept', 'kept_unassigned', 'dropped_claw', 'no_claw', 'no_claw_nan',
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
    monkeypatch.setattr(G, 'reading96', fake_reading)
    monkeypatch.setattr(G, 'wiring_wnorm', lambda w: wn)
    monkeypatch.setattr(G, 'wiring_counts', lambda w='real': (per_kc, per_claw))
    monkeypatch.setattr(G.E95, 'd_unit', lambda model, L=None: 3.0)
    return m, seen, (wn, per_kc, per_claw)


@pytest.mark.skipif(not _shared() or not (G.local_dir('W') / 'real_kc.parquet').exists(),
                    reason="needs the shared directory and part 0's per-KC tables")
def test_run_point_rows_on_saved_data(monkeypatch):
    m, seen, (wn, per_kc, per_claw) = _saved_model(monkeypatch)
    saved = {}
    monkeypatch.setattr(G, 'save_run96', lambda tag, o, sub='runs', **k: saved.update({tag: k}))
    vs = [v for a in G.ARMS for v in G.ARMS[a]] + ['U_uniform_d1_pool1', 'B1_offclaw_d1_permclaw1', 'D0', 'D1_uniform_d1']
    for v in vs:
        row, _ = G.run_point(m, v, part='test', cc={'queue_commit': 'x'})
        pv = G.parse_variant(v)
        tag = f'{v}_real'
        assert row['variant'] == v and row['tag'] == tag and row['component_unit'] == 3.0
        assert row['pool'] == pv['pool'] and row['dose'] == pv['dose'] and row['size'] == pv['size']
        assert all(k in row for k in ('rate_on', 'realism_flag', 'max_count_on', 'K_run', 'bouton_share', 'comp_share'))
        assert np.isclose(row['aplpn_multiple'], G.E94.aplpn_multiple('real', pv['dose']))
        rpv, wk = seen[tag]
        assert rpv['pool'] == pv['pool'] and rpv['dose'] == pv['dose']
        comp = saved[tag]['component']
        if pv['size'] > 0:
            add = wk.to_numpy() - m.wAPLKC.to_numpy()
            assert np.isclose(add.mean(), pv['size'] * 3.0 * wn.mean())          # d x unit x the connectome's mean
            assert row['comp_share'] > 0 and np.isclose(comp.component_sum.sum(), add.sum())
            if pv['mode'] in G.MODES:
                assert np.isclose(row['synapse_share'], G.synapse_share(pv['mode'], per_kc))
                a = G.offclaw_counts(pv['mode'], per_kc, pv['permute'], pv['seed'] or G.PERM_SEED)
                assert np.allclose(comp.set_index(H.KC_ID_LEVEL).a_k.sort_index(), a.sort_index())
        else:
            assert wk is None and comp is None and row['comp_share'] == 0
    # B4's permutations keep each claw count's values; B2 at d = 1 against B4 within a claw count
    c2 = saved['B2_offclawped_d1_real']['component'].set_index(H.KC_ID_LEVEL)
    c4 = saved['B4_offclawped_d1_permclaw1_real']['component'].set_index(H.KC_ID_LEVEL)
    for k, g in c2.groupby('n_claws'):
        assert sorted(g.a_k) == sorted(c4.loc[g.index, 'a_k'])
    assert not np.array_equal(c2.a_k.to_numpy(), c4.a_k.to_numpy())
