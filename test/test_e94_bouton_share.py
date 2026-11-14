#!/usr/bin/env python3
"""Tests for E94 (the PN boutons' share of APL's inhibition; experiments/E94_bouton_share.py) and E10_joint_fit.py's
--aplpn-hold: the variant names and their screens and kinds of equal; the dose reaches only the APL->bouton scale and
dose 0 runs a simulator at the fit's parameters with that scale 0 (a simulator built without the path); the form, the
dose, the weights and the static inhibition reset after a reading, also one that fails; the share of the claws' loss made
at the boutons on a toy with two claws, the boutons' readings, the KCs' drives and the control's matrix in its two forms;
the choice, the gate, the three kinds of equal, the anatomical dose, D1's rule and D2's plan on made-up rows; D1's
command and tag; the E10 option (without it the fit's call is unchanged); and, slow (skipped with QUICK=1), dose 0 on the
real wiring against a simulator built without the path."""

import inspect
import json
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
G = pytest.importorskip('experiments.E94_bouton_share')

ODORS = list(H.APL16_ODORS)
QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists() and (H.shared_dir() / 'results' / 'E67_apl_effect_fits' / 'fits.csv').exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# names

def test_variant_names_screens_and_kinds():
    cases = [(('sub', 2, 1), 'B_sub_k2_d1', 'ref', None), (('sub', 1, 1), 'B_sub_k1_d1', 'ref', None),
             (('gain', 2, 1), 'B_gain_k2_d1', 'E92', None), (('sub', 2, 0), 'B_sub_k2_d0', 'B1', 'lower_share'),
             (('sub', 2, 0.5), 'B_sub_k2_d0.5', 'B1', 'lower_share'), (('sub', 2, 2), 'B_sub_k2_d2', 'B1', None),
             (('sub', 1, 4), 'B_sub_k1_d4', 'B2', 'weaker_weighting'), (('gain', 2, 0.05), 'B_gain_k2_d0.05', 'B3', 'gain')]
    for (form, k, d), name, screen, kind in cases:
        assert G.variant_name(form, k, d) == name
        p = G.parse_variant(name)
        assert (p['form'], p['exponent'], p['dose'], p['screen'], p['equal_kind']) == (form, k, d, screen, kind)
        assert p['bouton_inh_divisive'] == (form == 'gain') and not p['control']
    c = G.parse_variant('B_sub_k2_d4_ctrl')
    assert c['control'] and c['screen'] == 'control' and c['control_of'] == 'B_sub_k2_d4' and c['equal_kind'] is None
    assert G.control_of('B_gain_k2_d0.15') == 'B_gain_k2_d0.15_ctrl'
    for bad in ('B_F_k2_bdiv', 'B_sub_k2', 'B_div_k2_d1', 'B_sub_k2_d1_ctrl_ctrl', 'C1_F_k1'):
        with pytest.raises(ValueError):
            G.parse_variant(bad)
    with pytest.raises(ValueError):
        G.variant_name('bdiv', 2, 1)
    assert G.REF_VARIANTS == ('B_sub_k2_d1', 'B_sub_k1_d1') and G.e92_variant() == 'B_gain_k2_d1'


def test_the_arms_sizes():
    n = len(G.S_WIRINGS)
    assert len(G.arm_variants('B1')) * n == 24 and len(G.arm_variants('B2')) * n == 8 and len(G.arm_variants('B3')) * n == 12
    assert 1.0 not in G.B1_DOSES + G.B2_DOSES + G.B3_DOSES                 # dose 1 is the reference of each arm
    assert set(G.ANATOMICAL_DOSES) <= set(G.B1_DOSES) and G.APL_BRACKET == (0.01, 10.0)
    assert all(G.parse_variant(v)['screen'] == arm for arm in G.ARMS for v in G.arm_variants(arm))


# ---------------------------------------------------------------------------
# the dose and the form on a simulator

def _mp():
    olfsysm = pytest.importorskip('olfsysm')
    mp = olfsysm.ModelParams()
    mp.pn.preset_wAPLPN = True
    return mp


P_FIT = dict(fixed_thr=200.0, wAPLKC_scale=0.004, wKCAPL_scale=0.002, wAPLPN_scale=0.376, wPNAPL_scale=0.03,
             wKCKC_scale=60.0)


def _fake_model(mp=None, fail: bool = False):
    """A WiringModel94 without a simulation: its simulator records the parameters, the form and the static inhibition
    of every run and returns fixed outputs."""
    seen = []
    idx = pd.MultiIndex.from_arrays([np.arange(3), ['g'] * 3], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    counts = pd.DataFrame(np.array([[3.0, 0.0], [0.0, 2.0], [0.0, 0.0]]), index=idx, columns=ODORS[:2])
    m = object.__new__(G.WiringModel94)
    m.label, m.target_on, m.size, m.p, m.dose = 'fake', 2, 6, dict(P_FIT), 1.0
    m.order = np.array([0, 1])
    m.saved = {'off': ((counts > 1.5).astype(float), counts), 'on': ((counts > 1.5).astype(float), counts)}
    m._kc0, m._kc = np.ones((3, 1)), np.ones((3, 1))
    m.W = None
    if mp is not None:
        mp.kc.N = 3
        mp.kc.ext_inh_scale = 1.0
        m.mp = mp
        m.rv = SimpleNamespace(kc=SimpleNamespace(ext_inh=np.zeros((0, 0))))

    def sim(params):
        seen.append(dict(params=dict(params), form=G.get_form94(m.mp) if mp is not None else None,
                         ext=np.asarray(m.rv.kc.ext_inh).shape if mp is not None else None))
        if fail:
            raise RuntimeError('a C++ check')
        return (counts > 1.5).astype(float), counts
    m.sim = sim
    return m, seen


def test_the_dose_reaches_only_the_apl_bouton_scale():
    m, seen = _fake_model()
    assert m.params_on(1.0) == P_FIT                                          # dose 1, factor 1: the fit exactly
    for dose in (0.0, 0.1, 0.25, 0.5, 2.0, 4.0, 0.05):
        m.dose = dose
        p = m.params_on(1.0)
        assert {k for k in P_FIT if p[k] != P_FIT[k]} <= {'wAPLPN_scale'}
        assert p['wAPLPN_scale'] == P_FIT['wAPLPN_scale'] * dose
        q = m.params_on(0.5)                                                  # the tied factor: E83's three scales
        assert {k for k in P_FIT if q[k] != p[k]} == set(G.E83.TIED_ON)
        assert all(q[k] == P_FIT[k] * 0.5 for k in G.E83.TIED_ON) and q['wAPLPN_scale'] == p['wAPLPN_scale']
        m.run_on(0.5)
        assert seen[-1]['params'] == q
    assert m.params_on(1.0, dose=3.0)['wAPLPN_scale'] == P_FIT['wAPLPN_scale'] * 3.0   # an explicit dose


def test_dose_0_is_a_simulator_built_without_the_path():
    """At dose 0 the parameters of every APL-on run are those a simulator built at the fit with wAPLPN_scale 0 runs at
    (part 0 checks the counts on the real wiring; the slow test below too)."""
    m, seen = _fake_model()
    m.dose = 0.0
    built_without = dict(P_FIT, wAPLPN_scale=0.0)
    assert m.params_on(1.0) == built_without
    m.run_on(1.0)
    assert seen[-1]['params'] == built_without
    m.dose = 1.0
    assert m.params_on(1.0) != built_without


def test_b_reading_sets_and_resets_the_form_and_the_dose(monkeypatch):
    E69 = pytest.importorskip('experiments.E69_kc_firing_ceiling')
    E90 = pytest.importorskip('experiments.E90_model_stages')
    monkeypatch.setattr(E69, 'bisect_factor', lambda run, n, size, br, step, label: (0.5, run(0.5), [{}], False))
    got = {}

    def fake_stages(mp, rv, W, order, odors, window):
        got['stage_form'] = G.get_form94(mp)
        got['save_inh_sims'] = bool(mp.kc.save_inh_sims)
        return {'odors': list(odors)}
    monkeypatch.setattr(E90, 'read_stages', fake_stages)
    monkeypatch.setattr(G.E92G, 'gain_stats', lambda rv, mp, odors, resp, form: dict(L=np.zeros((3, 2)), D=np.zeros((3, 2)),
                                                                                   row={}, form92=form))
    mp = _mp()
    m, seen = _fake_model(mp)
    res = G.b_reading94(m, None, 0.05, 'gain', ext=np.array([1.0, 2.0, 0.0]), label='t')
    assert res['apl_factor'] == 0.5 and res['gain']['form92'] == 'bdiv' and got['stage_form'] == 'gain'
    assert all(s['form'] == 'gain' and s['ext'] == (3, 2) for s in seen)       # set (with the static inhibition) in the runs
    assert all(s['params']['wAPLPN_scale'] == P_FIT['wAPLPN_scale'] * 0.05 for s in seen)
    assert res['params']['wAPLPN_scale'] == P_FIT['wAPLPN_scale'] * 0.05 and res['params']['wAPLKC_scale'] == P_FIT['wAPLKC_scale'] * 0.5
    assert len(seen) == 2                                                      # the bisection's run and the run again
    # everything reset after
    assert G.get_form94(mp) == 'sub' and not mp.pn.bouton_inh_divisive and m.dose == 1.0
    assert m.rv.kc.ext_inh.size == 0 and not mp.kc.save_inh_sims
    # the subtraction at a dose: the form stays the subtraction
    seen.clear()
    G.b_reading94(m, None, 4.0, 'sub', label='t')
    assert all(s['form'] == 'sub' and s['params']['wAPLPN_scale'] == P_FIT['wAPLPN_scale'] * 4.0 for s in seen)
    # reset also when the reading fails
    m2, _ = _fake_model(mp, fail=True)
    with pytest.raises(RuntimeError):
        G.b_reading94(m2, None, 0.15, 'gain', label='t')
    assert G.get_form94(mp) == 'sub' and m2.dose == 1.0 and m2.rv.kc.ext_inh.size == 0
    with pytest.raises(ValueError):
        G.set_form94(mp, 'bdiv')


# ---------------------------------------------------------------------------
# the claws' loss, the boutons, the drives, the control's matrix

def test_claw_loss_share_on_two_claws():
    # two claws, one bouton each (weights 1 and 2), one odor: APL off the boutons fire 10 and 5 Hz, APL on 8 and 5;
    # the claws keep 6 and 9 after APL's inhibition at the claws
    W = np.array([[1.0, 0.0], [0.0, 2.0]])
    s = G.claw_loss_shares(W, np.array([[10.0], [5.0]]), np.array([[8.0], [5.0]]), np.array([[6.0], [9.0]]),
                           kc_of=np.array([7, 8]))
    # drives off 10 and 10, on 8 and 10; loss 4 and 1 (5 of 20), of which 2 at the boutons
    assert np.isclose(s['drive_removed'], 5 / 20) and np.isclose(s['bouton_share'], 2 / 5) and np.isclose(s['claw_part'], 3 / 5)
    assert np.allclose(s['by_odor'], [0.4]) and np.isclose(s['bouton_share_1-4'], 0.4) and np.isnan(s['bouton_share_9+'])
    assert s['n_claws_1-4'] == 2 and s['n_claws_9+'] == 0
    # two odors; no subtraction at the boutons for the second: its share 0
    s2 = G.claw_loss_shares(W, np.array([[10.0, 4.0], [5.0, 2.0]]), np.array([[8.0, 4.0], [5.0, 2.0]]),
                            np.array([[6.0, 3.0], [9.0, 3.0]]))
    assert np.allclose(s2['by_odor'], [0.4, 0.0]) and np.isclose(s2['bouton_share'], 2 / (5 + 2))
    with pytest.raises(ValueError):
        G.claw_loss_shares(W, np.ones((2, 1)), np.ones((2, 1)), np.ones((3, 1)))


def test_claw_groups_by_the_kcs_claws():
    kc = np.array([1] * 3 + [2] * 6 + [3] * 10)
    n = G.claw_counts(kc)
    assert list(n[:3]) == [3] * 3 and list(n[3:9]) == [6] * 6 and list(n[9:]) == [10] * 10
    W = np.eye(len(kc))
    b_off = np.full((len(kc), 1), 10.0)
    b_on = np.r_[np.full(3, 5.0), np.full(6, 10.0), np.full(10, 10.0)][:, None]
    c_on = np.r_[np.full(3, 5.0), np.full(6, 5.0), np.full(10, 10.0)][:, None]
    s = G.claw_loss_shares(W, b_off, b_on, c_on, kc_of=kc)
    assert s['bouton_share_1-4'] == 1.0 and s['bouton_share_5-8'] == 0.0 and np.isnan(s['bouton_share_9+'])
    assert s['drive_removed_9+'] == 0.0 and np.isclose(s['bouton_share'], 15 / 45)


def test_bouton_realism():
    off = np.array([[10.0, 4.0], [6.0, 0.0]])
    on = np.array([[8.0, 0.0], [6.0, 0.0]])
    r = G.bouton_realism(off, on)
    assert np.isclose(r['bouton_share_removed'], 6 / 20) and r['bouton_zero_share'] == 0.5 and r['bouton_zero_share_off'] == 0.25
    assert np.isclose(r['bouton_sub_p50'], np.median([2, 4, 0, 0]))
    # the A3 transforms: dose 1 is the saved run, dose 0 APL off, dose 2 clipped
    assert np.array_equal(G.a3_transform(off, on, 'sub', 1.0), on) and np.array_equal(G.a3_transform(off, on, 'sub', 0.0), off)
    assert np.array_equal(G.a3_transform(off, on, 'sub', 2.0), np.array([[6.0, 0.0], [6.0, 0.0]]))
    assert np.allclose(G.a3_transform(off, on, 'gain', 0.5, np.full((2, 2), 2.0)), off / 2.0)


def _wiring_toy():
    idx = pd.MultiIndex.from_arrays([[11, 11, 12, 13], [0, 1, 0, 0]], names=[H.KC_ID_LEVEL, 'claw_id'])
    cols = pd.MultiIndex.from_arrays([['DA1', 'VA2'], [0, 1]], names=['glomerulus', 'bouton_id'])
    W = pd.DataFrame([[1.0, 0.0], [0.0, 2.0], [0.5, 0.0], [0.0, 1.0]], index=idx, columns=cols)
    kc_index = pd.MultiIndex.from_arrays([[13, 11, 12], ['g', 'ab', 'ab']], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    return W, kc_index


def test_kc_drive_sums_the_claws_in_the_outputs_order():
    W, kc_index = _wiring_toy()
    b = np.array([[10.0, 0.0], [4.0, 1.0]])                    # boutons x odors
    d = G.kc_drive(W, b, kc_index, ['2h', 'IaA'])
    assert d.index.equals(kc_index) and list(d.columns) == ['2h', 'IaA']
    assert np.allclose(d.to_numpy(), [[4.0, 1.0], [10.0 + 8.0, 0.0 + 2.0], [5.0, 0.0]])
    with pytest.raises(RuntimeError):
        G.kc_drive(W, b, pd.MultiIndex.from_arrays([[99], ['g']], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL]), ['2h', 'IaA'])


def test_the_control_matrix_in_its_two_forms():
    W, kc_index = _wiring_toy()
    fr = lambda a: pd.DataFrame(np.asarray(a, dtype=float), index=kc_index, columns=['2h', 'IaA'])
    d_off, d_ref = fr([[6, 4], [20, 10], [8, 2]]), fr([[5, 3], [16, 8], [6, 2]])
    # a dose above 1: every KC loses more than with the reference -> the reference's path plus the mean change
    d_hi = fr([[4, 2], [12, 7], [6, 1]])
    m, form, dose = G.control_matrix94(d_ref, d_hi, d_off)
    assert form == 'ref_plus_change' and dose == 1.0 and m.index.equals(kc_index)
    assert np.allclose(m.iloc[:, 0], [1.0, 2.5, 0.5]) and np.allclose(m.to_numpy(), m.to_numpy()[:, [0]])
    # a dose below 1 (or the gain): some KC gains drive -> the path removed plus the mean loss at the boutons
    d_lo = fr([[6, 4], [18, 9], [6, 2]])
    m2, form2, dose2 = G.control_matrix94(d_ref, d_lo, d_off)
    assert form2 == 'none_plus_removal' and dose2 == 0.0 and np.allclose(m2.iloc[:, 0], [0.0, 1.5, 1.0])
    assert (m2.to_numpy() >= 0).all() and (m.to_numpy() >= 0).all()          # olfsysm refuses a negative ext_inh_kc
    with pytest.raises(ValueError):
        G.control_matrix94(d_ref, d_hi.iloc[::-1], d_off)


# ---------------------------------------------------------------------------
# the choice, the gate, the equal rule, the anatomical dose, D1's rule, D2's plan

TARGETS = ('halfA', 'halfB', 'halfC', 'halfD', 'z8', 'z6')


def _fake_readings(patterns, flagged=(), at_bound=(), T=None, raw=None, mean=None, per_target=None):
    """Readings of made-up settings ({variant: pattern part on S, the same against every target unless `per_target`
    gives {variant: {target: pattern}}}), with T's ({variant: pattern}); raw-count losses and mean effects default to the
    reference's."""
    rows = []
    for v, p in patterns.items():
        pv = G.parse_variant(v)
        for unit, pat in (('S', p),) + ((('T', T[v]),) if T and v in T else ()):
            base = dict(arm='B', variant=v, screen=pv['screen'], form=pv['form'], exponent=pv['exponent'], dose=pv['dose'],
                        control=pv['control'], equal_kind=pv['equal_kind'], unit=unit, n_wirings=4,
                        n_flagged=int(v in flagged), n_at_bound=int(v in at_bound), n_above_bound=int(pv['dose'] > 1.08),
                        mean_effect=(mean or {}).get(v, 0.177))
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
    ref, k1 = G.ref_variant(), G.ref_k1_variant()
    rd = _fake_readings({ref: 0.050, 'B_sub_k2_d0.5': 0.0505, 'B_gain_k2_d0.05': 0.052})
    c = G.choose94(rd, 'A')
    assert c['chosen'] == ref and not c['confirm'] and c['gain'] == 0.0
    rd = _fake_readings({ref: 0.050, 'B_sub_k2_d2': 0.049})
    c = G.choose94(rd, 'B')
    assert c['chosen'] == 'B_sub_k2_d2' and np.isclose(c['gain'], 0.001) and not c['confirm']        # below 0.0015
    rd = _fake_readings({ref: 0.050, k1: 0.045, 'B_sub_k2_d4': 0.040, 'B_sub_k2_d0.1': 0.041, 'B_gain_k2_d0.15': 0.0484,
                         'B_sub_k1_d2': 0.0490, 'B_sub_k2_d4_ctrl': 0.030, 'B_gain_k2_d1': 0.030},
                        flagged=('B_sub_k2_d4',), at_bound=('B_sub_k2_d0.1',))
    c = G.choose94(rd, 'C')
    assert c['chosen'] == k1 and not c['confirm']                       # E89's point: read, never confirmed here
    assert {'B_sub_k2_d4', 'B_sub_k2_d0.1'} <= set(json.loads(c['excluded']))   # outside the bands, at a bound
    rd = rd[rd.variant != k1]
    c = G.choose94(rd, 'D')
    assert c['chosen'] == 'B_gain_k2_d0.15' and c['confirm'] and np.isclose(c['gain'], 0.0016) and c['screen'] == 'B3'
    assert np.isclose(c['chosen_shape'], 0.8 * 0.0484) and np.isclose(c['ref_amplitude'], 0.6 * 0.050)
    assert G.choose94(_fake_readings({'B_sub_k2_d2': 0.04}), 'A')['chosen'] is None                # no reference
    # the bouton entries at zero of the chosen setting, named above a fifth
    runs = pd.DataFrame(dict(variant=['B_gain_k2_d0.15'] * 4, wiring=list(G.S_WIRINGS), bouton_zero_share=[0.3] * 4,
                             aplpn_multiple=[0.4] * 4))
    c2 = G.choose94(rd, 'D', runs)
    assert c2['zero_share_named'] and np.isclose(c2['zero_share'], 0.3)


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the flies' mean APL effect)")
def test_the_gate_on_four_directions():
    v, ref = 'B_sub_k2_d2', G.ref_variant()
    ctl = G.control_of(v)
    runs = pd.DataFrame(dict(variant=[v] * 4, wiring=list(G.T_WIRINGS), flagged=[False] * 4))
    rd = _fake_readings({ref: 0.050, v: 0.046, ctl: 0.049}, T={ref: 0.0538, v: 0.0490, ctl: 0.0530})
    rows = [G.E92G.gate_direction92(rd, runs, h, v, ref, [ctl]) for h in G.HALVES]
    g = G.E92G.gate92(rows)
    assert g['i'] and g['ii'] and g['iii'] and g['iv'] and g['control'] and g['passes']
    rd5 = _fake_readings({ref: 0.050, v: 0.046, ctl: 0.049}, T={ref: 0.0538, v: 0.0490, ctl: 0.0510})
    g5 = G.E92G.gate92([G.E92G.gate_direction92(rd5, runs, h, v, ref, [ctl]) for h in G.HALVES])
    assert not g5['control'] and g5['passes_i_to_iv'] and not g5['passes']   # the control gains half as much


def test_the_three_kinds_of_equal():
    ref = G.ref_variant()
    rd = _fake_readings({ref: 0.050, 'B_sub_k2_d0.25': 0.052, 'B_sub_k2_d0': 0.0545, 'B_sub_k2_d2': 0.0505,
                         'B_sub_k1_d2': 0.0495, 'B_sub_k1_d4': 0.049, 'B_gain_k2_d0.05': 0.051, 'B_gain_k2_d0.15': 0.0501,
                         'B_sub_k2_d4_ctrl': 0.05, G.ref_k1_variant(): 0.0505},
                        raw={'B_gain_k2_d0.05': 0.090}, mean={'B_sub_k1_d4': 0.150},
                        per_target={'B_sub_k2_d0.25': {'halfA': 0.0525}})
    eq = {e['variant']: e for e in G.equal_rule94(rd)}
    assert set(eq) == {'B_sub_k2_d0.25', 'B_sub_k2_d0', 'B_sub_k2_d2', 'B_sub_k1_d2', 'B_sub_k1_d4', 'B_gain_k2_d0.05',
                       'B_gain_k2_d0.15'}                              # the B settings; no references or controls
    assert eq['B_sub_k2_d0.25']['equal'] and eq['B_sub_k2_d0.25']['kind'] == 'lower_share'   # 0.0025 off on half A
    assert not eq['B_sub_k2_d0']['equal'] and eq['B_sub_k2_d0']['kind'] == 'lower_share'      # 0.0045 worse
    assert not eq['B_sub_k2_d2']['equal'] and eq['B_sub_k2_d2']['kind'] is None               # above dose 1: no kind
    assert eq['B_sub_k1_d2']['equal'] and eq['B_sub_k1_d2']['better'] and eq['B_sub_k1_d2']['kind'] == 'weaker_weighting'
    assert not eq['B_sub_k1_d4']['equal'] and not eq['B_sub_k1_d4']['better']                # the mean APL effect 0.027 lower
    assert not eq['B_gain_k2_d0.05']['equal'] and eq['B_gain_k2_d0.05']['kind'] == 'gain'     # raw-count loss 0.012 above
    assert eq['B_gain_k2_d0.15']['equal'] and not eq['B_gain_k2_d0.15']['better']
    # flagged or at a bound: not an equal
    rd2 = _fake_readings({ref: 0.050, 'B_gain_k2_d0.15': 0.0501}, flagged=('B_gain_k2_d0.15',))
    assert not G.equal_rule94(rd2)[0]['equal']


def test_equal_window_is_three_thousandths():
    ref = G.ref_variant()
    rd = _fake_readings({ref: 0.050, 'B_sub_k2_d0.5': 0.050}, per_target={'B_sub_k2_d0.5': {'halfA': 0.0531}})
    assert not G.equal_rule94(rd)[0]['equal']                         # 0.0031 off against half A
    rd = _fake_readings({ref: 0.050, 'B_sub_k2_d0.5': 0.050}, per_target={'B_sub_k2_d0.5': {'halfB': 0.0529}})
    assert G.equal_rule94(rd)[0]['equal']


def test_the_anatomical_dose():
    rows = []
    for d, sh in ((0.1, 0.06), (0.25, 0.12), (0.5, 0.30)):
        for i, w in enumerate(G.S_WIRINGS):
            rows.append(dict(variant=G.variant_name('sub', 2.0, d), wiring=w, bouton_share=sh + 0.01 * (i - 1.5)))
    runs = pd.DataFrame(rows)
    if not _shared():
        pytest.skip("needs the fits' APL->bouton multiples (the shared directory)")
    a = G.anatomical_dose(runs)
    assert a['dose'] == 0.25 and np.isclose(a['share'], 0.12) and set(a['shares']) == {'0.1', '0.25', '0.5'}
    assert np.isclose(a['holds']['real'], 0.25 * G.aplpn_fit_multiple('real'))
    rows2 = [dict(r, bouton_share=0.104) if r['variant'] == 'B_sub_k2_d0.1' else r for r in rows]
    assert G.anatomical_dose(pd.DataFrame(rows2))['dose'] == 0.1
    with pytest.raises(RuntimeError):
        G.anatomical_dose(runs[runs.wiring != 'M1'])


def test_d1_rule():
    base = dict(pattern_all=0.0510, halfA=0.0570, halfB=0.0560, n_flagged=0, raw_loss=0.080, raw_loss_unchanged=0.078,
                mean_effect=0.180)
    assert G.d1_rule(base) == (True, False)
    assert G.d1_rule(dict(base, pattern_all=0.049, halfA=0.056, halfB=0.055)) == (True, True)
    assert G.d1_rule(dict(base, pattern_all=0.0535)) == (False, False)
    assert G.d1_rule(dict(base, n_flagged=1)) == (False, False)
    assert G.d1_rule(dict(base, raw_loss=0.084)) == (False, False)
    assert G.d1_rule(dict(base, mean_effect=0.165)) == (False, False)


def test_d2_plan():
    no = dict(confirm=False)
    p = G.d2_plan({'A': dict(confirm=True, screen='B1', dose=2.0, chosen='B_sub_k2_d2'), 'B': no}, [])
    assert p['due'] and p['arms'] == ['B1']
    assert not G.d2_plan({'A': dict(confirm=True, screen='B1', dose=0.25, chosen='B_sub_k2_d0.25')}, [])['due']
    p = G.d2_plan({'A': no}, [dict(equal=True, screen='B2', kind='weaker_weighting', variant='B_sub_k1_d2'),
                              dict(equal=True, screen='B1', kind='lower_share', variant='B_sub_k2_d0.5')])
    assert p['arms'] == ['B2']
    assert G.d2_plan({'A': dict(confirm=True, screen='B3', dose=0.05, chosen='B_gain_k2_d0.05')}, [])['arms'] == ['B3']
    with pytest.raises(NotImplementedError):
        G.d2_command('real', 'B3')


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the fits of C2's family)")
def test_d1_commands_and_tags():
    J = pytest.importorskip('experiments.E10_joint_fit')
    c = G.d1_command('real', 0.25)
    assert c[c.index('-k') + 1] == '2' and '--aplpn-fixed' not in c and '--aplpn-max' not in c
    assert np.isclose(float(c[c.index('--aplpn-hold') + 1]), 0.25 * G.aplpn_fit_multiple('real'))
    assert c[c.index('--init-params') + 1] == f'results/{G.EXPERIMENT}/d_start_params.json:real_d1'
    assert c[c.index('--experiment') + 1] == G.EXPERIMENT and c[c.index('--maxfev') + 1] == '150'
    for a, v in (('--rate-target', 'ctrl_frac'), ('--loss', 'apl'), ('--apl-readout', 'both'), ('--tag', '_d1')):
        assert c[c.index(a) + 1] == v
    assert G.d1_tag('real') == 'rate0.068_aplboth_aplpnhold_k2_d1'
    cp = G.d1_command('P1', 0.25)
    assert cp[cp.index('--wiring-seed') + 1] == '1'
    assert G.d1_tag('P1') == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_aplpnhold_k2_d1'
    assert G.d1_tag('M1').startswith('wiring-curveball_within_type_matched-seed1_')
    mults = [G.aplpn_fit_multiple(w) for w in G.ALL_WIRINGS]
    assert 2.4 < min(mults) and max(mults) < 2.81                                   # the entry's 2.4 to 2.8
    d2 = G.d2_command('P2', 'B2')
    assert d2[d2.index('-k') + 1] == '1' and d2[d2.index('--aplpn-max') + 1] == '12' and '--init-from-exact' not in d2
    assert '--init-from-exact' in G.d2_command('real', 'B1')
    assert J.aplpn_hold_tag(None) == ''


# ---------------------------------------------------------------------------
# E10_joint_fit.py's --aplpn-hold

def test_e10_aplpn_hold_helpers():
    J = pytest.importorskip('experiments.E10_joint_fit')
    assert J.aplpn_hold_kws(None) == {} and J.aplpn_hold_tag(None) == ''
    assert J.aplpn_hold_kws(0.7) == {'aplpn_hold': 0.7} and J.aplpn_hold_tag(0.7) == 'aplpnhold_'
    for bad in (-1.0, float('nan'), float('inf')):
        with pytest.raises(ValueError):
            J.aplpn_hold_kws(bad)
    assert inspect.signature(J.fit_one).parameters['aplpn_hold'].default is None
    src = inspect.getsource(J.main)
    assert "'--aplpn-hold'" in src and 'aplpn_hold_tag(args.aplpn_hold)' in src
    assert "ap.error('--aplpn-hold with --aplpn-fixed or --aplpn-max" in src
    assert "**({} if args.aplpn_hold is None else dict(aplpn_hold=args.aplpn_hold))" in src


def _fit_one_call(monkeypatch, **kw):
    """fit_one's call of optimize_mb_model_params (captured, not run) on made-up inputs."""
    J = pytest.importorskip('experiments.E10_joint_fit')
    import pyolfsysm.mb_model as M
    got = {}

    def fake(orn, **k):
        got.update(k)
        return dict(history=pd.DataFrame(), params=dict(k['initial_params']))
    monkeypatch.setattr(M, 'optimize_mb_model_params', fake)
    idx = pd.MultiIndex.from_arrays([[1, 1, 2], [0, 1, 0]], names=[H.KC_ID_LEVEL, 'claw_id'])
    bidx = pd.Index([0, 1], name='bouton_id')
    wnorm = dict(wAPLKC=pd.Series([1.0, 2.0, 1.0], index=idx), wKCAPL=pd.Series([1.0, 1.0, 1.0], index=idx),
                 wAPLPN=pd.Series([1.0, 3.0], index=bidx), wPNAPL=pd.Series([1.0, 1.0], index=bidx))
    breadth = pd.Series([3.0, 5.0], index=pd.Index([1, 2], name=H.KC_ID_LEVEL))
    initial = dict(fixed_thr=100.0, wAPLKC_scale=0.01, wKCAPL_scale=0.005, wAPLPN_scale=0.13, wPNAPL_scale=0.02)
    res, init, bounds = J.fit_one(2.0, pd.DataFrame(), initial, wnorm, breadth, lambda *a, **k: 0.0, kckc0=60.0, maxfev=5,
                                  method='Nelder-Mead', thr0_factor=0.7, verbose=False, **kw)
    return got, init, bounds, initial


def test_e10_fit_one_with_and_without_the_hold(monkeypatch):
    got0, init0, b0, initial = _fit_one_call(monkeypatch)
    assert got0['params_to_optimize'] == ['fixed_thr', 'wAPLKC_scale', 'wAPLPN_scale', 'wKCKC_scale']
    assert got0['tie_apl_weights'] is True and init0['wAPLPN_scale'] == initial['wAPLPN_scale'] and 'wAPLPN_scale' in b0
    # the hold: the APL->bouton scale fixed at F x the base's, the bouton->APL scale tied as without it
    got1, init1, b1, _ = _fit_one_call(monkeypatch, aplpn_hold=0.7)
    assert got1['params_to_optimize'] == ['fixed_thr', 'wAPLKC_scale', 'wKCKC_scale'] and 'wAPLPN_scale' not in b1
    assert tuple(got1['tie_apl_weights']) == ('wKCAPL_scale', 'wPNAPL_scale')
    assert init1['wAPLPN_scale'] == initial['wAPLPN_scale'] * 0.7 and init1['wPNAPL_scale'] == init0['wPNAPL_scale']
    assert {k: v for k, v in init1.items() if k != 'wAPLPN_scale'} == {k: v for k, v in init0.items() if k != 'wAPLPN_scale'}
    # the call is otherwise the same
    assert set(got0) == set(got1)
    assert all(got0[k] == got1[k] for k in ('condition_param_overrides', 'maxfev', 'method', 'kc_kc_inh'))
    # E31's --aplpn-fixed keeps its own behaviour (the bouton->APL scale fixed at F x the base's, only KC->APL tied)
    got2, init2, _, _ = _fit_one_call(monkeypatch, aplpn_fixed=0.7)
    assert tuple(got2['tie_apl_weights']) == ('wKCAPL_scale',) and init2['wPNAPL_scale'] == initial['wPNAPL_scale'] * 0.7
    with pytest.raises(ValueError):
        _fit_one_call(monkeypatch, aplpn_hold=0.7, aplpn_max=3.0)


# ---------------------------------------------------------------------------
# dose 0 on the real wiring (slow)

SLOW = pytest.mark.skipif(QUICK or not _shared(), reason="simulations of C2's family on the real wiring (minutes)")


@SLOW
def test_dose_0_on_the_real_wiring():
    E83 = pytest.importorskip('experiments.E83_cross_wiring_fit')
    m = G.WiringModel94('real', check=False)
    p0 = dict(m.p, wAPLPN_scale=0.0)
    sim0, _ = E83.build_simulator('real', params=p0)
    a = sim0(p0)
    m.dose = 0.0
    b = m.run_on(1.0)
    m.dose = 1.0
    assert np.array_equal(a[1].to_numpy(), b[1].to_numpy()) and np.array_equal(a[0].to_numpy(), b[0].to_numpy())
    on = m.run_on(1.0)
    assert np.array_equal(on[1].to_numpy(), m.saved['on'][1].to_numpy())          # dose 1: the fit again
