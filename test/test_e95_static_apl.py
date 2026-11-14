#!/usr/bin/env python3
"""Tests for E95 (APL's inhibition as a static threshold raise per KC; experiments/E95_static_apl.py): the variant names,
the arms' sizes and the kinds of equal; the losses per KC on a toy with two claws; the static matrices (the exponents
keep the mean, the permutations keep the values within a claw count, the per-KC forms the same for every odor); the
KC-by-odor setter (the odor order, refusals in Python); the inhibition absent from a run without APL and reset after a
run, also one that fails; the component's size in units of the bouton path and the weights with a component of 0 the
fit's; the realism readings; the two-way split; the choice, the gate, the three kinds of equal and part F's condition on
made-up rows; E10_joint_fit.py's --aplkc-flat."""

import json
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
G = pytest.importorskip('experiments.E95_static_apl')

ODORS = list(H.APL16_ODORS)
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
    per = {p: sum(len(v) for a, v in G.ARMS.items() if a[0] == p) * n for p in 'BCD'}
    assert per == {'B': 36, 'C': 24, 'D': 28}                                # the entry's 88 points
    for arm, vs in G.ARMS.items():
        for v in vs:
            p = G.parse_variant(v)
            assert p['arm'] == arm and p['screen'] == arm and not p['control'] and p['part'] == arm[0]
    assert G.APL_BRACKET == (0.01, 10.0) and G.STATIC_BRACKET == (0.05, 20.0)
    assert set(G.CHOICE_ARMS) == set(G.ARMS) - {'B4', 'C5'}


def test_variant_fields():
    p = G.parse_variant('B1_bk')
    assert (p['mode'], p['static'], p['dose'], p['equal_kind']) == ('bouton', 'b', 0.0, 'bouton_static')
    p = G.parse_variant('B3_nk0.5')
    assert p['static'] == 'nke_b' and p['exponent'] == 0.5 and p['equal_kind'] is None
    p = G.parse_variant('B4_bk_permclaw1')
    assert p['permute'] == 'claw' and p['seed'] == 1 and p['equal_kind'] is None and p['static'] == 'b'
    p = G.parse_variant('C2_bcko')
    assert (p['mode'], p['static'], p['equal_kind']) == ('whole', 'bc_o', 'whole_static')
    p = G.parse_variant('C6_ck')
    assert (p['mode'], p['static'], p['dose'], p['equal_kind']) == ('claw', 'c', 1.0, 'whole_static')
    p = G.parse_variant('C3_nk')
    assert p['static'] == 'nk_bc' and p['equal_kind'] is None
    p = G.parse_variant('D2_connectome_d0.5')
    assert (p['mode'], p['component'], p['size'], p['dose'], p['equal_kind']) == ('kcs', 'connectome', 0.5, 0.1, 'at_kcs')
    p = G.parse_variant('D3_uniform_d1_dose0')
    assert (p['component'], p['size'], p['dose'], p['equal_kind']) == ('uniform', 1.0, 0.0, 'at_kcs')
    p = G.parse_variant('D0')
    assert p['control'] and p['dose'] == 0.1 and p['size'] == 0.0 and p['screen'] == 'control'
    c = G.parse_variant('C1_bck_ctrl')
    assert c['control'] and c['control_of'] == 'C1_bck' and c['permute'] == 'all' and c['seed'] == G.PERM_SEED
    assert c['static'] == 'bc' and c['equal_kind'] is None and c['screen'] == 'control'
    assert G.control_of('B5_bko') == 'B5_bko_ctrl' and G.control_of('D1_uniform_d2') == 'D0'
    for bad in ('B4_bk_permall1', 'B1_bk_ctrl_ctrl', 'D0'):
        with pytest.raises(ValueError):
            G.control_of(bad)
    for bad in ('B3_nk1', 'B3_nk0.50', 'C5_bck_permclaw1', 'D1_connectome_d1', 'D2_uniform_d1', 'D1_uniform_d1.0',
                'D1_uniform_d1_ctrl', 'B4_bk_permall1_ctrl', 'B_sub_k2_d1', 'E1_bk'):
        with pytest.raises(ValueError):
            G.parse_variant(bad)
    r = G.parse_variant('ref')
    assert r['screen'] == 'ref' and r['mode'] == 'ref'


# ---------------------------------------------------------------------------
# the losses per KC on a toy

def _toy_wiring():
    """Three KCs: KC 11 with two claws (boutons DA1 and VA2), KC 12 and KC 13 with one claw each; the outputs' KC order
    13, 11, 12."""
    idx = pd.MultiIndex.from_arrays([[11, 11, 12, 13], [0, 1, 0, 0]], names=[H.KC_ID_LEVEL, 'claw_id'])
    cols = pd.MultiIndex.from_arrays([['DA1', 'VA2'], [0, 1]], names=['glomerulus', 'bouton_id'])
    W = pd.DataFrame([[1.0, 0.0], [0.0, 2.0], [0.5, 0.0], [0.0, 1.0]], index=idx, columns=cols)
    kc_index = pd.MultiIndex.from_arrays([[13, 11, 12], ['g', 'ab', 'ab']], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    return W, kc_index


def _toy_losses(monkeypatch):
    W, kc_index = _toy_wiring()
    odors = ['2h', 'IaA']
    b_off = pd.DataFrame([[10.0, 4.0], [5.0, 2.0]], index=W.columns, columns=odors)
    b_on = pd.DataFrame([[8.0, 4.0], [5.0, 1.0]], index=W.columns, columns=odors)
    # claw states with APL on: KC 11's claws 6 / 9 and 3 / 1; KC 12's 3 / 1; KC 13's 4 / 1
    c_on = pd.DataFrame([[6.0, 3.0], [9.0, 1.0], [3.0, 1.0], [4.0, 1.0]], index=W.index, columns=odors)
    c_off = pd.DataFrame(W.to_numpy() @ b_off.to_numpy(), index=W.index, columns=odors)
    st = {'on': dict(bouton_win=b_on, claw_win=c_on, apl=pd.DataFrame(dict(win=[2.0, 4.0], pre=[1.0, 1.0]), index=odors)),
          'off': dict(bouton_win=b_off, claw_win=c_off, apl=pd.DataFrame(dict(win=[0.0, 0.0], pre=[0.0, 0.0]), index=odors))}
    saved = pd.DataFrame(np.zeros((3, 2)), index=kc_index, columns=[f'{o} @ -3' for o in odors])
    monkeypatch.setattr(G.E94, 'wiring_W', lambda w: (W, np.ones(2)))
    monkeypatch.setattr(G.E94, 'fit_stage', lambda w, c: st[c])
    monkeypatch.setattr(G.E94, 'check_stage_shapes', lambda *a, **k: None)
    monkeypatch.setattr(G.E83, 'saved_fit', lambda w: {'saved': {'on': (saved, saved), 'off': (saved, saved)}})
    G._CACHE.pop(('losses', 'toy'), None)
    return G.kc_losses('toy'), W, kc_index


def test_the_losses_per_kc_on_a_toy_with_two_claws(monkeypatch):
    L, W, kc_index = _toy_losses(monkeypatch)
    assert L['Lb'].index.equals(kc_index) and L['odors'] == ['2h', 'IaA']
    # KC 11: drive off 10 + 10 = 20 / 4 + 4 = 8, from the boutons with APL on 8 + 10 = 18 / 4 + 2 = 6, states 15 / 4
    # KC 12: off 5 / 2, on 4 / 2, state 3 / 1; KC 13: off 5 / 2, on 5 / 1, state 4 / 1
    assert np.allclose(L['Lb'].to_numpy(), [[0.0, 1.0], [2.0, 2.0], [1.0, 0.0]])
    assert np.allclose(L['Lc'].to_numpy(), [[1.0, 0.0], [3.0, 2.0], [1.0, 1.0]])
    assert np.allclose(L['Doff'].to_numpy(), [[5.0, 2.0], [20.0, 8.0], [5.0, 2.0]])
    assert np.allclose(L['b'], [0.5, 2.0, 0.5]) and np.allclose(L['c'], [0.5, 2.5, 1.0])
    assert list(L['n']) == [1.0, 2.0, 1.0]
    # the mean loss a claw at the boutons: claws' losses (2, 0, 1, 0) and (0, 2, 0, 1) over 8 entries
    assert np.isclose(L['bclaw_mean'], 6 / 8) and np.isclose(L['apl_mean'], 3.0)


def test_the_static_matrices(monkeypatch):
    L, W, kc_index = _toy_losses(monkeypatch)
    n, b, c = np.array([1.0, 2.0, 1.0]), np.array([0.5, 2.0, 0.5]), np.array([0.5, 2.5, 1.0])
    rec = {}
    M = G.static_matrix(G.parse_variant('B1_bk'), L, rec)
    assert M.index.equals(kc_index) and list(M.columns) == ['2h', 'IaA']
    assert np.allclose(M.to_numpy(), np.repeat(b[:, None], 2, axis=1)) and rec['n_clipped'] == 0
    beta = (n @ b) / (n @ n)
    rec = {}
    assert np.allclose(G.static_matrix(G.parse_variant('B2_nk'), L, rec).iloc[:, 1], beta * n) and np.isclose(rec['slope'], beta)
    for e in G.B3_EXPONENTS:                                                     # the exponents keep b's mean
        Me = G.static_matrix(G.parse_variant(f'B3_nk{e:g}'), L)
        assert np.isclose(Me.iloc[:, 0].mean(), b.mean()) and np.allclose(Me.iloc[:, 0] / Me.iloc[:, 0].mean(), n ** e / (n ** e).mean())
    assert np.allclose(G.static_matrix(G.parse_variant('B3_nk0'), L).to_numpy(), b.mean())   # the same raise for every KC
    assert np.allclose(G.static_matrix(G.parse_variant('B5_bko'), L).to_numpy(), L['Lb'].to_numpy())
    assert np.allclose(G.static_matrix(G.parse_variant('C1_bck'), L).iloc[:, 0], b + c)
    assert np.allclose(G.static_matrix(G.parse_variant('C2_bcko'), L).to_numpy(), (L['Lb'] + L['Lc']).to_numpy())
    assert np.allclose(G.static_matrix(G.parse_variant('C3_nk'), L).iloc[:, 0], (n @ (b + c)) / (n @ n) * n)
    assert np.allclose(G.static_matrix(G.parse_variant('C4_nkck'), L).iloc[:, 0], beta * n + c)
    assert np.allclose(G.static_matrix(G.parse_variant('C6_ck'), L).iloc[:, 0], c)
    # a negative loss is clipped and recorded
    L2 = dict(L, Lb=L['Lb'] - 1.5)
    rec = {}
    M2 = G.static_matrix(G.parse_variant('B5_bko'), L2, rec)
    assert (M2.to_numpy() >= 0).all() and rec['n_clipped'] == 4 and rec['min_before_clip'] == -1.5


def test_the_permutations_keep_the_values_within_a_claw_count():
    rng = np.random.default_rng(0)
    n = rng.integers(1, 12, 500).astype(float)
    v = rng.gamma(2.0, 3.0, 500)
    p = G.permute_values(v, 1, n)
    for k in np.unique(n):
        assert np.array_equal(np.sort(p[n == k]), np.sort(v[n == k]))
    assert not np.array_equal(p, v)
    pa = G.permute_values(v, 1)
    assert np.array_equal(np.sort(pa), np.sort(v)) and not np.array_equal(pa, v)
    assert np.array_equal(G.permute_values(v, 1), pa) and not np.array_equal(G.permute_values(v, 2), pa)  # one seed
    A = np.repeat(v[:, None], 3, axis=1)
    PA = G.permute_values(A, 1, n)                                        # rows of a KCs x odors matrix
    assert np.array_equal(PA[:, 0], p) and np.array_equal(PA[:, 0], PA[:, 2])
    with pytest.raises(ValueError):
        G.permute_values(v, 1, n[:-1])


def test_the_permuted_settings_keep_their_values(monkeypatch):
    L, _, _ = _toy_losses(monkeypatch)
    for v in ('B4_bk_permall1', 'B4_bk_permclaw1', 'C5_bck_permall1', 'B1_bk_ctrl', 'C2_bcko_ctrl'):
        pv = G.parse_variant(v)
        base = G.static_matrix(dict(pv, permute=None), L).to_numpy()
        M = G.static_matrix(pv, L).to_numpy()
        assert np.array_equal(np.sort(M, axis=0), np.sort(base, axis=0))
    # within a claw count: KCs 13 and 12 (one claw each) may swap, KC 11 (two claws) keeps its own
    Mc = G.static_matrix(G.parse_variant('B4_bk_permclaw1'), L).to_numpy()
    assert Mc[1, 0] == 2.0


# ---------------------------------------------------------------------------
# the setter, the runs

def _mp(n_kc: int = 3):
    olfsysm = pytest.importorskip('olfsysm')
    mp = olfsysm.ModelParams()
    mp.kc.N = n_kc
    mp.kc.ext_inh_scale = 1.0
    return mp


def test_the_setter_maps_the_odor_order():
    mp = _mp(3)
    rv = SimpleNamespace(kc=SimpleNamespace(ext_inh=np.zeros((0, 0))))
    M = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 5.0], [3.0, 0.0, 0.0]])
    order = np.array([2, 0, 1])                                  # output column j is olfsysm's odor order[j]
    G.set_ext_inh_ko(mp, rv, M, order)
    A = np.asarray(rv.kc.ext_inh)
    assert np.array_equal(A[:, 2], M[:, 0]) and np.array_equal(A[:, 0], M[:, 1]) and np.array_equal(A[:, 1], M[:, 2])
    # E91's per-KC setter and the KC-by-odor one with the same value for every odor
    v = np.array([0.5, 1.5, 2.5])
    G.E91.set_ext_inh(mp, rv, v, 3)
    a91 = np.array(rv.kc.ext_inh)
    G.set_ext_inh_ko(mp, rv, np.repeat(v[:, None], 3, axis=1), order)
    assert np.array_equal(a91, np.asarray(rv.kc.ext_inh))
    G.set_ext_inh_ko(mp, rv, None)
    assert np.asarray(rv.kc.ext_inh).size == 0
    # refusals in Python (olfsysm's check would abort the process inside its parallel region)
    for bad in (np.full((3, 3), -1.0), np.full((3, 3), np.nan), np.ones((2, 3)), np.ones(3)):
        with pytest.raises(ValueError):
            G.set_ext_inh_ko(mp, rv, bad, order)
    with pytest.raises(ValueError):
        G.set_ext_inh_ko(mp, rv, np.ones((3, 3)), np.array([0, 0, 1]))
    rec = G.set_ext_inh_ko(mp, rv, np.array([[1.0, -1e-12, 0.0]] * 3), order)    # rounding below 0: clipped
    assert rec['n_clipped'] == 3 and (np.asarray(rv.kc.ext_inh) >= 0).all()
    mp.kc.ext_inh_scale = 2.0
    with pytest.raises(RuntimeError):
        G.set_ext_inh_ko(mp, rv, np.ones((3, 3)), order)


P_FIT = dict(fixed_thr=100.0, wAPLKC_scale=0.12, wKCAPL_scale=7e-5, wAPLPN_scale=0.376, wPNAPL_scale=9e-4, wKCKC_scale=26.0)


def _fake_model(fail: bool = False):
    """A StaticModel without a simulation: its simulator records the parameters and the static inhibition of every run."""
    seen = []
    idx = pd.MultiIndex.from_arrays([np.arange(3), ['g'] * 3], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    counts = pd.DataFrame(np.array([[3.0, 0.0], [0.0, 2.0], [0.0, 0.0]]), index=idx, columns=['2h @ -3', 'IaA @ -3'])
    m = object.__new__(G.StaticModel)
    m.label, m.target_on, m.size, m.p, m.dose, m.ext = 'fake', 2, 6, dict(P_FIT), 1.0, None
    m.order = np.array([1, 0])
    m.saved = {'off': ((counts > 1.5).astype(float), counts), 'on': ((counts > 1.5).astype(float), counts)}
    m._kc0, m._kc = np.ones((3, 1)), np.ones((3, 1))
    m.wAPLKC = pd.Series(np.ones(3), index=pd.Index(range(3), name='claw'))
    m.W = None
    m.mp = _mp(3)
    m.rv = SimpleNamespace(kc=SimpleNamespace(ext_inh=np.zeros((0, 0))))

    def sim(params):
        seen.append(dict(params=dict(params), ext=np.array(m.rv.kc.ext_inh, copy=True)))
        if fail:
            raise RuntimeError('a C++ check')
        return (counts > 1.5).astype(float), counts
    m.sim = sim
    return m, seen


def test_the_inhibition_is_absent_without_apl_and_reset_after_a_run():
    m, seen = _fake_model()
    M = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    m.ext = M
    m.run_on(0.5)
    assert np.array_equal(seen[-1]['ext'][:, m.order], M)                   # set in the APL-on run, odor order mapped
    assert seen[-1]['params']['wAPLKC_scale'] == P_FIT['wAPLKC_scale'] * 0.5
    assert np.asarray(m.rv.kc.ext_inh).size == 0                             # reset after
    m.run_static(2.0)
    assert np.array_equal(seen[-1]['ext'][:, m.order], 2.0 * M)
    assert all(seen[-1]['params'][k] == 0.0 for k in G.APL_SCALES)          # every APL scale 0
    assert {k: v for k, v in seen[-1]['params'].items() if k not in G.APL_SCALES} == {k: v for k, v in P_FIT.items() if k not in G.APL_SCALES}
    m.run_off()
    assert seen[-1]['ext'].size == 0 and seen[-1]['params']['wAPLKC_scale'] == 0.0 and seen[-1]['params']['wAPLPN_scale'] == 0.0
    assert np.array_equal(m.ext, M)                                         # run_off leaves the setting alone
    m.ext = None
    m.run_on(1.0)
    assert seen[-1]['ext'].size == 0 and seen[-1]['params'] == m.params_on(1.0)
    with pytest.raises(RuntimeError):
        m.run_static(1.0)
    # reset also when a run fails
    m2, seen2 = _fake_model(fail=True)
    m2.ext = M
    with pytest.raises(RuntimeError):
        m2.run_on(1.0)
    assert np.asarray(m2.rv.kc.ext_inh).size == 0 and seen2[-1]['ext'].size == 6


def test_the_reading_sets_and_resets_everything(monkeypatch):
    E69 = pytest.importorskip('experiments.E69_kc_firing_ceiling')
    E90 = pytest.importorskip('experiments.E90_model_stages')
    calls = {}

    def fake_bisect(run, n, size, br, step, label):
        calls['bracket'] = br
        return 1.25, run(1.25), [{}], False
    monkeypatch.setattr(E69, 'bisect_factor', fake_bisect)
    monkeypatch.setattr(E90, 'read_stages', lambda mp, rv, W, order, odors, window: {'odors': list(odors)})
    monkeypatch.setattr(G.E92G, 'gain_stats', lambda rv, mp, odors, resp, form: dict(L=np.zeros((3, 2)), row={}, form=form))
    monkeypatch.setattr(G.E94, 'set_form94', lambda mp, form: None)
    m, seen = _fake_model()
    M = pd.DataFrame(np.ones((3, 2)), index=m.saved['on'][1].index, columns=['2h', 'IaA'])
    r = G.reading95(m, G.parse_variant('C1_bck'), M, None, 't')
    assert calls['bracket'] == G.STATIC_BRACKET and r['whole'] and np.allclose(r['ext'], 1.25)
    assert all(np.array_equal(s['ext'][:, m.order], np.full((3, 2), 1.25)) for s in seen)
    assert all(all(s['params'][k] == 0.0 for k in G.APL_SCALES) for s in seen)
    assert m.ext is None and np.asarray(m.rv.kc.ext_inh).size == 0 and m.dose == 1.0 and not m.mp.kc.save_inh_sims
    seen.clear()
    wk = pd.Series([0.0, 0.0, 0.0], index=m.wAPLKC.index)
    r = G.reading95(m, G.parse_variant('C6_ck'), M, wk, 't')
    assert calls['bracket'] == G.APL_BRACKET and not r['whole'] and np.allclose(r['ext'], 1.0)
    assert all(s['params']['wAPLKC_scale'] == P_FIT['wAPLKC_scale'] * 1.25 and s['params']['wAPLPN_scale'] == P_FIT['wAPLPN_scale']
               for s in seen)                                                # dose 1: the bouton path as fitted
    assert np.array_equal(m._kc, m._kc0)                                    # the weights reset
    seen.clear()
    G.reading95(m, G.parse_variant('B1_bk'), M, None, 't')
    assert all(s['params']['wAPLPN_scale'] == 0.0 for s in seen)              # dose 0: the bouton path removed
    G.reading95(m, G.parse_variant('D1_uniform_d1'), None, None, 't')
    assert seen[-1]['params']['wAPLPN_scale'] == P_FIT['wAPLPN_scale'] * 0.1 and seen[-1]['ext'].size == 0
    with pytest.raises(ValueError):
        G.reading95(m, G.parse_variant('C1_bck'), None, None, 't')
    with pytest.raises(ValueError):                                         # the odor order is checked
        G.reading95(m, G.parse_variant('B1_bk'), M[['IaA', '2h']], None, 't')
    assert m.ext is None and m.dose == 1.0


# ---------------------------------------------------------------------------
# the component

def test_the_component_size_in_units_of_the_bouton_path():
    # at the fit's scale 0.12 and APL's mean activity 20, a claw of mean connectome weight 0.4 takes 0.12 x w x 20; d = 1
    # takes the mean loss a claw at the boutons, 9.5
    u = G.component_unit(9.5, 0.12, 20.0, 0.4)
    assert np.isclose(0.12 * (u * 0.4) * 20.0, 9.5)
    idx = pd.Index(range(4), name='claw')
    w_fam = pd.Series([0.1, 0.0, 2.0, 0.5], index=idx)
    wn = pd.Series([0.2, 0.0, 0.6, 0.8], index=idx)
    assert G.d_weights(w_fam, wn, 'uniform', 0.0, u).equals(w_fam)          # a component of 0: the fit's weights
    assert G.d_weights(w_fam, wn, 'connectome', 0.0, u).equals(w_fam)
    for d in (0.5, 1.0, 2.0):
        du = G.d_weights(w_fam, wn, 'uniform', d, u) - w_fam
        dc = G.d_weights(w_fam, wn, 'connectome', d, u) - w_fam
        assert np.allclose(du, d * u * wn.mean())                           # the same on every claw
        assert np.allclose(dc, d * u * wn) and np.isclose(dc.mean(), du.mean())   # the same mean
        assert np.isclose(0.12 * du.mean() * 20.0, d * 9.5)                 # d units of the bouton path's take
    assert G.d_weights(w_fam, wn, 'connectome', 1.0, u)[1] == 0.0            # no connectome weight, no component
    with pytest.raises(ValueError):
        G.d_weights(w_fam, wn, 'flat', 1.0, u)
    with pytest.raises(ValueError):
        G.d_weights(w_fam, wn.iloc[::-1], 'uniform', 1.0, u)
    with pytest.raises(ValueError):
        G.component_unit(9.5, 0.0, 20.0)


# ---------------------------------------------------------------------------
# realism, the split

def test_static_realism():
    thr = np.array([100.0, 50.0, 200.0])
    n = np.array([2.0, 6.0, 12.0])
    D = np.array([[140.0, 80.0], [60.0, 40.0], [150.0, 100.0]])
    M = np.array([[30.0, 10.0], [70.0, 70.0], [20.0, 20.0]])
    on = np.array([[3.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    off = np.array([[3.0, 2.0], [2.0, 0.0], [0.0, 0.0]])
    r = G.static_realism(M, thr, n, D, on, off)
    assert np.isclose(r['static_over_thr_max'], 70 / 50) and np.isclose(r['static_over_thr_p50'], 0.2)
    assert np.isclose(r['static_over_thr_p50_1-4'], 0.2) and np.isclose(r['static_over_thr_p50_9+'], 0.1)
    assert np.isclose(r['static_over_max_drive'], 1 / 3)                    # KC 2: 70 > 60
    assert np.isclose(r['silent_by_mean_drive'], 2 / 3)                     # KC 2 (-10 < 50) and KC 3 (130 < 200)
    assert np.isclose(r['silent_share_on'], 2 / 3) and np.isclose(r['silent_share_off'], 1 / 3)
    r0 = G.static_realism(None, thr, n, D)
    assert np.isnan(r0['static_over_thr_p50']) and np.isnan(r0['static_over_max_drive'])


def test_d_realism_splits_the_claws_loss():
    # two claws of two KCs, one odor: drives off 10 and 10, from the boutons with APL on 8 and 10, states 6 and 9
    W = np.eye(2)
    r = G.d_realism(W, np.array([[10.0], [10.0]]), np.array([[8.0], [10.0]]), np.array([[6.0], [9.0]]),
                    w_fam=np.array([1.0, 0.0]), w_comp=np.array([1.0, 2.0]), kc_of=np.array([7, 8]))
    # loss 5: 2 at the boutons, 3 at the claws (2 at claw 1, half of it the component's; 1 at claw 2, all the component's)
    assert np.isclose(r['bouton_share_d'], 2 / 5) and np.isclose(r['comp_share'], 2 / 5) and np.isclose(r['weighted_share'], 1 / 5)
    assert np.isclose(r['bouton_share_d'] + r['comp_share'] + r['weighted_share'], 1.0)
    r2 = G.d_realism(W, np.array([[10.0], [10.0]]), np.array([[8.0], [10.0]]), np.array([[6.0], [9.0]]),
                     w_fam=np.array([1.0, 4.0]), w_comp=np.zeros(2), kc_of=np.array([7, 8]), responders=np.array([7, 8]))
    assert r2['comp_share'] == 0.0 and np.isclose(r2['wkc_fam_resp_p90_over_p10'], r2['wkc_all_resp_p90_over_p10'])


def test_two_way_split():
    kc = np.array([1.0, 2.0, 3.0])[:, None]
    od = np.array([0.0, 0.5, -0.5, 0.0])[None, :]
    s = G.two_way_split(10 + kc + od)
    assert np.isclose(s['share_kc'] + s['share_odor'], 1.0) and np.isclose(s['share_rest'], 0.0)
    assert np.isclose(s['share_kc'], (4 * 2) / (4 * 2 + 3 * 0.5))
    s2 = G.two_way_split(np.repeat(np.array([1.0, 5.0])[:, None], 4, axis=1))
    assert np.isclose(s2['share_kc'], 1.0) and s2['cv_odors_median_kc'] == 0.0


# ---------------------------------------------------------------------------
# the identities on the real wiring (slow)

SLOW = pytest.mark.skipif(QUICK or not _shared(), reason="simulations of C2's family on the real wiring (minutes)")


@SLOW
def test_a_zero_static_inhibition_is_the_fit_on_the_real_wiring():
    m = G.StaticModel('real', check=False)
    m.ext = np.zeros(m.saved['on'][1].shape)
    on = m.run_on(1.0)
    off = m.run_static(1.0)
    assert np.array_equal(on[1].to_numpy(), m.saved['on'][1].to_numpy())
    assert np.array_equal(off[1].to_numpy(), m.saved['off'][1].to_numpy())
    assert np.asarray(m.rv.kc.ext_inh).size == 0


# ---------------------------------------------------------------------------
# E10_joint_fit.py's --aplkc-flat

def test_e10_aplkc_flat_helpers():
    import inspect
    J = pytest.importorskip('experiments.E10_joint_fit')
    assert J.aplkc_flat_spec(None) is None and J.aplkc_flat_tag(None) == ''
    s = J.aplkc_flat_spec('uniform', 3.8)
    assert s == dict(kind='uniform', unit=3.8, initial=1.0, bounds=(0.0, 8.0)) and J.aplkc_flat_tag(s) == 'flatu_'
    assert J.aplkc_flat_tag(J.aplkc_flat_spec('connectome', 3.8, 0.5)) == 'flatc_'
    for bad in (dict(kind='flat', unit=1.0), dict(kind='uniform'), dict(kind='uniform', unit=-1.0),
                dict(kind='uniform', unit=1.0, init=9.0), dict(kind=None, unit=1.0)):
        with pytest.raises(ValueError):
            J.aplkc_flat_spec(bad.pop('kind'), **bad)
    assert inspect.signature(J.fit_one).parameters['aplkc_flat'].default is None
    src = inspect.getsource(J.main)
    assert "'--aplkc-flat'" in src and 'aplkc_flat_tag(flat_spec)' in src
    assert "**({} if flat_spec is None else dict(aplkc_flat=flat_spec))" in src


def test_e10_aplkc_flat_sets_the_simulators_weights():
    J = pytest.importorskip('experiments.E10_joint_fit')
    idx = pd.MultiIndex.from_arrays([[1, 1, 2, 3], [0, 1, 0, 0]], names=[H.KC_ID_LEVEL, 'claw_id'])
    w0 = pd.Series([1.0, 3.0, 0.0, 4.0], index=idx)                   # the connectome's normalized weights
    wk = pd.Series([0.5, 1.5, 0.0, 9.0], index=idx)                   # the family's breadth-weighted weights

    def make_sim(w):
        wAPLKC_unscaled = w

        def simulate(params):
            return wAPLKC_unscaled
        return simulate
    held = wk.to_numpy().reshape(-1, 1).copy()
    sim = make_sim(held)
    p = J.aplkc_flat_param(w0, wk, J.aplkc_flat_spec('uniform', 2.0, 0.0))
    assert p['initial'] == 0.0 and p['bounds'] == (0.0, 8.0)
    p['set'](sim, 0.0)
    assert np.array_equal(held.ravel(), wk.to_numpy())                # size 0: the family's weights exactly
    p['set'](sim, 1.5)
    assert np.allclose(held.ravel(), wk.to_numpy() + 1.5 * 2.0 * w0.mean())
    p['set'](sim, 0.0)
    assert np.array_equal(held.ravel(), wk.to_numpy())
    pc = J.aplkc_flat_param(w0, wk, J.aplkc_flat_spec('connectome', 2.0, 1.0))
    held2 = wk.to_numpy().reshape(-1, 1).copy()
    pc['set'](make_sim(held2), 0.5)
    assert np.allclose(held2.ravel(), wk.to_numpy() + 0.5 * 2.0 * w0.to_numpy())
    # a simulator that does not hold the family's weights is refused (first call)
    with pytest.raises(RuntimeError):
        J.aplkc_flat_param(w0, wk, J.aplkc_flat_spec('uniform', 2.0))['set'](make_sim(np.zeros((4, 1))), 1.0)
    # and the components have the connectome's mean
    assert np.isclose(J.aplkc_flat_component(w0, 'uniform').mean(), J.aplkc_flat_component(w0, 'connectome').mean())


def _fit_one_call(monkeypatch, **kw):
    """fit_one's call of optimize_mb_model_params (captured, not run) on made-up inputs (as E94's test)."""
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
    J.fit_one(2.0, pd.DataFrame(), initial, wnorm, breadth, lambda *a, **k: 0.0, kckc0=60.0, maxfev=5,
              method='Nelder-Mead', thr0_factor=0.7, verbose=False, **kw)
    return got


def test_e10_fit_one_with_and_without_the_component(monkeypatch):
    J = pytest.importorskip('experiments.E10_joint_fit')
    got0 = _fit_one_call(monkeypatch, aplpn_hold=0.27)
    assert 'extra_params' not in got0                                    # off: the call as before
    got1 = _fit_one_call(monkeypatch, aplpn_hold=0.27, aplkc_flat=J.aplkc_flat_spec('uniform', 3.8, 1.2))
    assert set(got1) == set(got0) | {'extra_params'} and set(got1['extra_params']) == {'aplkc_flat'}
    e = got1['extra_params']['aplkc_flat']
    assert e['initial'] == 1.2 and e['bounds'] == (0.0, 8.0) and callable(e['set'])
    assert all(got0[k] == got1[k] for k in ('params_to_optimize', 'tie_apl_weights', 'initial_params', 'maxfev'))
    assert got1['params_to_optimize'] == ['fixed_thr', 'wAPLKC_scale', 'wKCKC_scale']


# ---------------------------------------------------------------------------
# the choice, the gate, the three kinds of equal, part F's condition (made-up rows)

TARGETS = ('halfA', 'halfB', 'halfC', 'halfD', 'z8', 'z6')


def _fake_readings(patterns, flagged=(), at_bound=(), T=None, raw=None, mean=None, per_target=None):
    """Readings of made-up settings ({variant: pattern part on S, the same against every target unless `per_target`
    gives {variant: {target: pattern}}}), with T's ({variant: pattern}); raw-count losses and mean effects default to the
    reference's (as E94's test)."""
    rows = []
    for v, p in patterns.items():
        pv = G.parse_variant(v)
        for unit, pat in (('S', p),) + ((('T', T[v]),) if T and v in T else ()):
            base = dict(variant=v, arm=pv['arm'], screen=pv['screen'], mode=pv['mode'], control=pv['control'],
                        equal_kind=pv['equal_kind'], unit=unit, n_wirings=4, n_flagged=int(v in flagged),
                        n_at_bound=int(v in at_bound), mean_effect=(mean or {}).get(v, 0.177))
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
    rd = _fake_readings({ref: 0.050, 'B1_bk': 0.0505, 'C1_bck': 0.052})
    c = G.choose95(rd, 'A')
    assert c['chosen'] == ref and not c['confirm'] and c['gain'] == 0.0 and c['control'] == ''
    rd = _fake_readings({ref: 0.050, 'B2_nk': 0.049})
    c = G.choose95(rd, 'B')
    assert c['chosen'] == 'B2_nk' and np.isclose(c['gain'], 0.001) and not c['confirm']           # below 0.0015
    rd = _fake_readings({ref: 0.050, 'B4_bk_permall1': 0.030, 'C5_bck_permall1': 0.031, 'C1_bck_ctrl': 0.029, 'D0': 0.03,
                         'B3_nk1.5': 0.040, 'C2_bcko': 0.041, 'D1_uniform_d2': 0.0484, 'C3_nk': 0.047},
                        flagged=('B3_nk1.5',), at_bound=('C2_bcko',))
    c = G.choose95(rd, 'C')
    assert c['chosen'] == 'C3_nk' and c['confirm'] and np.isclose(c['gain'], 0.003) and c['control'] == 'C3_nk_ctrl'
    assert {'B3_nk1.5', 'C2_bcko'} <= set(json.loads(c['excluded']))       # outside the bands, at a bound
    assert np.isclose(c['chosen_shape'], 0.8 * 0.047) and np.isclose(c['ref_amplitude'], 0.6 * 0.050)
    rd = rd[rd.variant != 'C3_nk']
    c = G.choose95(rd, 'D')
    assert c['chosen'] == 'D1_uniform_d2' and c['confirm'] and c['control'] == 'D0'          # B4, C5, controls excepted
    assert G.choose95(_fake_readings({'B1_bk': 0.04}), 'A')['chosen'] is None                # no reference


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the flies' mean APL effect)")
def test_the_gate_on_four_directions():
    v, ref = 'B1_bk', G.REF
    ctl = G.control_of(v)
    runs = pd.DataFrame(dict(variant=[v] * 4, wiring=list(G.T_WIRINGS), flagged=[False] * 4))
    rd = _fake_readings({ref: 0.050, v: 0.046, ctl: 0.049}, T={ref: 0.0538, v: 0.0490, ctl: 0.0530})
    g = G.E92G.gate92([G.E92G.gate_direction92(rd, runs, h, v, ref, [ctl]) for h in G.HALVES])
    assert g['i'] and g['ii'] and g['iii'] and g['iv'] and g['control'] and g['passes']
    rd5 = _fake_readings({ref: 0.050, v: 0.046, ctl: 0.049}, T={ref: 0.0538, v: 0.0490, ctl: 0.0510})
    g5 = G.E92G.gate92([G.E92G.gate_direction92(rd5, runs, h, v, ref, [ctl]) for h in G.HALVES])
    assert not g5['control'] and g5['passes_i_to_iv'] and not g5['passes']   # the control gains half as much
    d = 'D1_uniform_d1'                                                     # part D's control: D0 on T
    runs_d = pd.DataFrame(dict(variant=[d] * 4, wiring=list(G.T_WIRINGS), flagged=[False] * 4))
    rdd = _fake_readings({ref: 0.050, d: 0.046, 'D0': 0.055}, T={ref: 0.0538, d: 0.0490, 'D0': 0.0560})
    gd = G.E92G.gate92([G.E92G.gate_direction92(rdd, runs_d, h, d, ref, [G.control_of(d)]) for h in G.HALVES])
    assert gd['passes']


def test_the_three_kinds_of_equal():
    ref = G.REF
    rd = _fake_readings({ref: 0.050, 'B1_bk': 0.052, 'B2_nk': 0.0545, 'B5_bko': 0.0495, 'B3_nk0.5': 0.0501,
                         'C1_bck': 0.0505, 'C2_bcko': 0.051, 'C3_nk': 0.050, 'C4_nkck': 0.049, 'C6_ck': 0.0502,
                         'D1_uniform_d1': 0.0535, 'D2_connectome_d1': 0.0501, 'D3_uniform_d1_dose0': 0.0502,
                         'B4_bk_permall1': 0.050, 'C1_bck_ctrl': 0.050, 'D0': 0.050},
                        raw={'C2_bcko': 0.090}, mean={'C4_nkck': 0.150, 'D2_connectome_d1': 0.160},
                        per_target={'B1_bk': {'halfA': 0.0525}}, at_bound=('C6_ck',), flagged=('D3_uniform_d1_dose0',))
    eq = {e['variant']: e for e in G.equal_rule95(rd)}
    assert set(eq) == {'B1_bk', 'B2_nk', 'B5_bko', 'B3_nk0.5', 'C1_bck', 'C2_bcko', 'C3_nk', 'C4_nkck', 'C6_ck',
                       'D1_uniform_d1', 'D2_connectome_d1', 'D3_uniform_d1_dose0'}   # no reference, permutation or control
    assert eq['B1_bk']['equal'] and eq['B1_bk']['kind'] == 'bouton_static'          # 0.0025 off on half A
    assert not eq['B2_nk']['equal'] and eq['B2_nk']['kind'] == 'bouton_static'      # 0.0045 worse
    assert eq['B5_bko']['equal'] and eq['B5_bko']['better']
    assert not eq['B3_nk0.5']['equal'] and eq['B3_nk0.5']['kind'] is None           # no kind: chosen, never an equal
    assert eq['C1_bck']['equal'] and eq['C1_bck']['kind'] == 'whole_static' and not eq['C1_bck']['better']
    assert not eq['C2_bcko']['equal']                                               # raw-count loss 0.012 above
    assert eq['C3_nk']['kind'] is None and not eq['C3_nk']['equal']
    assert not eq['C4_nkck']['equal'] and not eq['C4_nkck']['better']               # the mean APL effect 0.027 lower
    assert not eq['C6_ck']['equal']                                                 # at a bound
    assert not eq['D1_uniform_d1']['equal'] and eq['D1_uniform_d1']['kind'] == 'at_kcs'
    assert not eq['D2_connectome_d1']['equal']                                      # the mean APL effect lost
    assert not eq['D3_uniform_d1_dose0']['equal']                                   # outside the bands
    rd2 = _fake_readings({ref: 0.050, 'B5_bko': 0.050}, per_target={'B5_bko': {'halfB': 0.0531}})
    assert not G.equal_rule95(rd2)[0]['equal']                                      # 0.0031 off against half B


def test_part_f_condition():
    ref = G.REF
    rd = _fake_readings({ref: 0.050, 'D1_uniform_d0.5': 0.0585, 'D1_uniform_d1': 0.0590, 'D3_uniform_d1_dose0': 0.0579,
                         'D2_connectome_d1': 0.0501})
    c = G.f_condition(rd)
    assert c['runs'] and c['best'] == 'D3_uniform_d1_dose0' and np.isclose(c['above'], 0.0079)
    assert np.isclose(c['best_pattern'], 0.0579)
    assert c['connectome'] and c['d2_equals'] == ['D2_connectome_d1']
    rd = _fake_readings({ref: 0.050, 'D1_uniform_d1': 0.0585, 'D2_connectome_d1': 0.0560})
    c = G.f_condition(rd)
    assert not c['runs'] and not c['connectome']                                    # 0.0085 above; D2 not an equal
    assert not G.f_condition(_fake_readings({ref: 0.050}))['runs']                  # D1 / D3 not run


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the fits of C2's family)")
def test_part_f_commands_and_tags():
    J = pytest.importorskip('experiments.E10_joint_fit')
    c = G.f_command('real', 'uniform', unit=3.84)
    assert c[c.index('-k') + 1] == '2' and c[c.index('--aplkc-flat') + 1] == 'uniform'
    assert np.isclose(float(c[c.index('--aplpn-hold') + 1]), 0.1 * G.E94.aplpn_fit_multiple('real'))
    assert float(c[c.index('--aplkc-flat-unit') + 1]) == 3.84 and c[c.index('--aplkc-flat-max') + 1] == '8'
    assert c[c.index('--init-params') + 1] == f'results/{G.EXPERIMENT}/f_start_params.json:real_f1_uniform'
    for a, v in (('--rate-target', 'ctrl_frac'), ('--loss', 'apl'), ('--apl-readout', 'both'), ('--tag', '_f1'),
                 ('--maxfev', '200'), ('--experiment', G.EXPERIMENT)):
        assert c[c.index(a) + 1] == v
    assert G.f_tag('real') == 'rate0.068_aplboth_aplpnhold_flatu_k2_f1'
    assert G.f_tag('P1', 'connectome') == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_aplpnhold_flatc_k2_f1'
    cp = G.f_command('M1', 'uniform', unit=3.2)
    assert '--wiring-file' in cp and '--wiring-tag' in cp
    assert J.aplkc_flat_tag(dict(kind='uniform')) == 'flatu_'


# ---------------------------------------------------------------------------
# run_point's row and the readings on saved runs (no simulation; the shared directory)

def _saved_model(monkeypatch):
    """A stand-in for a StaticModel on the real wiring at its fit, from the saved fit and E94's wiring cache (no
    simulator); reading95 replaced by a reading of the saved APL-on run with the fit's stage matrices."""
    E83 = G.E83
    f = E83.saved_fit('real')
    W, _ = G.E94.wiring_W('real')
    m = object.__new__(G.StaticModel)
    m.label, m.folder, m.tag, m.p, m.saved, m.K = 'real', f['folder'], f['tag'], dict(f['p']), f['saved'], f['K']
    m.target_on, m.target_off, m.W = f['target_on'], f['target_off'], W
    m.wAPLKC = pd.Series(np.ones(W.shape[0]), index=W.index)
    kc, _ = G.load_kc_tables('real')
    m.rv = SimpleNamespace(kc=SimpleNamespace(thr=kc.thr.to_numpy()))

    def fake_reading(model, pv, M, wk, label=''):
        whole = pv['mode'] == 'whole'
        params = model.params_static() if whole else dict(model.p, wAPLPN_scale=model.p['wAPLPN_scale'] * pv['dose'])
        ext = None if M is None else M.to_numpy()
        return dict(on=f['saved']['on'], off=f['saved']['off'], factor=1.2, at_bound=False, n_sims=3, params=params,
                    seconds=0.0, history=[], stages=dict(G.E94.fit_stage('real', 'on')),
                    gain=dict(L=np.zeros(f['saved']['on'][1].shape), row={}), bracket=G.STATIC_BRACKET if whole else G.APL_BRACKET,
                    ext=ext, whole=whole)
    monkeypatch.setattr(G, 'reading95', fake_reading)
    monkeypatch.setattr(G, 'wiring_wnorm', lambda w: pd.Series(np.ones(W.shape[0]), index=W.index))
    monkeypatch.setattr(G, 'd_unit', lambda model, L=None: 3.0)
    return m


@pytest.mark.skipif(not _shared() or not (G.local_dir('W') / 'real_kc.parquet').exists(),
                    reason="needs the shared directory and part 0's per-KC tables")
def test_run_point_rows_on_saved_data(monkeypatch, tmp_path):
    m = _saved_model(monkeypatch)
    saved = {}
    monkeypatch.setattr(G, 'save_run95', lambda tag, o, sub='runs', **k: saved.update({tag: k}))
    for v in ('B1_bk', 'C1_bck', 'C6_ck', 'D1_uniform_d1', 'B4_bk_permclaw1', 'C2_bcko_ctrl', 'D0'):
        row, _ = G.run_point(m, v, part='test', cc={'queue_commit': 'x'})
        pv = G.parse_variant(v)
        assert row['variant'] == v and row['tag'] == f'{v}_real' and row['mode'] == pv['mode']
        assert row['factor_kind'] == ('static' if pv['mode'] == 'whole' else 'tied_apl')
        assert all(k in row for k in ('rate_on', 'realism_flag', 'max_count_on', 'K_run', 'bouton_share'))
        if pv['static']:
            assert 0 < row['static_over_thr_p50'] < 1 and np.isfinite(row['static_over_fixed_thr_p50'])
            assert saved[row['tag']]['static'] is not None
        else:
            assert np.isnan(row['static_over_thr_p50'])
        if pv['part'] == 'D':
            assert row['component_unit'] == 3.0 and 'comp_share' in row
            assert (row['comp_share'] > 0) == (pv['size'] > 0)
        assert np.isclose(row['aplpn_multiple'], G.E94.aplpn_multiple('real', row['dose']))
    rB1 = G.run_point(m, 'B1_bk', cc={})[0]
    L = G.kc_losses('real')
    assert np.isclose(rB1['static_mean'], L['b'].mean())


@pytest.mark.skipif(not _shared() or not (G.H.outputs_dir('E94_bouton_share') / 'B').exists(),
                    reason="needs E94's saved runs (viki's outputs) and the shared directory")
def test_the_readings_and_the_table_on_saved_runs(monkeypatch, tmp_path):
    """The scoring path on E94's saved runs standing in for this experiment's points on S (B1_bk: E94's dose 0, C1_bck:
    dose 0.5, D1_uniform_d1: dose 0.25)."""
    stand = {'B1_bk': 'B_sub_k2_d0', 'C1_bck': 'B_sub_k2_d0.5', 'D1_uniform_d1': 'B_sub_k2_d0.25'}
    rows = []
    for v, e in stand.items():
        for w in G.S_WIRINGS:
            rows.append(dict(tag=f'{e}_{w}', variant=v, wiring=w, source='E94', at_bound=False, realism_flag='',
                             factor=1.0, factor_kind='tied_apl'))
    pd.DataFrame(rows).to_csv(tmp_path / 'runs_B.csv', index=False)
    rd = G.readings95(tmp_path, readouts=('cal', 'raw'))
    s = rd[(rd.unit == 'S') & (rd.readout == 'cal')].set_index('variant')
    assert set(s.index) >= {G.REF, 'B1_bk', 'C1_bck', 'D1_uniform_d1', G.D0}
    assert np.isclose(s.loc[G.REF, 'z8_pattern'], G.UNCHANGED_S['all'], atol=5e-5)
    e94 = pd.read_csv(G.H.results_dir('E94_bouton_share') / 'readings.csv')
    e94 = e94[(e94.unit == 'S') & (e94.readout == 'cal')].set_index('variant')
    assert np.isclose(s.loc['B1_bk', 'z8_pattern'], e94.loc['B_sub_k2_d0', 'z8_pattern'], rtol=1e-9)
    assert np.isclose(s.loc[G.D0, 'z8_pattern'], e94.loc['B_sub_k2_d0.1', 'z8_pattern'], rtol=1e-9)
    G.c1_groups(tmp_path)
    v = G.verdict95(tmp_path, rd)
    assert set(v['choice']) == set(G.HALVES) and 'F_condition' in v
    text = G.comparison_table(tmp_path)
    assert '| C1_bck |' in text and 'Verdict' in text
