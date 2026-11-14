#!/usr/bin/env python3
"""Tests for E95 (APL's inhibition as a static threshold raise per KC; experiments/E95_static_apl.py): the variant names,
the arms' sizes and the kinds of equal; the losses per KC on a toy with two claws; the static matrices (the exponents
keep the mean, the permutations keep the values within a claw count, the per-KC forms the same for every odor); the
KC-by-odor setter (the odor order, refusals in Python); the inhibition absent from a run without APL and reset after a
run, also one that fails; the component's size in units of the bouton path and the weights with a component of 0 the
fit's; the realism readings; the two-way split; the choice, the gate, the three kinds of equal and part F's condition on
made-up rows; E10_joint_fit.py's --aplkc-flat."""

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
