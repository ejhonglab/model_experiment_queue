#!/usr/bin/env python3
"""Tests for E91 (APL's inhibition on the KC's summed drive instead of claw by claw;
experiments/E91_summed_drive_inhibition.py) and E10_joint_fit.py's --claw-inh-pool: the variant names; the unabsorbed
inhibition, the claws carrying half of a drive and the leaky integration on hand-built claws; the claw statistics of a
made-up run; the control's matrix (one value per KC); E88's shared-responder classes as E88 counts them; the choice,
the gate and the equal rule on made-up tables; the E10 option's parsing; and part 0's identities (slow, skipped with
QUICK=1 and without the parameter: s = 0 through the new path reproduces the saved fit, at s = 1 the APL-off run is the
saved one, the parameter set on a simulator equals the simulator built with it, E89's C1 point F at exponent 1 on P1)."""

import inspect
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E91_summed_drive_inhibition')

ODORS = list(H.APL16_ODORS)
QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


def _has_pool() -> bool:
    try:
        import olfsysm
        return hasattr(olfsysm.ModelParams().kc, 'claw_inh_pool')
    except Exception:
        return False


# ---------------------------------------------------------------------------
# names

def test_variant_names_round_trip():
    cases = [(('B1', ('F', 2.0)), {}, 'B1_F_k2'), (('B1', ('F', 0.5)), dict(s=0.5), 'B1_F_k0.5_s0.5'),
             (('B1', ('F', 0.0)), dict(s=1.0), 'B1_none_s1'), (('B3', ('F', 2.0)), dict(net=True), 'B3_F_k2_net'),
             (('B1', ('F', 2.0)), dict(s=1.0, fixed=True), 'B1_F_k2_s1_fixed'),
             (('B2', ('F', 1.0)), dict(s=1.0, control=True), 'B2_F_k1_s1_ctrl')]
    for (arm, w), kw, name in cases:
        assert X.variant_name(arm, w, **kw) == name
        p = X.parse_variant(name)
        assert p['arm'] == arm and p['exponent'] == w[1] and p['net'] == kw.get('net', False)
        assert p['fixed'] == kw.get('fixed', False) and p['control'] == kw.get('control', False)
        assert p['s_of'] == kw.get('s', 0.0)
        assert p['s'] == (0.0 if kw.get('control') else kw.get('s', 0.0))     # a control runs at s = 0
    assert X.ref_variant('B1') == 'B1_F_k2' and X.ref_variant('B2') == 'B2_F_k2'
    assert X.parse_variant('B1_none')['measure'] == 'none'
    with pytest.raises(ValueError):
        X.parse_variant('B1_F_k2_bogus')


def test_settings_of_the_screen():
    st = X.b1_settings()
    assert len([s for s in st if not s.get('fixed')]) == 8 and len([s for s in st if s.get('fixed')]) == 1
    assert {X.variant_name('B1', s['w'], s['s']) for s in st if not s.get('fixed')} == {
        'B1_F_k2_s0.5', 'B1_F_k2_s1', 'B1_F_k1_s0.5', 'B1_F_k1_s1', 'B1_F_k0.5_s0.5', 'B1_F_k0.5_s1', 'B1_none_s0.5', 'B1_none_s1'}


# ---------------------------------------------------------------------------
# the claws' net drives

def test_unabsorbed_inhibition_of_a_made_up_kc():
    # KC 0: claws with net drives 5, -3, 0 (one claw takes more inhibition than its input); KC 1: 2 and -1, then 4, 1
    x = np.array([[5.0, 5.0], [-3.0, -3.0], [0.0, 0.0], [2.0, 4.0], [-1.0, 1.0]])
    kc = np.array([0, 0, 0, 1, 1])
    D, U = X.unabsorbed_of(x, kc, 2)
    assert np.allclose(D, [[5, 5], [2, 5]]) and np.allclose(U, [[3, 3], [1, 0]])
    # the summed drive at s: D - s U; at s = 1 the plain sum of the claws' net drives (no clipping)
    for s in (0.0, 0.5, 1.0):
        summed = D - s * U
        exp = np.array([[x[kc == i, t].clip(min=0).sum() + s * x[kc == i, t].clip(max=0).sum() for t in range(2)] for i in range(2)])
        assert np.allclose(summed, exp)
    assert np.allclose(D - U, [[2, 2], [1, 5]])


def test_claws_carrying_half_of_a_drive():
    assert X.n_claws_half([5, 3, 1, 1]) == 1
    assert X.n_claws_half([3, 3, 3, 1]) == 2
    assert X.n_claws_half([1, 1, 1, 1]) == 2
    assert X.n_claws_half([0, 0]) == 0
    assert X.n_claws_half([0.2, 4.0, 0.1]) == 1


def test_leaky_integration_is_olfsysms_euler_step():
    d = np.full((2, 50), 3.0)
    d[1] = 1.0
    V = X.leaky(d, 0.0005, 0.01, v0=np.zeros(2))
    a = 0.0005 / 0.01
    n = np.arange(1, 51)
    assert np.allclose(V[0], 3.0 * (1 - (1 - a) ** n)) and np.allclose(V[1], 1.0 * (1 - (1 - a) ** n))
    assert np.allclose(X.leaky(d, 0.0005, 0.01), d)          # from the drive's level it stays there


def _fake_run(monkeypatch, apl: bool = True):
    """A made-up run: 2 KCs (3 and 2 claws), 2 boutons, 2 odors, 5 time steps at t = -0.001 .. 0.003 with the pulse and
    the drive window patched to cover t >= 0."""
    import pyolfsysm.mb_model as mbm
    t = np.array([-0.001, 0.0, 0.001, 0.002, 0.003])
    monkeypatch.setattr(mbm, 'get_time_index', lambda mp, rv: pd.Index(t))
    monkeypatch.setattr(X, 'PULSE', (0.0, 0.0035))
    monkeypatch.setattr(X, 'DRIVE_WINDOW', (0.0, 0.003))
    Wp = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.0, 2.0], [1.0, 1.0]])   # claws x boutons
    c2k = np.array([0, 0, 0, 1, 1])
    wk = np.array([1.0, 1.0, 1.0, 0.5, 0.5])                                       # APL>claw weights (scaled)
    bouton = [np.array([[4.0] * 5, [1.0] * 5]), np.array([[0.0] * 5, [3.0] * 5])]  # odors: boutons x time
    inh = [np.array([[2.0] * 5]), np.array([[1.0] * 5])]                           # odors: APL x time
    rv = SimpleNamespace(kc=SimpleNamespace(wPNKC=Wp, wAPLKC=wk[:, None], claw_to_kc=c2k, inh_sims=inh if apl else [],
                                            thr=np.array([[1.0], [100.0]])),
                         pn=SimpleNamespace(bouton_sims=bouton))
    mp = SimpleNamespace(kc=SimpleNamespace(N=2, taum=0.01), time_dt=0.001)
    return rv, mp, Wp, c2k, wk, bouton, inh


def test_claw_stats_of_a_made_up_run(monkeypatch):
    rv, mp, Wp, c2k, wk, bouton, inh = _fake_run(monkeypatch)
    st = X.claw_stats(rv, mp, [0, 1], np.array([[True, False], [True, True]]))
    for j in range(2):
        x = Wp @ bouton[j][:, 0] - wk * inh[j][0, 0]       # constant over time: x(t) = W b - w APL(t - 1)
        D = np.array([x[c2k == k].clip(min=0).sum() for k in range(2)])
        U = np.array([(-x[c2k == k]).clip(min=0).sum() for k in range(2)])
        assert np.allclose(st['D'][:, j], D) and np.allclose(st['U'][:, j], U)
    # odor 0: KC 0's claws 4 - 2, 1 - 2, 2 - 2 -> drives 2, 0, 0 (one claw net-inhibited by 1); KC 1's 2 - 1, 5 - 1
    assert np.allclose(st['D'][:, 0], [2.0, 5.0]) and np.allclose(st['U'][:, 0], [1.0, 0.0])
    assert st['n_half'][0, 0] == 1 and st['n_half'][1, 0] == 1 and np.isnan(st['n_half'][0, 1])
    assert np.array_equal(st['n_driven'][:, 0], [1, 2])
    # without APL (the APL-off run) there is nothing to absorb
    rv0, mp0, *_ = _fake_run(monkeypatch, apl=False)
    st0 = X.claw_stats(rv0, mp0, [0, 1], np.ones((2, 2), dtype=bool), apl=False)
    assert np.allclose(st0['U'], 0.0)
    # the prediction's leaky integration: from the drive's level, so the peaks are D and D - U
    stp = X.claw_stats(rv, mp, [0, 1], np.ones((2, 2), dtype=bool), predict=True)
    assert np.allclose(stp['peak0'], stp['D']) and np.allclose(stp['peak1'], stp['D'] - stp['U'])
    assert np.allclose(stp['area0'][0, 0], (2.0 - 1.0) * 4 * 0.001) and np.allclose(stp['area1'][0, 0], 0.0)
    assert np.allclose(stp['area0'][1], 0.0)                           # KC 1 below its threshold of 100


def test_claw_stats_refuses_a_run_without_apl_activity(monkeypatch):
    rv, mp, *_ = _fake_run(monkeypatch, apl=False)
    with pytest.raises(RuntimeError):
        X.claw_stats(rv, mp, [0, 1], np.ones((2, 2), dtype=bool), apl=True)


def test_stats_summary_reads_the_responders():
    st = dict(n_half=np.array([[1.0, np.nan], [3.0, 2.0]]), U=np.array([[1.0, 0.0], [2.0, 4.0]]), D=np.array([[2.0, 0.0], [4.0, 8.0]]))
    counts = np.array([[3.0, 0.0], [2.0, 5.0]])
    r = X.stats_summary(st, counts, 1.0)
    assert r['claws_half_on_median'] == 2.0 and np.isclose(r['claws_half_on_mean'], 2.0)
    assert np.isclose(r['unabsorbed_resp_mean'], (1 + 2 + 4) / 3) and np.isclose(r['unabsorbed_over_drive_resp'], 0.5)
    r2 = X.stats_summary(dict(st), counts, 0.5)
    assert np.isclose(r2['pooled_over_drive_resp'], 0.25)
    off = X.stats_summary(st, counts, 0.0, 'off')
    assert set(off) == {'claws_half_off_median', 'claws_half_off_mean'}


# ---------------------------------------------------------------------------
# the control and E88's classes

def test_the_control_has_one_value_per_kc():
    rng = np.random.default_rng(0)
    idx = pd.MultiIndex.from_arrays([np.arange(50), ['g'] * 50], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    U = pd.DataFrame(rng.exponential(1.0, (50, 16)), index=idx, columns=ODORS)
    for s in (0.5, 1.0):
        cm = X.control_matrix(U, s)
        assert cm.shape == U.shape and cm.index.equals(U.index) and list(cm.columns) == ODORS
        assert np.allclose(cm.to_numpy(), cm.to_numpy()[:, [0]])                   # constant over the odors
        assert np.allclose(cm.iloc[:, 0].to_numpy(), s * U.mean(axis=1).to_numpy())
        assert np.isclose(cm.to_numpy().mean(), s * U.to_numpy().mean())          # the same total inhibition


def test_shared_classes_are_e88s():
    E88 = pytest.importorskip('experiments.E88_responder_sets')
    rng = np.random.default_rng(5)
    n = 300
    off = rng.poisson(np.where(rng.random((n, 16)) < 0.25, 6.0, 0.2)).astype(float)
    on = np.where(rng.random((n, 16)) < 0.55, off, 0.0)
    idx = pd.MultiIndex.from_arrays([np.arange(n), ['g'] * n], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    o = {c: ((pd.DataFrame(m, index=idx, columns=ODORS) > 1.5).astype(float), pd.DataFrame(m, index=idx, columns=ODORS))
         for c, m in (('on', on), ('off', off))}
    kc = pd.DataFrame(dict(breadth=1.0, breadth_factor=rng.uniform(0.5, 2, n), apl_weight=1.0, n_claws=6.0,
                           max_claws_one_glomerulus=2.0, kc_type='g'), index=pd.Index(np.arange(n), name=H.KC_ID_LEVEL))
    kc['glom:DM1'] = 1.0
    rows, _, _ = E88.shared_kcs(o, kc)
    want = {(r['group'], r['kc_class']): r['n_entries'] for r in rows}
    got = X.shared_classes(on, off)
    for (g, cl), nn in want.items():
        assert sum(1 for e in got[g] if e[4] == cl) == nn, (g, cl)


# ---------------------------------------------------------------------------
# the choice, the gate, the equal rule

def _fake_readings(patterns, s_of=None, flagged=(), T=None, raw=None, mean=None, arm='B1'):
    """Readings of made-up settings ({variant: pattern part on S, the same against every half and all recordings}),
    with T's ({variant: pattern}) when given; raw-count losses and mean effects default to the reference's."""
    rows = []
    for v, p in patterns.items():
        pv = X.parse_variant(v)
        for unit, pat in (('S', p), ) + ((('T', T[v]),) if T and v in T else ()):
            base = dict(arm=pv['arm'], variant=v, weighting=pv['wtag'], exponent=pv['exponent'], s=pv['s'], s_of=pv['s_of'],
                        net=pv['net'], fixed=pv['fixed'], control=pv['control'], permute_seed=np.nan, unit=unit,
                        n_wirings=4, n_flagged=int(v in flagged), n_at_bound=0, mean_effect=(mean or {}).get(v, 0.177))
            for t in ('halfA', 'halfB', 'halfC', 'halfD', 'z8', 'z6'):
                base[f'{t}_pattern'] = pat
                base[f'{t}_loss'] = pat + 0.01
            cal, rw = dict(base, readout='cal'), dict(base, readout='raw', z8_loss=(raw or {}).get(v, 0.078))
            rows += [cal, rw]
            if unit == 'T':
                for w in ('P3', 'P4', 'M2', 'M3'):
                    rows += [dict(cal, unit=w), dict(rw, unit=w)]
    return pd.DataFrame(rows)


def test_choice_needs_s_above_0_and_the_gain():
    rd = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1': 0.048, 'B1_F_k2_s1': 0.0495})
    c = X.choose91(rd, 'B1', 'A')
    assert c['chosen'] == 'B1_F_k1' and not c['confirm']                 # s = 0 wins: nothing to confirm
    rd = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1': 0.052, 'B1_F_k2_s1': 0.0490})
    c = X.choose91(rd, 'B1', 'B')
    assert c['chosen'] == 'B1_F_k2_s1' and np.isclose(c['gain'], 0.001) and not c['confirm']
    rd = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1_s1': 0.048, 'B1_none_s1': 0.040, 'B1_F_k2_s1_fixed': 0.030,
                         'B1_F_k2_s1_ctrl': 0.030}, flagged=('B1_none_s1',))
    c = X.choose91(rd, 'B1', 'C')
    assert c['chosen'] == 'B1_F_k1_s1' and c['confirm'] and 'B1_none_s1' in c['excluded']
    assert X.choose91(_fake_readings({'B1_F_k1_s1': 0.04}), 'B1', 'A')['chosen'] is None   # no reference: runs missing


@pytest.mark.skipif(not _shared(), reason="needs the shared directory (the flies' mean APL effect)")
def test_gate_on_four_directions():
    runs = pd.DataFrame(dict(variant=['B1_F_k1_s1'] * 4, wiring=['P3', 'P4', 'M2', 'M3'], flagged=[False] * 4))
    rd = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1_s1': 0.046, 'B1_F_k1_s1_ctrl': 0.049},
                        T={'B1_F_k2': 0.0538, 'B1_F_k1_s1': 0.0490, 'B1_F_k1_s1_ctrl': 0.0530})
    rows = [X.gate_direction(rd, runs, 'B1', h, 'B1_F_k1_s1', 'B1_F_k2', ['B1_F_k1_s1_ctrl']) for h in X.HALVES]
    assert all(np.isclose(r['g'], 0.0048) and r['n_singles_lower'] == 4 and r['v'] is True for r in rows)
    assert [r['split'] for r in rows] == [0, 0, 1, 1]
    g = X.gate91(rows)
    assert g['i'] and g['ii'] and g['n_singles'] == 16 and g['iii'] and g['control']
    assert np.isclose(g['sq_reduction_all'], 1 - 0.049 ** 2 / 0.0538 ** 2) and g['iv'] and g['passes']
    # one half chose nothing: g 0 there, its four cells not lower
    rows2 = rows[:3] + [X.gate_direction(rd, runs, 'B1', 'D', None, 'B1_F_k2', [])]
    g2 = X.gate91(rows2)
    assert np.isclose(g2['g_mean'], 0.75 * 0.0048) and g2['n_singles_lower'] == 12 and g2['ii'] and g2['i']
    assert not g2['iv'] and not g2['passes']                              # 0.75 x 17 % of the squared error
    # a smaller gain: (i) holds, (iv) does not (14 % < 15 %)
    rd1 = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1_s1': 0.046}, T={'B1_F_k2': 0.0538, 'B1_F_k1_s1': 0.0500})
    g1 = X.gate91([X.gate_direction(rd1, runs, 'B1', h, 'B1_F_k1_s1', 'B1_F_k2', []) for h in X.HALVES])
    assert g1['i'] and not g1['iv'] and not g1['control'] and not g1['passes']   # (and no control run)
    # a control that gains half as much fails the control rule
    rd3 = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1_s1': 0.046, 'B1_F_k1_s1_ctrl': 0.049},
                         T={'B1_F_k2': 0.0538, 'B1_F_k1_s1': 0.0490, 'B1_F_k1_s1_ctrl': 0.0510})
    rows3 = [X.gate_direction(rd3, runs, 'B1', h, 'B1_F_k1_s1', 'B1_F_k2', ['B1_F_k1_s1_ctrl']) for h in X.HALVES]
    assert not X.gate91(rows3)['control'] and X.gate91(rows3)['passes_i_to_iv']
    # a run on T outside the bands fails (iii); so does a raw-count loss 0.012 above the reference's
    rows4 = [X.gate_direction(rd, runs.assign(flagged=[True, False, False, False]), 'B1', h, 'B1_F_k1_s1', 'B1_F_k2',
                              ['B1_F_k1_s1_ctrl']) for h in X.HALVES]
    assert not X.gate91(rows4)['iii']
    rd5 = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1_s1': 0.046}, T={'B1_F_k2': 0.0538, 'B1_F_k1_s1': 0.0490},
                         raw={'B1_F_k1_s1': 0.090})
    assert not X.gate91([X.gate_direction(rd5, runs, 'B1', h, 'B1_F_k1_s1', 'B1_F_k2', []) for h in X.HALVES])['iii']


def test_equal_rule():
    rd = _fake_readings({'B1_F_k2': 0.050, 'B1_F_k1_s1': 0.052, 'B1_F_k0.5_s1': 0.0545, 'B1_none_s0.5': 0.049,
                         'B1_F_k2_s1': 0.0501, 'B1_F_k1': 0.050, 'B1_none_s1': 0.050},
                        raw={'B1_none_s1': 0.090}, mean={'B1_none_s0.5': 0.150})
    eq = {e['variant']: e for e in X.equal_rule91(rd, 'B1')}
    assert set(eq) == {'B1_F_k1_s1', 'B1_F_k0.5_s1', 'B1_none_s0.5', 'B1_none_s1'}   # s > 0 and exponent <= 1 only
    assert eq['B1_F_k1_s1']['equal'] and not eq['B1_F_k1_s1']['better']
    assert not eq['B1_F_k0.5_s1']['equal']                                # 0.0045 worse
    assert not eq['B1_none_s0.5']['equal'] and eq['B1_none_s0.5']['better']   # the mean APL effect 0.027 lower
    assert not eq['B1_none_s1']['equal']                                  # the raw-count loss 0.012 above


# ---------------------------------------------------------------------------
# E10_joint_fit.py's option

def test_e10_claw_inh_pool_option():
    J = pytest.importorskip('experiments.E10_joint_fit')
    assert J.claw_inh_pool_kws(None) == {} and J.claw_inh_pool_tag(None) == ''
    assert J.claw_inh_pool_kws(1) == {'claw_inh_pool': 1.0} and J.claw_inh_pool_tag(1.0) == 'pool1_'
    assert J.claw_inh_pool_kws(0.0) == {'claw_inh_pool': 0.0} and J.claw_inh_pool_tag(0.5) == 'pool0.5_'
    for bad in (-0.1, 1.5, float('nan'), float('inf')):
        with pytest.raises(ValueError):
            J.claw_inh_pool_kws(bad)
    assert inspect.signature(J.fit_one).parameters['claw_inh_pool'].default is None
    src = inspect.getsource(J.main)
    assert "'--claw-inh-pool'" in src and 'claw_inh_pool_tag(args.claw_inh_pool)' in src


@pytest.mark.skipif(not _shared(), reason='needs the shared directory (the fits of C2\'s family)')
def test_d_commands():
    J = pytest.importorskip('experiments.E10_joint_fit')
    c = X.d_command('real', 2.0)
    assert c[c.index('-k') + 1] == '2' and c[c.index('--claw-inh-pool') + 1] == '1' and '--init-from-exact' in c
    assert c[c.index('--experiment') + 1] == X.EXPERIMENT and '--wiring' not in c and '--wiring-file' not in c
    assert X.d_tag('real', 2.0) == 'rate0.068_aplboth_pool1_k2_cap3'
    c1 = X.d_command('P1', 1.0)
    assert '--init-from-exact' not in c1 and c1[c1.index('--wiring-seed') + 1] == '1'
    assert X.d_tag('P1', 2.0) == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_pool1_k2_cap3'
    cm = X.d_command('M1', 2.0)
    assert cm[cm.index('--wiring-tag') + 1] == 'curveball_within_type_matched-seed1'
    assert X.d_tag('M1', 2.0) == 'wiring-curveball_within_type_matched-seed1_rate0.068_aplboth_pool1_k2_cap3'
    assert J.claw_inh_pool_tag(1.0) in X.d_tag('M1', 2.0)


# ---------------------------------------------------------------------------
# part 0's identities (slow)

SLOW = pytest.mark.skipif(QUICK or not _shared() or not _has_pool(),
                          reason="part 0's identities: simulations of C2's family (minutes; needs olfsysm's claw_inh_pool)")


@SLOW
def test_s0_through_the_new_path_and_s1_apl_off_reproduce_the_saved_fit():
    from experiments import E83_cross_wiring_fit as E
    m = X.WiringModel91('real', check=True)
    X.set_pool(m.mp, 0.0)
    ns, nr = E._same({'on': m.run_on(1.0), 'off': m.sim(m.p_off())}, m.saved)
    assert ns == 0 and nr == 0
    X.set_pool(m.mp, 1.0)
    try:
        off1, on1 = m.sim(m.p_off()), m.run_on(1.0)
    finally:
        X.set_pool(m.mp, 0.0)
    assert np.array_equal(off1[1].to_numpy(), m.saved['off'][1].to_numpy())
    assert (on1[1].to_numpy() != m.saved['on'][1].to_numpy()).sum() > 0        # s = 1 acts with APL on
    assert (on1[1].to_numpy() > 1.5).sum() <= (m.saved['on'][1].to_numpy() > 1.5).sum()
    # the simulator built with the kwarg: s = 0 reproduces the saved fit, s = 1 equals the parameter set on a simulator
    for s, ref_on in ((0.0, m.saved['on']), (1.0, on1)):
        mk = X.WiringModel91('real', check=False, claw_inh_pool=s)
        assert mk.mp.kc.claw_inh_pool == s
        ns, nr = E._same({'on': mk.run_on(1.0), 'off': mk.sim(mk.p_off())}, {'on': ref_on, 'off': m.saved['off']})
        assert ns == 0 and nr == 0


@SLOW
def test_e89s_c1_point_reproduced():
    E89 = pytest.importorskip('experiments.E89_apl_breadth_weighting')
    w, wt = 'P1', ('F', 1.0)
    ref = E89.load_run(E89.run_tag('C1', wt, w))
    if ref is None:
        pytest.skip("E89's C1 run F at exponent 1 on P1 is not on this host or in the archive")
    m = X.WiringModel91(w, check=False)
    wk, _ = X._setting_weights(m, wt)
    res = X.b_reading(m, wk, s=0.0, label='test')
    r89 = pd.read_csv(H.results_dir(E89.EXPERIMENT) / 'runs_C1.csv').set_index('tag').loc[E89.run_tag('C1', wt, w)]
    assert res['apl_factor'] == float(r89['apl_factor'])
    assert np.array_equal(res['on'][1].to_numpy(), ref['on'][1].to_numpy())
    assert np.isfinite(res['stats']['U']).all() and (res['stats']['U'] >= 0).all()
