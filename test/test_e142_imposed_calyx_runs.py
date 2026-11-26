"""Tests of experiments/E142_imposed_calyx_runs.py (E142's step 1, the recorded PN dendrites imposed under A1's calyx
factors): the 56 points and their tags, the factors on the two calyx scales, the guards, and one point's run through
E135's runner functions (all faked here). No model is run."""
import json
import types

import numpy as np
import pandas as pd
import pytest

from experiments import E135_antennal_lobe_dendrites as E135
from experiments import E142_imposed_calyx_runs as M


@pytest.fixture(autouse=True)
def no_simulation(monkeypatch):
    """No test here runs the model: E83's simulator, E77's hold and E135's family stage raise unless a test plants its
    own."""
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83

    def refuse(*a, **k):
        raise RuntimeError('a test reached a model simulation')
    for mod, name in ((E83, 'build_simulator'), (S, 'hold_rates'), (E135, 'family_stage'), (E135, '_build'),
                      (E135, '_run')):
        monkeypatch.setattr(mod, name, refuse)


A1F = dict(aplpn=0.155, kckc=3.24, z=dict(aplpn=np.log(0.155), kckc=np.log(3.24)), source='x', sha256=None,
           measure='strength')


def test_planned_points_and_tags():
    pts = M.planned_step1()
    assert len(pts) == 56 == len(M.STEP1_ARMS) * len(E135.WIRINGS)
    tags = [M.tag_of(a, w) for a, w in pts]
    assert len(set(tags)) == 56
    assert sum(a.panel == 'apl16' for a, _ in pts) == 8
    assert M.tag_of(M.STEP1_ARMS[4], 'P1') == 'megamat_imposed_a1_glom_perm_P1'
    assert M.planned_step1(reverse=True) == pts[::-1]
    assert {a.table for a in M.STEP1_ARMS} <= {'imposed', 'imposed_glom_perm', 'imposed_odor_perm'}


def test_factors_on_the_two_calyx_scales():
    p = dict(fixed_thr=100.0, wAPLKC_scale=2.0, wKCAPL_scale=1.0, wAPLPN_scale=4.0, wPNAPL_scale=1.0, wKCKC_scale=10.0)
    by = {a.name + '_' + a.panel: a for a in M.STEP1_ARMS}
    got = {k: M.arm_params(p, M.arm_factors(a, A1F)) for k, a in by.items()}
    assert got['imposed_a1_megamat']['wAPLPN_scale'] == pytest.approx(4.0 * 0.155)
    assert got['imposed_a1_megamat']['wKCKC_scale'] == pytest.approx(32.4)
    assert got['imposed_aplpn_megamat']['wKCKC_scale'] == 10.0
    assert got['imposed_aplpn_megamat']['wAPLPN_scale'] == pytest.approx(0.62)
    assert got['imposed_kckc_megamat']['wAPLPN_scale'] == 4.0
    assert got['imposed_kckc_aplpn_none_megamat']['wAPLPN_scale'] == 0.0
    assert got['imposed_kckc_aplpn_none_megamat']['wKCKC_scale'] == pytest.approx(32.4)
    for k, q in got.items():
        assert {x: q[x] for x in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wPNAPL_scale')} == \
            dict(fixed_thr=100.0, wAPLKC_scale=2.0, wKCAPL_scale=1.0, wPNAPL_scale=1.0), k
    with pytest.raises(ValueError):
        M.arm_factors(M.Step1Arm('megamat', 'x', 'imposed', 'half', 'A1'), A1F)
    with pytest.raises(KeyError):
        M.arm_params(dict(wAPLPN_scale=1.0), dict(aplpn=1.0, kckc=1.0))


def test_a1_factors_read_e136s_classed_point(monkeypatch, tmp_path):
    from experiments import E136_stage_population_fit as F
    seen = {}
    monkeypatch.setattr(F, 'out_dir', lambda: tmp_path)

    def fz(arm):
        seen['measure'] = F.MEASURE
        return dict(aplpn=np.log(0.155), kckc=np.log(3.24), inhadd=3.48, inhsc=0.09)
    monkeypatch.setattr(F, 'fitted_z', fz)
    old = F.MEASURE
    a = M.a1_factors()
    assert seen['measure'] == 'strength' and F.MEASURE == old
    assert a['aplpn'] == pytest.approx(0.155) and a['kckc'] == pytest.approx(3.24) and set(a['z']) == {'aplpn', 'kckc'}


def test_runs_wait_for_the_authorization(tmp_path):
    with pytest.raises(PermissionError):
        M.require_step1_allowed(tmp_path)
    rec = M.write_step1_authorization('Yang', '2026-11-26 06:00', 'yes', 'viki', root=tmp_path)
    assert M.require_step1_allowed(tmp_path)['who'] == 'Yang' and rec['allowed']


@pytest.fixture
def env(tmp_path, monkeypatch):
    """One point through `run_point1` with E135's functions faked: a stage of 4 glomeruli and 3 odors, a table, a hold
    that returns the scales it was given, and a run of 6 KCs."""
    from experiments import E134_population_baseline as B
    res = tmp_path / 'res'
    res.mkdir()
    monkeypatch.setattr(M, 'out_dir', lambda: res)
    monkeypatch.setattr(M, 'step1_dirs', lambda: [tmp_path / 'local'])
    gl, od = ['g0', 'g1', 'g2', 'g3'], ['a', 'b', 'c']
    table = pd.DataFrame(np.arange(12.0).reshape(4, 3), index=gl, columns=od)
    p = dict(fixed_thr=100.0, wAPLKC_scale=2.0, wKCAPL_scale=1.0, wAPLPN_scale=4.0, wPNAPL_scale=1.0, wKCKC_scale=10.0)
    calls = dict(hold=[], save=[], stage=0, identity=0)

    def point_table(arm, wiring, ctx):
        calls['stage'] += 1
        return table, gl[:2], dict(seed=E135.arm_seed(wiring), set_odors=od, clip=dict(n_set=6, n_clipped=1),
                                   n_mapped=6), object(), dict(port=None, lane_identity=dict(exact=True))

    def hold(wiring, orn, fam, params, extra, targets, label):
        calls['hold'].append(dict(params=dict(params), extra=extra, targets=targets))
        idx = pd.MultiIndex.from_tuples([(i, 'ab') for i in range(6)], names=['kc_id', 'kc_type'])
        s = pd.DataFrame(np.ones((6, 3)), index=idx, columns=od)
        return dict(params=dict(params, fixed_thr=90.0, wAPLKC_scale=3.0),
                    params_off=dict(params, fixed_thr=90.0, wAPLKC_scale=0.0, wAPLPN_scale=0.0),
                    on=(s, s), off=(s, s), thr_factor=0.9, apl_factor=1.5, at_bound_thr=False, at_bound_apl=False,
                    n_sims=5)

    built = types.SimpleNamespace(stages=('pn', 'bouton', 'claw'), info=dict(model='C2'), orn=pd.DataFrame([[1.0]]),
                                  meta={})
    s6 = pd.DataFrame(np.ones((6, 3)), columns=od)
    monkeypatch.setattr(M, 'point_table', point_table)
    monkeypatch.setattr(M, 'hold', hold)
    monkeypatch.setattr(E135, 'arm_input', lambda st, t, g, route, set_odors=None: (
        dict(extra_kws=dict(pn_rates=t)), dict(route=1)))
    monkeypatch.setattr(E135, 'panel_orn', lambda panel: pd.DataFrame(np.zeros((4, 3)), index=gl, columns=od))
    monkeypatch.setattr(E135, 'input_columns', lambda t, orn: t)
    monkeypatch.setattr(E135, 'family_parts', lambda w: dict(p=dict(p), breadth=None, n_kc=6, size=18, target_on=2,
                                                             target_off=4))
    monkeypatch.setattr(B, 'held_targets', lambda pt, size, lock=None: (8, 4, 8 / size, 4 / size))
    monkeypatch.setattr(E135, '_build', lambda panel, wiring, params, orn, extra: built)
    monkeypatch.setattr(E135, '_run', lambda panel, wiring, b, params, params_off: dict(
        on=dict(spike_counts=s6, responses=s6), off=dict(spike_counts=s6, responses=s6)))
    monkeypatch.setattr(E135, '_hold_identity', lambda panel, wiring, out, h: dict(exact=True))
    monkeypatch.setattr(E135, 'run_pn_check', lambda *a, **k: dict(max_abs_diff_set=0.0, max_abs_diff_unset=0.0))
    monkeypatch.setattr(E135, 'imposed_run_check', lambda *a, **k: dict(max_abs_diff_hz=0.1, named=False))
    monkeypatch.setattr(E135, 'e134_pn_rates', lambda cond: table)
    monkeypatch.setattr(E135, 'stage_rates', lambda out: dict(bouton_share_above_300=0.0, bouton_cap_excluded=False))
    monkeypatch.setattr(E135, 'save_point', lambda tag, out, b, info, t, dirs=None: calls['save'].append(
        dict(tag=tag, info=info, dirs=dirs)) or [])
    monkeypatch.setattr(B, '_counts_row', lambda out: dict(rate_on=0.09, rate_off=0.2, rate_ratio=2.2))
    monkeypatch.setattr(M.H, 'model_provenance', lambda orn_deltas=None: dict(planted=True))
    monkeypatch.setattr(M.H, 'code_commits', lambda: dict(queue_commit='abc', queue_dirty=False))
    ctx = types.SimpleNamespace(lane='lane-t', a1=dict(A1F), lock134={}, dendrites=lambda: table)
    return types.SimpleNamespace(res=res, ctx=ctx, calls=calls, p=p, tmp=tmp_path)


def test_run_point1_applies_the_factors_and_saves(env):
    arm = M.STEP1_ARMS[3]    # KC-KC at A1's, APL->bouton at none
    tp, pp = env.res / 'runs_1_h_lane-t.csv', env.res / 'provenance_1_h_lane-t.json'
    row = M.run_point1(arm, 'P1', env.ctx, tp, pp)
    h = env.calls['hold'][0]['params']
    assert h['wAPLPN_scale'] == 0.0 and h['wKCKC_scale'] == pytest.approx(32.4) and h['fixed_thr'] == 100.0
    assert row['wAPLPN_scale'] == 0.0 and row['wKCKC_scale'] == pytest.approx(32.4) and row['wAPLKC_scale'] == 3.0
    assert row['aplpn'] == 'none' and row['kckc_factor'] == pytest.approx(3.24) and row['seed'] == E135.arm_seed('P1')
    assert env.calls['hold'][0]['targets']['target_on'] == 4
    s = env.calls['save'][0]
    assert s['tag'] == 'megamat_imposed_kckc_aplpn_none_P1' and s['dirs'] == [env.tmp / 'local']
    assert s['info']['factors'] == dict(aplpn=0.0, kckc=pytest.approx(3.24)) and s['info']['table_arm'] == 'imposed'
    assert M.runs1_table(env.res).tag.tolist() == ['megamat_imposed_kckc_aplpn_none_P1']
    prov = json.loads(pp.read_text())['points']['megamat_imposed_kckc_aplpn_none_P1']
    assert prov['model'] == dict(planted=True) and prov['info']['a1']['aplpn'] == 0.155
    assert (env.res / 'claims' / 'megamat_imposed_kckc_aplpn_none_P1').is_dir()
    # done or claimed points are skipped
    assert M.run_point1(arm, 'P1', env.ctx, tp, pp, done={'megamat_imposed_kckc_aplpn_none_P1'}) is None
    assert M.run_point1(arm, 'P1', env.ctx, tp, pp) is None


def test_run_point1_releases_its_claim_on_failure(env, monkeypatch):
    def moving_hold(wiring, orn, fam, params, extra, targets, label):
        raise RuntimeError('planted failure')
    monkeypatch.setattr(M, 'hold', moving_hold)
    with pytest.raises(RuntimeError, match='planted failure'):
        M.run_point1(M.STEP1_ARMS[0], 'real', env.ctx, env.res / 't.csv', env.res / 'p.json')
    assert not (env.res / 'claims' / 'megamat_imposed_a1_real').exists()


def test_run_point1_refuses_a_hold_that_moves_a_calyx_scale(env, monkeypatch):
    old = M.hold

    def moved(*a, **k):
        h = old(*a, **k)
        h['params']['wKCKC_scale'] *= 2
        return h
    monkeypatch.setattr(M, 'hold', moved)
    with pytest.raises(RuntimeError, match='moved wKCKC_scale'):
        M.run_point1(M.STEP1_ARMS[0], 'real', env.ctx, env.res / 't.csv', env.res / 'p.json')
    assert not (env.res / 'claims' / 'megamat_imposed_a1_real').exists()
