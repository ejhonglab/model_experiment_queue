"""Tests of experiments/E142_imposed_calyx_runs.py (E142's step 1, the recorded PN dendrites imposed under A1's calyx
factors): the lock's 56 points and the amendment's 40 and their tags, the factors on the two calyx scales, A1's PN stage
(its kwargs, its identity, its hold), the guards, one point's run through E135's runner functions (all faked here), the
classes as amended and the lock's reading of them, outcome (g), and route 1 on the PN layer alone (olfsysm's ORN, LN
and PN layers on a planted model, no mushroom body). No model of the mushroom body is run."""
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
    lock = M.planned_step1(arms=M.LOCK_ARMS)
    assert len(lock) == 56 == M.N_LOCK_RUNS == len(M.LOCK_ARMS) * len(E135.WIRINGS)
    pts = M.planned_step1()
    assert len(pts) == M.N_LOCK_RUNS + M.N_AMENDED_RUNS == 96 and len(M.AMENDED_ARMS) * 8 == M.N_AMENDED_RUNS
    assert pts[:56] == lock          # the lock's points first, in their order: existing tags and lanes stay valid
    tags = [M.tag_of(a, w) for a, w in pts]
    assert len(set(tags)) == 96
    assert sum(a.panel == 'apl16' for a, _ in pts) == 16
    assert M.tag_of(M.STEP1_ARMS[4], 'P1') == 'megamat_imposed_a1_glom_perm_P1'
    assert [M.tag_of(a, 'P3') for a in M.AMENDED_ARMS] == [
        'megamat_a1stage_imposed_P3', 'megamat_a1stage_imposed_glom_perm_P3', 'megamat_a1stage_imposed_odor_perm_P3',
        'apl16_a1stage_imposed_P3', 'megamat_imposed_aplpn_none_P3']
    assert {a.pn_stage for a in M.LOCK_ARMS} == {'family'}
    assert [a.pn_stage for a in M.AMENDED_ARMS] == ['A1', 'A1', 'A1', 'A1', 'family']
    assert all((a.aplpn, a.kckc) == ('A1', 'A1') for a in M.AMENDED_ARMS if a.pn_stage == 'A1')
    none = M.AMENDED_ARMS[-1]
    assert (none.aplpn, none.kckc, none.table) == ('none', 'family', 'imposed')
    assert M.planned_step1(reverse=True) == pts[::-1]
    assert {a.table for a in M.STEP1_ARMS} <= {'imposed', 'imposed_glom_perm', 'imposed_odor_perm'}
    assert set(M.ARM_KINDS) == {a.name for a in M.STEP1_ARMS}


def test_select_tags_and_the_cli(monkeypatch):
    pts = M.planned_step1()
    one = M.select_tags(pts, ['megamat_a1stage_imposed_P3'])
    assert [M.tag_of(a, w) for a, w in one] == ['megamat_a1stage_imposed_P3'] and one[0][0].pn_stage == 'A1'
    assert M.select_tags(pts, None) == pts and M.select_tags(pts, []) == pts
    with pytest.raises(ValueError, match='not planned'):
        M.select_tags(pts, ['megamat_a1stage_imposed_X9'])
    seen = {}
    monkeypatch.setattr(M, 'part_step1', lambda lane, rev, reclaim, tags: seen.update(lane=lane, tags=tags))
    monkeypatch.setattr(M, 'part_step1_dry', lambda wiring, tags: seen.update(dry=(wiring, tags)))
    M.main(['--part', 'step1', '--lane', 'hal-a', '--tags', 'megamat_a1stage_imposed_P3'])
    assert seen == dict(lane='hal-a', tags=['megamat_a1stage_imposed_P3'])
    M.main(['--part', 'step1_dry', '--wiring', 'P3', '--tags', 'apl16_a1stage_imposed_P3'])
    assert seen['dry'] == ('P3', ['apl16_a1stage_imposed_P3'])


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
    assert got['imposed_aplpn_none_megamat']['wAPLPN_scale'] == 0.0
    assert got['imposed_aplpn_none_megamat']['wKCKC_scale'] == 10.0
    assert got['a1stage_imposed_megamat']['wKCKC_scale'] == pytest.approx(32.4)
    for k, q in got.items():
        assert {x: q[x] for x in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wPNAPL_scale')} == \
            dict(fixed_thr=100.0, wAPLKC_scale=2.0, wKCAPL_scale=1.0, wPNAPL_scale=1.0), k
    with pytest.raises(ValueError):
        M.arm_factors(M.Step1Arm('megamat', 'x', 'imposed', 'half', 'A1'), A1F)
    with pytest.raises(ValueError, match="A1's classed point whole"):
        M.arm_factors(M.Step1Arm('megamat', 'x', 'imposed', 'A1', 'family', 'A1'), A1F)
    with pytest.raises(ValueError, match='pn_stage'):
        M.arm_factors(M.Step1Arm('megamat', 'x', 'imposed', 'A1', 'A1', 'E135'), A1F)
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


Z_A1 = dict(aplpn=np.log(0.155), kckc=np.log(3.24), inhadd=3.477822174445807, inhsc=0.08950048164176685)
A1PT = dict(z=dict(Z_A1), pn_kws=dict(pn_inhadd=1017.3, pn_inhsc=598.85), aplpn=0.155, kckc=3.24, source='x',
            sha256=None, measure='strength', e135_mode='search', best_eval=122, classed_is_best=True)


def test_a1_point_is_a1s_classed_point_whole(monkeypatch, tmp_path):
    from experiments import E136_stage_population_fit as F
    monkeypatch.setattr(F, 'out_dir', lambda: tmp_path)
    (tmp_path / 'fit_A1.json').write_text(json.dumps(dict(best_eval=122, best_z=Z_A1, classed_z=Z_A1)))
    seen = {}

    def fz(arm):
        seen['measure'] = F.MEASURE
        return dict(Z_A1)
    monkeypatch.setattr(F, 'fitted_z', fz)
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='search'))
    old = F.MEASURE
    pt = M.a1_point()
    assert seen['measure'] == 'strength' and F.MEASURE == old
    assert pt['z'] == pytest.approx(Z_A1) and pt['classed_is_best'] and pt['best_eval'] == 122
    # E136's pn_norm_kws: pn_inhadd olfsysm's default x exp(inhadd), pn_inhsc E77's inhsc_for x exp(inhsc)
    assert pt['pn_kws']['pn_inhadd'] == pytest.approx(1017.30, abs=0.01)
    assert pt['pn_kws']['pn_inhsc'] == pytest.approx(598.85, abs=0.01)
    assert pt['aplpn'] == pytest.approx(0.155) and pt['kckc'] == pytest.approx(3.24)
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='table'))
    with pytest.raises(RuntimeError, match="mode 'table'"):
        M.a1_point()
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='search'))
    monkeypatch.setattr(F, 'fitted_z', lambda arm: dict(aplpn=0.1, kckc=0.2))
    with pytest.raises(RuntimeError, match='lacks'):
        M.a1_point()


def _fake_stage(inhadd, inhsc):
    return types.SimpleNamespace(sp=types.SimpleNamespace(pn_inhadd=inhadd, pn_inhsc=inhsc))


def test_a1_stage_builds_e135s_stage_with_a1s_normalization(monkeypatch):
    seen = {}

    def fs(panel, wiring, purpose='step2', extra_kws=None):
        seen.update(panel=panel, wiring=wiring, extra_kws=extra_kws)
        return _fake_stage(extra_kws['pn_inhadd'], extra_kws['pn_inhsc']), dict(sim=None, identity=dict(ok=True))
    monkeypatch.setattr(E135, 'family_stage', fs)
    st, finfo = M.a1_stage('apl16', 'P3', A1PT)
    assert seen == dict(panel='apl16', wiring='P3', extra_kws=A1PT['pn_kws']) and st.sp.pn_inhadd == 1017.3
    monkeypatch.setattr(E135, 'family_stage', lambda panel, wiring, purpose='step2', extra_kws=None: (
        _fake_stage(31.4088, 368.6631), {}))
    with pytest.raises(RuntimeError, match="not A1's"):
        M.a1_stage('megamat', 'P3', A1PT)


def _saved_run(z, kws, params, params_off, counts):
    idx = pd.MultiIndex.from_tuples([(i, 'ab') for i in range(4)], names=['kc_id', 'kc_type'])
    cond = {c: dict(spike_counts=pd.DataFrame(counts[c], index=idx, columns=['2h @ -3', 'eb @ -3']),
                    responses=pd.DataFrame((counts[c] >= 2).astype(float), index=idx, columns=['2h @ -3', 'eb @ -3']))
            for c in ('on', 'off')}
    return types.SimpleNamespace(tag='A1_megamat_P3', conditions=cond,
                                 info=dict(z=z, params=params, params_off=params_off, **kws))


def test_a1_identity_reproduces_a1s_saved_read(monkeypatch):
    on = np.array([[0, 2], [3, 1], [0, 0], [5, 2]], float)
    off = on + 1
    run = _saved_run(dict(Z_A1), A1PT['pn_kws'], dict(fixed_thr=50.0), dict(fixed_thr=50.0, wAPLKC_scale=0.0),
                     dict(on=on, off=off))
    monkeypatch.setattr(M, 'a1_saved', lambda panel, wiring, point: run)
    calls = []

    def sim(p):
        calls.append(p)
        c = off if p.get('wAPLKC_scale') == 0.0 else on
        idx = pd.RangeIndex(4)
        return (pd.DataFrame((c >= 2).astype(float), index=idx, columns=['2h', 'eb']),
                pd.DataFrame(c, index=idx, columns=['2h', 'eb']))
    rec = M.a1_identity('megamat', 'P3', sim, A1PT)
    assert rec['exact'] and rec['z_equal'] and rec['pn_kws_equal'] and len(rec['checks']) == 4
    assert calls == [dict(fixed_thr=50.0), dict(fixed_thr=50.0, wAPLKC_scale=0.0)]
    # a count that differs, or a run at another point
    bad = run.conditions['on']['spike_counts'].copy()
    bad.iloc[0, 0] = 9
    run.conditions['on']['spike_counts'] = bad
    assert M.a1_identity('megamat', 'P3', sim, A1PT)['exact'] is False
    run.conditions['on']['spike_counts'].iloc[0, 0] = 0
    assert M.a1_identity('megamat', 'P3', sim, A1PT)['exact']
    assert not M.a1_identity('megamat', 'P3', sim, dict(A1PT, z=dict(Z_A1, kckc=1.0)))['exact']
    monkeypatch.setattr(M, 'a1_saved', lambda panel, wiring, point: None)
    assert M.a1_identity('megamat', 'P1', sim, A1PT)['exact'] is None


def test_a1_saved_reads_and_fit_evaluations(monkeypatch):
    from experiments import E136_stage_population_fit as F
    seen = []
    monkeypatch.setattr(F, 'load_read', lambda pt, claw_after=False: seen.append(('read', tuple(pt), F.MEASURE)) or 'r')
    monkeypatch.setattr(F, 'megamat_runs', lambda arm, ws, half: seen.append(('fit', arm, tuple(ws), half)) or
                        {ws[0]: 'f'})
    assert M.a1_saved('megamat', 'P3', A1PT) == 'r' and M.a1_saved('apl16', 'real', A1PT) == 'r'
    assert M.a1_saved('megamat', 'P1', A1PT) == 'f'
    assert M.a1_saved('megamat', 'P1', dict(A1PT, classed_is_best=False)) is None
    assert seen == [('read', ('A1', 'megamat', 'P3', 'held'), 'strength'), ('read', ('A1', 'apl16', 'real', None),
                                                                          'strength'),
                    ('fit', 'A1', ('P1',), 'held')]


def test_a1_lane_identity_is_cached_and_refuses(monkeypatch, tmp_path):
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(M.H, 'code_commits', lambda: dict(queue_commit='abc'))
    n = dict(calls=0)

    def ident(panel, wiring, sim, point):
        n['calls'] += 1
        return dict(exact=True, panel=panel, wiring=wiring)
    monkeypatch.setattr(M, 'a1_identity', ident)
    assert M.a1_lane_identity('megamat', 'P3', 'l', dict(sim=None), A1PT)['exact']
    assert M.a1_lane_identity('megamat', 'P3', 'l', dict(sim=None), A1PT)['cached'] and n['calls'] == 1
    ids = json.loads(M.identity_path('l').read_text())
    assert 'A1:megamat:P3' in ids
    monkeypatch.setattr(M, 'a1_identity', lambda panel, wiring, sim, point: dict(exact=False))
    with pytest.raises(RuntimeError, match='does not reproduce'):
        M.a1_lane_identity('apl16', 'P3', 'l', dict(sim=None), A1PT)
    monkeypatch.setattr(M, 'a1_identity', lambda panel, wiring, sim, point: dict(exact=None, note='none saved'))
    assert M.a1_lane_identity('megamat', 'real', 'l', dict(sim=None), A1PT)['exact'] is None


def test_step1_context_builds_each_pn_stage_once(monkeypatch):
    built = []
    monkeypatch.setattr(E135, 'read_lock', lambda: dict(glomeruli=['g0', 'g1']))
    monkeypatch.setattr(E135, 'family_stage', lambda panel, wiring, purpose='step2', extra_kws=None: (
        built.append((panel, wiring, extra_kws)) or _fake_stage(*(
            (extra_kws['pn_inhadd'], extra_kws['pn_inhsc']) if extra_kws else (31.4088, 368.6631))),
        dict(identity='port', sim=None, mp='mp')))
    ctx = M.Step1('l', a1=A1F, identity=False, a1pt=A1PT, keep_vars=True)
    s1, r1 = ctx.stage('megamat', 'P3', 'A1')
    s2, r2 = ctx.stage('megamat', 'P3', 'A1')
    s3, r3 = ctx.stage('megamat', 'P3')
    assert s1 is s2 and s3 is not s1 and len(built) == 2
    assert built == [('megamat', 'P3', A1PT['pn_kws']), ('megamat', 'P3', None)]
    assert r1['pn_stage'] == 'A1' and r1['pn_inhadd'] == 1017.3 and r3['pn_inhadd'] == 31.4088
    assert set(ctx._vars) == {('megamat', 'P3', 'A1'), ('megamat', 'P3', 'family')}
    ctx.free_vars()
    assert not ctx._vars
    with pytest.raises(ValueError):
        ctx.stage('megamat', 'P3', 'E135')


def test_hold_a1_is_e136s_hold_at_with_the_table(monkeypatch):
    from experiments import E136_stage_population_fit as F
    seen = {}

    def hold_at(arm, z, wiring, panel, label='E136', lock=None, extra_kws=None):
        seen.update(arm=arm, z=z, wiring=wiring, panel=panel, extra=extra_kws, measure=F.MEASURE)
        return dict(hold='h')
    monkeypatch.setattr(F, 'hold_at', hold_at)
    old = F.MEASURE
    assert M.hold_a1('P3', 'apl16', dict(pn_rates='T'), A1PT, 'lab') == dict(hold='h')
    assert seen == dict(arm='A1', z=A1PT['z'], wiring='P3', panel='apl16', extra=dict(pn_rates='T'),
                        measure='strength') and F.MEASURE == old


def test_e136_hold_at_adds_the_extra_kwargs_last(monkeypatch):
    from experiments import E136_stage_population_fit as F
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E134_population_baseline as B
    monkeypatch.setattr(F, 'start_params', lambda w, panel='megamat': (dict(wAPLPN_scale=1.0, wKCKC_scale=10.0,
                                                                            fixed_thr=50.0), {}))
    monkeypatch.setattr(F, 'a1_kws', lambda w, panel, orn, z: dict(pn_inhadd=1017.3, pn_inhsc=598.85))
    monkeypatch.setattr(F, 'read_targets', lambda w, panel, n, lock=None: (10, 5, 0.2, 0.1))
    monkeypatch.setattr(B, 'panel_input', lambda panel: pd.DataFrame(np.zeros((2, 3))))
    monkeypatch.setattr(E83, 'wiring_parts', lambda w: dict(breadth='b'))
    seen = {}

    def hold_rates(build, **k):
        seen.update(k)
        return dict(params=k['params'])
    monkeypatch.setattr(S, 'hold_rates', hold_rates)
    h = F.hold_at('A1', dict(aplpn=np.log(0.5), kckc=np.log(2.0), inhadd=3.0, inhsc=0.1), 'P3', 'megamat', lock={},
                  extra_kws=dict(pn_rates='T'))
    assert h['kws'] == dict(pn_inhadd=1017.3, pn_inhsc=598.85, pn_rates='T')
    assert h['p0']['wAPLPN_scale'] == pytest.approx(0.5) and h['p0']['wKCKC_scale'] == pytest.approx(20.0)
    assert h['targets'] == (10, 5, 0.2, 0.1) and seen['target_off'] == 10 and seen['apl_bracket'] == F.HOLD_RANGE
    assert F.hold_at('A1', {}, 'P3', 'megamat', lock={})['kws'] == dict(pn_inhadd=1017.3, pn_inhsc=598.85)


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
    builds = []
    monkeypatch.setattr(E135, '_build', lambda panel, wiring, params, orn, extra: builds.append(
        dict(panel=panel, params=dict(params), extra=dict(extra))) or built)
    ctx = types.SimpleNamespace(lane='lane-t', a1=dict(A1F), lock134={}, dendrites=lambda: table, point=lambda: A1PT)
    return types.SimpleNamespace(res=res, ctx=ctx, calls=calls, p=p, tmp=tmp_path, builds=builds, table=table)


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


def test_run_point1_on_a1s_stage_holds_as_a1s_reads(env, monkeypatch):
    """An A1-stage point: held by E136's hold_at of A1 (faked) with the table beside A1's normalization, built with
    both, its factors A1's, its targets the hold's."""
    seen = {}
    p0 = dict(env.p, wAPLPN_scale=4.0 * 0.155, wKCKC_scale=32.4, fixed_thr=55.0)

    def hold_a1(wiring, panel, extra, point, label):
        seen.update(wiring=wiring, panel=panel, extra=extra, point=point)
        idx = pd.MultiIndex.from_tuples([(i, 'ab') for i in range(6)], names=['kc_id', 'kc_type'])
        s = pd.DataFrame(np.ones((6, 3)), columns=['a', 'b', 'c'], index=idx)
        h = dict(params=dict(p0, fixed_thr=50.0, wAPLKC_scale=3.0),
                 params_off=dict(p0, fixed_thr=50.0, wAPLKC_scale=0.0, wAPLPN_scale=0.0), on=(s, s), off=(s, s),
                 thr_factor=0.91, apl_factor=1.2, at_bound_thr=False, at_bound_apl=False, n_sims=7)
        return dict(hold=h, p0=dict(p0), kws=dict(A1PT['pn_kws'], **extra), targets=(9, 4, 0.5, 0.0915))
    monkeypatch.setattr(M, 'hold_a1', hold_a1)
    monkeypatch.setattr(M, 'hold', lambda *a, **k: pytest.fail("the lock's hold on an A1-stage point"))
    arm = M.AMENDED_ARMS[0]
    tp, pp = env.res / 'runs_1_h_lane-t.csv', env.res / 'provenance_1_h_lane-t.json'
    row = M.run_point1(arm, 'P3', env.ctx, tp, pp)
    assert seen['wiring'] == 'P3' and seen['panel'] == 'megamat' and set(seen['extra']) == {'pn_rates'}
    b = env.builds[-1]
    assert b['extra']['pn_inhadd'] == 1017.3 and b['extra']['pn_inhsc'] == 598.85 and 'pn_rates' in b['extra']
    assert b['params']['fixed_thr'] == 50.0 and b['params']['wKCKC_scale'] == pytest.approx(32.4)
    assert row['tag'] == 'megamat_a1stage_imposed_P3' and row['pn_stage'] == 'A1'
    assert row['aplpn_factor'] == 0.155 and row['kckc_factor'] == 3.24 and row['start_fixed_thr'] == 55.0
    assert row['target_on'] == 4 and row['frac_on_target'] == 0.0915 and row['thr_factor'] == 0.91
    assert row['start'].startswith("the wiring's held megamat point")
    s = env.calls['save'][-1]
    assert s['info']['pn_stage'] == 'A1' and s['info']['a1_point'] == A1PT
    assert s['info']['targets']['hold_rule'] == M.A1_HOLD_RULE
    # a hold whose kwargs are not A1's normalization is refused, and the claim released
    monkeypatch.setattr(M, 'hold_a1', lambda *a, **k: dict(hold_a1(*a, **k), kws=dict(pn_rates='T')))
    with pytest.raises(RuntimeError, match="not A1's normalization"):
        M.run_point1(arm, 'P4', env.ctx, tp, pp)
    assert not (env.res / 'claims' / 'megamat_a1stage_imposed_P4').exists()


def test_run_point1_none_arm_takes_the_locks_hold(env):
    arm = M.AMENDED_ARMS[-1]
    row = M.run_point1(arm, 'M3', env.ctx, env.res / 't.csv', env.res / 'p.json')
    h = env.calls['hold'][-1]['params']
    assert h['wAPLPN_scale'] == 0.0 and h['wKCKC_scale'] == 10.0 and h['fixed_thr'] == 100.0
    assert row['tag'] == 'megamat_imposed_aplpn_none_M3' and row['pn_stage'] == 'family'
    assert row['start'] == "the wiring's fit (E83's wiring_parts)" and set(env.builds[-1]['extra']) == {'pn_rates'}


# ---------------------------------------------------------------------------
# the score's classes on planted values

def _F():
    from experiments import E136_stage_population_fit as F
    return F


def _vals(bouton, claw=0.9, soma=0.95, den=0.8, wirings=None, jitter=0.0):
    F = _F()
    ws = list(wirings or F.WIRINGS)
    out = {}
    for i, w in enumerate(ws):
        b = bouton[w] if isinstance(bouton, dict) else bouton
        d = dict(pn_dendrite=den, pn_bouton=b + jitter * i, kc_claw=claw + jitter * i, kc_soma=soma + jitter * i)
        d['downstream'] = float(np.sqrt(np.mean([d[s] ** 2 for s in M.DOWNSTREAM])))
        d['stage'] = float(np.sqrt(np.mean([d[s] ** 2 for s in ('pn_dendrite',) + M.DOWNSTREAM])))
        out[w] = d
    return out


def _mega(carry_bouton=0.8, a1_bouton=1.0, ref_bouton=1.10, apl_bouton=1.0, ctrl_bouton=1.2):
    F = _F()
    mega = {F.FAMILY: _vals(1.19, jitter=0.01), 'A1': _vals(a1_bouton, jitter=0.001),
            M.APL_REFERENCE: _vals(ref_bouton)}
    for a in M.APL_ARMS + M.LOCK_APL_ARMS:
        mega[M.variant_of(a)] = _vals(apl_bouton)
    for arm in (M.CARRY_ARM, M.LOCK_CARRY_ARM):
        mega[M.variant_of(arm)] = _vals(carry_bouton, jitter=0.001)
    for c in M.CONTROL_ARMS + M.LOCK_CONTROL_ARMS:
        mega[M.variant_of(c)] = _vals(ctrl_bouton, claw=1.0, soma=1.05)
    fam = mega[F.FAMILY]
    spreads = dict(downstream=F._sd(v['downstream'] for v in fam.values()),
                   **{s: F._sd(v[s] for v in fam.values()) for s in ('pn_dendrite',) + M.DOWNSTREAM})
    return mega, spreads


def test_apl_hides_step_class():
    F = _F()
    assert M.APL_ARMS == ('imposed_aplpn', 'imposed_aplpn_none')
    mega, spreads = _mega(carry_bouton=1.0, apl_bouton=1.0, ref_bouton=1.10)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['apl_hides_step']
    assert c['holds'] and c['cls'] == 'holds' and not c['reference_meets'] and len(c['arms_meeting']) == 2
    # the family's calyx gets below the yardstick too: not hidden by APL
    mega, spreads = _mega(carry_bouton=1.0, apl_bouton=1.0, ref_bouton=1.0)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['apl_hides_step']
    assert not c['holds'] and c['reference_meets'] and c['cls'] == 'does not hold'
    # below on only 5 of 8 wirings
    b = {w: (1.0 if i < 5 else 1.2) for i, w in enumerate(F.WIRINGS)}
    mega, spreads = _mega(carry_bouton=1.2, apl_bouton=1.2)
    for a in M.APL_ARMS:
        mega[M.variant_of(a)] = _vals(b)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['apl_hides_step']
    assert not c['holds'] and all(c['arms'][M.variant_of(a)]['n_below'] == 5 for a in M.APL_ARMS)
    # a failed hold takes the arm out, named; the other arm is read
    mega, spreads = _mega(carry_bouton=1.2, apl_bouton=1.2)
    mega[M.variant_of('imposed_aplpn')] = _vals(1.0)
    bad = {M.variant_of('imposed_aplpn'): ['megamat_imposed_aplpn_P3']}
    c = M.classify_step1(mega, spreads, 1.084, bad, {}, F)['apl_hides_step']
    assert not c['holds'] and c['arms'][M.variant_of('imposed_aplpn')]['reasons']
    assert c['left_out'] == [M.variant_of('imposed_aplpn')] and c['kept'] == [M.variant_of('imposed_aplpn_none')]
    assert c['cls'] == 'does not hold'
    mega[M.variant_of('imposed_aplpn_none')] = _vals(1.0)
    c = M.classify_step1(mega, spreads, 1.084, bad, {}, F)['apl_hides_step']
    assert c['holds'] and c['arms_meeting'] == [M.variant_of('imposed_aplpn_none')]
    # both left out: no class
    bad[M.variant_of('imposed_aplpn_none')] = ['megamat_imposed_aplpn_none_M3']
    c = M.classify_step1(mega, spreads, 1.084, bad, {}, F)['apl_hides_step']
    assert c['cls'] == 'no class' and not c['kept'] and len(c['left_out']) == 2
    # the reference's runs missing: no class
    mega, spreads = _mega(carry_bouton=1.2, apl_bouton=1.0)
    del mega[M.APL_REFERENCE]
    assert M.classify_step1(mega, spreads, 1.084, {}, {}, F)['apl_hides_step']['cls'] == 'no class'


def test_the_locks_reading_kept_as_the_record():
    """The lock's own reading: imposed_a1's failed holds give its classes no class (no arm left out)."""
    F = _F()
    mega, spreads = _mega(carry_bouton=0.8, a1_bouton=1.0, apl_bouton=1.0, ref_bouton=1.10)
    bad = {M.variant_of('imposed_a1'): [f'megamat_imposed_a1_{w}' for w in F.WIRINGS]}
    kw = dict(apl_arms=M.LOCK_APL_ARMS, carry_arm=M.LOCK_CARRY_ARM, control_arms=M.LOCK_CONTROL_ARMS,
              leave_out_failed=False)
    lk = M.classify_step1(mega, spreads, 1.084, bad, {}, F, **kw)
    assert lk['dendrites_carry']['cls'] == 'no class' and lk['dendrites_carry']['arm'] == 'E142_imposed_a1'
    assert any('holds fail' in r for r in lk['dendrites_carry']['reasons'])
    assert not lk['apl_hides_step']['left_out'] and lk['apl_hides_step']['holds']   # two other arms meet
    assert set(lk['apl_hides_step']['arms']) == {M.variant_of(a) for a in M.LOCK_APL_ARMS} | {M.APL_REFERENCE}
    # the amended reading of the same values: the main arm on A1's stage is read
    am = M.classify_step1(mega, spreads, 1.084, bad, {}, F)
    assert am['dendrites_carry']['arm'] == 'E142_a1stage_imposed' and am['dendrites_carry']['cls'] == 'holds'


def test_step1_outcome_and_not_readable():
    F = _F()
    mega, spreads = _mega(carry_bouton=0.80, a1_bouton=1.0, ctrl_bouton=1.2)
    main = M.variant_of(M.CARRY_ARM)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)
    o = M.step1_outcome(c, {}, {})
    assert o['letters'] == ['d', 'c'] and not o['not_readable']
    # the main arm's hold fails (here on an S wiring, which the class itself does not read): (g)
    bad = {main: ['megamat_a1stage_imposed_P1']}
    c = M.classify_step1(mega, spreads, 1.084, bad, {}, F)
    o = M.step1_outcome(c, bad, {})
    assert o['letters'][0] == 'g' and o['not_readable'] and o['holds_failed'] == ['megamat_a1stage_imposed_P1']
    # it leaves a band (its APL16 runs): (g)
    o = M.step1_outcome(M.classify_step1(mega, spreads, 1.084, {}, {main: ['real']}, F), {}, {main: ['real']})
    assert o['letters'][0] == 'g' and o['bands_left'] == ['real']
    # no better than A1: (e); a control's hold fails: no letter for the carry class, pending
    mega, spreads = _mega(carry_bouton=1.0, a1_bouton=1.0, apl_bouton=1.2)
    o = M.step1_outcome(M.classify_step1(mega, spreads, 1.084, {}, {}, F), {}, {})
    assert o['letters'] == ['e']
    bad = {M.variant_of(M.CONTROL_ARMS[0]): ['megamat_a1stage_imposed_glom_perm_P3']}
    o = M.step1_outcome(M.classify_step1(mega, spreads, 1.084, bad, {}, F), bad, {})
    assert o['letters'] == [] and o['pending'] and o['dendrites_carry'] == 'no class'


def test_holds_table_and_stopped_claims(tmp_path):
    F = _F()
    (tmp_path / 'claims' / 'megamat_imposed_kckc_P1').mkdir(parents=True)
    (tmp_path / 'claims' / 'megamat_imposed_kckc_P1' / 'claim.json').write_text(json.dumps(dict(
        tag='megamat_imposed_kckc_P1', lane='stopped', reason='repeats imposed_a1')))
    (tmp_path / 'claims' / 'megamat_imposed_aplpn_P1').mkdir(parents=True)
    (tmp_path / 'claims' / 'megamat_imposed_aplpn_P1' / 'claim.json').write_text(json.dumps(dict(
        tag='megamat_imposed_aplpn_P1', lane='viki-e142')))
    stopped = M.stopped_claims(tmp_path)
    assert stopped == {'megamat_imposed_kckc_P1': 'repeats imposed_a1'}
    for w in F.WIRINGS:
        if w != 'P1':
            stopped[f'megamat_imposed_kckc_{w}'] = 'x'
        if w != 'M3':
            stopped[f'megamat_imposed_a1_glom_perm_{w}'] = 'x'
    row = dict(rate_on=0.0915, frac_on_target=0.0915, rate_off=0.2, frac_off_target=0.2, at_bound_thr=False,
               at_bound_apl=False, thr_factor=0.5, apl_factor=1.0)
    rows = [dict(row, tag=f'megamat_imposed_a1_{w}', rate_off=0.18, at_bound_thr=True, thr_factor=0.1)
            for w in F.WIRINGS]
    rows += [dict(row, tag=f'megamat_imposed_aplpn_{w}') for w in F.WIRINGS]
    rows += [dict(row, tag='megamat_imposed_a1_glom_perm_M3', rate_off=0.19, at_bound_thr=True)]
    rows += [dict(row, tag=f'megamat_a1stage_imposed_{w}') for w in ('P3', 'P4')]
    ht = M.holds_table(pd.DataFrame(rows), 0.005, stopped, F.WIRINGS).set_index(['panel', 'arm'])
    assert len(ht) == len(M.STEP1_ARMS)
    assert ht.loc[('megamat', 'imposed_a1'), 'record'] == 'no class: a hold fails'
    assert ht.loc[('megamat', 'imposed_a1'), 'n_at_bound_thr'] == 8
    assert ht.loc[('megamat', 'imposed_a1'), 'thr_factor_min'] == 0.1
    assert ht.loc[('megamat', 'imposed_kckc'), 'record'].startswith('no class: a hold fails (its APL-off run is')
    assert ht.loc[('megamat', 'imposed_kckc'), 'n_stopped'] == 8
    assert ht.loc[('megamat', 'imposed_aplpn'), 'record'] == 'held on every wiring'
    assert ht.loc[('megamat', 'imposed_a1_glom_perm'), 'record'] == 'no class: a hold fails'
    assert ht.loc[('megamat', 'a1stage_imposed'), 'record'] == '2 of 8 read'
    assert ht.loc[('megamat', 'a1stage_imposed'), 'pn_stage'] == 'A1'
    assert not ht.loc[('megamat', 'a1stage_imposed'), 'lock'] and ht.loc[('apl16', 'imposed_a1'), 'lock']


def test_bouton_beside():
    F = _F()
    mega, _ = _mega(carry_bouton=1.0, a1_bouton=1.0)
    mega[M.variant_of(M.CARRY_ARM)]['P3']['pn_bouton'] = 1.2
    b = M.bouton_beside(mega, 1.084, F)
    v = b['variants'][M.variant_of(M.CARRY_ARM)]
    assert v['n_below_T'] == 3 and v['n_below_all'] == 7 and v['n_wirings'] == 8
    assert b['variants']['A1']['n_below_T'] == 4


def test_dendrites_carry_class():
    F = _F()
    mega, spreads = _mega(carry_bouton=0.80, a1_bouton=1.0, ctrl_bouton=1.2)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['dendrites_carry']
    assert c['holds'] and c['cls'] == 'holds' and c['downstream']['improves_2se'] and c['below_controls']
    assert set(c['controls']) == {M.variant_of(x) for x in M.CONTROL_ARMS}
    # a control below the arm at the boutons
    mega, spreads = _mega(carry_bouton=0.80, a1_bouton=1.0, ctrl_bouton=0.7)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['dendrites_carry']
    assert not c['holds'] and not c['below_controls'] and c['cls'] == 'does not hold'
    # no better than A1
    mega, spreads = _mega(carry_bouton=1.0, a1_bouton=1.0)
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['dendrites_carry']
    assert not c['holds'] and not c['downstream']['improves_2se']
    # a stage worse than A1's by more than its spread
    mega, spreads = _mega(carry_bouton=0.6, a1_bouton=1.0)
    for w in F.WIRINGS:
        mega[M.variant_of(M.CARRY_ARM)][w]['kc_soma'] = 2.0
    c = M.classify_step1(mega, spreads, 1.084, {}, {}, F)['dendrites_carry']
    assert 'kc_soma' in c['stage_worse'] and not c['holds']
    # a band left on a T wiring: no class
    mega, spreads = _mega(carry_bouton=0.80, a1_bouton=1.0)
    c = M.classify_step1(mega, spreads, 1.084, {}, {M.variant_of(M.CARRY_ARM): ['P3']}, F)['dendrites_carry']
    assert c['cls'] == 'no class' and not c['holds']


def test_variant_runs_and_holds(monkeypatch, tmp_path):
    from kcq import stage_readings as R
    seen = []

    def load(tag, dirs=None):
        seen.append((tag, [str(d) for d in dirs]))
        raise FileNotFoundError(tag)
    monkeypatch.setattr(R, 'load_step3_run', load)
    monkeypatch.setattr(M, 'step1_dirs', lambda: [tmp_path / 's1'])
    monkeypatch.setattr(E135, 'step2_dirs', lambda: [tmp_path / 's2'])
    F = _F()
    assert M.load_variant_run(F.ReadPoint('E142_imposed_a1', 'apl16', 'M2', None)) is None
    assert M.load_variant_run(F.ReadPoint('E135_imposed_odor_perm', 'megamat', 'P1', 'held')) is None
    assert seen == [('apl16_imposed_a1_M2', [str(tmp_path / 's1')]),
                    ('megamat_imposed_odor_perm_P1', [str(tmp_path / 's2')])]
    with pytest.raises(ValueError):
        M.load_variant_run(F.ReadPoint('A1', 'megamat', 'P1', 'held'))
    assert M.step1_variants() == ['E142_' + n for n in ('imposed_a1', 'imposed_aplpn', 'imposed_kckc',
                                                         'imposed_kckc_aplpn_none', 'imposed_a1_glom_perm',
                                                         'imposed_a1_odor_perm', 'a1stage_imposed',
                                                         'a1stage_imposed_glom_perm', 'a1stage_imposed_odor_perm',
                                                         'imposed_aplpn_none')]
    seen.clear()
    assert M.load_variant_run(F.ReadPoint('E142_a1stage_imposed', 'apl16', 'P3', None)) is None
    assert seen == [('apl16_a1stage_imposed_P3', [str(tmp_path / 's1')])]
    rows = pd.DataFrame([dict(tag='a', rate_on=0.0915, rate_off=0.2, frac_on_target=0.0915, frac_off_target=0.2,
                              at_bound_thr=False, at_bound_apl=False),
                         dict(tag='b', rate_on=0.0915, rate_off=0.2, frac_on_target=0.0915, frac_off_target=0.2,
                              at_bound_thr='True', at_bound_apl=False),
                         dict(tag='c', rate_on=0.0960, rate_off=0.2, frac_on_target=0.0915, frac_off_target=0.2,
                              at_bound_thr=False, at_bound_apl=False)])
    assert M.holds_failed(rows, 0.005) == ['b', 'c']


# ---------------------------------------------------------------------------
# the routes' check (route 2, E135's solved input), on E135's planted PN stage (olfsysm's parameters, no fit_mb_model)

def _toy_stage():
    import olfsysm as osm
    rng = np.random.default_rng(1)
    mp = osm.ModelParams()
    for k, v in dict(time_pre_start=-0.6, time_start=-0.3, time_end=0.7).items():
        setattr(mp, k, v)
    G, O = 10, 4
    delta = rng.uniform(-20, 150, (G, O))
    spont = rng.uniform(5, 30, (G, 1))
    delta[-2:] = 0
    mp.orn.data.spont = spont
    mp.orn.data.delta = delta
    mp.orn.n_physical_gloms = 12
    mp.kc.N = 10
    return E135.PNStage.from_mp(mp, [f'g{i}' for i in range(G)], [f'o{j}' for j in range(O)])


@pytest.fixture(scope='module')
def toy():
    st = _toy_stage()
    own = st.run(reading='preset')
    gl = [g for g in st.glomeruli if np.any(st.delta.loc[g].to_numpy() != 0)]
    return st, own, gl


def test_route2_tables_are_the_four_step1_uses():
    got = [(a.panel, a.table) for a in M.route2_tables()]
    assert got == [('megamat', 'imposed'), ('megamat', 'imposed_glom_perm'), ('megamat', 'imposed_odor_perm'),
                   ('apl16', 'imposed')]


def test_route2_target_sets_only_the_set_entries(toy):
    st, own, gl = toy
    table = own * 0 + 7.0
    has, od, tgt = M.route2_target(st, table, ['g0', 'g1', 'g9'], ['o1', 'o2'])
    assert has == gl and od == ['o1', 'o2']
    assert (tgt.loc[['g0', 'g1'], ['o1', 'o2']] == 7.0).all().all()
    pd.testing.assert_frame_equal(tgt.loc[gl, ['o0', 'o3']], own.loc[gl, ['o0', 'o3']])
    pd.testing.assert_frame_equal(tgt.loc[gl[2:], ['o1', 'o2']], own.loc[gl[2:], ['o1', 'o2']])


def test_route2_mean_map_finds_the_stages_own_fixed_point(toy):
    st, own, gl = toy
    m = M.route2_mean_map(st, own.loc[gl], gl, st.odors)
    assert list(m.odor) == st.odors
    assert m.root_in_grid.all()
    assert (m.F_at_own_mean.abs() < M.ROUTE2_ROOT_TOL).all()
    np.testing.assert_allclose(m.own_total_hz, m.target_total_hz)


def test_route2_mean_map_names_a_target_out_of_reach(toy):
    st, own, gl = toy
    tgt = own.loc[gl].copy()
    tgt['o0'] = tgt['o0'] * 1.6       # under the stage's normalization no input returns this much for o0
    tgt['o1'] = tgt['o1'] * 0.6       # less is within reach
    m = M.route2_mean_map(st, tgt, gl, st.odors).set_index('odor')
    assert not m.loc['o0', 'root_in_grid'] and m.loc['o0', 'min_F'] > M.ROUTE2_ROOT_TOL
    assert m.loc['o1', 'root_in_grid']
    assert m.loc['o0', 'target_total_hz'] == pytest.approx(1.6 * m.loc['o0', 'own_total_hz'])


def test_routes_verdict():
    row = dict(wiring='real', panel='megamat', table='imposed', converged=True, max_abs_diff_hz=0.1, at_edge=0,
               identity_ok=True)
    ok = pd.DataFrame([dict(row, wiring=w) for w in ('real', 'P1')])
    v = M.routes_verdict(ok, 2, 0.5)
    assert v['check1']['passes'] and v['check2'] == 'to run' and v['route'].startswith('2')
    assert not M.routes_verdict(ok, 3, 0.5)['check1']['passes']          # a table not read
    bad = ok.copy()
    bad.loc[1, 'max_abs_diff_hz'] = 20.9
    v = M.routes_verdict(bad, 2, 0.5)
    assert not v['check1']['passes'] and v['check1']['failing'] == ['megamat imposed P1']
    assert v['check2'].startswith('not run') and v['route'].startswith('1')
    bad = ok.copy()
    bad.loc[0, 'identity_ok'] = False
    assert not M.routes_verdict(bad, 2, 0.5)['check1']['passes']


def test_part_routes_check_and_summary_on_the_planted_stage(toy, monkeypatch, tmp_path):
    st, own, gl = toy
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(E135, 'require_e134_gate', lambda: None)
    monkeypatch.setattr(M, 'a1_factors', lambda: A1F)
    reach = own.copy()
    reach.loc[['g0', 'g1'], ['o0', 'o1']] = own.loc[['g0', 'g1'], ['o0', 'o1']].to_numpy() * 0.9

    def point_table(arm, wiring, ctx):
        return reach, ['g0', 'g1'], dict(set_odors=['o0', 'o1']), st, {}
    monkeypatch.setattr(M, 'point_table', point_table)
    rows = M.part_routes_check(['real'], 'test')
    assert len(rows) == 4 and rows.converged.all() and (rows.max_abs_diff_hz <= 0.5).all()
    assert rows.identity_ok.all() and (rows.n_odors_without_root == 0).all()
    assert (tmp_path / f'routes_check_1_{M._host()}_test.csv').exists()
    odors = pd.read_csv(tmp_path / f'routes_check_1_odors_{M._host()}_test.csv')
    assert len(odors) == 8 and (odors.solve_max_abs_diff_hz <= 0.5).all()
    v = M.part_routes_summary()
    assert v['check1']['n_tables'] == 4 and not v['check1']['passes']   # 4 of the 32 tables read
    assert json.loads((tmp_path / M.ROUTES_FILE).read_text())['check2'].startswith('not run')
    tabs = pd.read_csv(tmp_path / 'routes_check_tables.csv')
    assert [(r.panel, r.table) for r in tabs.itertuples()] == [(a.panel, a.table) for a in M.route2_tables()]
    assert len(pd.read_csv(tmp_path / 'routes_check_odors.csv')) == 8


# ---------------------------------------------------------------------------
# route 1 on the PN layer alone (olfsysm's ORN, LN and PN layers on a planted model; no mushroom body), and the dry run

def _toy_model(pn_kws=None, seed=1):
    """A planted model's PN stage and `family_stage`-like record (mp, rv, a wiring naming the glomeruli), at olfsysm's
    default normalization or at `pn_kws`."""
    import olfsysm as osm
    rng = np.random.default_rng(seed)
    mp = osm.ModelParams()
    for k, v in dict(time_pre_start=-0.6, time_start=-0.3, time_end=0.7).items():
        setattr(mp, k, v)
    G, O = 10, 4
    delta = rng.uniform(-20, 150, (G, O))
    delta[-2:] = 0
    mp.orn.data.spont = rng.uniform(5, 30, (G, 1))
    mp.orn.data.delta = delta
    mp.orn.n_physical_gloms = 12
    mp.kc.N = 10
    if pn_kws:
        mp.pn.inhadd, mp.pn.inhsc = float(pn_kws['pn_inhadd']), float(pn_kws['pn_inhsc'])
    rv = osm.RunVars(mp)
    gl = [f'g{i}' for i in range(G)]
    st = E135.PNStage.from_mp(mp, gl, [f'o{j}' for j in range(O)])
    return st, dict(mp=mp, rv=rv, info=dict(W=pd.DataFrame(np.zeros((3, G)), columns=gl)), identity=dict(ok=True),
                    sim=None)


def _toy_dendrites(st, n=6, seed=3):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.uniform(0, 40, (n, len(st.odors))), index=st.glomeruli[:n], columns=st.odors)


@pytest.mark.parametrize('pn_kws', [None, A1PT['pn_kws']])
def test_pn_only_route1_sets_the_table_at_either_normalization(pn_kws):
    st, finfo = _toy_model(pn_kws)
    if pn_kws:
        assert (st.sp.pn_inhadd, st.sp.pn_inhsc) == (1017.3, 598.85)
    own = st.run(reading='preset')
    gl = st.glomeruli[:6]
    table, rec = E135.imposed_rates(own, _toy_dendrites(st), gl, bounds=st.rate_bounds(), out_of_bounds='clip')
    assert rec['n_mapped'] == 24 and not np.allclose(table.loc[gl].to_numpy(), own.loc[gl].to_numpy())
    chk = M.pn_only_route1(finfo, st, table, gl, st.odors)
    assert chk['ok'] and chk['n_set'] == 24 and chk['n_unset'] == 16
    assert chk['max_abs_diff_set'] < 1e-8 and chk['max_abs_diff_unset'] < 1e-9
    assert chk['pn_inhadd'] == (1017.3 if pn_kws else st.sp.pn_inhadd)
    # the preset rates are cleared after: the layers then give the stage's own rates again
    import olfsysm as osm
    osm.run_ORN_LN_sims(finfo['mp'], finfo['rv'])
    osm.run_PN_sims(finfo['mp'], finfo['rv'])
    np.testing.assert_allclose(np.asarray(osm.pn_evoked_rates(finfo['mp'], finfo['rv'])), own.to_numpy(), atol=1e-9)


def test_dry_run_on_the_amended_tables(monkeypatch, tmp_path):
    from experiments import E136_stage_population_fit as F
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(M, 'a1_factors', lambda: A1F)
    monkeypatch.setattr(M, 'a1_point', lambda: A1PT)
    monkeypatch.setattr(M.H, 'code_commits', lambda: dict(queue_commit='abc'))
    st0, _ = _toy_model()
    monkeypatch.setattr(E135, 'read_lock', lambda: dict(glomeruli=st0.glomeruli[:6]))
    monkeypatch.setattr(E135, 'dendrite_means', lambda units, glomeruli: _toy_dendrites(st0).loc[list(glomeruli)])
    builds = []

    def fs(panel, wiring, purpose='step2', extra_kws=None):
        builds.append((panel, wiring, extra_kws))
        return _toy_model(extra_kws)
    monkeypatch.setattr(E135, 'family_stage', fs)
    monkeypatch.setattr(E135, 'family_parts', lambda w: dict(p=dict(wAPLPN_scale=4.0, wKCKC_scale=10.0,
                                                                    fixed_thr=100.0)))
    monkeypatch.setattr(F, 'start_params', lambda w, panel='megamat': (dict(wAPLPN_scale=2.0, wKCKC_scale=20.0,
                                                                            fixed_thr=50.0), {}))
    out = M.part_step1_dry('P3')
    assert len(out['arms']) == len(M.STEP1_ARMS) and len(builds) == 4      # one build per (panel, PN stage)
    assert sorted((p, bool(k)) for p, _, k in builds) == [('apl16', False), ('apl16', True), ('megamat', False),
                                                          ('megamat', True)]
    assert all(a['pn_only_route1']['ok'] for a in out['arms'].values())
    a1 = out['arms']['megamat_a1stage_imposed_P3']
    assert a1['pn_stage'] == 'A1' and a1['pn_inhadd'] == 1017.3 and a1['n_set_glomeruli'] == 6
    assert a1['wAPLPN_scale'] == pytest.approx(2.0 * 0.155) and a1['wKCKC_scale'] == pytest.approx(20.0 * 3.24)
    assert a1['own_min_hz'] <= a1['table_min_hz'] and a1['table_max_hz'] <= a1['own_max_hz'] + 1e-9
    none = out['arms']['megamat_imposed_aplpn_none_P3']
    assert none['pn_stage'] == 'family' and none['wAPLPN_scale'] == 0.0 and none['wKCKC_scale'] == 10.0
    assert out['a1_point'] == A1PT and set(out['stages']) == {'megamat/family', 'megamat/A1', 'apl16/family',
                                                             'apl16/A1'}
    saved = json.loads((tmp_path / 'step1_dry_P3.json').read_text())
    assert saved['arms']['apl16_a1stage_imposed_P3']['pn_stage'] == 'A1'
    builds.clear()
    one = M.part_step1_dry('P3', ['megamat_a1stage_imposed_P3'])
    assert list(one['arms']) == ['megamat_a1stage_imposed_P3'] and len(builds) == 1


def test_e143_beside_on_planted_views(monkeypatch, tmp_path):
    """E143's interval on the amended main arm's change from A1, through E143's own sides and errors, on planted views:
    the arm near the flies' mean, A1 further, the controls furthest."""
    from experiments import E143_megamat_fly_sampling as S
    rng = np.random.default_rng(0)
    flies = [f'f{i}' for i in range(5)]
    entries = [f'o{j}' for j in range(6)]
    keys = [S.megamat_key(s) for s in S.LOSS_STAGES]
    paths = {k: tmp_path / f"flies_{k.replace('/', '_')}.parquet" for k in keys}
    model = {M.E143_A1: [], **{M.variant_of(a): [] for a in (M.CARRY_ARM, *M.CONTROL_ARMS)}}
    noise = {M.E143_A1: 0.4, M.variant_of(M.CARRY_ARM): 0.05, **{M.variant_of(c): 1.0 for c in M.CONTROL_ARMS}}
    for k in keys:
        rows = []
        for q in ('strength/rate', 'strength/mean'):
            X = rng.normal(1.0, 0.3, (len(flies), len(entries)))
            rows += [dict(key=k, part=q, fly=f, entry=e, value=X[i, j], unit=0.3)
                     for i, f in enumerate(flies) for j, e in enumerate(entries)]
            mean = X.mean(axis=0)
            for nm in model:
                for w in S.T_WIRINGS:
                    v = mean + rng.normal(0, noise[nm], len(entries))
                    model[nm] += [dict(setting=nm, input='megamat', key=k, wiring=w, part=q, entry=e, value=v[j])
                                  for j, e in enumerate(entries)]
        pd.DataFrame(rows).to_parquet(paths[k], index=False)
    pd.DataFrame(model[M.E143_A1]).to_parquet(tmp_path / 'a1.parquet', index=False)
    # (E143's resamples draw the held-out panels after the stages, from one generator, so the lock lists them too)
    monkeypatch.setattr(S, 'require_lock', lambda: dict(recordings={k: dict(flies=flies) for k in
                                                                    keys + list(S.HOLDS_PANELS + S.BESIDE_PANELS)},
                                                        sha256='abc'))
    monkeypatch.setattr(S, 'views_path', lambda name, panel: tmp_path / 'a1.parquet')
    monkeypatch.setattr(S, 'flies_path', lambda key: paths[key])
    monkeypatch.setattr(S, '_ENTRIES', {})
    monkeypatch.setattr(M, 'load_point', lambda tag: tag)
    monkeypatch.setattr(S, 'read_runs', lambda name, panel, runs: pd.DataFrame(model[name]))
    spreads = dict(downstream=0.01, pn_bouton=0.01, kc_claw=0.01, kc_soma=0.01)
    b = M.e143_beside(spreads, 'holds')
    assert b['status'] == 'given' and b['n_boot'] == S.N_BOOT and b['ci'] == list(S.CI)
    ch = b['changes']
    assert set(ch) == {'stage', 'downstream', *S.LOSS_STAGES}
    assert ch['downstream']['all_flies'] < 0 and ch['downstream']['hi'] < 0 and ch['downstream']['excludes_zero']
    assert b['class_all_flies'] == 'holds' and b['share_holds'] > 0.9 and b['share_scored'] == b['share_holds']
    # a missing run or A1's views not on this host: not given
    monkeypatch.setattr(M, 'load_point', lambda tag: (_ for _ in ()).throw(FileNotFoundError(tag)))
    assert M.e143_beside(spreads, 'holds')['status'].startswith('not given')
    monkeypatch.setattr(S, 'views_path', lambda name, panel: tmp_path / 'none.parquet')
    assert M.e143_beside(spreads, 'holds')['status'].startswith('not given')


def test_write_variants_fills_the_family_stage_normalization(tmp_path):
    """step1_variants.csv, the archive's params table: a row written before the amendment (no pn_stage, no
    normalization: the family's PN stage) gets olfsysm's defaults; the amendment's rows keep what they record."""
    rows = [dict(tag='megamat_imposed_a1_P1', panel='megamat', arm='imposed_a1', table_arm='imposed', wiring='P1',
                 seed=1, aplpn_factor=0.155, kckc_factor=3.24, fixed_thr=8.3, wKCKC_scale=133.3),
            dict(tag='megamat_a1stage_imposed_P1', panel='megamat', arm='a1stage_imposed', table_arm='imposed',
                 wiring='P1', seed=1, aplpn_factor=0.155, kckc_factor=3.24, pn_stage='A1', pn_inhadd=1017.303,
                 pn_inhsc=598.8538, fixed_thr=46.1, wKCKC_scale=133.3),
            dict(tag='megamat_imposed_aplpn_none_P2', panel='megamat', arm='imposed_aplpn_none', table_arm='imposed',
                 wiring='P2', seed=2, aplpn_factor=0.0, kckc_factor=1.0, pn_stage='family', pn_inhadd=31.4088,
                 pn_inhsc=368.6631, fixed_thr=30.3, wKCKC_scale=22.9)]
    pd.DataFrame(rows[:1]).to_csv(tmp_path / 'runs_1_hal_hal-s1.csv', index=False)
    pd.DataFrame(rows[1:]).to_csv(tmp_path / 'runs_1_hal_hal-s3.csv', index=False)
    v = M.write_variants(tmp_path).set_index('tag')
    assert (tmp_path / 'step1_variants.csv').exists()
    assert list(v.pn_stage) == ['family', 'A1', 'family']
    assert v.loc['megamat_imposed_a1_P1', 'pn_inhadd'] == E135.OLFSYSM_DEFAULTS['pn_inhadd']
    assert v.loc['megamat_imposed_a1_P1', 'pn_inhsc'] == E135.OLFSYSM_DEFAULTS['pn_inhsc']
    assert (v.loc['megamat_a1stage_imposed_P1', 'pn_inhadd'], v.loc['megamat_a1stage_imposed_P1', 'pn_inhsc']) == \
        (1017.303, 598.8538)
    assert v.loc['megamat_imposed_aplpn_none_P2', 'pn_inhadd'] == 31.4088
    assert v.loc['megamat_imposed_a1_P1', 'wKCKC_scale'] == 133.3 and (v.route == 1).all()
    assert "A1's APL->bouton and KC-KC factors" in v.loc['megamat_imposed_a1_P1', 'kind']
