"""Tests of experiments/E147_reads.py (E147's runner: the reads of steps 1 and 2 on model flies with each ORN fly's own
input): the lanes and the claims (two lanes never take one read, a reclaim retakes only the lane's own), the read list's
order and counts, the held flags on planted counts, a read's row, run and provenance on a planted model, a lane's
resume and two lanes sharing a step, the input checks (an input fly with no deviation is the lab's input), the identity
on planted runs, and E136's `hold_at` with an input and a start. No model is built or run."""
import json
import os
import time
import types

import numpy as np
import pandas as pd
import pytest

from experiments import E136_stage_population_fit as F
from experiments import E147_megamat_input_flies as M
from experiments import E147_reads as E
from kcq import harness as H
from kcq import stage_readings as R

from test_e136_stage_population_fit import _fake_read_pipeline, _out_and_built

COMMITS = dict(queue_commit='q', queue_dirty=False, pyolfsysm_commit='p', pyolfsysm_dirty=False)
P = dict(fixed_thr=40.0, wAPLKC_scale=0.2, wKCAPL_scale=1e-4, wAPLPN_scale=0.4, wPNAPL_scale=1e-3, wKCKC_scale=26.0)


def _planted_sim(free=False, miss=0):
    out, built = _out_and_built()
    built.info = {'model': 'C2'}
    n_on = E.n_responses(out['on']['responses'])
    n_off = E.n_responses(out['off']['responses'])
    hold = None if free else dict(thr_factor=0.9, apl_factor=1.1, at_bound_thr=False, at_bound_apl=False, n_sims=7,
                                  history=[])
    return dict(out=out, built=built, params=dict(P), params_off=dict(P, wAPLKC_scale=0.0, wAPLPN_scale=0.0),
                hold=hold, kws={}, targets=(n_off, n_on + miss, 0.2, 0.1), arm='family', z={}, start='fit',
                p_start=dict(P), free=free)


@pytest.fixture
def planted(monkeypatch, tmp_path):
    """atlas's folders under tmp_path, the host 'hal', the locks, the commits, the provenance and the lab-input
    readings planted; `E.simulate` replaced by a planted model fly that records its calls."""
    from experiments import E134_population_baseline as B
    res, loc = tmp_path / 'res', tmp_path / 'loc'
    monkeypatch.setattr(E, '_host', lambda: 'hal')
    monkeypatch.setattr(B, '_host', lambda: 'hal')
    monkeypatch.setattr(M, 'results_dir', lambda: (res.mkdir(parents=True, exist_ok=True), res)[1])
    monkeypatch.setattr(M, 'local_dir', lambda part='': loc / part)
    monkeypatch.setattr(M, 'run_dirs', lambda step, setting, wiring, fly: [loc / f'step{step}' / setting / wiring / fly])
    monkeypatch.setattr(M, 'require_lock', lambda: {'planted': True})
    monkeypatch.setattr(R, 'require_lock', lambda: {})
    monkeypatch.setattr(E, 'commits', lambda: dict(COMMITS))
    monkeypatch.setattr(H, 'model_provenance', lambda orn_deltas=None: {'planted': True})
    monkeypatch.setattr(E, 'own_params', lambda base, panel, wiring: dict(params=dict(P, fixed_thr=50.0),
                                                                           params_off=dict(P), source='planted'))
    calls = []

    def sim(setting, panel, wiring, orn, lock, lock134, label='E147'):
        calls.append((setting, panel, wiring, orn.attrs.get('fly')))
        return _planted_sim(free=setting.endswith('free'))
    monkeypatch.setattr(E, 'simulate', sim)

    def orn_of(fly, panel):
        _, built = _out_and_built()
        o = built.orn.copy()
        o.attrs['fly'] = fly
        return o
    monkeypatch.setattr(E, 'input_of', orn_of)
    monkeypatch.setattr(E, 'input_checks', lambda panel: {'panel': panel})
    return types.SimpleNamespace(res=res, loc=loc, calls=calls)


def test_lanes(monkeypatch):
    monkeypatch.setattr(E, '_host', lambda: 'hal')
    assert E.lane_name() == 'hal' and E.lane_name('hal-2') == 'hal-2'
    for bad in ('atlas', '1', 'hal_2', 'hal-', 'halx'):
        with pytest.raises(ValueError):
            E.lane_name(bad)


def test_settings():
    assert E.split_setting('family') == ('family', False) and E.split_setting('A1free') == ('A1', True)
    with pytest.raises(ValueError):
        E.split_setting('A0')
    assert set(E.STARTS) == set(M.HELD_SETTINGS)
    assert E.HOLD_TOL == F.HOLD_TOL == M.HOLD_TOL


def test_read_list_order_and_counts():
    s1, s2 = E.planned(1), E.planned(2)
    assert len(s1) == 72 and {(p.setting, p.panel) for p in s1} == {('family', 'megamat')}
    assert [p.wiring for p in s1[:36]] == [w for w in M.T_WIRINGS for _ in range(9)]
    assert [p.fly for p in s1[:9]] == list(M.FLY_LABELS)
    assert {p.wiring for p in s1[36:]} == set(M.S_WIRINGS)
    assert len(s2) == 108 and {p.wiring for p in s2} == set(M.T_WIRINGS)
    assert [(p.setting, p.panel) for p in s2[::36]] == [('A1', 'megamat'), ('A1', 'apl16'), ('family', 'apl16')]
    for k in range(3):
        assert len({(p.setting, p.panel) for p in s2[36 * k:36 * (k + 1)]}) == 1
    tags = [M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in s1 + s2]
    assert len(set(tags)) == 180
    assert E.planned(1, reverse=True) == s1[::-1]
    f2 = E.planned(2, free=True)
    assert len(f2) == 108 and {p.setting for p in f2} == {'A1free', 'familyfree'}


def test_claims_two_lanes_never_take_one_read(monkeypatch, tmp_path):
    from experiments import E134_population_baseline as B
    monkeypatch.setattr(E, '_host', lambda: 'hal')
    monkeypatch.setattr(B, '_host', lambda: 'hal')
    assert E.claim('t', 'hal', tmp_path)
    assert not E.claim('t', 'hal-2', tmp_path)
    assert not E.claim('t', 'hal-2', tmp_path, reclaim=True)          # another lane's claim stays its own
    assert E.claim('t', 'hal', tmp_path, reclaim=True)                # this lane's stale claim
    d = tmp_path / 'claims' / 'u'
    d.mkdir(parents=True)                                              # an orphan: no claim.json
    assert not E.claim('u', 'hal', tmp_path, reclaim=True)
    old = time.time() - 2 * E.ORPHAN_CLAIM_S
    os.utime(d, (old, old))
    assert E.claim('u', 'hal', tmp_path, reclaim=True)
    assert json.loads((d / 'claim.json').read_text())['lane'] == 'hal'


def test_held_flags():
    """(held_exact, held) on planted counts; held is atlas's `is_held` (the lock's, 0.5 %)."""
    cases = [((2694, 5916, 2694, 5916, False), (True, True)),
             ((2694, 5917, 2694, 5916, False), (False, True)),       # E134's held run on real: one count off
             ((2696, 5916, 2694, 5916, False), (False, True)),
             ((2694, 5960, 2694, 5916, False), (False, False)),      # 0.74 % off
             ((2694, 5916, 2694, 5916, True), (False, False))]
    for (n_on, n_off, t_on, t_off, b), want in cases:
        assert E.held_flags(n_on, n_off, t_on, t_off, b) == want
        row = dict(n_on=n_on, n_off=n_off, target_on=t_on, target_off=t_off, at_bound_thr=b, at_bound_apl=False)
        assert M.is_held(row) == want[1]


def test_read_one_row_run_and_provenance(planted):
    root, d = planted.res, [planted.loc / 'r']
    tag = 'family_megamat_f3_P4'
    row = E.read_one(tag, 1, 'family', 'megamat', 'f3', 'P4', lambda: E.input_of('f3', 'megamat'), 'hal', root, d,
                     lock={})
    assert list(row)[:len(M.READ_COLUMNS)] == list(M.READ_COLUMNS)
    sim = _planted_sim()
    assert row['n_on'] == E.n_responses(sim['out']['on']['responses']) == row['target_on']
    assert row['held'] and row['held_exact'] and row['reading'] == 'held' and row['start'] == 'fit'
    assert row['thr_factor'] == 0.9 and row['apl_factor'] == 1.1 and row['hold_n_sims'] == 7
    assert row['thr_vs_own'] == pytest.approx(40.0 / 50.0) and row['own_source'] == 'planted'
    assert row['fixed_thr'] == 40.0 and row['off_wAPLKC_scale'] == 0.0 and row['queue_commit'] == 'q'
    tab = pd.read_csv(root / 'reads_hal.csv')
    assert list(tab.tag) == [tag] and list(tab.columns)[:len(M.READ_COLUMNS)] == list(M.READ_COLUMNS)
    assert list(M.reads_table().tag) == [tag]
    run = R.load_step3_run(tag, dirs=d)
    info = run.info
    assert (info['step'], info['setting'], info['fly'], info['start'], info['wiring']) == (1, 'family', 'f3', 'fit', 'P4')
    assert info['input_sha256'] == row['input_sha256']
    prov = json.loads((root / 'provenance_reads_hal.json').read_text())
    assert set(prov['points']) == {tag} and prov['points'][tag]['model'] == {'planted': True}
    # done, then claimed elsewhere: skipped without a simulation
    n = len(planted.calls)
    assert E.read_one(tag, 1, 'family', 'megamat', 'f3', 'P4', lambda: None, 'hal', root, d, done={tag}, lock={}) is None
    assert E.claim('family_megamat_f4_P4', 'hal-2', root)
    assert E.read_one('family_megamat_f4_P4', 1, 'family', 'megamat', 'f4', 'P4', lambda: None, 'hal', root, d,
                      lock={}) is None
    assert len(planted.calls) == n


def test_unheld_read_is_kept(planted, monkeypatch):
    monkeypatch.setattr(E, 'simulate', lambda *a, **k: _planted_sim(miss=1))
    row = E.read_one('family_megamat_f0_P3', 1, 'family', 'megamat', 'f0', 'P3', lambda: E.input_of('f0', 'megamat'),
                     'hal', planted.res, [planted.loc / 'r'], lock={})
    assert not row['held'] and not row['held_exact']          # one count of the planted 36 is beyond HOLD_TOL
    assert R.load_step3_run('family_megamat_f0_P3', dirs=[planted.loc / 'r']) is not None
    assert list(M.reads_table().tag) == ['family_megamat_f0_P3']


def test_failed_read_releases_its_claim(planted, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError('planted failure')
    monkeypatch.setattr(E, 'simulate', boom)
    with pytest.raises(RuntimeError, match='planted failure'):
        E.read_one('family_megamat_f0_P3', 1, 'family', 'megamat', 'f0', 'P3', lambda: E.input_of('f0', 'megamat'),
                   'hal', planted.res, [planted.loc / 'r'], lock={})
    assert not (planted.res / 'claims' / 'family_megamat_f0_P3').exists()


def test_free_read_row(planted, monkeypatch):
    row = E.read_one('familyfree_megamat_f0_P3', 1, 'familyfree', 'megamat', 'f0', 'P3',
                     lambda: E.input_of('f0', 'megamat'), 'hal', planted.res, [planted.loc / 'r'], lock={})
    assert row['reading'] == 'free' and not row['held'] and row['thr_factor'] == 1.0 and row['hold_n_sims'] == 2


def test_lanes_share_a_step_and_resume(planted, monkeypatch):
    monkeypatch.setattr(E, 'require_identity', lambda: {})
    s = E.part_reads(1, 'hal', max_reads=3)
    assert s['made_here'] == 3 and s['done_of_step'] == 3
    s = E.part_reads(1, 'hal-2', reverse=True, max_reads=2)
    assert s['made_here'] == 2
    pts = E.planned(1)
    first = [M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in pts]
    tab = M.reads_table()
    assert set(tab.tag) == set(first[:3]) | set(first[-2:])
    assert set(tab[tab.lane == 'hal-2'].tag) == set(first[-2:])
    # a read held by a lane still running is left to it
    assert E.claim(first[10], 'hal-3', planted.res)
    s = E.part_reads(1, 'hal')
    assert s['made_here'] == 72 - 5 - 1 and s['unheld_here'] == []
    tags = [M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in pts]
    sims = [(c[2], c[3]) for c in planted.calls]
    assert len(sims) == len(set(sims)) == 71                       # no read simulated twice
    assert first[10] not in set(M.reads_table().tag) and len(M.reads_table()) == 71
    assert [planted.calls[i][2:] for i in range(3)] == [(p.wiring, p.fly) for p in pts[:3]]
    assert set(M.reads_table().tag) == set(tags) - {first[10]}


def test_rows_and_runs_are_what_the_score_reads(planted, monkeypatch):
    """atlas's score readers on the rows and runs this runner writes: every READ_COLUMNS column in every row, the tags
    parsed, every run found where the score looks (`available`, `load_point`), held by `is_held`, the counts as
    written (`fill_counts` leaves them), the realism columns filled, and every input_sha256 that of `input_fly`
    (`check_inputs`)."""
    monkeypatch.setattr(E, 'require_identity', lambda: {})
    E.part_reads(2, 'hal', max_reads=4)
    E.part_reads(2, 'hal-2', reverse=True, max_reads=3)
    E.part_reads(1, 'hal', free=True, max_reads=2)
    rows = M.reads_table()
    assert len(rows) == 9 and set(M.READ_COLUMNS) <= set(rows.columns)
    assert not rows[list(M.READ_COLUMNS)].drop(columns=['threads']).isna().any().any()
    pts = M.available(M.points_from_rows(rows))
    assert len(pts) == 9 and set(pts.group) == {'A1/megamat', 'family/apl16', 'familyfree/megamat'}
    assert pts.held.all()
    filled = M.fill_counts(pts)
    assert (filled.n_on == pts.n_on).all() and (filled.n_off == pts.n_off).all()
    for r in pts.itertuples():
        assert M.load_point(r.tag, r.step) is not None
    real = M.megamat_realism_rows(pts[pts.panel == 'megamat'], {})
    assert not real[['rate_on', 'rate_off', 'thr_factor', 'fixed_thr', 'max_count_on']].isna().any().any()
    monkeypatch.setattr(M, 'input_fly', lambda f, panel='megamat', *a, **k: E.input_of(f, panel))
    assert M.check_inputs(pts) == []


def test_reads_wait_for_the_identity(planted, monkeypatch):
    monkeypatch.setattr(E, 'identity_root', lambda: (planted.res / 'identity').mkdir(parents=True, exist_ok=True)
                        or planted.res / 'identity')
    with pytest.raises(RuntimeError, match='identity'):
        E.part_reads(1, 'hal', max_reads=1)
    assert planted.calls == []


def _planted_inputs(monkeypatch, bad=None):
    from test_e147_megamat_input_flies import planted_long
    Fp = M.orn_flies(planted_long())
    monkeypatch.setattr(E, 'input_of', lambda f, p: M.input_fly(f, p, F=Fp))
    if bad == 'zero':
        monkeypatch.setattr(E, 'lab_of', lambda p: M.lab_input(p) * 1.001)
    else:
        monkeypatch.setattr(E, 'lab_of', lambda p: M.input_fly('f0', p, 'zero', F=Fp))
    return Fp


def test_input_with_no_deviation_is_the_lab_input(monkeypatch):
    _planted_inputs(monkeypatch)
    for panel in M.PANELS:
        c = E.input_checks(panel)
        assert c['lab_sha256'] == c['zero_sha256'] and c['max_mean_diff'] <= E.INPUT_MEAN_TOL
        assert all(v['n_entries_changed'] > 0 for v in c['flies'].values())
    assert set(E.input_checks('apl16')['kept_odors']) == set(M.APL16_KEPT_AT_LAB)


def test_input_checks_fail(monkeypatch):
    _planted_inputs(monkeypatch, bad='zero')
    with pytest.raises(RuntimeError, match='zero-deviation'):
        E.input_checks('megamat')
    Fp = _planted_inputs(monkeypatch)
    monkeypatch.setattr(E, 'input_of', lambda f, p: M.input_fly(f, p, F=Fp) + (1.0 if f == 'f2' else 0.0))
    with pytest.raises(RuntimeError, match='mean over the input flies'):
        E.input_checks('megamat')


def test_identity_on_planted_runs(planted, monkeypatch):
    """The identity read reproduces a planted reference exactly; a one-count difference fails and blocks the reads."""
    root = planted.res / 'identity'
    monkeypatch.setattr(E, 'identity_root', lambda: (root.mkdir(parents=True, exist_ok=True), root)[1])
    monkeypatch.setattr(E, 'lab_of', lambda panel: E.input_of('zero', panel))
    out, built = _out_and_built()
    R.save_step3_run('ref', out, built.meta, built.orn, dict(model='C2', wiring='P3', panel='megamat', reading='held'),
                     [planted.loc / 'ref'])
    ref_run = R.load_step3_run('ref', dirs=[planted.loc / 'ref'])       # as E134's held run is loaded
    ref = dict(params=dict(P), params_off=dict(P), run=ref_run, source='planted reference',
               outputs={c: (out[c]['responses'], out[c]['spike_counts']) for c in R.CONDITIONS})
    monkeypatch.setattr(E, 'lab_reading', lambda setting, panel, wiring: ref)
    from experiments import E134_population_baseline as B
    monkeypatch.setattr(B, 'runs_table', lambda: pd.DataFrame(columns=['tag', 'input_sha256']))
    rec = E.part_identity('hal', wirings=('P3',))
    r = rec['identity_family_megamat_P3']
    assert r['exact'] and set(r['counts'].values()) == {0} and r['frames'] and set(r['frames'].values()) == {0}
    assert r['params_max_rel_diff'] == 0.0
    assert E.require_identity()['identity_family_megamat_P3']['exact']
    assert len(planted.calls) == 1
    assert E.part_identity('hal', wirings=('P3',))['identity_family_megamat_P3']['exact']
    assert len(planted.calls) == 1                                      # recorded: not read again
    # a one-count difference on another wiring
    bad = out['on']['spike_counts'].copy()
    bad.iloc[0, 0] += 1
    ref2 = dict(ref, outputs=dict(ref['outputs'], on=(out['on']['responses'], bad)))
    monkeypatch.setattr(E, 'lab_reading', lambda setting, panel, wiring: ref2)
    with pytest.raises(RuntimeError, match='identity fails'):
        E.part_identity('hal', wirings=('M2',))
    with pytest.raises(RuntimeError, match='identity fails'):
        E.require_identity()


def _recording_hold(monkeypatch):
    """`_fake_read_pipeline` with E77's hold and E83's simulator recorded: the hold's kwargs and what its build()
    passes to E83's `build_simulator`."""
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83
    _fake_read_pipeline(monkeypatch)
    out, _ = _out_and_built()
    seen = []
    monkeypatch.setattr(E83, 'build_simulator', lambda w, **kw: ('sim', dict(kw, wiring=w)))

    def hold(build, **kw):
        seen.append(dict(kw, built=build()))
        return dict(params=dict(kw['params']), params_off=dict(kw['params']),
                    on=(out['on']['responses'], out['on']['spike_counts']),
                    off=(out['off']['responses'], out['off']['spike_counts']), thr_factor=1.0, apl_factor=1.0,
                    at_bound_thr=False, at_bound_apl=False, n_sims=2)
    monkeypatch.setattr(S, 'hold_rates', hold)
    return seen


def test_hold_with_input_is_hold_at(monkeypatch):
    """On the panel's input from the held megamat point, `hold_with_input` hands E77's hold and E83's simulator exactly
    what E136's `hold_at` does, for every arm of E136 (E136 is not changed); with an input and a start it holds on that
    input from those parameters, the caller's start left as it was."""
    from experiments import E134_population_baseline as B
    seen = _recording_hold(monkeypatch)
    z = dict(aplpn=np.log(0.5), kckc=np.log(2.0), inhadd=np.log(314.0 / 31.4088), inhsc=0.1, claw_thr=0.8)
    for arm in (F.FAMILY, 'A0', 'A1', 'A2', 'A3', 'A4', F.CONTROL):
        za = {k: v for k, v in z.items() if k in F.arm_searched(arm)}
        for panel in ('megamat', 'apl16'):
            seen.clear()
            a = F.hold_at(arm, za, 'P3', panel, label='x', lock={})
            b = E.hold_with_input(arm, za, 'P3', panel, B.panel_input(panel), None, {}, label='x')
            assert len(seen) == 2
            ka, kb = seen
            assert ka['built'][1].pop('orn') is kb['built'][1].pop('orn')
            assert ka == kb, (arm, panel)
            for k in ('kws', 'breadth', 'claw_thr', 'targets', 'p0', 'off_extra'):
                assert a[k] == b[k], (arm, panel, k)
    _, built = _out_and_built()
    other = built.orn * 2.0
    p0 = dict(P, fixed_thr=100.0, wKCKC_scale=10.0)
    seen.clear()
    h = E.hold_with_input('A0', dict(kckc=np.log(2.0)), 'P3', 'megamat', other, p0, {})
    assert h['orn'] is other and seen[0]['built'][1]['orn'] is other
    assert seen[0]['params']['fixed_thr'] == pytest.approx(100.0) and h['p0']['wKCKC_scale'] == pytest.approx(20.0)
    assert seen[0]['built'][1]['params']['wKCKC_scale'] == pytest.approx(20.0)
    assert p0['wKCKC_scale'] == 10.0                                    # the caller's start is not changed
    assert F.MEASURE == 'all'                                           # E136's measure restored


def test_simulate_starts(monkeypatch):
    """The family's hold starts from the wiring's fit, A1's from E136's start with its classed point (the lock's, which
    must be E136's), the free reading at the lab-input reading's parameters with no hold."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E77_input_screen as S
    got = _fake_read_pipeline(monkeypatch)
    _, built = _out_and_built()
    fit = dict(P, fixed_thr=90.0, wAPLKC_scale=0.05)
    monkeypatch.setattr(E83, 'saved_fit', lambda w: {'p': dict(fit)})
    z = dict(aplpn=np.log(0.5), kckc=0.0, inhadd=np.log(1000.0 / 31.4088), inhsc=0.0)
    monkeypatch.setattr(F, 'fitted_z', lambda arm: dict(z))
    lock = {'settings': {'A1_classed_z': dict(z)}}
    s = E.simulate('family', 'megamat', 'P3', built.orn, lock, {})
    assert got['hold_params'] == pytest.approx(fit) and s['start'] == 'fit' and s['p_start'] == pytest.approx(fit)
    assert s['kws'] == {} and s['hold'] is not None and s['arm'] == F.FAMILY
    s = E.simulate('A1', 'apl16', 'P3', built.orn, lock, {})
    assert got['hold_params']['fixed_thr'] == pytest.approx(40.0) and got['hold_params']['wAPLPN_scale'] == pytest.approx(0.2)
    assert s['kws']['pn_inhadd'] == pytest.approx(1000.0) and s['start'] == 'held_megamat_point'
    monkeypatch.setattr(E, 'own_params', lambda b, p, w: dict(params=dict(P, fixed_thr=55.0), params_off=dict(P),
                                                               source='own'))
    monkeypatch.setattr(S, 'hold_rates', lambda *a, **k: pytest.fail('the free reading holds nothing'))
    s = E.simulate('A1free', 'megamat', 'P3', built.orn, lock, {})
    assert s['hold'] is None and s['params']['fixed_thr'] == 55.0 and s['kws']['pn_inhadd'] == pytest.approx(1000.0)
    assert s['targets'] == (20, 10, 0.2, 0.1) and s['start'] == 'own'
    with pytest.raises(RuntimeError, match='A1_classed_z'):
        E.simulate('A1', 'megamat', 'P3', built.orn, {'settings': {'A1_classed_z': dict(z, kckc=0.1)}}, {})
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='table', rung='R1', outcome='a'))
    with pytest.raises(RuntimeError, match='table'):
        E.simulate('A1', 'megamat', 'P3', built.orn, lock, {})
