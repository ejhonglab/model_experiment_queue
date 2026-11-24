"""Tests of experiments/E136_stage_population_fit.py (E136, the best family fitted to all megamat stages at once): the
arms and wirings, the log factors on each wiring's parameters, the hold's secant, the bounds, the permutations of the
control, the halves' references in all flies' units, the divisive claw reader, the in-memory run against a saved and
reloaded one, and the loss on planted runs and planted flies. No model is built or run; the claw identity against
olfsysm runs only with E136_REAL=1."""
import os
import types

import numpy as np
import pandas as pd
import pytest

from experiments import E136_stage_population_fit as F
from kcq import population_fit as P
from kcq import stage_readings as R

from test_e134_stage_readings import ODORS, planted_flies, planted_run

LOCK = {'remy_megamat_criterion_rate': {'mean': 0.0915}}
GL = ['G0', 'G1', 'G2', 'G3']                  # the planted run's bouton glomeruli read (atlas's tests' subset)


def planted_full(odors=ODORS):
    return {st: planted_flies(st, odors=odors, named=R.NAMED_UNITS[st], n_units=6 if R.NAMED_UNITS[st] else 60,
                              seed=hash(st) % 1000) for st in F.LOSS_STAGES}


def test_arms_and_wirings():
    F.check_wirings()
    assert set(F.WIRINGS) == set(F.S_WIRINGS) | set(F.T_WIRINGS) and not set(F.S_WIRINGS) & set(F.T_WIRINGS)
    assert tuple(F.ARMS) == ('A0', 'A1', 'A2', 'A3', 'A4')
    n = {name: len(arm.searched) + len(F.HELD) for name, arm in F.ARMS.items()}
    assert n == {'A0': 4, 'A1': 4, 'A2': 4, 'A3': 5, 'A4': 3}
    for name, arm in F.ARMS.items():
        assert arm.name == name and n[name] <= F.MAX_FREE and set(arm.searched) <= set(F.BOUNDS)
        assert name in F.MAXFEV
    assert F.ARMS['A4'].searched == ('aplpn',) and F.ARMS['A4'].kws == {'kckc_fixed': 0.0}
    assert F.CONTROL not in F.ARMS and F.base_arm(F.CONTROL) == 'A0' and F.base_arm(F.FAMILY) == 'A0'
    assert F.arm_kws('A2') == {'apl_inh_divisive': True} and F.arm_kws('A4') == {} and F.arm_kws('A0') == {}
    assert F.claw_reader_of('A2') is F.claw_inputs_divisive and F.claw_reader_of('A0') is None
    assert [F.seed_index(w) for w in F.WIRINGS] == list(range(8))
    with pytest.raises(RuntimeError, match='E135'):
        F.arm_kws('A1')


def _model(arm='A0'):
    m = F.ArmModel.__new__(F.ArmModel)
    m.arm, m.searched, m.wirings = arm, F.ARMS[F.base_arm(arm)].searched, ('real',)
    p = dict(fixed_thr=40.0, wAPLKC_scale=0.2, wKCAPL_scale=1e-4, wAPLPN_scale=0.4, wPNAPL_scale=1e-3, wKCKC_scale=26.0)
    m.base = {'real': p}
    m.base_off = {'real': dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0, wKCAPL_scale=7e-5, wPNAPL_scale=9e-4)}
    return m


def test_params_for_factors_and_ties():
    m = _model()
    on, off, ct = m.params_for('real', {})
    assert on == m.base['real'] and off == m.base_off['real'] and ct is None
    on, off, _ = m.params_for('real', dict(aplpn=np.log(0.5), kckc=np.log(0.25)), x_thr=np.log(2.0), x_apl=np.log(3.0))
    assert on['fixed_thr'] == pytest.approx(80.0) and off['fixed_thr'] == on['fixed_thr']
    assert on['wAPLKC_scale'] == pytest.approx(0.6)
    # the KC->APL and bouton->APL scales keep the wiring's own ratios to the APL->KC scale
    assert on['wKCAPL_scale'] / on['wAPLKC_scale'] == pytest.approx(1e-4 / 0.2)
    assert on['wPNAPL_scale'] / on['wAPLKC_scale'] == pytest.approx(1e-3 / 0.2)
    assert on['wAPLPN_scale'] == pytest.approx(0.2) and on['wKCKC_scale'] == pytest.approx(6.5)
    assert off['wAPLKC_scale'] == 0.0 and off['wAPLPN_scale'] == 0.0 and off['wKCKC_scale'] == on['wKCKC_scale']
    with pytest.raises(ValueError):
        m.params_for('real', dict(thr=0.1))         # set by the hold, not searched
    with pytest.raises(ValueError):
        _model('A4').params_for('real', dict(kckc=0.1))
    _, _, ct = _model('A3').params_for('real', dict(claw_thr=0.6))
    assert ct == pytest.approx(30.0)
    assert _model('A3').params_for('real', {})[2] == 0.0


def _counts_run(k, n0=10000.0):
    """A planted run whose number of responses falls with the log factor as n0 * exp(k x) (rounded)."""
    calls = []

    def run(x):
        calls.append(x)
        n = int(round(min(n0 * np.exp(k * x), 40000)))
        return pd.DataFrame(np.r_[np.ones(n), np.zeros(40000 - n)].reshape(-1, 1)), None
    return run, calls


def test_secant_hold():
    run, calls = _counts_run(-2.5)
    x, out, info = F.secant_hold(run, 6000, 0.0, np.log(0.1), np.log(10.0), -1.0)
    assert info['ok'] and abs(info['n'] - 6000) <= 0.005 * 6000 and info['n_sims'] <= F.HOLD_MAX_SIMS
    assert F.n_responses(out[0]) == info['n'] and x == info['x']
    # at the start already: one simulation, no step
    run, calls = _counts_run(-2.5)
    x, _, info = F.secant_hold(run, 10000, 0.0, np.log(0.1), np.log(10.0), -3.0)
    assert x == 0.0 and info['n_sims'] == 1 and info['ok']
    # out of reach: the hold ends at an end of its range and says so
    run, calls = _counts_run(-0.1)
    x, _, info = F.secant_hold(run, 100, 0.0, np.log(0.1), np.log(10.0), -1.0)
    assert not info['ok'] and info['at_end'] and x == pytest.approx(np.log(10.0)) and info['miss'] > 0.5


def test_bounds_and_penalty():
    zc, pen = F.clip_coords(('aplpn', 'kckc'), [0.5, -0.1])
    assert zc == (F.BOUNDS['aplpn'][1], -0.1) and pen == pytest.approx((0.5 - F.BOUNDS['aplpn'][1]) / 0.1)
    assert F.clip_coords(('aplpn',), [-1.0]) == ((-1.0,), 0.0)
    assert F.at_bound(('aplpn', 'claw_thr'), [F.BOUNDS['aplpn'][1] - 0.01, 0.01]) == ['aplpn high', 'claw_thr low']
    assert F.at_bound(('aplpn', 'claw_thr'), [-1.0, 0.6]) == []


def test_permutations():
    od = list(ODORS)
    ps = F.permutations(od)
    assert len(ps) == F.N_PERMS and ps == F.permutations(od)
    for p in ps:
        assert sorted(p) == sorted(od) and sorted(p.values()) == sorted(od) and all(p[o] != o for o in od)
    pop = P.Population(pd.DataFrame(np.arange(12.0).reshape(2, 6), columns=od), pd.DataFrame(np.ones((2, 6)), columns=od))
    q = F.permute_population(pop, ps[0])
    assert list(q.values.columns) == od
    for o in od:
        assert np.array_equal(q.values[ps[0][o]].to_numpy(), pop.values[o].to_numpy())


def test_halves_and_references():
    full = planted_full()
    fl = F.Flies('fit', F.LOSS_STAGES, full=full)
    for st in F.LOSS_STAGES:
        names = fl.names[st]
        assert len(names['fit']) == 2 and len(names['held']) == 2 and not set(names['fit']) & set(names['held'])
        assert list(fl.cache.get((st, F.PANEL, 'on')).flies) == names['fit']
        ref = fl.reference(st, list(ODORS))
        allref = F.Flies(None, F.LOSS_STAGES, full=full).reference(st, list(ODORS))
        assert ref.flies == names['fit'] and allref.flies == [str(f) for f in full[st].flies]
        assert ref.one_fly == allref.one_fly              # the half's targets in all flies' one-fly units
        assert fl.reference(st, list(ODORS)) is ref       # cached
    assert fl.cache.get(('kc_soma', F.PANEL, 'off')) is None


def test_claw_inputs_divisive():
    rng = np.random.default_rng(0)
    W = (rng.random((5, 3)) < 0.4).astype(float)
    b = rng.gamma(1.0, 10.0, (3, 7))
    inh = rng.normal(1.0, 1.0, 7)
    wk = rng.random(5)
    x = F.claw_inputs_divisive(W, b, inh, wk)
    expect = np.maximum((W @ b) / (1.0 + np.maximum(np.outer(wk, inh), 0.0)), 0.0)
    assert np.array_equal(x, expect)
    assert np.array_equal(F.claw_inputs_divisive(W, b, np.zeros(7), wk), W @ b)
    after = F.claw_reader_after_threshold(5.0, divisive=True)(W, b, inh, wk)
    assert np.array_equal(after, np.maximum(expect - 5.0, 0.0))


def _out_and_built(seed=0):
    run = planted_run(seed)
    out = {}
    for c in R.CONDITIONS:
        st = dict(run.conditions[c])
        for k in R.STAGE_FILES:
            st[f'{k}_pre_settled'] = st[f'{k}_pre'] * 0.9
        st['apl'] = pd.DataFrame({'win': np.linspace(1, 2, len(ODORS)), 'pre': 0.1, 'pre_settled': 0.05},
                                 index=list(ODORS))
        out[c] = st
    built = types.SimpleNamespace(meta=run.meta, orn=run.orn)
    return out, built


def test_in_memory_run_reads_as_a_saved_run(tmp_path):
    out, built = _out_and_built()
    info = dict(model='C2', wiring='real', panel='megamat', reading='fit', tag='t_real')
    mem = F.in_memory_run('t_real', out, built, info)
    R.save_step3_run('t_real', out, built.meta, built.orn, info, [tmp_path])
    saved = R.load_step3_run('t_real', dirs=[tmp_path])
    assert mem.stages() == saved.stages()
    for c in R.CONDITIONS:
        assert set(mem.conditions[c]) == set(saved.conditions[c])
        for k, x in mem.conditions[c].items():
            y = saved.conditions[c][k]
            assert np.array_equal(x.to_numpy(dtype=float), y.to_numpy(dtype=float)), (c, k)
            assert list(map(str, x.columns)) == list(map(str, y.columns)), (c, k)
    flies = F.Flies(None, F.LOSS_STAGES, full=planted_full()).cache
    a = R.read_run(mem, stages=F.LOSS_STAGES, lock=LOCK, flies=flies, seed_index=0, subset_glomeruli=GL)
    b = R.read_run(saved, stages=F.LOSS_STAGES, lock=LOCK, flies=flies, seed_index=0, subset_glomeruli=GL)
    assert set(a.populations) == set(b.populations) and a.populations
    for k, p in a.populations.items():
        pd.testing.assert_frame_equal(p.values, b.populations[k].values)
        pd.testing.assert_frame_equal(p.responds, b.populations[k].responds)


def _planted_runs():
    runs = {}
    for i, w in enumerate(('real', 'P1')):
        out, built = _out_and_built(seed=i)
        runs[w] = F.in_memory_run(w, out, built, dict(model='C2', wiring=w, panel='megamat'))
    return runs


def test_stage_loss_on_planted_runs():
    full = planted_full()
    runs = _planted_runs()
    loss = F.StageLoss(F.LOSS_STAGES, F.Flies('fit', F.LOSS_STAGES, full=full), LOCK, subset_glomeruli=GL, check_pn=True)
    res = loss(runs)
    errs = [res[f'{s}_error'] for s in F.LOSS_STAGES]
    assert all(np.isfinite(errs))
    assert res['stage_error'] == pytest.approx(float(np.sqrt(np.mean(np.square(errs)))))
    assert res['stage_error_downstream'] == pytest.approx(float(np.sqrt(np.mean(np.square(errs[1:])))))
    assert res['loss_stages'] == res['stage_error'] and 'stage_error_perm' not in res
    assert 0 < res['rate_on'] < 1 and res['rate_off'] > res['rate_on']
    assert res == loss(runs)                       # deterministic, and the PN-dendrite check passes
    # the control: the mean over the permutations, the recorded labels beside
    ctl = F.StageLoss(F.LOSS_STAGES, F.Flies('fit', F.LOSS_STAGES, perms=F.permutations(list(ODORS), 3), full=full),
                      LOCK, subset_glomeruli=GL)(runs)
    assert ctl['stage_error'] == pytest.approx(res['stage_error'])
    assert ctl['loss_stages'] == ctl['stage_error_perm'] != ctl['stage_error']
    with pytest.raises(ValueError):
        F.StageLoss(('orn',), F.Flies(None, full=full), LOCK)
    # the hold penalty: one unit per 10 % of relative miss of a hold that did not reach its target
    runs['real'].info['hold'] = dict(on=dict(ok=False, miss=0.05), off=dict(ok=True, miss=0.001))
    assert F.hold_penalty(runs) == pytest.approx(0.5)


def test_readable_stages():
    c = {'stages': {'megamat_held': {'pn_dendrite': {'cls': 'between'}, 'pn_bouton': {'cls': 'not readable'},
                                     'kc_claw': {'cls': 'far'}, 'kc_soma': {'cls': 'inside'}}}}
    assert F.readable_stages(c) == ('pn_dendrite', 'kc_claw', 'kc_soma')
    with pytest.raises(RuntimeError):
        F.readable_stages({'stages': {}})


def test_fit_refuses_before_the_gate(monkeypatch, tmp_path):
    monkeypatch.setattr(F, 'e134_dir', lambda: tmp_path)
    with pytest.raises(RuntimeError, match='measure gate'):
        F.require_gate()
    (tmp_path / F.GATE_FILE).write_text('{"verdict": "not usable"}')
    with pytest.raises(RuntimeError, match='not .usable'):
        F.require_gate()


def test_planned_reads_and_tags():
    pts = F.planned_reads()
    assert len(pts) == len(set(pts)) == 8 + 28 * len(F.ARMS)
    fam = [pt for pt in pts if pt.arm == F.FAMILY]
    assert len(fam) == 8 and {pt.panel for pt in fam} == {'validation2'}
    a0 = [pt for pt in pts if pt.arm == 'A0']
    assert [pt.wiring for pt in a0 if pt.panel == F.PANEL] == list(F.T_WIRINGS)
    assert all(pt.half == 'held' for pt in a0 if pt.panel == F.PANEL)
    assert {pt.panel for pt in a0} == set(F.READ_PANELS)
    assert len({F.read_tag(pt) for pt in pts}) == len(pts)
    assert F.read_tag(F.ReadPoint('A2', 'apl16', 'M3', None)) == 'A2_apl16_M3'
    with pytest.raises(ValueError):
        F.read_tag(F.ReadPoint('A9', 'apl16', 'M3', None))
    assert F.fitted_z(F.FAMILY) == {}
    assert F.planned_reads(['A4']) == [pt for pt in pts if pt.arm == 'A4']


@pytest.mark.skipif(os.environ.get('E136_REAL') != '1', reason='builds the family on the real wiring (E136_REAL=1)')
def test_claw_identity_divisive():
    r = F.claw_identity('A2')
    assert r['exact'], r
