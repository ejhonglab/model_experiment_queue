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
    assert F.arm_kws('A1') == {}


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
    # far from the start on a curved count (a saturating response curve, as A4's first hold): the wider budget reaches it
    def curved(x):
        n = int(round(40000 / (1.0 + np.exp(2.0 * (x - 0.3)) ** 1.5)))
        return pd.DataFrame(np.r_[np.ones(n), np.zeros(40000 - n)].reshape(-1, 1)), None
    x, _, info = F.secant_hold(curved, 2500, 0.0, np.log(0.1), np.log(10.0), -3.0, F.HOLD_FIRST_MAX_SIMS)
    assert info['ok'], info
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


def test_strength_split():
    """The strength split adds back to each strength part's error (level^2 + spread^2 + order^2 = error^2), for the
    wirings' mean views and for each wiring, and changes nothing else of the loss."""
    full = planted_full()
    runs = _planted_runs()
    flies = F.Flies('fit', F.LOSS_STAGES, full=full)
    plain = F.StageLoss(F.LOSS_STAGES, flies, LOCK, subset_glomeruli=GL)(runs)
    res = F.StageLoss(F.LOSS_STAGES, flies, LOCK, subset_glomeruli=GL, split=True)(runs)
    assert {k: v for k, v in res.items() if '_split_' not in k} == plain
    n = 0
    for st in F.LOSS_STAGES:
        for part in ('strength_rate', 'strength_mean'):
            if f'{st}_part_{part}' not in res:
                continue
            for suffix, err in [('', res[f'{st}_part_{part}'])] + \
                    [(f'_{w}', None) for w in runs]:
                lv, sp, od = (res[f'{st}_split_{part}_{x}{suffix}'] for x in ('level', 'spread', 'order'))
                if err is not None:
                    assert lv ** 2 + sp ** 2 + od ** 2 == pytest.approx(err ** 2, rel=1e-9)
                assert min(lv, sp, od) >= 0 and res[f'{st}_split_{part}_sd_ratio{suffix}'] > 0
                n += 1
    assert n >= len(F.LOSS_STAGES) * (1 + len(runs))
    # each wiring's split against its own error, from hal's wirings_as_flies
    from kcq import population_fit as P
    reads = F.StageLoss(F.LOSS_STAGES, flies, LOCK, subset_glomeruli=GL).read(runs)
    st = 'kc_soma'
    pops = [reads[w].populations[(st, 'megamat', 'on')] for w in runs]
    ref = flies.reference(st, list(pops[0].values.columns))
    tab = P.wirings_as_flies(pops, ref)
    sp = F.strength_split(pops, ref)
    for name in [*[p.name for p in pops], 'common']:
        e = tab.loc[name, 'part_strength_rate']
        s = sp.loc[name]
        assert s.strength_rate_level ** 2 + s.strength_rate_spread ** 2 + s.strength_rate_order ** 2 == \
            pytest.approx(e ** 2, rel=1e-9)
    assert sp.loc['mean', 'strength_rate_order'] == pytest.approx(sp.loc[[p.name for p in pops], 'strength_rate_order'].mean())


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
    assert len(pts) == len(set(pts)) == 8 + 28 * len(F.ARMS) + 4
    ctl = [pt for pt in pts if pt.arm == F.CONTROL]
    assert [(pt.panel, pt.wiring, pt.half) for pt in ctl] == [(F.PANEL, w, 'held') for w in F.T_WIRINGS]
    with pytest.raises(ValueError):
        F.read_tag(F.ReadPoint(F.CONTROL, 'apl16', 'P3', None))
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


def _mega(stage_by_wiring, stages=('kc_claw', 'kc_soma')):
    return {w: dict(stage=v, downstream=v, **{s: v for s in stages}) for w, v in stage_by_wiring.items()}


def test_classes():
    stages = ('kc_claw', 'kc_soma')
    fam = {w: 1.0 for w in F.WIRINGS}
    spreads = dict(stage=0.04, downstream=0.04, kc_claw=0.04, kc_soma=0.04, **{'/'.join(k): 0.05 for k in F.HOLDS_PANELS})
    panels_fam = {k: {w: 1.0 for w in F.WIRINGS} for k in F.HOLDS_PANELS}
    base = dict(mega={F.FAMILY: _mega(fam, stages)}, panels={F.FAMILY: panels_fam}, spreads=spreads,
                realism={F.FAMILY: {}, 'A0': {}}, fit={}, control_fall=0.01)

    def run(arm_vals, panel_shift=0.0, **over):
        kw = {k: (dict(v) if isinstance(v, dict) else v) for k, v in base.items()}
        kw['mega'][ 'A0'] = _mega(arm_vals, stages)
        kw['panels']['A0'] = {k: {w: 1.0 + panel_shift for w in F.WIRINGS} for k in F.HOLDS_PANELS}
        kw.update(over)
        return F.classify('A0', stages, **kw)
    better = {w: 0.9 for w in F.T_WIRINGS}           # a fall of 0.1 against a standard error of 0.02
    assert run(better)['cls'] == 'holds'
    assert run(better, panel_shift=0.1)['cls'] == 'megamat only'
    assert run({w: 0.99 for w in F.T_WIRINGS})['cls'] == 'does not improve'
    assert run(better, fit=dict(at_bound=['aplpn high']))['cls'] == 'no class'
    assert run(better, control_fall=0.2)['cls'] == 'no class'
    r = run(better, realism={F.FAMILY: {}, 'A0': {'P3': 'APL-on response rate 0.200'}})
    assert r['cls'] == 'no class' and 'leaves a realism band on P3' in r['reasons']
    # a fall carried by one wiring does not pass: the SD of the changes over T sets the standard error
    assert run({'P3': 0.6, 'P4': 1.0, 'M2': 1.0, 'M3': 1.0})['cls'] == 'does not improve'
    c = F.compare({w: 0.9 for w in F.T_WIRINGS}, {w: 1.0 for w in F.T_WIRINGS}, 0.04)
    assert c['change'] == pytest.approx(-0.1) and c['se'] == pytest.approx(0.02) and c['improves']


def test_measures():
    assert F.MEASURE == 'all' and F.measure_parts(None) is None
    assert F.measure_parts(('strength',)) == ['strength/rate', 'strength/mean']
    assert F.error_column('pn_dendrite', None) == 'error_with_matrix'
    assert F.error_column('pn_dendrite', ('strength',)) == 'error' == F.error_column('kc_soma', None)
    full = planted_full()
    runs = _planted_runs()
    try:
        F.set_measure('strength')
        fl = F.Flies('fit', F.LOSS_STAGES, full=full)
        assert fl.views == ('strength',)
        assert set(fl.reference('kc_soma', list(ODORS)).parts) == {'strength/rate', 'strength/mean'}
        res = F.StageLoss(F.LOSS_STAGES, fl, LOCK, subset_glomeruli=GL)(runs)
        assert res['kc_soma_error'] == pytest.approx(res['kc_soma_view_strength'])
        assert 'kc_soma_view_sharing' not in res
        assert F._measure_sub() == 'measure_strength'
    finally:
        F.set_measure('all')
    with pytest.raises(ValueError):
        F.set_measure('sharing_only')


def test_gate_on_the_adopted_measure(monkeypatch, tmp_path):
    import json
    monkeypatch.setattr(F, 'e134_dir', lambda: tmp_path)
    (tmp_path / 'score_classes.json').write_text('{"stages": {}}')
    (tmp_path / 'score_stage_table.csv').write_text('x\n')
    (tmp_path / F.GATE_FILE).write_text(json.dumps(dict(verdict='usable', verdict_as_read='not usable',
                                                        measure=dict(views=['strength'], adopted=True))))
    try:
        with pytest.raises(RuntimeError, match='--measure strength'):
            F.require_gate()
        F.set_measure('strength_sharing')
        with pytest.raises(RuntimeError, match='not .strength_sharing'):
            F.require_gate()
        F.set_measure('strength')
        assert F.require_gate()['gate']['measure']['views'] == ['strength']
        (tmp_path / F.GATE_FILE).write_text(json.dumps(dict(verdict='usable',
                                                            measure=dict(views=['strength'], adopted=False))))
        with pytest.raises(RuntimeError, match='not adopted'):
            F.require_gate()
    finally:
        F.set_measure('all')
    assert set(F.START_CHECK['strength']) == {'stage_error', *(f'{s}_error' for s in F.LOSS_STAGES)}


def test_archive_best(monkeypatch, tmp_path):
    loc, arch = tmp_path / 'local', tmp_path / 'arch'
    monkeypatch.setattr(F, 'local_dir', lambda part='': (loc / part).mkdir(parents=True, exist_ok=True) or loc / part)
    monkeypatch.setattr(F, 'archive_dir', lambda: arch)
    for w in F.S_WIRINGS:
        for suf in ('on_spike_counts.parquet', 'run.json'):
            (loc / 'fit_A0').mkdir(parents=True, exist_ok=True)
            (loc / 'fit_A0' / f'A0_e3_{w}_{suf}').write_text(w)
            (loc / 'fit_A0' / f'A0_e4_{w}_{suf}').write_text('other')
    out = F.archive_best('A0', 3)
    assert len(out) == 8 and sorted(p.name for p in (arch / 'fit_A0').iterdir()) == sorted(p.name for p in out)
    assert (arch / 'fit_A0' / 'A0_e3_real_run.json').read_text() == 'real'
    assert F.fit_dirs('A0') == [loc / 'fit_A0', arch / 'fit_A0']
    with pytest.raises(FileNotFoundError):
        F.archive_best('A0', 9)


def test_archive_across_hosts(monkeypatch, tmp_path):
    """A fit resumed on another host: the first host archives its best (`archive_evals`), and the second's
    `archive_best` takes the archived evaluation it never ran."""
    loc, arch, res = tmp_path / 'local', tmp_path / 'arch', tmp_path / 'res'
    res.mkdir()
    monkeypatch.setattr(F, 'local_dir', lambda part='': (loc / part).mkdir(parents=True, exist_ok=True) or loc / part)
    monkeypatch.setattr(F, 'archive_dir', lambda: arch)
    monkeypatch.setattr(F, 'out_dir', lambda: res)
    pd.DataFrame(dict(eval=[0, 1, 2, 3], loss_sim=[0.97, 0.95, 0.96, 0.90],
                      host=['first', 'first', 'first', 'second'])).to_csv(res / 'fit_history_A3.csv', index=False)
    (loc / 'fit_A3').mkdir(parents=True)
    for i in (0, 1, 2):
        for w in F.S_WIRINGS:
            for suf in ('on_spike_counts.parquet', 'run.json'):
                (loc / 'fit_A3' / f'A3_e{i}_{w}_{suf}').write_text(f'{i}{w}')
    monkeypatch.setattr(F, '_host', lambda: 'third')
    with pytest.raises(ValueError):
        F.archive_evals('A3')                       # 'third' ran none
    monkeypatch.setattr(F, '_host', lambda: 'second')
    with pytest.raises(FileNotFoundError):
        F.archive_evals('A3')                       # 'second' ran eval 3, whose files are not here
    monkeypatch.setattr(F, '_host', lambda: 'first')
    assert F.archive_evals('A3') == [1]             # the best of the first host's evaluations
    assert sorted(p.name for p in (arch / 'fit_A3').iterdir()) == sorted(
        f'A3_e1_{w}_{suf}' for w in F.S_WIRINGS for suf in ('on_spike_counts.parquet', 'run.json'))
    for p in list((loc / 'fit_A3').iterdir()):      # the second host has none of the first's files
        p.unlink()
    out = F.archive_best('A3', 1)
    assert len(out) == 8 and (arch / 'fit_A3' / 'A3_e1_real_run.json').read_text() == '1real'
    with pytest.raises(FileNotFoundError, match='--part archive'):
        F.archive_best('A3', 2)


def test_guard_on_the_views_beside():
    stages = ('kc_claw', 'kc_soma')
    fam = {w: 1.0 for w in F.WIRINGS}
    spreads = dict(stage=0.04, downstream=0.04, kc_claw=0.04, kc_soma=0.04, **{'/'.join(k): 0.05 for k in F.HOLDS_PANELS})
    panels = {k: {w: 1.0 for w in F.WIRINGS} for k in F.HOLDS_PANELS}
    better = {w: 0.9 for w in F.T_WIRINGS}

    def run(guard_soma):
        g_mega = {F.FAMILY: _mega(fam, stages), 'A0': _mega({w: 1.0 for w in F.T_WIRINGS}, stages)}
        for w in F.T_WIRINGS:
            g_mega['A0'][w]['kc_soma'] = guard_soma
        guard = dict(mega=g_mega, panels={F.FAMILY: panels, 'A0': panels},
                     spreads=dict(stage=0.04, kc_claw=0.04, kc_soma=0.05, **{'/'.join(k): 0.05 for k in F.HOLDS_PANELS}))
        return F.classify('A0', stages, {F.FAMILY: _mega(fam, stages), 'A0': _mega(better, stages)},
                          {F.FAMILY: panels, 'A0': panels}, spreads, {F.FAMILY: {}, 'A0': {}}, {}, 0.01, guard)
    assert run(1.0)['cls'] == 'holds'
    r = run(1.2)                                     # the somata's error over all views up by 0.2, its spread 0.05
    assert r['cls'] == "strength at another view's cost" and r['guard_worse'] == ['kc_soma'] and not r['improves']
    # a fit's classed point's bounds stand in for the best point's
    c = F.classify('A0', stages, {F.FAMILY: _mega(fam, stages), 'A0': _mega(better, stages)},
                   {F.FAMILY: panels, 'A0': panels}, spreads, {F.FAMILY: {}, 'A0': {}},
                   dict(at_bound=['aplpn high'], at_bound_classed=[]), 0.01)
    assert c['cls'] == 'holds'


def test_class_ignores_kmix():
    # the kiwi mixtures are scored beside the class (Tom, 2026-11-25 10:42), the control mixtures stay in it
    assert ('mixtures', 'kmix', 'on') not in F.HOLDS_PANELS and ('mixtures', 'kmix', 'on') in F.BESIDE_PANELS
    assert ('mixtures', 'cmix', 'on') in F.HOLDS_PANELS
    stages = ('kc_claw', 'kc_soma')
    fam = {w: 1.0 for w in F.WIRINGS}
    spreads = dict(stage=0.04, downstream=0.04, kc_claw=0.04, kc_soma=0.04,
                   **{'/'.join(k): 0.05 for k in F.SCORED_PANELS})
    flat = {k: {w: 1.0 for w in F.WIRINGS} for k in F.SCORED_PANELS}
    better = {w: 0.9 for w in F.T_WIRINGS}

    def run(worse_key, in_guard=False):
        arm_p = {k: {w: 1.0 + (0.2 if k == worse_key else 0.0) for w in F.WIRINGS} for k in F.SCORED_PANELS}
        guard = None
        if in_guard:
            g_mega = {F.FAMILY: _mega(fam, stages), 'A0': _mega({w: 1.0 for w in F.T_WIRINGS}, stages)}
            guard = dict(mega=g_mega, panels={F.FAMILY: flat, 'A0': arm_p}, spreads=spreads)
        p = {F.FAMILY: flat, 'A0': flat if in_guard else arm_p}
        return F.classify('A0', stages, {F.FAMILY: _mega(fam, stages), 'A0': _mega(better, stages)}, p, spreads,
                          {F.FAMILY: {}, 'A0': {}}, {}, 0.01, guard)
    for in_guard in (False, True):
        r = run(('mixtures', 'kmix', 'on'), in_guard)
        assert r['cls'] == 'holds' and r['panel_worse'] == [] and 'mixtures/kmix/on' not in r['panels']
        if in_guard:
            assert r['guard']['mixtures/kmix/on']['worse_than_spread']
        else:
            assert r['panels_beside']['mixtures/kmix/on']['worse_than_spread']
        r = run(('mixtures', 'cmix', 'on'), in_guard)
        assert r['cls'] == 'megamat only' and r['panel_worse'] == ['mixtures/cmix/on']


def test_e137_beside_compare():
    rows = []
    for arm, v in ((F.FAMILY, 1.0), ('A0', 0.9)):
        for i, w in enumerate(F.WIRINGS):
            rows.append(dict(arm=arm, panel='kmix', reading='components_first_block', views='strength', wiring=w,
                             error=v + 0.01 * (i % 2)))
        rows.append(dict(arm=arm, panel='kmix', reading='components_first_block', views='strength', wiring='mean',
                         error=v))
    out = F.e137_beside(pd.DataFrame(rows), ['A0'])
    c = out['A0']['kmix/components_first_block/strength']
    assert c['complete'] and c['change'] == pytest.approx(-0.1)
    assert c['spread'] == pytest.approx(np.std([0.01 * (i % 2) for i in range(8)], ddof=1))
    assert F.e137_beside(pd.DataFrame(), ['A0']) == {}


def test_save_eval_tags(monkeypatch, tmp_path):
    monkeypatch.setattr(F, 'local_dir', lambda part='': (tmp_path / part).mkdir(parents=True, exist_ok=True) or tmp_path / part)
    out, built = _out_and_built()
    run = F.in_memory_run('x', out, built, dict(model='C2', wiring='real', panel='megamat'))
    run.info['_out'] = out
    F.save_eval('A0', 3, {'real': run})
    F.save_eval('A0', 'det_aplpn0', {'real': run})
    names = {p.name for p in (tmp_path / 'fit_A0').iterdir()}
    assert 'A0_e3_real_run.json' in names and 'A0_det_aplpn0_real_run.json' in names


def test_consistency(monkeypatch, tmp_path):
    monkeypatch.setattr(F, 'out_dir', lambda: tmp_path)
    rng = np.random.default_rng(0)
    x = rng.random(30)
    rows = []
    for i, v in enumerate(x):
        for half, e in (('fit', v), ('held', v + 0.01 * rng.standard_normal())):
            rows.append(dict(arm='A0', eval=i, measure='strength', half=half, stage_error=e,
                             **{f'{s}_error': e for s in F.LOSS_STAGES}))
    pd.DataFrame(rows).to_csv(tmp_path / 'rescore_A0.csv', index=False)
    c = F.consistency(['A0'], 'strength')
    assert len(c) == 5 and (c.rho > 0.9).all() and (c.n_evals == 30).all()


def test_rescore_merges_hosts(monkeypatch, tmp_path):
    monkeypatch.setattr(F, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(R, 'require_lock', lambda: {})
    monkeypatch.setattr(F, 'Flies', lambda half, stages, views=None: half)
    monkeypatch.setattr(F, 'StageLoss', lambda stages, flies, lock: (lambda runs: dict(stage_error=next(iter(runs.values())) + (
        0.5 if flies == 'held' else 0.0))))
    saved = {}

    def load(tag, dirs=None):
        i = int(tag.split('_e')[1].split('_')[0])
        if i not in saved:
            raise FileNotFoundError(tag)
        return saved[i]
    monkeypatch.setattr(R, 'load_step3_run', load)
    saved.update({i: float(i) for i in (1, 2, 3)})                 # host 1 has evaluations 1 to 3
    F.rescore('A0', measures=['strength'], evals=range(6))
    saved.clear()
    saved.update({i: 10.0 + i for i in (3, 4, 5)})                 # host 2 has 3 (again) to 5
    df = F.rescore('A0', measures=['strength'], evals=range(6))
    on_disk = pd.read_csv(tmp_path / 'rescore_A0.csv')
    assert len(df) == len(on_disk) == 10
    fit = on_disk[on_disk.half == 'fit'].set_index('eval').stage_error
    assert fit.to_dict() == {1: 1.0, 2: 2.0, 3: 13.0, 4: 14.0, 5: 15.0}
    held = on_disk[on_disk.half == 'held'].set_index('eval').stage_error
    assert held[2] == 2.5 and held[5] == 15.5


def _verdict(monkeypatch, tmp_path, outcome, **kw):
    import json
    from experiments import E135_antennal_lobe_dendrites as AL
    monkeypatch.setattr(AL, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(AL, 'read_step1_fits', lambda root=None: {'step2_rung': 'R3', 'highest_kept': 'R3',
                                                                   'all_flies': {}})
    if outcome is not None:
        (tmp_path / AL.SCORE_VERDICT_FILE).write_text(json.dumps({'cls': 'x', 'outcome': dict(outcome=outcome, **kw)}))


def test_require_e135(monkeypatch, tmp_path):
    _verdict(monkeypatch, tmp_path, None)
    with pytest.raises(RuntimeError, match='waits on E135'):
        F.require_e135()
    with pytest.raises(RuntimeError, match='waits on E135'):
        F.arm_searched('A1')
    assert F.arm_searched('A0') == ('aplpn', 'kckc') and F.arm_searched('A3') == ('aplpn', 'kckc', 'claw_thr')
    for o in ('incomplete', 'needs the dF/F repeat', 'unit-dependent', 'reading-dependent'):
        _verdict(monkeypatch, tmp_path, o)
        with pytest.raises(RuntimeError, match='waits on E135'):
            F.require_e135()
    _verdict(monkeypatch, tmp_path, 'a', improvement=True)
    h = F.require_e135()
    assert (h['mode'], h['rung'], h['outcome'], h['improvement']) == ('table', 'R1', 'a', True)
    assert F.arm_searched('A1') == ('aplpn', 'kckc')
    _verdict(monkeypatch, tmp_path, 'b', improvement=False, at_a_cost_beside=['kc_soma breadth'])
    h = F.require_e135()
    assert (h['mode'], h['rung'], h['at_a_cost_beside']) == ('table', 'R3', ['kc_soma breadth'])
    for o in ('c', 'd', 'e', 'f'):
        _verdict(monkeypatch, tmp_path, o)
        h = F.require_e135()
        assert (h['mode'], h['rung']) == ('search', None)
        searched = F.arm_searched('A1')
        assert searched == ('aplpn', 'kckc', 'inhadd', 'inhsc')
        assert len(searched) + len(F.HELD) <= F.MAX_FREE and set(searched) <= set(F.BOUNDS) <= set(F.FIRST_STEP)


def test_restart_point(monkeypatch, tmp_path):
    """The restart starts from the best when the first run stopped, and a resumed fit starts it there again although
    the restart has since found a better point; a resume with no record takes the history's best."""
    monkeypatch.setattr(F, 'out_dir', lambda: tmp_path)
    hist = [dict(eval=i, loss_sim=v) for i, v in enumerate((0.97, 0.96, 0.965))]
    assert F.restart_point('A3', hist, resumed=False)['eval'] == 1
    rec = F._read_json(tmp_path / 'fit_restart_A3.json')
    assert (rec['from_eval'], rec['n_evals']) == (1, 3)
    hist += [dict(eval=3, loss_sim=0.95)]              # found by the restart
    assert F.restart_point('A3', hist, resumed=True)['eval'] == 1
    (tmp_path / 'fit_restart_A3.json').unlink()
    assert F.restart_point('A3', hist, resumed=True)['eval'] == 3
    F._write_json(tmp_path / 'fit_restart_A3.json', dict(from_eval=7))
    with pytest.raises(RuntimeError, match='not in the history'):
        F.restart_point('A3', hist, resumed=True)


def test_pn_norm_kws():
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    from experiments.E77_input_screen import inhsc_for
    assert F.pn_norm_kws({}) == {} and F.pn_norm_kws(dict(aplpn=-1.0, inhadd=0.0, inhsc=0.0)) == {}
    k = F.pn_norm_kws(dict(inhadd=F.FIRST_STEP['inhadd']))
    assert k['pn_inhadd'] == pytest.approx(100.0) and k['pn_inhsc'] == pytest.approx(inhsc_for(100.0))
    k = F.pn_norm_kws(dict(inhsc=np.log(2.0)))
    assert k['pn_inhadd'] == pytest.approx(PN_INHADD_DEFAULT) and k['pn_inhsc'] == pytest.approx(2 * inhsc_for(PN_INHADD_DEFAULT))
    lo, hi = F.BOUNDS['inhadd']
    assert PN_INHADD_DEFAULT * np.exp(hi) > 10000.0 and lo < 0.0 < F.FIRST_STEP['inhadd'] < hi
    assert F.coord_value('claw_thr', 1.0) == F.CLAW_THR_UNIT and F.coord_value('inhadd', np.log(3.0)) == pytest.approx(3.0)
    assert F.coord_value('aplpn', np.log(0.5)) == pytest.approx(0.5)


@pytest.mark.skipif(os.environ.get('E136_REAL') != '1', reason='builds the family on the real wiring (E136_REAL=1)')
@pytest.mark.parametrize('outcome,panel', [('a', 'megamat'), ('b', 'validation2')])
def test_a1_table_real(monkeypatch, tmp_path, outcome, panel):
    """A1's kept table on the real wiring (E135's step-1 fits from the shared directory) entering by route 1: the built
    model's PN rates equal the table on every glomerulus and odor of the input."""
    import json
    from experiments import E134_population_baseline as B
    from experiments import E135_antennal_lobe_dendrites as AL
    fits = AL.read_step1_fits()
    monkeypatch.setattr(AL, 'read_step1_fits', lambda root=None: fits)
    real_out = AL.out_dir()
    monkeypatch.setattr(F, 'e135_verdict_path', lambda: tmp_path / 'v.json')
    (tmp_path / 'v.json').write_text(json.dumps({'outcome': dict(outcome=outcome)}))
    assert real_out.exists()
    F._A1_STAGES.clear()
    orn = B.panel_input(panel)
    kws = F.a1_kws('real', panel, orn)
    table, stage, rec = F.a1_stage('real', panel, orn)
    assert rec['rung'] == ('R1' if outcome == 'a' else AL.step2_rung_of(fits).replace('m', '').replace('f', ''))
    assert list(kws['pn_rates'].columns) == list(orn.columns) and np.isfinite(kws['pn_rates'].to_numpy()).all()
    p, _ = F.start_params('real')
    built = B.build_family('real', panel, dict(p), orn=orn, extra_kws=kws)
    built.sim(dict(p))
    chk = F.a1_pn_check(built, 'real', panel)
    assert chk['ok'] and chk['n_set'] == table.size



def test_write_variants(monkeypatch, tmp_path):
    """variants.csv: one row per run the archive holds, from its run.json, both conditions' parameters and A1's PN
    terms (a read's at the top level, a fit's in pn_kws)."""
    import json
    arch, res = tmp_path / 'arch', tmp_path / 'res'
    res.mkdir()
    monkeypatch.setattr(F, 'out_dir', lambda: res)
    on = dict(fixed_thr=49.8, wAPLKC_scale=1.08, wKCAPL_scale=6e-4, wAPLPN_scale=0.056, wPNAPL_scale=3e-3,
              wKCKC_scale=138.3)
    off = dict(on, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    runs = {'read/A1_megamat_P3_run.json': dict(model='C2', wiring='P3', panel='megamat', reading='read', arm='A1',
                                                z=dict(aplpn=-1.86, inhadd=3.48), params=on, params_off=off,
                                                claw_thr=None, tag='A1_megamat_P3', host='viki', lane='viki',
                                                pn_inhadd=1017.3, pn_inhsc=598.9),
            'fit_A1/A1_e122_real_run.json': dict(model='C2', wiring='real', panel='megamat', reading='fit', arm='A1',
                                                 z=dict(aplpn=-1.86, inhadd=3.48), params=on, params_off=off,
                                                 claw_thr=None, tag='A1_e122_real', eval=122, host='viki',
                                                 pn_kws=dict(pn_inhadd=1017.3, pn_inhsc=598.9))}
    for rel, d in runs.items():
        (arch / rel).parent.mkdir(parents=True, exist_ok=True)
        (arch / rel).write_text(json.dumps(d))
    df = F.write_variants(arch).set_index('tag')
    assert list(df.index) == ['A1_e122_real', 'A1_megamat_P3']
    assert df.loc['A1_e122_real', 'subdir'] == 'fit_A1' and df.loc['A1_e122_real', 'eval'] == 122
    assert (df.pn_inhadd == 1017.3).all() and (df.off_wAPLKC_scale == 0).all() and (df.wKCKC_scale == 138.3).all()
    assert df.loc['A1_megamat_P3', 'z_inhadd'] == 3.48 and (res / 'variants.csv').exists()
