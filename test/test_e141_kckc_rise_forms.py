"""Tests of experiments/E141_kckc_rise_forms.py (the rise of the inhibition among KCs beside the weaker antennal-lobe
normalization, along its dose and in other forms): the arms, their names and coordinates, their E136 arms, the planned
reads and the reads' dispatch. No model is run."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments import E136_stage_population_fit as F
from experiments import E141_kckc_rise_forms as L


def test_arms_and_names():
    assert L.STEP1 == ('N', 'P', 'D1.8', 'D5.8', 'G') and set(L.STEP1) | {'M', 'T', 'S'} == set(L.ARMS)
    assert L.ARMS['D1.8'].kckc == pytest.approx(1.8) and L.ARMS['D5.8'].kckc == pytest.approx(5.832)
    assert L.ARMS['N'].kckc == 3.24 and L.ARMS['N'].kckc_off == 1.0 and L.ARMS['G'].kckc == 1.0
    assert L.variant('D1.8') == 'D1.8_inhadd1000' and L.parse_variant('N_inhadd1000') == 'N'
    for bad in ('B_inhadd1000', 'N_inhadd314', 'N'):
        with pytest.raises(ValueError):
            L.parse_variant(bad)
    assert L.point_tag(L.Point('P', 'apl16', 'M3')) == 'P_inhadd1000_apl16_M3'
    assert all(c.isascii() and (c.isalnum() or c in '._-') for c in ''.join(L.variant(a) for a in L.ARMS))


def test_coordinates_and_terms():
    from experiments.E77_input_screen import inhsc_for
    t = L.pn_terms()
    assert t['pn_inhadd'] == pytest.approx(1000.0) and t['pn_inhsc'] == pytest.approx(inhsc_for(1000.0))
    d = L.arm_z('D5.8')
    assert d['aplpn'] == pytest.approx(np.log(0.155)) and d['kckc'] == pytest.approx(1.5 * np.log(3.24))
    assert np.exp(d['inhadd']) * 31.4088 == pytest.approx(1000.0, rel=1e-3) and d['inhsc'] == 0.0
    # a form's normalization is in its arm's kwargs, not in its coordinates
    assert set(L.arm_z('N')) == {'aplpn', 'kckc'} and L.arm_z('G')['kckc'] == 0.0
    # D1.8 differs from E139's K at 1000 in the KC-KC factor alone
    from experiments import E139_normalization_ladder as E139
    k = E139.point_z('K', 1000.0)
    assert {x: v for x, v in L.arm_z('D1.8').items() if x != 'kckc'} == pytest.approx(
        {x: v for x, v in k.items() if x != 'aplpn' and x != 'kckc'} | {'aplpn': np.log(0.155)})


def test_arm_defs(monkeypatch):
    from experiments import E122_feedback_neuron as E122
    M = pd.DataFrame(np.ones((3, 3)))
    monkeypatch.setattr(E122, 'load_matrix', lambda arm: M * {'U': 1.0, 'TY': 2.0, 'SH': 3.0}[arm])
    n = L.arm_def('N')
    assert n.name == 'N_inhadd1000' and n.searched == ('aplpn', 'kckc')
    assert n.kws['kckc_off_factor'] == pytest.approx(1 / 3.24) and n.kws['pn_inhadd'] == pytest.approx(1000.0)
    assert L.arm_def('G').kws['apl_inh_divisive'] is True
    assert float(L.arm_def('P').kws['_wKCKC'].to_numpy().sum()) == 9.0
    assert float(L.arm_def('T').kws['_wKCKC'].to_numpy().sum()) == 18.0
    with pytest.raises(ValueError):
        L.arm_def('D1.8')
    with pytest.raises(NotImplementedError):
        L.arm_def('M')
    with F.extra_arms(L.arm_def('N'), L.arm_def('G')):
        assert F.arm_kws('N_inhadd1000') == L.pn_terms()
        assert F.off_extra_of('N_inhadd1000', dict(wKCKC_scale=32.4)) == pytest.approx(dict(wKCKC_scale=10.0))
        assert F.claw_reader_of('G_inhadd1000') is F.claw_inputs_divisive and F.claw_reader_of('N_inhadd1000') is None
    assert 'N_inhadd1000' not in F.ARMS


def test_planned_points_and_step2(monkeypatch, tmp_path):
    monkeypatch.setattr(L, 'out_dir', lambda: tmp_path)
    pts = L.planned_points(['1'])
    assert len(pts) == 160 and pts[0] == L.Point('N', 'megamat', F.WIRINGS[0])
    assert [p.arm for p in pts[::32]] == list(L.STEP1) and pts[-1] == L.Point('G', 'apl16', F.WIRINGS[-1])
    assert len(L.planned_points(arms=['P'])) == 32
    with pytest.raises(RuntimeError):
        L.step2_arms()
    (tmp_path / 'classes.json').write_text(json.dumps(dict(classes={'P_inhadd1000': dict(cls='holds')})))
    assert L.step2_arms() == ('M',)
    (tmp_path / 'classes.json').write_text(json.dumps(dict(classes={'P_inhadd1000': dict(cls='megamat only')})))
    assert L.step2_arms() == ('T', 'S') and len(L.planned_points(['2'])) == 64


def test_read_point_dispatch(monkeypatch, tmp_path):
    """A dose is read through E136's A1 at its coordinates; a form through its own arm inside `extra_arms`; both with
    the fields of `read_extra`; a done tag is skipped without a read."""
    calls = []
    monkeypatch.setattr(L, 'out_dir', lambda: tmp_path / 'res')
    monkeypatch.setattr(L, 'read_dirs', lambda: [tmp_path / 'loc'])
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='search'))

    def read_at(tag, arm, z, wiring, panel, half, lane, root, dirs, **kw):
        calls.append(dict(tag=tag, arm=arm, z=z, wiring=wiring, panel=panel, kws=dict(F.arm_kws(arm)) if arm in F.ARMS
                          else None, off=F.off_extra_of(arm, dict(wKCKC_scale=3.24)), **kw))
        return dict(tag=tag)
    monkeypatch.setattr(F, 'read_at', read_at)
    assert L.read_point(L.Point('D1.8', 'validation2', 'P3'), 'x')['tag'] == 'D1.8_inhadd1000_validation2_P3'
    d = calls[-1]
    assert d['arm'] == 'A1' and d['z']['kckc'] == pytest.approx(np.log(1.8)) and d['label'] == 'E141'
    assert d['extra']['variant'] == 'D1.8_inhadd1000' and d['extra']['form_arm'] == 'D1.8' and d['extra']['rung'] == 1000.0
    L.read_point(L.Point('N', 'apl16', 'M2'), 'x')
    d = calls[-1]
    assert d['arm'] == 'N_inhadd1000' and d['kws'] == L.pn_terms() and d['off'] == pytest.approx(dict(wKCKC_scale=1.0))
    assert d['extra']['kckc_off_factor'] == 1.0 and d['extra']['pn_inhadd'] == pytest.approx(1000.0)
    assert 'N_inhadd1000' not in F.ARMS
    n = len(calls)
    assert L.read_point(L.Point('N', 'apl16', 'M2'), 'x', done={'N_inhadd1000_apl16_M2'}) is None and len(calls) == n
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='table'))
    with pytest.raises(RuntimeError):
        L.read_point(L.Point('D5.8', 'megamat', 'real'), 'x')
    with pytest.raises(ValueError):
        L.read_point(L.Point('N', 'extended', 'real'), 'x')


def _mega(stage_by_wiring):
    return {w: dict(stage=v) for w, v in stage_by_wiring.items()}


def test_against_a1_and_entry_classes():
    """The comparisons with A1 on T (megamat stage error, each held-out panel against the family's spread, the pattern
    part against each target) and the entry's classes read off them."""
    T = F.T_WIRINGS
    key = F.HOLDS_PANELS[0]
    k = '/'.join(key)
    spreads = dict(stage=0.01, **{'/'.join(x): 0.05 for x in F.HOLDS_PANELS})
    mega = {'A1': _mega({w: 0.93 for w in T}), 'D1.8_inhadd1000': _mega({w: 0.90 for w in T}),
            'P_inhadd1000': _mega({w: 0.935 for w in T})}
    pv = {'A1': {x: {w: 1.0 for w in T} for x in F.HOLDS_PANELS},
          'D1.8_inhadd1000': {x: {w: (0.9 if x == key else 1.0) for w in T} for x in F.HOLDS_PANELS},
          'P_inhadd1000': {x: {w: (1.1 if x == key else 1.0) for w in T} for x in F.HOLDS_PANELS}}
    parts = {v: {w: {f'{t}_pattern': base for t in ('z8', 'halfA', 'halfB')} for w in T}
             for v, base in (('A1', 0.070), ('D1.8_inhadd1000', 0.065), ('P_inhadd1000', 0.069))}
    d = L.against_a1(F, 'D1.8_inhadd1000', mega, pv, spreads, parts)
    assert d['megamat']['change'] == pytest.approx(-0.03) and d['megamat']['improves']
    assert d['panels'][k]['better'] and not d['panels'][k]['worse'] and d['price']['lower']
    assert d['price']['by_target']['z8'] == pytest.approx(-0.005)
    p = L.against_a1(F, 'P_inhadd1000', mega, pv, spreads, parts)
    assert not p['megamat']['improves'] and p['panels'][k]['worse'] and p['price']['lower'] is False
    assert L.entry_classes('D1.8_inhadd1000', dict(cls='holds'), d) == dict(stands_in=None, does_more_than_a1=True,
                                                                            lower_price=True)
    assert L.entry_classes('D1.8_inhadd1000', dict(cls='megamat only'), d)['does_more_than_a1'] is False
    assert L.entry_classes('P_inhadd1000', dict(cls='holds'), p) == dict(stands_in=True, does_more_than_a1=None,
                                                                         lower_price=False)
    # a dose no better than A1 on megamat, better on one panel and worse on another, does not do more
    pv['D1.8_inhadd1000'][F.HOLDS_PANELS[1]] = {w: 1.2 for w in T}
    mega['D1.8_inhadd1000'] = _mega({w: 0.93 for w in T})
    d2 = L.against_a1(F, 'D1.8_inhadd1000', mega, pv, spreads, parts)
    assert L.entry_classes('D1.8_inhadd1000', dict(cls='holds'), d2)['does_more_than_a1'] is False
    # missing APL16 parts leave the price comparison undecided
    assert L.against_a1(F, 'D1.8_inhadd1000', mega, pv, spreads, {})['price']['lower'] is None


def test_labels():
    assert L.label_of('N_inhadd1000') == dict(setting='N', rung=1000.0)
    assert L.label_of('B_inhadd1000') == dict(setting='B', rung=1000.0)
    lab = L.label_of('A1')
    assert lab['setting'] == 'A1' and np.isnan(lab['rung'])


def test_beside_dispatch_and_lanes(monkeypatch, tmp_path):
    """The beside step: E139's `read_beside` in this experiment's context (a dose through A1, a form inside its arm);
    its arms the doses plus every arm that holds; a lane's job list by steps and arms."""
    from experiments import E139_normalization_ladder as E139
    monkeypatch.setattr(L, 'out_dir', lambda: tmp_path)
    got = []

    def rb(name, wiring, lane, reclaim, done, arm=None, z=None, ctx=None):
        got.append(dict(name=name, wiring=wiring, arm=arm, z=z, ctx=ctx, in_arms=arm in F.ARMS,
                        kws=F.arm_kws(arm) if arm in F.ARMS else None))
    monkeypatch.setattr(E139, 'read_beside', rb)
    L.read_beside('D5.8', 'P3', 'x')
    g = got[-1]
    assert g['name'] == 'D5.8_inhadd1000' and g['arm'] == 'A1' and g['z']['kckc'] == pytest.approx(np.log(3.24 ** 1.5))
    c = g['ctx']
    assert c.root == tmp_path and c.label == 'E141' and c.label_of('D5.8_inhadd1000') == ('D5.8', 1000.0)
    assert c.ext_dir('D5.8_inhadd1000', 1000.0, 'local').parts[-4:] == ('extended', 'D5.8_inhadd1000', 'B',
                                                                         'inhadd1000')
    L.read_beside('N', 'M1', 'x')
    assert got[-1]['arm'] == 'N_inhadd1000' and got[-1]['in_arms'] and got[-1]['kws'] == L.pn_terms()
    assert 'N_inhadd1000' not in F.ARMS
    assert L.beside_arms() == ('D1.8', 'D5.8')
    (tmp_path / 'classes.json').write_text(json.dumps(dict(classes={'N_inhadd1000': dict(cls='holds'),
                                                                    'P_inhadd1000': dict(cls='megamat only')})))
    assert L.beside_arms() == ('D1.8', 'D5.8', 'N') and len(L.beside_points()) == 24
    jobs = []
    monkeypatch.setattr(L, 'read_point', lambda pt, lane, reclaim, done: jobs.append(('point', pt.arm)))
    monkeypatch.setattr(L, 'read_beside', lambda a, w, lane, reclaim, done: jobs.append(('beside', a)))
    monkeypatch.setattr(L, 'reads_table', lambda: pd.DataFrame(columns=['tag']))
    L.part_read('x', ['beside'], arms=['G'])
    assert jobs == [('beside', 'G')] * 8
    jobs.clear()
    L.part_read('x', ['1'], arms=['P'])
    assert jobs == [('point', 'P')] * 32
    jobs.clear()
    L.part_read('x', ['1', 'beside'], reverse=True)
    assert len(jobs) == 160 + 24 and jobs[0] == ('beside', 'N') and jobs[-1] == ('point', 'N')


def test_write_variants(monkeypatch, tmp_path):
    arch, res = tmp_path / 'arch', tmp_path / 'res'
    (arch / 'read').mkdir(parents=True)
    (arch / 'extended/N_inhadd1000/B/inhadd1000').mkdir(parents=True)
    res.mkdir()
    monkeypatch.setattr(L, 'out_dir', lambda: res)
    on = dict(fixed_thr=66.3, wAPLKC_scale=0.37, wKCAPL_scale=2e-4, wAPLPN_scale=0.058, wPNAPL_scale=1e-3,
              wKCKC_scale=85.6)
    off = dict(on, wAPLKC_scale=0.0, wAPLPN_scale=0.0, wKCKC_scale=26.4)
    (arch / 'read' / 'N_inhadd1000_megamat_P3_run.json').write_text(json.dumps(dict(
        model='C2', wiring='P3', panel='megamat', reading='read', arm='N_inhadd1000', z=L.arm_z('N'), params=on,
        params_off=off, tag='N_inhadd1000_megamat_P3', variant='N_inhadd1000', form_arm='N')))
    (arch / 'extended/N_inhadd1000/B/inhadd1000/P3_point.json').write_text(json.dumps(dict(
        tag='N_inhadd1000_extended_P3', variant='N_inhadd1000', setting='N', arm='N_inhadd1000', wiring='P3',
        z=L.arm_z('N'), params=on, params_off=off)))
    df = L.write_variants(arch).set_index(['subdir', 'tag'])
    assert len(df) == 2 and (res / 'variants.csv').exists()
    r = df.loc[('read', 'N_inhadd1000_megamat_P3')]
    assert r.arm == 'N' and r.form == 'on_only' and r.off_wKCKC_scale == 26.4 and r.pn_inhadd == pytest.approx(1000.0)
    e = df.loc[('extended/N_inhadd1000/B/inhadd1000', 'P3')]
    assert e.panel == 'extended' and e.part == 'extended' and "extended input" in e.kind and e.kind.isascii()


def test_read_extra_does_not_repeat_a1s_record():
    """A dose is read through A1, whose `read_at` puts its PN record (pn_mode, pn_rung, e135_outcome, pn_inhadd,
    pn_inhsc) in the row beside `extra`: a dose's extra holds none of them (viki's crash, 2026-11-26 05:37); a form's
    records the normalization's terms itself."""
    pn_rec = {'pn_mode', 'pn_rung', 'e135_outcome', 'pn_inhadd', 'pn_inhsc', 'pn_check_max_abs_diff'}
    row_keys = {'tag', 'arm', 'panel', 'wiring', 'half', 'claw_thr', 'target_on', 'target_off', 'seconds', 'host',
                'lane', 'time', 'threads', 'input_sha256'}
    for a in L.DOSES:
        assert not set(L.read_extra(a)) & (pn_rec | row_keys), a
    for a in ('N', 'P', 'G'):
        x = L.read_extra(a)
        assert x['pn_inhadd'] == pytest.approx(1000.0) and not set(x) & row_keys


def test_stopped_class():
    c = L.stopped_class('P_inhadd1000', ['P_inhadd1000_megamat_real'], 3)
    assert c['cls'] == 'no class' and c['stands_in'] is False and c['does_more_than_a1'] is None
    assert c['n_reads'] == 3 and '3 of 32' in c['reasons'][1] and c['setting'] == 'P'
    assert L.stopped_class('D5.8_inhadd1000', ['x'], 1)['does_more_than_a1'] is False
