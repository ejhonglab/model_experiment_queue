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
