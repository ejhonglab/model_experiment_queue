"""Tests of experiments/E149_two_view_picks.py on planted values (no simulation; shared data only where marked)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from experiments import E149_two_view_picks as P

SHARED = Path('/mnt/tensor/Tom/claude_kc_remapping/results/E148_megamat_other_views/losses_picks.csv')


def _pick(cell, i, loss, steps, searched, read, fit_half, read_half):
    z = dict(aplpn=-1.0, kckc=1.0, inhadd=3.5, inhsc=0.1)
    return dict(name=f'{cell}e{i}', cell=cell, eval=i, loss=loss, z=z, factors={k: float(np.exp(v)) for k, v in z.items()},
                pn_kws=dict(pn_inhadd=1000.0, pn_inhsc=500.0), searched=list(searched), read=list(read),
                fit_half=fit_half, read_half=read_half, steps=list(steps), arm='A1')


S, T = list(P.S_WIRINGS), list(P.T_WIRINGS)
LOCK = dict(picks=[_pick('A1', 46, 'L1', (1, 2), S, T, 'a', 'b'), _pick('X', 21, 'L1', (1, 2), T, S, 'b', 'a'),
                   _pick('F', 86, 'L1', (1, 2), S, T, 'b', 'a'), _pick('W', 24, 'L1', (1, 2), T, S, 'a', 'b'),
                   _pick('F', 116, 'L0', (2,), S, T, 'b', 'a'),
                   _pick('A1', 45, 'L0_stopped', (3,), S, T, 'a', 'b'), _pick('X', 14, 'L0_stopped', (3,), T, S, 'b', 'a'),
                   _pick('F', 80, 'L0_stopped', (3,), S, T, 'b', 'a'), _pick('W', 15, 'L0_stopped', (3,), T, S, 'a', 'b')])


def test_tags_round_trip():
    assert P.point_tag('Xe21', 'megamat', 'real') == 'Xe21_megamat_real'
    assert P.parse_tag('A1e46_validation2_P3') == ('A1e46', 'validation2', 'P3')
    assert P.parse_tag('Fe116_extended_M1') == ('Fe116', 'extended', 'M1')
    assert P.parse_pick('We24') == ('W', 24)
    for bad in ('Ye1_megamat_real', 'Xe21_apl99_real', 'Xe21_megamat_Q1', 'X21_megamat_real'):
        with pytest.raises(ValueError):
            P.parse_tag(bad)


def test_planned_points():
    s1 = P.planned_points(1, lock=LOCK)
    assert len(s1) == 112 and len({p.tag for p in s1}) == 112
    assert [p.panel for p in s1[:16]] == ['megamat'] * 16
    assert [p.wiring for p in s1[:4]] == T and [p.wiring for p in s1[4:8]] == S   # A1 read on T, X on S
    assert [p.panel for p in s1[16:48]] == ['apl16'] * 32 and s1[-1].panel == 'mixtures'
    s2 = P.planned_points(2, lock=LOCK)
    assert len(s2) == 60 and [p.panel for p in s2[:32]] == ['extended'] * 32
    assert all(p.pick == 'Fe116' for p in s2[32:]) and [p.wiring for p in s2[32:36]] == T
    s3 = P.planned_points(3, cells=['A1', 'W'], lock=LOCK)
    assert len(s3) == 56 and {p.pick for p in s3} == {'A1e45', 'We15'}
    assert not {p.tag for p in s1} & {p.tag for p in s2} and not {p.tag for p in s3} & {p.tag for p in s1 + s2}
    pt = s1[0]
    assert (pt.cell, pt.eval, pt.half) == ('A1', 46, 'held') and s1[20].half is None
    ex = P.read_extra(pt, lock=LOCK)
    assert ex == dict(step=1, pick='A1e46', cell='A1', eval=46, loss='L1', searched='S', fit_half='a', read_half='b')
    assert P.read_extra(s1[20], lock=LOCK)['read_half'] is None
    assert P.point_step('Fe116_megamat_P3', lock=LOCK) == 2 and P.point_step('A1e46_extended_P3', lock=LOCK) == 2
    assert P.point_step('We15_apl16_real', lock=LOCK) == 3 and P.point_step('We24_apl16_real', lock=LOCK) == 1


def test_folders(tmp_path, monkeypatch):
    monkeypatch.setitem(P.DIRS, 'results', tmp_path / 'r')
    monkeypatch.setitem(P.DIRS, 'local', tmp_path / 'l')
    monkeypatch.setitem(P.DIRS, 'archive', tmp_path / 'a')
    d = P.run_dirs(P.Point(1, 'Xe21', 'megamat', 'real'))
    assert d[0] == tmp_path / 'l' / 'step1' / 'Xe21' / 'megamat' / 'real' and d[0].is_dir()
    assert d[1] == tmp_path / 'a' / 'step1' / 'Xe21' / 'megamat' / 'real' and not d[1].exists()
    e = P.extended_dir('Xe21', 'local', lock=LOCK)
    assert e.parent.name == 'B' and e.parent.parent == tmp_path / 'l' / 'step2' / 'Xe21' / 'extended'
    assert P.extended_dir('Xe21', 'archive', lock=LOCK).parts[-4:-2] == ('Xe21', 'extended')
    assert P.reads_table().empty
    pd.DataFrame([dict(tag='Xe21_megamat_real', rate_on=0.09)]).to_csv(tmp_path / 'r' / 'reads_h_h.csv', index=False)
    pd.DataFrame([dict(tag='Xe21_megamat_real', rate_on=0.1)]).to_csv(tmp_path / 'r' / 'reads_v_v.csv', index=False)
    assert len(P.reads_table()) == 1


def test_is_held():
    r = dict(rate_on=0.0915, frac_on_target=0.0915, rate_off=0.2, frac_off_target=0.2004, at_bound_thr=False,
             at_bound_apl='False')
    assert P.is_held(r)
    assert not P.is_held(dict(r, rate_off=0.19))
    assert not P.is_held(dict(r, at_bound_apl=True))
    assert not P.is_held(dict(r, frac_on_target=np.nan))


@pytest.mark.skipif(not SHARED.exists(), reason='the shared directory is not mounted')
def test_pick_table_from_the_saved_tables():
    t = P.pick_table()
    assert [p.eval for p in P.l1_picks(t)] == [46, 21, 86, 24]
    assert sorted(p.eval for p in t.values() if p.loss == 'L0_stopped') == [14, 15, 45, 80]
    assert t['Fe116'].loss == 'L0'
    x = t['Xe21']
    assert x.read == P.S_WIRINGS and x.read_half == 'a' and np.isclose(x.pn_kws['pn_inhadd'], 2272.5, atol=0.1)
