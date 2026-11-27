"""E145 (the fit of all megamat stages made again in the other cells): the cells, E136's `cell` and its swapped halves,
the reads' plan and tags, and the lock's amendment rule. No data and no simulation."""
from types import SimpleNamespace

import pytest

from experiments import E136_stage_population_fit as F
from experiments import E145_stage_fit_replication as E


def test_cells_partition_the_wirings_and_halves():
    assert set(E.CELLS) == {'A1', 'X', 'F', 'W'} and E.FIT_CELLS == ('X', 'F', 'W')
    seen = set()
    for k, c in E.CELLS.items():
        s, r = c.searched_wirings(), c.read_wirings()
        assert len(s) == len(r) == 4 and set(s) | set(r) == set(F.WIRINGS) and not set(s) & set(r)
        assert {c.fit_half, c.read_half} == {'a', 'b'}
        seen.add((c.searched, c.fit_half))
    assert seen == {('S', 'a'), ('T', 'b'), ('S', 'b'), ('T', 'a')}
    assert E.CELLS['A1'].searched_wirings() == F.S_WIRINGS and E.CELLS['A1'].fit_part == 0
    x = E.CELLS['X']
    assert x.searched_wirings() == F.T_WIRINGS and x.fit_half == 'b' and x.read_on == 'S' and x.read_half == 'a'


def test_e136_cell_swaps_and_restores():
    s0, t0, p0, c0 = F.S_WIRINGS, F.T_WIRINGS, F.FIT_PART, F.START_CHECK
    with F.cell(F.T_WIRINGS, 1, dict(stage_error=1.5)):
        assert F.S_WIRINGS == t0 and F.T_WIRINGS == s0 and F.FIT_PART == 1
        assert F.START_CHECK == {F.MEASURE: dict(stage_error=1.5)}
        assert F.WIRINGS == s0 + t0 and F.seed_index('P3') == 4
        assert F.ArmModel.__init__.__defaults__[0] is None
    assert (F.S_WIRINGS, F.T_WIRINGS, F.FIT_PART, F.START_CHECK) == (s0, t0, p0, c0)
    with F.cell(F.S_WIRINGS, 1):
        assert F.START_CHECK == {}
    with F.cell(F.S_WIRINGS, 0):
        assert F.START_CHECK is c0
    with pytest.raises(RuntimeError):
        with F.cell(('P3', 'P4', 'real', 'P1'), 1):
            raise RuntimeError('inside')
    assert (F.S_WIRINGS, F.T_WIRINGS, F.FIT_PART, F.START_CHECK) == (s0, t0, p0, c0)
    for bad in ((('P3', 'P3', 'M2', 'M3'), 0), (('P3', 'P4', 'M2'), 0), (('P3', 'P4', 'M2', 'Q9'), 0),
                (F.T_WIRINGS, 2)):
        with pytest.raises(ValueError):
            with F.cell(*bad):
                pass


def test_split_names_exchanged_inside_a_cell():
    flies = SimpleNamespace(flies=[f'f{i}' for i in range(7)])
    a = F.split_names(flies)
    assert len(a['fit']) == 3 and len(a['held']) == 4
    with F.cell(F.T_WIRINGS, 1):
        b = F.split_names(flies)
    with F.cell(F.T_WIRINGS, 0):
        c = F.split_names(flies)
    assert b == dict(fit=a['held'], held=a['fit']) and c == a


def test_planned_points_and_tags():
    pts = E.planned_points()
    assert len(pts) == 84 and len({E.point_tag(p) for p in pts}) == 84
    for k in E.FIT_CELLS:
        mine = [p for p in pts if p.cell == k]
        assert len(mine) == 28
        mega = [p.wiring for p in mine if p.panel == 'megamat']
        assert tuple(mega) == E.CELLS[k].read_wirings()
        for pan in ('validation2', 'mixtures', 'apl16'):
            assert tuple(p.wiring for p in mine if p.panel == pan) == F.WIRINGS
    p = E.Point('W', 'apl16', 'M3')
    assert E.point_tag(p) == 'W_A1_apl16_M3' and E.parse_tag(E.point_tag(p)) == p
    for bad in ('A1_A1_megamat_real', 'X_A0_megamat_real', 'X_A1_kiwi_real'):
        with pytest.raises(ValueError):
            E.parse_tag(bad)


def test_lock_amendments_only():
    cur = dict(a=1, amendments=[dict(n=1)], written='t', host='h')
    assert E._same_lock(cur, dict(cur, written='u', host='g'))
    assert E._amended_only(cur, dict(cur, amendments=[dict(n=1), dict(n=2)]))
    assert not E._amended_only(cur, dict(cur, a=2, amendments=[dict(n=1), dict(n=2)]))
    assert not E._amended_only(cur, dict(cur, amendments=[dict(n=3), dict(n=2)]))
    assert not E._amended_only(cur, cur)


def test_repeat_ranges_are_the_entrys():
    assert (E.REPEAT['inhadd']['lo'], E.REPEAT['inhadd']['hi']) == (562.0, 1778.0)
    assert (E.REPEAT['aplpn']['lo'], E.REPEAT['aplpn']['hi']) == (0.078, 0.31) and E.REPEAT['aplpn']['determined']
    assert E.REPEAT['kckc']['lo'] == 1.8 and E.REPEAT['kckc']['hi'] is None and E.REPEAT['kckc']['determined']
