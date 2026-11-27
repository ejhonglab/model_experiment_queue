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


def test_the_four_factor_control():
    assert F.CONTROL_A1 == E.CONTROL == 'A1_perm' and F.CONTROLS == {F.CONTROL: 'A0', F.CONTROL_A1: 'A1'}
    assert F.base_arm(F.CONTROL_A1) == 'A1' and F.base_arm(F.CONTROL) == 'A0' and F.base_arm(F.FAMILY) == 'A0'
    assert F.base_arm('A3') == 'A3'
    pts = E.control_points()
    assert [p.wiring for p in pts] == list(E.CELLS[E.CONTROL_CELL].read_wirings()) and all(p.arm == E.CONTROL for p in pts)
    tags = [E.point_tag(p) for p in pts]
    assert tags[0] == 'X_A1perm_megamat_real' and [E.parse_tag(t) for t in tags] == pts
    assert not set(tags) & {E.point_tag(p) for p in E.planned_points()}
    for bad in ('F_A1perm_megamat_P3', 'X_A1perm_validation2_real'):
        with pytest.raises(ValueError):
            E.parse_tag(bad)
    with pytest.raises(ValueError):
        F.register_variant(F.CONTROL_A1, lambda pt, claw_after: None)


def _rec(cls, change, se, kckc=True, at_bound=None, inhadd=True, aplpn=True):
    return dict(cls=dict(cls=cls, megamat=dict(complete=True, change=change, se=se)),
                repeats=dict(inhadd=dict(repeats=inhadd), aplpn=dict(repeats=aplpn), kckc=dict(repeats=kckc)),
                fit=dict(at_bound_classed=at_bound or []))


def test_score_fall_and_experiment_classes():
    from experiments import E145_score as S
    assert S.fall_repeats_in(_rec('holds', -0.10, 0.01), 0.04) is True
    assert S.fall_repeats_in(_rec('holds', -0.03, 0.01), 0.04) is False     # beyond the bar, not the control's fall
    assert S.fall_repeats_in(_rec('holds', -0.05, 0.03), 0.04) is False     # beyond the control's fall, not the bar
    assert S.fall_repeats_in(dict(cls=dict(megamat=dict(complete=False))), 0.04) is None
    a1 = _rec('holds', -0.115, 0.01)
    allhold = dict(A1=a1, X=_rec('holds', -0.1, 0.01), F=_rec('holds', -0.09, 0.01), W=_rec('megamat only', -0.1, 0.01))
    c = S.experiment_classes(allhold, 0.04)
    assert c['complete'] and c['fall_repeats'] and c['hold_repeats'] and c['outcome'][0] == 'a'
    assert c['the_wirings'] and not c['the_flies'] and 'd' in c['outcome']
    b = dict(A1=a1, X=_rec('megamat only', -0.1, 0.01, kckc=False), F=_rec('megamat only', -0.1, 0.01, kckc=False),
             W=_rec('holds', -0.1, 0.01))
    c = S.experiment_classes(b, 0.04)
    assert 'b' in c['outcome'] and not c['hold_repeats'] and not c['megamat_fixes']['kckc']
    e = dict(A1=a1, X=_rec('holds', -0.1, 0.01), F=_rec('does not improve', 0.0, 0.01),
             W=_rec('no class', -0.1, 0.01, at_bound=['kckc high']))
    c = S.experiment_classes(e, 0.04)
    assert 'e' in c['outcome'] and 'f' in c['outcome'] and c['no_class_cells'] == ['W'] and c['at_bound_cells'] == ['W']
    assert S.experiment_classes(dict(A1=a1, X=_rec('holds', -0.1, 0.01)), 0.04)['complete'] is False
    # a cell classed before its held-out panels are all read does not count as scored
    w = _rec("strength at another view's cost", -0.1, 0.01)
    w['cls']['panel_missing'] = ['apl16/apl16/on']
    c = S.experiment_classes(dict(e, W=w), 0.04)
    assert c['complete'] is False and c['cells_scored'] == ['X', 'F']


def test_score_repeats_on_a_planted_fit():
    import numpy as np
    from experiments import E145_score as S
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    z = dict(aplpn=float(np.log(0.155)), kckc=float(np.log(3.24)), inhadd=float(np.log(1000.0 / PN_INHADD_DEFAULT)),
             inhsc=0.0)
    det = {k: dict(determined=True) for k in z}
    r = S.repeats('X', dict(classed_z=z, best_z=z, determined=det))
    assert r['complete'] and r['inhadd']['repeats'] and r['aplpn']['repeats'] and r['kckc']['repeats']
    assert abs(r['inhadd']['value'] - 1000.0) < 1e-6 and r['inhsc_gain']['read'] is False
    det['kckc'] = dict(determined=False)
    z2 = dict(z, inhadd=float(np.log(3141.0 / PN_INHADD_DEFAULT)), aplpn=float(np.log(0.05)))
    r = S.repeats('X', dict(classed_z=z2, best_z=z2, determined=det))
    assert not r['kckc']['repeats'] and not r['inhadd']['repeats'] and not r['aplpn']['repeats']
    assert S.repeats('X', dict(best_z=z))['complete'] is False
    assert S.variant('W') == 'E145_W'


def test_resampled_derived_classes():
    import numpy as np
    from experiments import E145_score as S
    K = 4
    cls = {'A1': np.array(['holds'] * K, dtype=object),
           'X': np.array(['holds', 'holds', 'megamat only', 'holds'], dtype=object),
           'F': np.array(['holds', 'megamat only', 'megamat only', 'does not improve'], dtype=object),
           'W': np.array(['holds', 'holds', 'megamat only', 'megamat only'], dtype=object)}
    comps = {k: dict(megamat=dict(improves=np.array([True, True, True, False]), change=np.array([-0.1, -0.1, -0.1, 0.0])))
             for k in cls}
    rep = dict(inhadd=dict(repeats=True), aplpn=dict(repeats=True), kckc=dict(repeats=True))
    recs = {k: dict(repeats=rep) for k in cls}
    d = S.derived_per_reference(cls, comps, recs, 0.04, K)
    assert list(d['hold_repeats']) == [True, True, False, False]
    assert list(d['fall_repeats']) == [True, True, True, False]
    assert list(d['outcome_a']) == [True, True, False, False]
    assert list(d['the_wirings']) == [False, False, False, False] and list(d['the_flies']) == [False, True, False, False]
    s = S._share_of(d['hold_repeats'], True)
    assert s['all_flies'] is False and abs(s['share_true'] - 2 / 3) < 1e-12 and abs(s['share_scored'] - 2 / 3) < 1e-12
    assert S.derived_per_reference({'X': cls['X']}, comps, recs, 0.04, K) == {}


def test_saved_evals_for_a_moved_fit(monkeypatch, tmp_path):
    import pandas as pd
    out, loc = tmp_path / 'out', tmp_path / 'loc'
    monkeypatch.setattr(E, 'out_dir', lambda cell=None: (out / cell if cell else out))
    monkeypatch.setattr(E, 'local_dir', lambda cell=None, part='': (loc / cell / part))
    (out / 'F').mkdir(parents=True)
    pd.DataFrame(dict(eval=[0, 1, 2], loss_sim=[1.0, 0.9, 0.95])).to_csv(out / 'F' / 'fit_history_A1.csv', index=False)
    d = loc / 'F' / 'fit_A1'
    d.mkdir(parents=True)
    for i in (0, 2):
        for w in E.CELLS['F'].searched_wirings():
            (d / f'A1_e{i}_{w}_run.json').write_text('{}')
    (d / 'A1_e1_real_run.json').write_text('{}')     # one wiring only: not a saved evaluation
    assert E.saved_evals('F') == [0, 2]
    assert E.saved_evals('W') == []
    with pytest.raises(ValueError):
        E.part_archive('F', control=True)


def test_write_variants_on_a_planted_archive(tmp_path, monkeypatch):
    import json
    import pandas as pd
    out = tmp_path / 'results'
    out.mkdir()
    monkeypatch.setattr(E, 'out_dir', lambda cell=None: out)
    arch = tmp_path / 'archive'
    params = dict(fixed_thr=50.0, wAPLKC_scale=1.0, wKCAPL_scale=2.0, wAPLPN_scale=0.05, wPNAPL_scale=3.0,
                  wKCKC_scale=100.0)
    runs = [('X', 'fit_A1', dict(tag='A1_e115_P3', arm='A1', panel='megamat', wiring='P3', reading='fit', eval=115,
                                 z=dict(aplpn=-1.7, kckc=1.99), params=params, params_off=params,
                                 pn_kws=dict(pn_inhadd=2144.0, pn_inhsc=273.4), host='viki', time='t')),
            ('W', 'fit_A1', dict(tag='A1_e115_P3', arm='A1', panel='megamat', wiring='P3', reading='fit', eval=115,
                                 z={}, params=params, params_off=params, host='hal', time='t')),
            ('F', 'read', dict(tag='F_A1_megamat_P3', arm='A1', panel='megamat', wiring='P3', reading='read',
                               params=params, params_off=params, pn_inhadd=1633.0, pn_inhsc=468.5, host='disco',
                               lane='disco', time='t')),
            ('X', 'fit_A1_perm', dict(tag='A1_perm_e0_M2', arm='A1_perm', panel='megamat', wiring='M2',
                                      reading='fit', eval=0, params=params, params_off=params, host='hal',
                                      time='t'))]
    for cell, part, d in runs:
        (arch / cell / part).mkdir(parents=True, exist_ok=True)
        (arch / cell / part / f"{d['tag']}_run.json").write_text(json.dumps(d))
    df = E.write_variants(arch).set_index('subdir')
    assert len(df) == 4 and (out / 'variants.csv').exists()
    assert df.loc['X/fit_A1', 'pn_inhadd'] == 2144.0 and df.loc['X/fit_A1', 'z_kckc'] == 1.99
    assert df.loc['F/read', 'pn_inhadd'] == 1633.0 and df.loc['F/read', 'searched'] == 'S'
    assert df.loc['F/read', 'fit_half'] == 'b' and df.loc['W/fit_A1', 'read_half'] == 'b'
    assert 'evaluation 115 on searched wiring P3' in df.loc['W/fit_A1', 'kind']
    assert 'four-factor control' in df.loc['X/fit_A1_perm', 'kind']
    assert 'classed point, read on the megamat input' in df.loc['F/read', 'kind']
    assert df.loc['X/fit_A1', 'off_wKCKC_scale'] == 100.0
    (arch / 'X' / 'fit_A1' / 'A1_e115_P3b_run.json').write_text(json.dumps(dict(runs[0][2], tag='A1_e115_P3')))
    with pytest.raises(ValueError):
        E.write_variants(arch)
