"""Tests of experiments/E148_megamat_other_views.py on planted values (no simulation, no shared data except where
marked)."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from experiments import E148_megamat_other_views as M

SHARED = Path('/mnt/tensor/Tom/claude_kc_remapping/results/E136_stage_population_fit/measure_strength/rescore_A1.csv')


@pytest.fixture
def planted(tmp_path, monkeypatch):
    """The module's folders in a temporary directory, its caches empty, and planted populations of the four stages
    (7, 5, 4 and 4 flies) for the functions that only need their flies' names."""
    monkeypatch.setitem(M.DIRS, 'results', tmp_path / 'results')
    monkeypatch.setitem(M.DIRS, 'local', tmp_path / 'local')
    for name in ('_VIEWS', '_SIDES', '_FLIES'):
        monkeypatch.setattr(M, name, {})
    sizes = dict(pn_dendrite=7, pn_bouton=5, kc_claw=4, kc_soma=4)
    for st, n in sizes.items():
        M._FLIES[st] = SimpleNamespace(flies=tuple(f'{st}_f{i}' for i in range(n)))
    return tmp_path


def test_sub_names():
    assert M.sub_name((2, 0)) == 'sub:0-2'
    assert M.sub_indices('sub:0-2', 4) == (0, 2)
    assert M.sub_indices('all', 3) == (0, 1, 2)


def test_readings_of_stage(planted):
    from experiments import E143_megamat_fly_sampling as E143m
    claw = M.readings_of_stage('kc_claw')
    assert list(claw)[0] == 'all' and len(claw) == 7          # every pair of four flies, and all
    assert {len(v) for k, v in claw.items() if k != 'all'} == {2}
    den = M.readings_of_stage('pn_dendrite')
    a, b = E143m.e136_split(7)
    assert set(den) == {'all', M.sub_name(a), M.sub_name(b)}
    groups = M.reading_groups()
    assert groups[0] == {st: 'all' for st in M.LOSS_STAGES}
    assert len(groups) == 7 and set(groups[1]) == set(M.LOSS_STAGES) and set(groups[3]) == set(M.CALYX_STAGES)


def test_joint_splits_fallback(planted):
    js = M.joint_splits_for(False)
    assert len(js) == 9 and int(js['e136'].sum()) == 1
    assert (js['pn_dendrite'] == M.e136_split_index('pn_dendrite')).all()
    js2 = M.joint_splits_for(True)
    assert len(js2) == M.N_JOINT + 1 and bool(js2['e136'].iloc[-1])


def _plant_flies(planted, stage='kc_soma'):
    """Flies' parts of one stage: 4 flies, strength/rate on 3 odors (one NaN), unit 0.5."""
    from experiments import E143_megamat_fly_sampling as E143m
    rows = []
    X = np.array([[1.0, 2.0, 3.0], [1.2, np.nan, 2.6], [0.8, 2.2, 3.4], [1.0, 1.8, 3.0]])
    for i in range(4):
        for j, o in enumerate(('o1', 'o2', 'o3')):
            rows.append(dict(key=f'{stage}/megamat/on', part='strength/rate', fly=f'{stage}_f{i}', entry=o,
                             value=X[i, j], unit=0.5))
    with E143m.views_at(M.views_dir()):
        path = E143m.flies_path(f'{stage}/megamat/on')
    pd.DataFrame(rows).to_parquet(path, index=False)
    return X


def _plant_views(cell='Z', stage='kc_soma'):
    """Two evaluations on two wirings, one NaN entry."""
    vals = {('e0', 'w1'): [1.0, 2.0, 3.0], ('e0', 'w2'): [3.0, np.nan, 1.0],
            ('e1', 'w1'): [1.1, 2.1, 2.9], ('e1', 'w2'): [0.9, 1.9, 3.1]}
    rows = []
    for (nm, w), v in vals.items():
        for o, x in zip(('o1', 'o2', 'o3'), v):
            rows.append(dict(cell=cell, name=nm, eval=int(nm[1:]), wiring=w, reading='all', stage=stage,
                             part='strength/rate', entry=o, value=x))
    M._VIEWS[cell] = pd.DataFrame(rows)
    return vals


def test_model_rows_common_is_nan_aware_mean(planted):
    _plant_views()
    m = M.model_rows('Z', 'kc_soma', 'all', common=True, each=False)
    e0 = m[m.setting == 'Z:e0'].set_index('entry')['value']
    assert e0['o1'] == pytest.approx(2.0) and e0['o2'] == pytest.approx(2.0) and e0['o3'] == pytest.approx(2.0)


def test_view_errors_against_halves(planted):
    """The error of the common row against a half: RMS over the entries both define of (model - half mean), over the
    unit."""
    X = _plant_flies(planted)
    _plant_views()
    sd = M.make_side('kc_soma', M.model_rows('Z', 'kc_soma', 'all'))
    W = M.single_weights(4, (0, 1))
    e = M.view_errors(sd, W, ('strength',))[0]
    half = np.nanmean(X[[0, 1]], axis=0)
    model = np.array([1.0, 2.0, 2.0])     # 'common' of e1 is [1.0, 2.0, 3.0]; check e0
    k = sd.rows.get_loc(('Z:e0', 'common'))
    want = np.sqrt(np.mean((model - half) ** 2)) / 0.5
    model0 = np.array([2.0, 2.0, 2.0])
    want0 = np.sqrt(np.mean((model0 - half) ** 2)) / 0.5
    assert e[k] == pytest.approx(want0)
    assert want != want0


def test_rms_and_losses():
    a, b = np.array([[3.0, 0.0]]), np.array([[4.0, 0.0]])
    assert M.rms([a, b])[0, 0] == pytest.approx(np.sqrt(12.5))
    assert [len(M.LOSSES[k][1]) for k in ('L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6')] == [4, 4, 4, 6, 3, 2, 2]
    assert all(s != 'pn_bouton' for s, _ in M.LOSSES['L4'][1])


def test_pick_and_regret():
    fit = np.array([1.0, 0.8, 0.7, 0.9])
    other = np.array([1.0, 0.6, 0.7, 0.5])
    p, r = M.pick_and_regret(fit, other)
    assert p == 2 and r == pytest.approx((0.7 - 0.5) / (1.0 - 0.5))


def test_decided():
    assert M.decided([0.9, 0.8, 0.85])
    assert not M.decided([0.9, 0.85, 0.4])        # one split under 0.5
    assert not M.decided([0.79, 0.7, 0.9])        # median under 0.8
    assert M.decided([0.9, np.nan, 0.85])         # a split without a value is left out
    assert not M.decided([np.nan])


def test_holds_part_counts():
    """Evaluation 1 falls by 0.05 and holds; 2 falls by 0.05 with its claws up beyond the spread; 3 falls by 0.03."""
    rows = []
    stage = {0: 1.0, 1: 0.95, 2: 0.95, 3: 0.97}
    claw = {0: 1.0, 1: 1.0, 2: 1.02, 3: 1.0}
    for i in range(4):
        for meas in ('strength', 'all'):
            r = dict(eval=i, measure=meas, half_ab='a', stage_error=stage[i])
            for st in M.LOSS_STAGES:
                r[f'{st}_error'] = claw[i] if st == 'kc_claw' and meas == 'strength' else 1.0
            rows.append(r)
    rs = pd.DataFrame(rows)
    cl = dict(spreads={s: 0.01 for s in M.LOSS_STAGES}, beside=dict(spreads_all_views={s: 0.01 for s in M.LOSS_STAGES}))
    out = M.holds_part_counts(rs, 'a', [0, 1, 2, 3], cl)
    assert out == dict(n=4, lowering=2, holds_part=1, claws_worse_among_lowering=1,
                       somata_all_views_worse_among_lowering=0)
    assert M.holds_part_counts(rs, 'a', [0, 1, 2, 3], cl, claw_spread=0.03)['claws_worse_among_lowering'] == 0


def test_claws_class():
    rows = []
    for c in M.CELLS:
        for ref, worse in (('e136_a', True), ('e136_b', False)):
            rows.append(dict(stage='kc_claw', reference=ref, reading="E136's own", cell=c, worse_than_spread=worse))
        for i in range(4):
            rows.append(dict(stage='kc_claw', reference=f'fly{i}', reading='read once', cell=c,
                             worse_than_spread=i == 2))
    out = M.claws_class(pd.DataFrame(rows))
    assert out['read_half'] and out['half'] == 'a' and out['kind'] == "one fly's"
    df = pd.DataFrame(rows)
    df.loc[(df.reference == 'fly1'), 'worse_than_spread'] = True
    assert M.claws_class(df)['kind'] == 'two against two'
    df.loc[(df.cell == 'W') & (df.reference == 'e136_b'), 'worse_than_spread'] = True
    assert not M.claws_class(df)['read_half']


def test_digits_and_close():
    assert M._digits(0.155) == 3 and M._digits(3.24) == 2 and M._digits(32) == 0
    assert M._close((5, 111), (5, 111), 0.005) and not M._close((5, 110), (5, 111), 0.005)
    assert M._close(0.3405, 0.34, 0.005)


def test_require_lock_refuses_untracked(planted):
    copy = planted / 'E148_lock.json'
    copy.write_text('{}\n')
    with pytest.raises(PermissionError):
        M.require_lock(copy=copy, shared=copy)


def test_half_to_half_symmetric(planted):
    X = _plant_flies(planted)
    a, b = (0, 2), (1, 3)
    h = M.half_to_half('kc_soma', a, b, ('strength',))
    d = np.nanmean(X[list(a)], axis=0) - np.nanmean(X[list(b)], axis=0)
    assert h == pytest.approx(np.sqrt(np.mean(d ** 2)) / 0.5)
    assert M.half_to_half('kc_soma', b, a, ('strength',)) == pytest.approx(h)


@pytest.mark.skipif(not SHARED.exists(), reason='the shared directory is not mounted')
def test_probe_reproduced_from_rescore_tables():
    """The planner's probe from E136's and E145's rescore tables: every printed value but the claws' 429 (430 at the
    spread as scored; 429 at the spread rounded to 0.017)."""
    t = M.probe_table()
    bad = t[~t['reproduced']]
    assert list(zip(bad['statistic'], bad['cell'])) == [('claws_worse_among_lowering', 'a')]
    assert bad['note'].iloc[0].startswith('429 ')
