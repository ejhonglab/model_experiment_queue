"""Tests of E120 (experiments/E120_bouton_stage_measured.py) and its gate in E102's part D; synthetic data only, no
shared-directory reads."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments import E102_antennal_lobe_normalization as E102
from experiments import E120_bouton_stage_measured as E


def test_indicator_change_forms():
    win = pd.DataFrame([[4.0, 1.0], [9.0, 0.0]], columns=['2h', 'IaA'])
    pre = pd.DataFrame([[1.0, 2.0], [4.0, 0.0]], columns=['2h', 'IaA'])
    lin = E.indicator_change(win, pre, 1.0, 'evoked')
    np.testing.assert_allclose(lin.to_numpy(), (win - pre).to_numpy())
    np.testing.assert_allclose(E.indicator_change(win, pre, 1.0, 'baseline').to_numpy(), lin.to_numpy())
    np.testing.assert_allclose(E.indicator_change(win, pre, 2.0, 'baseline').to_numpy(), [[15.0, -3.0], [65.0, 0.0]])
    # the evoked form keeps the sign of the change
    np.testing.assert_allclose(E.indicator_change(win, pre, 2.0, 'evoked').to_numpy(), [[9.0, -1.0], [25.0, 0.0]])
    with pytest.raises(ValueError):
        E.indicator_change(win, pre, 2.0, 'other')


def _summ(rows):
    return pd.DataFrame(rows)


def test_choose_level_rule():
    def row(step, reading, n_in, n_in_est=0):
        return dict(step=E.e102_step_tag(step), reading=reading, criterion=E.VERDICT_CRITERION, totals_cv_in=n_in,
                    totals_cv_on_est_in=n_in_est)
    base = [row(s, rd, 0) for s in E.LEVEL_STEPS for rd in E.READINGS]
    lv = E.choose_level(_summ(base))
    assert lv['level'] == 10000.0 and not lv['found_inside'] and lv['part_b_runs']
    # the smallest step inside on 6 of 8 wirings, 'rate' reading only
    rows = [row(s, rd, 6 if (s >= 3141 and rd == 'rate') else 0) for s in E.LEVEL_STEPS for rd in E.READINGS]
    lv = E.choose_level(_summ(rows))
    assert lv['level'] == 3141.0 and lv['found_inside'] and lv['part_b_runs']
    # 5 of 8 is not inside; the other readings do not decide
    rows = [row(s, rd, 5 if rd == 'rate' else 8) for s in E.LEVEL_STEPS for rd in E.READINGS]
    assert E.choose_level(_summ(rows))['level'] == 10000.0
    # at 1000 E102's part D is the answer
    rows = [row(s, rd, 7 if s >= 1000 else 0) for s in E.LEVEL_STEPS for rd in E.READINGS]
    lv = E.choose_level(_summ(rows))
    assert lv['level'] == 1000.0 and not lv['part_b_runs']


def test_a1_pass_rule():
    rows = []
    for n, ins in ((1.0, (8, 8, 8)), (2.0, (6, 7, 5)), (2.5, (6, 6, 6))):
        rows.append(dict(stage='bouton', n=n, form='evoked', reading='rate', criterion=E.VERDICT_CRITERION,
                         **{f'{s}_in': k for s, k in zip(E.THREE, ins)}, **{s: 0.0 for s in E.THREE}))
    rows.append(dict(stage='bouton', n=1.0, form='evoked', reading='rate', criterion=1.0,
                     **{f'{s}_in': 0 for s in E.THREE}, **{s: 0.0 for s in E.THREE}))
    ps = E.a1_pass(pd.DataFrame(rows))
    assert list(ps.passes) == [True, False, True]


def test_dose_outcome_classes():
    def rows(inside_by_reading):
        out = []
        for rd, ins in inside_by_reading.items():
            for d in (0.0, 0.1, 0.5, 1.0, 2.0):
                out.append(dict(variant=f'v{d}', form='sub', exponent=2.0, dose=d, reading=rd, criterion=E.VERDICT_CRITERION,
                                inhibited=d, inhibited_in=4 if d in ins else 1))
        return pd.DataFrame(out)
    assert E.dose_outcome(rows({'rate': (0.1, 0.5)}))['overall'] == 'below_fitted_only'
    assert E.dose_outcome(rows({'rate': (0.0, 0.1, 0.5, 1.0, 2.0)}))['overall'] == 'every_dose'
    assert E.dose_outcome(rows({'rate': (), 'session_z': (0.1,)}))['overall'] == 'reading_decides'
    assert E.dose_outcome(rows({'rate': ()}))['overall'] == 'none'


def test_stats_of_matches_e119_on_synthetic_units(monkeypatch):
    rng = np.random.default_rng(0)
    odors = list(E.ODORS)
    V = pd.DataFrame(rng.gamma(2.0, 1.0, (50, len(odors))), columns=odors)
    sp = {'rate': dict(lo=0.2, hi=0.4, mean=0.3)}
    for s in ('opr_mean', 'inhibited', 'resp_cv', 'top1', 'totals_cv'):
        sp[s] = dict(lo=0.0, hi=10.0, mean=1.0)
    sp['per_odor'] = dict(mean=np.full(len(odors), 0.3), band=1.0)
    sp['opr'] = dict(mean=np.full(5, 0.2), band=1.0)
    monkeypatch.setattr(E, 'fly_spreads', lambda stage: {c: sp for c in E.CRITERIA})
    out = E.stats_of(V, V, 'bouton')
    assert [r['criterion'] for r in out] == list(E.CRITERIA)
    r = out[0]
    assert r['rate'] == pytest.approx(0.3, abs=0.01)
    assert r['totals_cv'] == pytest.approx(E._cv(V.mean(axis=0)))
    assert r['totals_cv_side'] == 0


def test_part_b_commands():
    cmds = E.part_b_commands(10000.0)
    assert len(cmds) == 3
    assert '--part D --level 10000 --wirings real,P1 --e120-level' in cmds[0]
    assert '--part Bstage --level 10000 --wirings P2,M1' in cmds[1]
    assert 'OMP_THREAD_LIMIT=16' in cmds[0] and 'nice -n 10' in cmds[0]
    assert '--score-refits' in cmds[2] and '--part Bread' in cmds[2]


def test_e102_part_d_takes_e120_level(tmp_path, monkeypatch):
    path = tmp_path / 'level.json'
    monkeypatch.setattr(E102, 'e120_level_path', lambda: path)
    ok, why = E102.refit_due(10000.0, ['real'], e120=True)
    assert not ok and 'no E120 level' in why
    path.write_text(json.dumps(dict(level=10000.0, part_b_runs=True)))
    assert E102.refit_due(10000.0, ['real', 'P1'], e120=True)[0]
    assert not E102.refit_due(3141.0, ['real'], e120=True)[0]
    assert not E102.refit_due(10000.0, ['P3'], e120=True)[0]
    path.write_text(json.dumps(dict(level=1000.0, part_b_runs=False)))
    assert not E102.refit_due(1000.0, ['real'], e120=True)[0]
    # without the flag part C's levels decide, as before
    monkeypatch.setattr(E102, 'out_dir', lambda: tmp_path)
    (tmp_path / 'levels.json').write_text(json.dumps(dict(refit_levels=[1000.0])))
    assert E102.refit_due(1000.0, ['real'])[0]
    assert not E102.refit_due(10000.0, ['real'])[0]
