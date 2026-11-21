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


def test_outcome_letter():
    assert E.outcome_letter(True, dict(found_inside=True, level=100.0)) == 'a'
    assert E.outcome_letter(False, dict(found_inside=True, level=1000.0)) == 'b'
    assert E.outcome_letter(False, dict(found_inside=True, level=314.0)) == 'b'
    assert E.outcome_letter(False, dict(found_inside=True, level=3141.0)) == 'c'
    assert E.outcome_letter(False, dict(found_inside=False, level=10000.0)) == 'd'


def test_write_variants_and_part_b_reading(tmp_path, monkeypatch):
    """variants.csv from the stage records (absolute parameters, the refit's tag, response rates from the saved
    responses) and part B's reading from b_stage_summary and E102's class."""
    shared = tmp_path / 'shared'
    res = shared / 'results' / E.EXPERIMENT
    res.mkdir(parents=True)
    monkeypatch.setattr(E, 'out_dir', lambda: res)
    monkeypatch.setattr(E.H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(E.H, 'outputs_dir', lambda e: tmp_path / 'outputs' / e)
    d = shared / 'data' / 'model_outputs' / E.EXPERIMENT / 'B' / 'inhadd10000'
    d.mkdir(parents=True)
    (d / 'real_stage.json').write_text(json.dumps(dict(
        wiring='real', level=10000.0, tag='rate0.068_x', params=dict(fixed_thr=132.4, wAPLKC_scale=0.15),
        pn_kws=dict(pn_inhadd=10000.0, pn_inhsc=2177.8), exact=True)))
    pd.DataFrame([[1.0, 0.0], [0.0, 0.0]]).to_parquet(d / 'real_on_responses.parquet')
    pd.DataFrame([[1.0, 1.0], [0.0, 0.0]]).to_parquet(d / 'real_acute_off_responses.parquet')
    v = E.write_variants(10000.0, ('real', 'P1'))
    assert list(v.tag) == ['real'] and v.fixed_thr.iloc[0] == 132.4 and v.pn_inhsc.iloc[0] == 2177.8
    assert v.rate_on.iloc[0] == 0.25 and v.rate_acute_off.iloc[0] == 0.5 and v.ref_fit.iloc[0] == 'rate0.068_x'
    assert (res / 'variants.csv').exists()
    assert E.part_b_reading(10000.0) is None
    c = E.VERDICT_CRITERION
    pd.DataFrame(dict(model=['refit inhadd10000', 'default'], stage='bouton', cond='on', reading='rate', criterion=c,
                      n_wirings=4, totals_cv=[0.2, 0.1], totals_cv_in=[0, 0])).to_csv(
        res / 'b_stage_summary_inhadd10000.csv', index=False)
    e102 = shared / 'results' / 'E102_antennal_lobe_normalization'
    e102.mkdir(parents=True)
    (e102 / 'd_class_inhadd10000.json').write_text(json.dumps(dict(S=dict(
        **{'class': 'price'}, pattern=0.058, pattern_unchanged=0.050, d_pattern=dict(z8=0.008), d_raw_loss=-0.02,
        n_flagged=0))))
    pb = E.part_b_reading(10000.0)
    assert pb['totals_cv_on_refit'] == 0.2 and pb['totals_cv_on_default'] == 0.1 and pb['n_wirings'] == 4
    assert pb['refit_class'] == 'price' and pb['pattern_refit'] == 0.058
