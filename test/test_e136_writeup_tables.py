"""Tests of E136's write-up tables (analyses/E136_writeup_tables.py) on planted score and fit files; no shared-directory
reads and no simulation."""
import json

import numpy as np
import pandas as pd
import pytest

from analyses import E136_writeup_tables as W
from experiments import E136_stage_population_fit as F

STAGES = F.LOSS_STAGES


def _cmp(arm, ref, change, se, spread, improves=False, worse=False):
    return dict(n=4, complete=True, changes=[change] * 4, change=change, se=se, spread=spread, arm=arm, ref=ref,
                improves=improves, worse_than_spread=worse)


def _plant(d):
    rows = []
    for arm, base in (('family', 1.05), ('A0', 1.04)):
        for half in ('held', 'all', 'fit'):
            ws = list(F.WIRINGS) if (arm == 'family' and half == 'held') else list(F.T_WIRINGS)
            for i, w in enumerate(ws + ['common']):
                r = dict(arm=arm, half=half, wiring=w, stage=base + 0.001 * i, downstream=base)
                r.update({s: base + 0.01 * j for j, s in enumerate(STAGES)})
                if w == 'common':
                    for s in STAGES:
                        for p in ('rate', 'mean'):
                            r.update({f'{s}_split_strength_{p}_{k}': 0.5 for k in
                                      ('level', 'spread', 'order', 'sd_ratio', 'r')})
                rows.append(r)
    pd.DataFrame(rows).to_csv(d / 'score_megamat.csv', index=False)
    pd.DataFrame(columns=['arm', 'input', 'panel']).to_csv(d / 'score_panels.csv', index=False)
    apl = []
    for arm in ('family', 'A0'):
        for w in F.WIRINGS:
            apl.append(dict(arm=arm, wiring=w, rate_on=0.068, rate_off=0.15, rate_ratio=2.2, rate_rho_on=0.6,
                            rate_rho_off=0.7, opk_on_tv=0.03, opk_off_tv=0.08, max_count_on=40, max_count_off=70,
                            realism_flag='response rate ρ, APL on, 0.30 (below 0.37)' if (arm, w) == ('A0', 'real')
                            else '', cal_loss=0.06, cal_mean=0.16, cal_pattern=0.058, cal_rho=0.45))
    pd.DataFrame(apl).to_csv(d / 'score_apl16.csv', index=False)
    keys = ['/'.join(k) for k in F.HOLDS_PANELS]
    beside = ['/'.join(k) for k in F.BESIDE_PANELS]
    a0 = dict(arm='A0', cls='no class', reasons=['no better than the permuted control'], fall=0.010,
              megamat=_cmp(1.039, 1.049, -0.010, 0.010, 0.020),
              stages={s: _cmp(1.0, 1.0, -0.01, 0.01, 0.02, worse=(s == 'kc_soma')) for s in STAGES},
              guard={**{s: _cmp(1.0, 1.0, 0.02, 0.01, 0.05) for s in STAGES},
                     **{k: _cmp(1.0, 1.0, 0.0, 0.01, 0.05) for k in keys}},
              panels={k: _cmp(1.08, 1.10, -0.02, 0.02, 0.04) for k in keys},
              panels_beside={k: _cmp(2.50, 2.53, -0.03, 0.06, 0.117, worse=False) for k in beside},
              stage_worse=['kc_soma'], guard_worse=[], panel_worse=[], panel_missing=[])
    classes = dict(stages=list(STAGES), measure='strength', control_fall=0.040,
                   spreads=dict(stage=0.020, **{s: 0.02 for s in STAGES}, **{k: 0.04 for k in keys}),
                   beside=dict(yardsticks={s: dict(flat=1.5, half_to_half=0.8) for s in STAGES},
                               two_view={'A0': _cmp(1.23, 1.20, 0.03, 0.01, np.nan)},
                               spreads_all_views=dict(stage=0.03, **{s: 0.05 for s in STAGES}),
                               e137_readings={'A0': {f'{p}/{r}/{v}': _cmp(2.378, 2.425, -0.047, 0.07, 0.139)
                                                     for p, rs in F.E137_BESIDE.items() for r in rs
                                                     for v in ('strength', 'all')}}),
                   classes={'A0': a0})
    (d / 'classes.json').write_text(json.dumps(classes))
    (d / 'provenance_score.json').write_text(json.dumps(dict(time='t', host='h', queue_commit='abcdef0', arms=['A0'])))
    pd.DataFrame([dict(measure='strength', statistic='stage_error', n_evals=10, rho=0.76, yardstick=0.8)]).to_csv(
        d / 'consistency_strength.csv', index=False)
    fit = dict(arm='A0', searched=['aplpn', 'kckc'], n_evals=2, budget=80, best_eval=1, best=dict(loss_sim=0.9565),
               classed_z=dict(aplpn=float(np.log(0.154)), kckc=0.0), classed_loss=0.9572, at_bound=[],
               at_bound_classed=[], messages=['Optimization terminated successfully.'],
               determined=dict(aplpn=dict(rise=0.02, determined=True), kckc=dict(rise=0.0007, determined=False)))
    (d / 'fit_A0.json').write_text(json.dumps(fit))
    h = []
    for i, (loss, apn) in enumerate(((0.9747, 0.35), (0.9565, 0.054))):
        r = dict(eval=i, loss_stages=loss, loss_sim=loss)
        for w in F.S_WIRINGS:
            r.update({f'rate_on_{w}': 0.0915, f'rate_off_{w}': 0.2, f'wAPLPN_scale_{w}': apn, f'wAPLKC_scale_{w}': 0.3,
                      f'wKCKC_scale_{w}': 30.0, f'fixed_thr_{w}': 40.0})
        h.append(r)
    pd.DataFrame(h).to_csv(d / 'fit_history_A0.csv', index=False)


def test_sites_row():
    from kcq import harness as H
    kc = H.KC_ID_LEVEL
    claws = pd.MultiIndex.from_tuples([(0, 0), (0, 1), (1, 2), (1, 3)], names=[kc, 'claw'])
    boutons = pd.MultiIndex.from_tuples([('DA1', 0), ('DL1', 1)], names=['glomerulus', 'bouton'])
    Wm = pd.DataFrame([[1, 0], [0, 1], [1, 1], [0, 2]], index=claws, columns=boutons, dtype=float)
    odors = ['2h', 'IaA']
    b = lambda x: pd.DataFrame(x, index=boutons, columns=odors, dtype=float)
    off = dict(bouton_win=b([[10, 4], [6, 8]]))
    on = dict(bouton_win=b([[8, 4], [6, 6]]),
              claw_win=pd.DataFrame([[7, 4], [5, 6], [12, 10], [12, 1]], index=claws, columns=odors, dtype=float),
              responses=pd.DataFrame([[1, 0], [0, 1]], index=pd.Index([0, 1], name=kc), columns=odors))
    # x_off sums to 84, x_on to 72, the claw states to 57: 27 removed, 12 of it at the boutons
    r = W.sites_row(Wm, on, off, divisive=False)
    assert r == pytest.approx(dict(bouton_share=12 / 27, drive_removed=27 / 84))
    # the responders' claws: 2h's KC 0 (8/7, 6/5), IaA's KC 1 (10/10, 12/1)
    r = W.sites_row(Wm, on, off, divisive=True)
    assert r['divisor_p50'] == pytest.approx((8 / 7 + 6 / 5) / 2)
    assert r['divisor_above10'] == pytest.approx(0.25)


def test_tables_planted(tmp_path, monkeypatch):
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    monkeypatch.setattr(F, 'e134_dir', lambda: tmp_path / 'no_e134')
    monkeypatch.setattr(W, 'data_targets', lambda: {})
    monkeypatch.setattr(B, 'held_rate', lambda panel, lock=None: 0.0915)
    monkeypatch.setattr(R, 'require_lock', lambda *a, **k: {})
    prev = F.MEASURE
    F.set_measure('strength')
    try:
        _plant(tmp_path)
        out = W.main(score_dir=tmp_path, fit_dir=tmp_path, with_apl=False)
    finally:
        F.set_measure(prev)
    md = out.read_text()
    for s in ('## The stage table', '## The classes', 'no better than the permuted control', 'APL→bouton ×0.154',
              'KC–KC +0.0007 (not determined)', 'real: response rate ρ, APL on, 0.30', '0.010 vs 0.040',
              'half to half', 'kiwi mixtures (beside)', "kiwi mixtures, the components' first block",
              '2.378 (-0.047; -0.047)', '2.500 (-0.030; '):
        assert s in md, s
    v = pd.read_csv(tmp_path / 'writeup_values.csv')
    pick = lambda **kw: v[np.logical_and.reduce([v[k] == x for k, x in kw.items()])].value
    assert pick(arm='A0', stage='four_stages', statistic='change').iloc[0] == pytest.approx(-0.010)
    assert pick(arm='A0', stage='four_stages', statistic='stage_error_S_fit').iloc[0] == pytest.approx(0.9572)
    assert pick(arm='family', stage='four_stages', statistic='stage_error_S_fit').iloc[0] == pytest.approx(0.9747)
    assert pick(arm='A0', stage='factor', statistic='aplpn').iloc[0] == pytest.approx(0.154)
    assert pick(arm='A0', stage='kc_soma', statistic='n_flagged').iloc[0] == 1
    assert pick(arm='yardstick_half_to_half', stage='four_stages', statistic='stage_error').iloc[0] == \
        pytest.approx(0.8)
    assert pick(arm=F.FAMILY, panel='mixtures/kmix/on', half='components_second_block',
                statistic='strength_mean').iloc[0] == pytest.approx(2.425)
