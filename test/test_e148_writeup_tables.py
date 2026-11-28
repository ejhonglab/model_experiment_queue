"""Tests of E148's write-up tables (analyses/E148_writeup_tables.py) on planted records; no shared-directory reads and
no simulation."""
import json

import numpy as np
import pandas as pd

from analyses import E148_writeup_tables as T


def test_helpers():
    assert T._rng(pd.Series([0.16, 0.1601])) == '0.160'
    assert T._rng(pd.Series([0.1, 0.2]), '.2f') == '0.10 to 0.20'
    assert T._rng(pd.Series([np.nan])) == '–'
    assert T._truth(True) and T._truth('True') and not T._truth('False') and not T._truth(np.nan)
    assert T._views_of('pn_dendrite')[-1] == 'matrix' and 'matrix' not in T._views_of('kc_soma')


def _probe():
    return pd.DataFrame([
        dict(statistic='consistency_L0', cell='A1', printed='0.34', value='0.338', value_held='0.357', reproduced=True,
             note=np.nan),
        dict(statistic='holds_part_half_a', cell='X', printed='[1, 116]', value='[1, 116]', value_held='[1, 116]',
             reproduced=True, note=np.nan),
        dict(statistic='claws_worse_among_lowering', cell='a', printed='[429, 456]', value='[430, 456]',
             value_held='[426, 452]', reproduced=False, note="429 with the claws' spread rounded to 0.017")])


def test_probe_table():
    V = T.W.Values()
    L = T.probe_table(dict(step0_probe=_probe()), V)
    rows = [r for r in L if r.startswith('| ')]
    assert '| consistency_L0 | A1 | 0.34 | 0.338 | 0.357 | yes |' in rows[1]
    assert "120 held evaluations" in rows[1]
    assert '| holds_part_half_a | X | 1 / 116 | 1 / 116 | 1 / 116 | yes |' in rows[2]
    assert '| claws_worse_among_lowering | a | 429 / 456 | 430 / 456 | 426 / 452 | no |' in rows[3]
    assert any('2 of 3 printed values reproduced' in r for r in L)
    f = V.frame()
    assert set(f.statistic) >= {'probe_consistency_L0', 'probe_holds_part_half_a_1'}


def _ids():
    chk = dict(checks=[dict(half='a', eval='e0', own=0.975, once=0.915, diff=0.06),
                       dict(half='b', eval='e0', own=1.04, once=1.051, diff=0.011)],
               max_diff_every_eval={'a/stage_error_strength': 0.066, 'b/stage_error_strength': 0.023})
    return pd.DataFrame([
        dict(cell='F', identity="rescore reproduced from the views (E136's own reading)", value=5.5e-16,
             expected=1e-9, passed=True, needs_views=True, note='', kind='identity'),
        dict(cell='F', identity='consistency F from the views (own reading)', value=0.8016, expected=0.802,
             passed=True, needs_views=True, note='', kind='identity'),
        dict(cell='F', identity='consistency F from the views (all reading)', value=0.346, expected=0.802,
             passed=np.nan, needs_views=True, note='beside', kind='beside'),
        dict(cell='F', identity="E143's read-once check", value=0.06, expected=0.02, passed=False, needs_views=True,
             note=json.dumps(chk), kind='check')])


def test_identities_and_read_once():
    V = T.W.Values()
    d = dict(identities=_ids())
    L = T.identities_table(d, V)
    assert any('| F | rescore reproduced' in r and '5.50e-16' in r and '| yes |' in r for r in L)
    assert any('(all reading)' in r and '| beside |' in r for r in L)
    R = T.read_once_table(d, V)
    assert '| F | a | e0 | 0.975 | 0.915 | 0.060 | 0.066 | no |' in R
    assert '| F | b | e0 | 1.040 | 1.051 | 0.011 | 0.023 | yes |' in R
    assert '| F | 0.802 | 0.346 |' in R


def test_decided_table():
    rows = []
    for k in T.SETS:
        rows.append(dict(stage='kc_soma', view='sharing', set=k, read_once=k.startswith(('settings', 'inhadd', 'family')),
                         n_models=9, n_splits=3, of_splits=3, median=0.97, min=0.9, e136=0.97, decided=True,
                         range_all_flies=0.5, spread=0.08, no_reading=False))
        rows.append(dict(stage='pn_dendrite', view='strength', set=k, read_once=False, n_models=9, n_splits=1,
                         of_splits=35, median=0.42, min=0.42, e136=0.42, decided=False, range_all_flies=0.01,
                         spread=0.04, no_reading=True))
    V = T.W.Values()
    L = T.decided_table(dict(views_decided_summary=pd.DataFrame(rows)), V)
    soma = next(r for r in L if r.startswith('| KC somata | sharing'))
    assert soma.count('**0.97 / 0.90**') == len(T.SETS)
    den = next(r for r in L if r.startswith('| PN dendrites | strength'))
    assert '0.42 (one split) nr' in den
    V.frame()


def test_losses_picks_and_calyx():
    summ, picks = [], []
    for c in T.CELLS:
        for L_ in ('L0', 'L1'):
            summ.append(dict(set=c, loss=L_, read_once=False, n_models=120, n_joint=9, median=0.9 if L_ == 'L1' else 0.6,
                             min=0.8 if L_ == 'L1' else 0.1, e136=0.85, decided=L_ == 'L1', pick_e136=46,
                             regret_e136=0.04, transfers_e136=True, holds_part_e136=False, regret_median=0.05,
                             share_transfers=1.0, share_holds_part=0.0))
            picks.append(dict(cell=c, loss=L_, split='e136', pick=46 if L_ == 'L1' else 122, regret=0.04,
                              transfers=True, holds_part_other_half=c == 'A1', fit_loss=0.9, other_loss=0.95,
                              other_min=0.93, start_other=1.04, aplpn=0.225, kckc=3.44, inhadd=43.4, pn_inhadd=1362.0,
                              inhsc=1.17, loss_sim=0.9))
    summ += [dict(set=p, loss=L_, read_once=False, n_models=250, n_joint=9, median=0.97, min=0.93, e136=0.97,
                  decided=True) for p in ('pool_S', 'pool_T') for L_ in ('L0', 'L1')]
    d = dict(losses_summary=pd.DataFrame(summ), losses_picks=pd.DataFrame(picks),
             classes=dict(losses=dict(L0=dict(decided=False), L1=dict(decided=True))))
    V = T.W.Values()
    L = T.losses_table(d, V)
    l1 = next(r for r in L if r.startswith('| L1 |'))
    assert l1.count('**0.90 / 0.80 (0.85)**') == 4 and l1.endswith('| yes |')
    P = T.picks_table(d, V)
    a1 = next(r for r in P if r.startswith('| A1 (S, half a) | L1 |'))
    assert '| 46 | ×0.225 | ×3.44 | ×43.4 | 1362 | ×1.17 | 0.040 | yes | yes |' in a1
    ws = T.M.get_cell('A1').wirings
    hist = {'eval': [0, 46, 122], 'z_aplpn': [0.0, np.log(0.225), np.log(0.155)], 'z_kckc': [0.0, np.log(3.44), 0.0],
            'z_inhadd': [0.0, np.log(43.4), 0.0], 'z_inhsc': [0.0, 0.0, 0.0], 'rate_on': [0.0915] * 3,
            'rate_off': [0.205] * 3, 'rate_ratio': [2.24] * 3}
    for w, k in zip(ws, range(4)):
        hist[f'fixed_thr_{w}'] = [30.0 + k, 50.0 + k, 40.0]
        hist[f'wKCKC_scale_{w}'] = [20.0 + k, 80.0 + k, 70.0]
        hist[f'wAPLKC_scale_{w}'] = [0.2, 0.6, 0.7]
        hist[f'wAPLPN_scale_{w}'] = [0.33, 0.075, 0.05]
    d.update(hist={'A1': pd.DataFrame(hist)}, flies_rate=0.0915)
    C = T.calyx_table(d, V)
    row = next(r for r in C if r.startswith('| A1 (S, half a) | the L1 pick |'))
    assert '| 46 | ×0.225 | ×3.44 | ×43.4 | 50.0 to 53.0 | 80.0 to 83.0 | 0.600 | 0.0750 | 0.0915 | 0.2050 | 2.24 |' in row
    fam = next(r for r in C if r.startswith('| A1 (S, half a) | evaluation 0'))
    assert '×1.000 | ×1.00 | ×1.0 |' in fam
    V.frame()


def test_claws_and_prediction():
    rows = []
    for c in T.CELLS:
        for ref, ch in (('e136_a', 0.05), ('e136_b', -0.06)):
            rows.append(dict(stage='kc_claw', reference=ref, flies='', reading="E136's own", cell=c, classed='x',
                             change=ch, se=0.01, spread=0.017, worse_than_spread=ch > 0.017, improves=False))
        for i in range(4):
            rows.append(dict(stage='kc_claw', reference=f'fly{i}', flies='', reading='read once', cell=c, classed='x',
                             change=0.06 if i == 2 else -0.02, se=0.01, spread=0.017, worse_than_spread=i == 2,
                             improves=False))
    cls = dict(claws=dict(read_half=True, half='a', kind="one fly's",
                          fly_sign=dict(fly0=False, fly1=False, fly2=True, fly3=False)),
               reads_calyx=dict(holds=True, decided_over_the_fits=['kc_soma/sharing']),
               loss_with_sharing_decided=['L1'], no_boutons=dict(orders=False))
    V = T.W.Values()
    L = T.claws_table(dict(claws=pd.DataFrame(rows), classes=cls), V)
    a1 = next(r for r in L if r.startswith("| A1 (S, half a) | E136's own |"))
    assert '+0.050 worse | -0.060 | – | – | – | – |' in a1
    assert any("the read half's: yes (the half: a); one fly's" in r for r in L)
    pred = dict(per_stage_view={'kc_soma/sharing': dict(each={c: 0.97 for c in T.CELLS},
                                                        predicted_decided_over_every_cell=True)},
                written='2026-11-27 20:37', host='disco', claws_cost_is_the_read_half='about 9 in 10',
                without_boutons_orders_the_settings='about 4 in 5')
    summ = pd.DataFrame([dict(stage='kc_soma', view='sharing', set=c, decided=True) for c in T.CELLS])
    P = T.prediction_table(dict(prediction=pred, views_decided_summary=summ, classes=cls), V)
    assert '| KC somata | sharing between odors | 0.97 | 0.97 | 0.97 | 0.97 | yes | yes |' in P
    assert "| the claws' cost is the read half's | about 9 in 10 | yes |" in P
    assert '| without the boutons the megamat stage error orders the settings | about 4 in 5 | no |' in P
    V.frame()


def test_main_on_a_planted_folder(tmp_path):
    _probe().to_csv(tmp_path / 'step0_probe.csv', index=False)
    _ids().to_csv(tmp_path / 'identities.csv', index=False)
    (tmp_path / 'classes.json').write_text(json.dumps(dict(cells={c: dict(identities=True, read_once=False,
                                                                           read_by="E136's own") for c in T.CELLS},
                                                           outcomes=['a', 'd', 'f'], time='t', host='h')))
    (tmp_path / 'provenance_score.json').write_text(json.dumps(dict(queue_commit='abc1234', time='t', host='h')))
    out = T.main(tmp_path, histories=False)
    txt = out.read_text()
    assert txt.startswith('# E148 write-up tables') and 'queue abc1234' in txt
    assert '## The outcomes' in txt and 'Outcomes: a, d, f.' in txt
    v = pd.read_csv(tmp_path / 'writeup_values.csv')
    assert {'arm', 'stage', 'panel', 'condition', 'half', 'statistic', 'value'} <= set(v.columns)
    assert ((v.arm == 'experiment') & (v.statistic == 'outcome_a')).any()
