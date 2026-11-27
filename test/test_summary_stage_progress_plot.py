"""Tests of the Summary's stage-target progress plot (reports/summary_stage_progress_plot.py) on planted points: the
stage error's root mean square, each experiment's point (lowest filled, else lowest hollow), the staircase over filled
points only, no point for undrawn variants, and the figure's data; each adapter from E139 on, on small planted tables
in its experiment's layout; the table route against E139's own classes on the shared tables (skipped where they are not
mounted)."""
import html
import json
import re
from datetime import datetime

import numpy as np
import pandas as pd
import pytest

from reports import summary_progress_plot as SP
from reports import summary_stage_progress_plot as S


def test_per_wiring():
    rows = pd.DataFrame([dict(wiring=w, stage=s, view_strength=v) for w, vals in (('real', (1, 1, 1, 1)),
                                                                                   ('P1', (2, 2, 2, 2)))
                         for s, v in zip(S.STAGES, vals)])
    p = S._per_wiring(rows)
    assert p.loc['P1', 'stage_error'] == pytest.approx(2.0) and p.kc_soma.mean() == pytest.approx(1.5)
    rows.loc[(rows.wiring == 'real') & (rows.stage == 'pn_dendrite'), 'view_strength'] = 3.0
    assert S._per_wiring(rows).loc['real', 'stage_error'] == pytest.approx(np.sqrt((9 + 3) / 4))
    assert S._per_wiring(rows[rows.stage != 'kc_soma']).empty


def test_collect_and_figure(monkeypatch):
    idx = pd.DataFrame([dict(exp='E134', name='the models as populations (design)', completed=datetime(2026, 11, 24, 15)),
                        dict(exp='E135', name='the antennal lobe (design)', completed=datetime(2026, 11, 25, 7))])
    monkeypatch.setattr(SP, 'experiment_index', lambda: idx)
    P = S._point
    monkeypatch.setattr(S, 'ADAPTERS', {
        'E134': lambda: [P('E134', 'family', 'stage_error', 0.92, 'all'), P('E134', 'family', 'kc_soma', 0.84, 'all'),
                         P('E134', 'reference', 'stage_error', np.nan, 'one wiring', drawn=False)],
        'E135': lambda: [P('E135', 'imposed', 'stage_error', 0.97, 'all', filled=False, why='dendrites set'),
                         P('E135', 'route 2', 'stage_error', 0.90, 'real', drawn=False, why='one wiring'),
                         P('E135', 'imposed', 'kc_soma', 1.06, 'all'), P('E135', 'R1', 'kc_soma', 1.56, 'all')],
        'E136': lambda: [P('E136', 'A0', 'stage_error', 0.908, 'T', filled=False, why='no class', arm=1.039,
                           family=1.049, change=-0.010, se=0.010, bar=-0.020, cls='no class', control_fall=0.04,
                           where='T'),
                         P('E136', 'A3', 'stage_error', 0.95, 'T', arm=1.08, family=1.049, change=0.03, se=0.01,
                           bar=-0.02, cls='does not improve', control_fall=0.04, where='T')]})
    pts = S.collect()
    assert set(pts.exp) == {'E134', 'E135'}                 # E136 is not completed
    assert pts[pts.exp == 'E134'].name.iloc[0] == 'the models as populations'
    b = pts[pts.best_of_exp].set_index(['exp', 'series'])
    assert ('E134', 'stage_error') in b.index and b.loc[('E134', 'stage_error'), 'new_best']
    # route 2 has no point; the hollow imposed arm is E135's point, and the staircase does not step on it
    assert b.loc[('E135', 'stage_error'), 'variant'] == 'imposed' and not b.loc[('E135', 'stage_error'), 'filled']
    assert not b.loc[('E135', 'stage_error'), 'new_best']
    assert b.loc[('E135', 'stage_error'), 'running_best'] == pytest.approx(0.92)
    assert b.loc[('E135', 'kc_soma'), 'variant'] == 'imposed'
    pts = S.collect(include=['E136'])
    e = pts[(pts.exp == 'E136') & pts.best_of_exp].iloc[0]
    # the filled A3 is E136's point although the hollow A0 is lower; it is not a new best (0.95 > 0.92)
    assert e.variant == 'A3' and e['name'].endswith('(running)') and not e.new_best
    fig = S.plot_html(pts, dict(stage_error=0.32, stage_min=0.30, stage_max=0.34, apl_loss=0.039, apl_pattern=0.029),
                      dict(stage_error=0.918, kc_soma=0.84))
    panels = json.loads(html.unescape(re.search(r"data-panels='(.*?)'>", fig, re.S).group(1)))
    assert [p['key'] for p in panels] == ['stage']
    st = {s['key']: s for s in panels[0]['series']}
    pe = st['stage_error']['points']
    assert [p['value'] for p in pe] == [0.92, 0.97, 0.95] and [p['filled'] for p in pe] == [True, False, True]
    assert pe[2]['fit']['cls'] == 'does not improve' and pe[2]['level'] == pytest.approx(0.918)
    assert [r['y'] for r in panels[0]['refs']] == [1.0, 0.32]


# ---------------------------------------------------------------------------
# the adapters from E139 on, on small planted tables in the layout of each experiment's results folder

LV = dict(stage_error=0.9, pn_dendrite=0.9, pn_bouton=1.0, kc_claw=0.8, kc_soma=0.8, validation2=1.0, cmix=0.7,
          apl16_on=1.0, apl16_off=1.1, apl_loss=0.06, apl_pattern=0.05)
SPREADS = dict(stage=0.02, pn_dendrite=0.04, pn_bouton=0.015, kc_claw=0.017, kc_soma=0.057,
               **{'validation2/validation2/on': 0.04, 'mixtures/cmix/on': 0.03, 'apl16/apl16/on': 0.05,
                  'apl16/apl16/off': 0.09})


@pytest.fixture
def results(monkeypatch, tmp_path):
    monkeypatch.setattr(S, 'RESULTS', tmp_path)
    monkeypatch.setattr(S, 'levels', lambda: dict(LV))
    monkeypatch.setattr(S, '_e136_dir', lambda: tmp_path / 'E136_stage_population_fit' / 'measure_strength')
    return tmp_path


def _cmp(change, se=0.01, arm=1.0, ref=1.1):
    return dict(n=4, complete=True, changes=[change] * 4, change=change, se=se, spread=2 * se, arm=arm, ref=ref)


def _class_rec(change, cls='megamat only', **extra):
    return dict(megamat=_cmp(change), stages={s: _cmp(change / 2) for s in S.STAGES},
                panels={k: _cmp(0.01) for k in S.E136_PANELS.values()}, cls=cls, reasons=[], control_fall=None,
                **extra)


def _mega(rows, key='variant', halves=('all', 'held')):
    """score_megamat rows: per variant a value per wiring ('stage' and every stage), plus the summary rows."""
    out = []
    for name, val in rows.items():
        for half in halves:
            for w in [*S.WIRINGS, 'mean_T', 'common']:
                v = val(w) if callable(val) else val
                out.append({key: name, 'half': half, 'wiring': w, 'stage': v, **{s: v for s in S.STAGES}})
    return pd.DataFrame(out)


def _panels(rows, key='arm', wirings=S.WIRINGS):
    return pd.DataFrame([{key: name, 'stage': 'kc_soma', 'panel': pan, 'condition': cond, 'wiring': w,
                          'error': v(w) if callable(v) else v}
                         for name, v in rows.items() for pan, cond in S.PANEL_SERIES.values()
                         for w in [*wirings, 'mean']])


def _apl(rows, key='arm', wirings=S.WIRINGS, flags=None):
    flags = flags or {}
    return pd.DataFrame([{key: name, 'wiring': w, 'cal_loss': v(w) if callable(v) else v,
                          'cal_pattern': (v(w) if callable(v) else v) - 0.01,
                          'realism_flag': flags.get((name, w), np.nan)}
                         for name, v in rows.items() for w in [*wirings, 'mean']])


def _write(folder, files):
    folder.mkdir(parents=True, exist_ok=True)
    for name, obj in files.items():
        if isinstance(obj, pd.DataFrame):
            obj.to_csv(folder / name, index=False)
        else:
            (folder / name).write_text(json.dumps(obj))


def _pt(pts, variant, series):
    x = [p for p in pts if p['variant'] == variant and p['series'] == series]
    assert len(x) == 1, (variant, series, len(x))
    return x[0]


def test_e139_unfitted_rungs_on_all_wirings_and_fitted_settings_on_t(results):
    cl = dict(variants=['F_inhadd1000', 'B_inhadd1000', 'K_inhadd1000'], fitted=['A0_inhadd314'],
              reference='family', beside_arm='A1', spreads=SPREADS,
              classes={'F_inhadd1000': _class_rec(-0.09, setting='F'), 'B_inhadd1000': _class_rec(-0.11, setting='B'),
                       'K_inhadd1000': _class_rec(-0.10, setting='K'),
                       'A0_inhadd314': _class_rec(-0.06, cls='no class', setting='A0')})
    cl['classes']['A0_inhadd314']['reasons'] = ['a hold fails']
    # F is 0.1 below the family on every wiring against all flies (0.2 on the held-out half: not read for F)
    mega = _mega({'family': 1.0, 'F_inhadd1000': lambda w: 0.9}, halves=('all',))
    mega = pd.concat([mega, _mega({'family': 1.1, 'F_inhadd1000': 0.9}, halves=('held',))])
    apl = _apl({'family': 0.06, 'F_inhadd1000': 0.07,
                # B: +0.01 on T, +0.05 on S (S is not read for a fitted setting)
                'B_inhadd1000': lambda w: 0.07 if w in S.T_WIRINGS else 0.11, 'K_inhadd1000': 0.06},
               flags={('K_inhadd1000', 'P3'): 'response rate ρ, APL on, 0.36 (below 0.37)'})
    _write(results / 'E139_normalization_ladder',
           {'classes.json': cl, 'score_megamat.csv': mega, 'score_apl16.csv': apl,
            'score_panels.csv': _panels({'family': 1.0, 'F_inhadd1000': 1.03})})
    pts = S.e139_points()
    f = _pt(pts, "the family's calyx, pn_inhadd 1000", 'stage_error')
    assert f['value'] == pytest.approx(0.9 - 0.1) and f['drawn'] and f['filled'] and f['reading'] == S.ALL8
    fp = json.loads(f['popup'])
    assert fp['bar'] is None and fp['family'] == pytest.approx(1.0) and 'read on T' in fp['cls']
    assert _pt(pts, "the family's calyx, pn_inhadd 1000", 'validation2')['value'] == pytest.approx(1.03)
    assert _pt(pts, "the family's calyx, pn_inhadd 1000", 'apl_loss')['value'] == pytest.approx(0.07)
    b = _pt(pts, 'APL→bouton ×0.155, pn_inhadd 1000', 'stage_error')
    assert b['value'] == pytest.approx(0.9 - 0.11) and b['filled'] and b['reading'] == S.T_HELD
    assert json.loads(b['popup'])['bar'] == pytest.approx(-0.02)
    assert _pt(pts, 'APL→bouton ×0.155, pn_inhadd 1000', 'kc_soma')['value'] == pytest.approx(0.8 - 0.055)
    assert _pt(pts, 'APL→bouton ×0.155, pn_inhadd 1000', 'apl_loss')['value'] == pytest.approx(0.07)
    a0 = [p for p in pts if p['variant'] == 'A0 refitted at pn_inhadd 314']
    assert a0 and all(p['drawn'] and not p['filled'] and 'a hold fails' in p['why'] for p in a0)
    assert not any(p['drawn'] for p in pts if p['variant'] == 'KC–KC ×3.24, pn_inhadd 1000')
    assert not _pt(pts, 'A1 (beside)', 'stage_error')['drawn']


def test_e141_arms_on_t_and_the_arms_without_a_point(results):
    cl = dict(variants=['D1.8_inhadd1000', 'S_inhadd1000', 'T_inhadd1000'], reference='family',
              beside_e139=['B_inhadd1000'], beside_arm='A1',
              classes={'D1.8_inhadd1000': _class_rec(-0.12, cls='holds', flagged_wirings={}),
                       'S_inhadd1000': _class_rec(-0.10, cls='no class',
                                                  flagged_wirings={'P3': 'response rate ρ, APL on, 0.36'}),
                       'T_inhadd1000': dict(cls='no class', reasons=['a hold fails: T_inhadd1000_megamat_real'])})
    _write(results / 'E141_kckc_rise_forms',
           {'classes.json': cl, 'score_apl16.csv': _apl({'family': 0.06, 'D1.8_inhadd1000': 0.075})})
    pts = S.e141_points()
    d = _pt(pts, "D1.8: KC–KC ×1.8, on E139's B at 1000", 'stage_error')
    assert d['drawn'] and d['filled'] and d['value'] == pytest.approx(0.9 - 0.12)
    assert json.loads(d['popup'])['cls'] == 'holds'
    assert _pt(pts, "D1.8: KC–KC ×1.8, on E139's B at 1000", 'apl_pattern')['value'] == pytest.approx(0.05 + 0.015)
    s = [p for p in pts if p['variant'].startswith('S:')]
    assert s and not any(p['drawn'] for p in s) and 'realism' in s[0]['why']
    t = _pt(pts, "T: one neuron per KC type, on E139's B at 1000", 'stage_error')
    assert not t['drawn'] and 'a hold fails' in t['why']
    assert not _pt(pts, 'B_inhadd1000 (beside)', 'stage_error')['drawn']


def test_e142_imposed_arms_hollow_at_the_dendrites_and_without_a_class(results):
    holds = [dict(panel='megamat', variant=v, n_read=n) for v, n in
             (('E142_a1stage_imposed', 8), ('E142_imposed_a1', 8), ('E142_imposed_aplpn_none', 1),
              ('E142_a1stage_imposed_glom_perm', 8))]
    sc = dict(variants=['family', 'A1', 'E135_imposed', 'E142_imposed_a1', 'E142_a1stage_imposed',
                        'E142_a1stage_imposed_glom_perm', 'E142_imposed_aplpn_none'],
              classes=dict(dendrites_carry=dict(arm='E142_a1stage_imposed', cls='does not hold'),
                           apl_hides_step=dict(kept=[], cls='does not hold')),
              holds_failed={'E142_imposed_a1': [f'megamat_imposed_a1_{w}' for w in S.WIRINGS]},
              holds_table=holds, spreads=dict(stage=0.02, **{s: 0.02 for s in S.STAGES}))
    # the main arm: +0.2 at every stage on T against the held-out half (and -0.5 on S, not read)
    main = lambda w: 1.2 if w in S.T_WIRINGS else 0.5
    mega = _mega({'family': 1.0, 'E142_a1stage_imposed': main, 'E142_imposed_a1': 0.95}, halves=('held',))
    mega = pd.concat([mega, _mega({'family': 2.0, 'E142_a1stage_imposed': 0.0}, halves=('all',))])
    pan = _panels({'family': 1.0, 'E142_a1stage_imposed': 1.1, 'E142_imposed_a1': 1.2})
    pan = pan[pan.panel == 'apl16']
    apl = _apl({'family': 0.06, 'E142_a1stage_imposed': 0.08, 'E142_imposed_a1': 0.09}, key='variant')
    apl.loc[apl.wiring == 'mean', 'cal_loss'] = 9.0               # the mean row is not a wiring
    _write(results / 'E142_dendrite_bouton_step',
           {'step1_classes.json': sc, 'score_step1_megamat.csv': mega, 'score_step1_panels.csv': pan,
            'score_step1_apl16.csv': apl})
    _write(S._e136_dir(), {'classes.json': dict(spreads=SPREADS)})
    pts = S.e142_points()
    m = "the dendrites imposed on A1's own PN stage (the amended main arm)"
    bo = _pt(pts, m, 'pn_bouton')
    assert bo['filled'] and bo['value'] == pytest.approx(1.0 + 0.2)
    assert 'does not hold' in json.loads(bo['popup'])['cls']
    for k in ('pn_dendrite', 'stage_error'):
        p = _pt(pts, m, k)
        assert p['drawn'] and not p['filled'] and 'imposed' in p['why']
    assert _pt(pts, m, 'apl16_on')['value'] == pytest.approx(1.0 + 0.1) and _pt(pts, m, 'apl16_on')['filled']
    assert _pt(pts, m, 'apl_loss')['value'] == pytest.approx(0.06 + 0.02)
    assert not [p for p in pts if p['variant'] == m and p['series'] in ('validation2', 'cmix')]
    lock = [p for p in pts if p['variant'].endswith("(the lock's main arm)")]
    assert lock and all(p['drawn'] and not p['filled'] and p['why'].startswith('no class') for p in lock)
    assert not _pt(pts, 'E135_imposed', 'stage_error')['drawn']
    assert 'permuted' in _pt(pts, 'E142_a1stage_imposed_glom_perm', 'stage_error')['why']
    assert 'read on 1 of 8' in _pt(pts, 'E142_imposed_aplpn_none', 'stage_error')['why']


def test_e144_settings_against_their_calyx_family(results):
    ws = ['real', 'P1', 'P2', 'P3']
    rec = lambda ch, cls, cx, st: dict(_class_rec(ch, cls=cls), calyx=cx, setting=st, wirings=ws,
                                       reference=f'{cx}_family', flagged_wirings={})
    cl = dict(classes={'L_A1': rec(-0.15, "strength at another view's cost", 'L', 'A1'),
                       'L_D5.8': dict(rec(-0.16, 'no class', 'L', 'D5.8'), reasons=['a hold fails: L_D5.8_megamat_real']),
                       'R_D5.8': rec(-0.12, 'holds', 'R', 'D5.8'),
                       'H_A1': dict(rec(-0.11, 'megamat only', 'H', 'A1'), wirings=['real'])})
    # the male calyces' wirings are named as the hemibrain's are; each setting is read against its own calyx's family
    apl = _apl({'L_family': 0.10, 'L_A1': 0.09, 'R_family': 0.05, 'R_D5.8': 0.07, 'L_D5.8': 0.10}, key='variant',
               wirings=ws)
    _write(results / 'E144_a1_malecns_calyces', {'classes.json': cl, 'score_apl16.csv': apl})
    pts = S.e144_points()
    la = "the male CNS's left calyx: A1 (B with the calyx's KC–KC ×3.24)"
    p = _pt(pts, la, 'stage_error')
    assert p['filled'] and p['value'] == pytest.approx(0.9 - 0.15) and "left calyx's 4 wirings" in p['reading']
    assert _pt(pts, la, 'apl_loss')['value'] == pytest.approx(0.06 - 0.01)
    assert _pt(pts, "the male CNS's right calyx: D5.8 (B with the calyx's KC–KC ×5.83)", 'apl_loss')['value'] == \
        pytest.approx(0.06 + 0.02)
    ld = [p for p in pts if p['variant'].startswith("the male CNS's left calyx: D5.8")]
    assert ld and all(p['drawn'] and not p['filled'] for p in ld)
    h = _pt(pts, "the hemibrain with 'argmax' claws: A1 (B with the calyx's KC–KC ×3.24)", 'stage_error')
    assert not h['drawn'] and h['why'] == 'read on 1 wiring'
    assert not _pt(pts, "the male CNS's left calyx: the family", 'stage_error')['drawn']


def test_compare_is_e136s():
    a = {w: 1.0 + 0.01 * i for i, w in enumerate(S.T_WIRINGS)}
    b = {w: 1.1 for w in S.T_WIRINGS}
    c = S._compare(a, b, 0.02, S.T_WIRINGS)
    ch = [a[w] - b[w] for w in S.T_WIRINGS]
    assert c['change'] == pytest.approx(np.mean(ch)) and c['se'] == pytest.approx(max(0.01, np.std(ch, ddof=1) / 2))
    assert not S._compare(a, b, 0.02, S.WIRINGS)['complete']


@pytest.mark.skipif(not (S.RESULTS / 'E139_normalization_ladder' / 'classes.json').exists(),
                    reason='the shared results are not mounted')
def test_table_route_reproduces_e139s_classes_on_the_shared_tables():
    """On E139's own tables, the comparisons read from the score tables on T against the held-out half (E142's route)
    are those of its classes.json (E136's `compare`)."""
    d = S.RESULTS / 'E139_normalization_ladder'
    cl = json.loads((d / 'classes.json').read_text())
    mega, pan = pd.read_csv(d / 'score_megamat.csv'), pd.read_csv(d / 'score_panels.csv')
    for v in ('B_inhadd1000', 'F_inhadd314'):
        comps = S._table_comps(mega, pan, pd.DataFrame(), v, 'family', S.T_WIRINGS, 'held', cl['spreads'], 'T', 'T',
                               'T')
        for k, (c, _) in S._class_comps(cl['classes'][v]).items():
            assert comps[k][0]['change'] == pytest.approx(c['change'], abs=1e-12), (v, k)
            assert comps[k][0]['se'] == pytest.approx(c['se'], abs=1e-12), (v, k)
