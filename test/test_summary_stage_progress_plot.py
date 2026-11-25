"""Tests of the Summary's stage-target progress plot (reports/summary_stage_progress_plot.py) on planted points: the
stage error's root mean square, each experiment's point (lowest filled, else lowest hollow), the staircase over filled
points only, no point for undrawn variants, and the figure's data; no shared-directory reads."""
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
