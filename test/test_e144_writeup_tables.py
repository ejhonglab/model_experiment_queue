"""Tests of E144's write-up tables (analyses/E144_writeup_tables.py) on planted score files; no shared-directory reads
and no simulation."""
import pandas as pd

from analyses import E144_writeup_tables as T


def _cmp(change, spread=0.05, worse=False, ref=1.0):
    return dict(complete=True, change=change, spread=spread, worse_than_spread=worse, ref=ref)


def _classes():
    ws = ['real', 'P1', 'P2', 'P3']
    def rec(cls, ws_, **kw):
        return dict(dict(cls=cls, wirings=ws_, megamat=_cmp(-0.12, 0.02), bar=0.02, stage_worse=[], guard_worse=[],
                         panel_worse=[], reasons=[], stages={s: _cmp(-0.1) for s in T.STAGES},
                         panels={k: _cmp(0.01) for k in T.PANEL_KEYS[:4]},
                         panels_beside={T.PANEL_KEYS[4]: _cmp(0.2, 0.1)}, guard={s: _cmp(0.0) for s in T.STAGES}),
                    **kw)
    held = rec("strength at another view's cost", ws[1:], guard_worse=['kc_soma'])
    failed = rec('no class', ws, reasons=['a hold fails: L_D5.8_megamat_real'], on_wirings_held=held,
                 price_on_apl_effect=dict(wirings=ws, price=True, price_value=0.0058))
    worse = rec('megamat only', ws, panel_worse=['apl16/apl16/on', 'apl16/apl16/off'],
                price_on_apl_effect=dict(wirings=ws, price=False, price_value=0.001))
    lad = {s: dict(change=c) for s, c in zip(('B', 'D1.8', 'A1', 'D5.8'), (0.06, -0.04, -0.18, -0.24))}
    cal = {c: dict(dose_ladder=lad) for c in ('H', 'L', 'R')}
    return dict(classes={'L_D5.8': failed, 'R_A1': worse}, calyces=cal)


def test_helpers():
    assert T._rng(pd.Series([0.16, 0.1601])) == '0.160'
    assert T._rng(pd.Series([0.1, 0.2]), '.2f') == '0.10 to 0.20'
    assert T._names(['apl16/apl16/on', 'kc_soma']) == 'APL16, APL on; KC somata'
    assert T._names([]) == 'none'
    assert T.label('L_D5.8@held') == 'L, D5.8: B with KC–KC ×5.83'
    assert len(T.variants()) == 19 and 'L_Dabs' in T.variants() and 'R_Dabs' not in T.variants()


def test_classes_and_ladder_tables():
    d = dict(classes=_classes(), shares=pd.DataFrame([dict(arm='L_D5.8', share_scored=1.0),
                                                      dict(arm='L_D5.8@held', share_scored=0.327),
                                                      dict(arm='R_A1', share_scored=0.261)]))
    V = T.Values()
    rows = [r for r in T.classes_table(d, V) if r.startswith('| L, D5.8') or r.startswith('| R, A1')]
    assert len(rows) == 3
    assert 'no class (a hold fails: L_D5.8_megamat_real)' in rows[0] and 'yes (+0.0058)' in rows[0]
    assert 'on the wirings whose holds held' in rows[1] and '| P1, P2, P3 |' in rows[1] and '0.327' in rows[1]
    assert 'KC somata' in rows[1] and rows[1].rstrip().endswith('| – | 0.327 |')
    assert 'APL16, APL on; APL16, APL off' in rows[2] and 'no (+0.0010)' in rows[2]
    v = V.frame().set_index(['variant', 'statistic'])['value']
    assert v[('L_D5.8@held', 'resampled_share_scored')] == 0.327 and v[('R_A1', 'change')] == -0.12
    lad = T.ladder_table(d, V)
    assert any(r.startswith("| the hemibrain ('first' claws; E141) | four held out | +0.080") for r in lad)
    four = [r for r in lad if r.startswith('| L: the male CNS') and '| four |' in r]
    assert len(four) == 1 and four[0].endswith('| – | – | – | +0.010 |')     # a dose's four wirings, from its class
    assert any('| real | +0.060 | -0.040 | -0.180 | -0.240 |' in r for r in lad)
