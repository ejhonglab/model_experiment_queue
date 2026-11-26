"""Tests of E139's write-up tables (analyses/E139_writeup_tables.py) on planted score files; no shared-directory reads
and no simulation."""
import numpy as np
import pandas as pd
import pytest

from analyses import E136_writeup_tables as W
from analyses import E139_writeup_tables as T
from experiments import E136_stage_population_fit as F


def test_variant_order_and_labels():
    names = ['A1', 'K_inhadd1000', 'F_inhadd3141', 'family', 'B_inhadd314', 'F_inhadd100', 'B_inhadd1000']
    assert T.variant_order(names) == ['family', 'F_inhadd100', 'F_inhadd3141', 'B_inhadd314', 'B_inhadd1000',
                                      'K_inhadd1000', 'A1']
    assert T.variant_label('F_inhadd314') == "the family's calyx, pn_inhadd 314"
    assert T.variant_label('B_inhadd1000') == 'APL→bouton ×0.155, pn_inhadd 1000'
    assert T.variant_label('K_inhadd1000') == 'KC–KC ×3.24, pn_inhadd 1000'


def test_as_e136_restores():
    order, names = W.ARM_ORDER, dict(W.ARM_NAMES)
    with T.as_e136(['family', 'F_inhadd100']):
        assert W.ARM_ORDER == ('family', 'F_inhadd100')
        assert W.arm_label('F_inhadd100', {}) == "the family's calyx, pn_inhadd 100"
    assert W.ARM_ORDER == order and W.ARM_NAMES == names
    with pytest.raises(RuntimeError):
        with T.as_e136(['family']):
            raise RuntimeError
    assert W.ARM_ORDER == order and W.ARM_NAMES == names


def test_read_step():
    assert T.read_step(dict(setting='F', panel='megamat')) == '1'
    assert T.read_step(dict(setting='B', panel='apl16')) == '2'
    assert T.read_step(dict(setting='K', panel='mixtures')) == '2'
    assert T.read_step(dict(setting='B', panel='extended')) == 'beside'
    assert T.read_step(dict(setting='A0', panel='megamat')) == 'fit'


def _classes():
    mg = dict(complete=True, change=-0.09, se=0.01, improves=True)
    price = dict(complete=True, price=True, price_value=0.0078,
                 by_target={t: dict(change_T=0.007, change_8=0.008) for t in ('z8', 'halfA', 'halfB')})
    return dict(classes={'F_inhadd1000': dict(cls='megamat only', megamat=mg, stage_worse=[], guard_worse=[],
                                              panel_worse=['validation2/validation2/on'], holds_failed=[],
                                              price_on_apl_effect=price, serves_both=False, reasons=[])},
                step3_candidates=[], named_rungs={'F': dict(values={'100': 0.98, '1000': 0.92}, named=1000.0)})


def _mega():
    rows = []
    for v, base in (('family', 1.05), ('F_inhadd1000', 0.96)):
        for half in ('held', 'all', 'fit'):
            for w in F.WIRINGS:
                rows.append(dict(arm=v, half=half, wiring=w, stage=base))
    return pd.DataFrame(rows)


def test_classes_and_halves_tables():
    d = dict(classes=_classes(), mega=_mega())
    V = W.Values()
    lines = T.classes_table(d, V)
    row = [x for x in lines if x.startswith("| the family's calyx, pn_inhadd 1000")][0]
    assert '-0.090 (0.010)' in row and '+0.0078 (yes)' in row and 'megamat only' in row and 'validation2' in row
    lines = T.halves_table(d, V)
    assert any('1000 (100 0.980, 1000 0.920)' in x for x in lines)
    f = V.frame().set_index(['arm', 'stage', 'statistic', 'half']).value
    assert f[('F_inhadd1000', 'four_stages', 'change', 'held')] == pytest.approx(-0.09)
    assert f[('F_inhadd1000', 'kc_soma', 'price_value', 'T')] == pytest.approx(0.0078)
    assert f[('F_inhadd1000', 'four_stages', 'stage_error_mean', 'fit_S')] == pytest.approx(0.96)
    assert f[('setting_F', 'four_stages', 'named_rung', 'fit_S')] == 1000.0


def test_price_table():
    rows = []
    for v, pat in (('family', 0.056), ('F_inhadd1000', 0.064)):
        for w in F.WIRINGS:
            rows.append(dict(arm=v, wiring=w, price_mean_effect=0.18,
                             **{f'price_{t}_{k}': (pat if k == 'pattern' else 0.5) for t in ('z8', 'halfA', 'halfB')
                                for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'rho', 'pearson')}))
    d = dict(classes=_classes(), apl=pd.DataFrame(rows), mega=_mega())
    V = W.Values()
    lines = T.price_table(d, V)
    row = [x for x in lines if x.startswith("| the family's calyx, pn_inhadd 1000")][0]
    assert '0.064 / 0.064' in row and '+0.0070' in row
    f = V.frame().set_index(['arm', 'statistic']).value
    assert f[('F_inhadd1000', 'pattern_z8_change_8')] == pytest.approx(0.008)
    assert f[('family', 'pattern_z8_mean_T')] == pytest.approx(0.056)
    assert not np.isnan(f[('family', 'mean_effect_mean_8')])
