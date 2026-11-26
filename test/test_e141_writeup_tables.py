"""Tests of E141's write-up tables (analyses/E141_writeup_tables.py) on planted score files; no shared-directory reads
and no simulation."""
import pandas as pd
import pytest

from analyses import E136_writeup_tables as W
from analyses import E139_writeup_tables as W139
from analyses import E141_writeup_tables as T
from experiments import E136_stage_population_fit as F


def test_variant_order_and_labels():
    names = ['K_inhadd1000', 'G_inhadd1000', 'A1', 'D5.8_inhadd1000', 'family', 'B_inhadd1000', 'N_inhadd1000',
             'D1.8_inhadd1000', 'P_inhadd1000', 'F_inhadd1000', 'other']
    assert T.variant_order(names) == ['family', 'F_inhadd1000', 'B_inhadd1000', 'D1.8_inhadd1000', 'A1',
                                      'D5.8_inhadd1000', 'N_inhadd1000', 'G_inhadd1000', 'P_inhadd1000',
                                      'K_inhadd1000', 'other']
    assert T.variant_label('D1.8_inhadd1000') == 'D1.8: KC–KC ×1.8'
    assert T.variant_label('D5.8_inhadd1000') == 'D5.8: KC–KC ×5.83'
    assert T.variant_label('N_inhadd1000') == 'N: KC–KC ×3.24 with APL on, ×1 with APL off'
    assert T.variant_label('B_inhadd1000') == "KC–KC ×1 (E139's B at 1000)"
    assert T.variant_label('other') == 'other'


def test_as_e136_and_as_e139_restore():
    order, names = W.ARM_ORDER, dict(W.ARM_NAMES)
    old = W139.variant_order, W139.present, W139.variant_label
    with pytest.raises(RuntimeError):
        with T.as_e136(['family', 'N_inhadd1000']), T.as_e139():
            assert W.arm_label('N_inhadd1000', {}) == T.variant_label('N_inhadd1000')
            assert W139.variant_label('G_inhadd1000') == T.variant_label('G_inhadd1000')
            raise RuntimeError
    assert W.ARM_ORDER == order and W.ARM_NAMES == names
    assert (W139.variant_order, W139.present, W139.variant_label) == old


def test_read_step():
    assert T.read_step(dict(tag='N_inhadd1000_megamat_real', panel='megamat')) == '1'
    assert T.read_step(dict(tag='D5.8_inhadd1000_apl16_M1', panel='apl16')) == '1'
    assert T.read_step(dict(tag='S_inhadd1000_validation2_P2', panel='validation2')) == '2'
    assert T.read_step(dict(tag='D1.8_inhadd1000_extended_M3', panel='extended')) == 'beside'


def _panel(arm, change, spread, worse=False):
    return dict(complete=True, arm=arm, change=change, spread=spread, worse_than_spread=worse)


def _classes():
    keys = T._keys()
    mg = dict(complete=True, change=-0.118, se=0.010, improves=True, arm=0.931)
    held = dict(cls='holds', megamat=mg, guard_worse=[], panel_worse=[], holds_failed=[], reasons=[],
                panels={k: _panel(1.0, 0.0, 0.05) for k in keys}, stands_in=None, does_more_than_a1=False,
                lower_price=False,
                against_a1=dict(megamat=dict(change=-0.0033, se=0.010),
                                panels={k: dict(change=0.01, better=False, worse=False) for k in keys},
                                price=dict(by_target={'z8': -0.002, 'halfA': -0.001, 'halfB': -0.002}, lower=False)))
    held['panels']['apl16/apl16/on'] = _panel(1.087, 0.05129, 0.05135)
    stopped = dict(cls='no class', reasons=['a hold fails: P_inhadd1000_megamat_real'],
                   holds_failed=['P_inhadd1000_megamat_real'], n_reads=3, stands_in=False, does_more_than_a1=None,
                   lower_price=None)
    return dict(classes={'D1.8_inhadd1000': held, 'P_inhadd1000': stopped}, arms_read=2,
                arms_holding=['D1.8_inhadd1000'], spreads={'stage': 0.020, **{k: 0.05 for k in keys}})


def test_classes_and_against_a1_tables():
    d = dict(classes=_classes())
    V = W.Values()
    lines = T.classes_table(d, V)
    row = [x for x in lines if x.startswith('| D1.8: KC–KC ×1.8')][0]
    assert '-0.118 (0.010)' in row and 'APL16, APL on: +0.00006' in row and row.endswith('| holds |')
    assert len(row.split(' | ')) == len(lines[2].split(' | '))
    row = [x for x in lines if x.startswith('| P: one uniform neuron')][0]
    assert '| 1 |' in row and row.endswith('| no class |')
    assert len(row.split(' | ')) == len(lines[2].split(' | '))
    lines = T.against_a1_table(d, V)
    row = [x for x in lines if x.startswith('| D1.8: KC–KC ×1.8')][0]
    assert '-0.0033 (0.010)' in row and '-0.0020 / -0.0010 / -0.0020' in row and row.endswith('| no |')
    f = V.frame().set_index(['arm', 'stage', 'panel', 'statistic']).value
    assert f[('D1.8_inhadd1000', 'kc_soma', 'apl16/apl16/on', 'margin_to_spread')] == pytest.approx(0.00006)
    assert f[('D1.8_inhadd1000', 'four_stages', 'megamat', 'change_from_a1')] == pytest.approx(-0.0033)
    assert f[('arms', 'all', 'megamat', 'n_arms_holding')] == 1


def test_dose_table():
    c = _classes()
    reads = pd.DataFrame([dict(variant='D1.8_inhadd1000', panel='apl16', wKCKC_scale=s) for s in (41.2, 82.9)])
    holds = pd.DataFrame([dict(variant='A1', input='apl16', condition='on', wKCKC_scale=s) for s in (74.2, 149.1)])
    d = dict(classes=c, reads=reads, reads_e139=pd.DataFrame(columns=['variant', 'panel', 'wKCKC_scale']),
             step0_holds=holds)
    V = W.Values()
    lines = T.dose_table(d, V)
    row = [x for x in lines if '| D1.8: KC–KC ×1.8 |' in x][0]
    assert row.startswith('| ×1.8 |') and '0.931 (-0.118)' in row and '41.2 to 82.9' in row and row.endswith('| holds |')
    row = [x for x in lines if x.startswith('| ×3.24 |')][0]
    assert '74.2 to 149.1' in row
    f = V.frame().set_index(['arm', 'statistic']).value
    assert f[('D5.8_inhadd1000', 'kckc_factor')] == pytest.approx(3.24 ** 1.5)


def test_spread_table_keys_the_input():
    rows = []
    for inp, pan, cond, _ in T.SOMA_PANELS:
        for w in F.T_WIRINGS:
            rows.append(dict(arm='family', input=inp, stage='kc_soma', panel=pan, condition=cond, wiring=w,
                             split_strength_rate_sd_ratio=1.5 if pan == 'megamat' else 1.0,
                             split_strength_rate_r=0.5))
    rows.append(dict(arm='family', input='megamat', stage='kc_soma', panel='megamat', condition='on', wiring='P3',
                     split_strength_rate_sd_ratio=9.0, split_strength_rate_r=0.0))
    V = W.Values()
    lines = T.spread_table(dict(panels=pd.DataFrame(rows)), V)
    row = [x for x in lines if x.startswith('| the family')][0]
    assert row.endswith('| 1.50 (0.50) |')
    f = V.frame().set_index(['arm', 'panel', 'condition', 'statistic']).value
    assert f[('family', 'apl16/megamat', 'on', 'sd_ratio_mean_T')] == pytest.approx(1.5)
    assert f[('family', 'apl16/apl16', 'off', 'sd_ratio_mean_T')] == pytest.approx(1.0)
