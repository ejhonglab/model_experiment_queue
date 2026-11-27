"""Tests of E147's input flies, points and score helpers on planted values (no simulation)."""
import numpy as np
import pandas as pd
import pytest

from experiments import E147_megamat_input_flies as M
from kcq import harness as H


def planted_long(seed=0, missing=(('2023-04-22', 2, 'aa'),), shift5=0.0):
    """A long table like olfdb dataset 1: nine flies, megamat's 17 odors at -3 on the lab input's 36 glomeruli, three
    trials each (each fly a gain times the lab's megamat input plus trial noise), and one diagnostics row per fly that
    sets its largest |trial mean|."""
    rng = np.random.default_rng(seed)
    meg = H.abbrev_columns(M.lab_input('megamat'))
    keys = [('2023-04-22', 2), ('2023-04-22', 3), ('2023-04-26', 2), ('2023-04-26', 3), ('2023-05-08', 1),
            ('2023-05-08', 3), ('2023-05-09', 1), ('2023-05-10', 1), ('2023-06-22', 1)]
    rows = []
    for i, (d, f) in enumerate(keys):
        gain = 1.0 + 0.1 * i
        for o in meg.columns:
            if (d, f, o) in missing:
                continue
            for g in meg.index:
                base = gain * float(meg.loc[g, o]) / 100.0 + (shift5 if o not in H.MEGAMAT12 else 0.0)
                for t in range(3):
                    rows.append((d, f, str(g), 'megamat', f'{o} @ -3', t, base + 0.01 * rng.standard_normal()))
        rows.append((d, f, 'D', 'diagnostics', 'x @ -5', 0, 10.0))
    return pd.DataFrame(rows, columns=['date', 'fly_num', 'roi', 'panel', 'odor1', 'repeat', 'value'])


@pytest.fixture(scope='module')
def planted():
    return M.orn_flies(planted_long())


def test_mean_over_flies_is_the_lab_input(planted):
    for panel in M.PANELS:
        lab = M.lab_input(panel).to_numpy(dtype=float)
        ins = [M.input_fly(f, panel, F=planted).to_numpy(dtype=float) for f in planted['labels']]
        assert np.allclose(np.mean(ins, axis=0), lab, atol=1e-9)


def test_unmeasured_and_kept_entries_are_the_lab_input(planted):
    meg = M.input_fly('f0', 'megamat', F=planted)
    lab = M.lab_input('megamat')
    j = [H.odor_abbrev(c) for c in lab.columns].index('aa')
    assert (meg.iloc[:, j].to_numpy() == lab.iloc[:, j].to_numpy()).all()
    a16 = M.input_fly('f3', 'apl16', F=planted)
    lab16 = M.lab_input('apl16')
    for k, c in enumerate(lab16.columns):
        same = (a16.iloc[:, k].to_numpy() == lab16.iloc[:, k].to_numpy()).all()
        assert same == (H.odor_abbrev(c) in M.APL16_KEPT_AT_LAB)
    assert a16.index.equals(lab16.index) and a16.columns.equals(lab16.columns)


def test_deviation_is_reliability_times_trial_mean_minus_mean(planted):
    from experiments import E85_orn_fly_population as E85
    x = planted['x']
    m = np.nanmean(x, axis=3)
    j = planted['odors'].index('eb')
    k = 5
    M_ = np.nanmean(m[:, j, k])
    r = M.reliabilities(planted)[k]
    dev = M.megamat_deviations(F=planted)
    assert np.allclose(dev[:, j, k], r * (m[:, j, k] - M_))
    i12 = [planted['odors'].index(o) for o in H.MEGAMAT12]
    assert np.allclose(M.reliabilities(planted),
                       np.clip(E85.trial_stats(x[:, i12])['r_g'] * M.e85_part_a()['near_check']['c'], 0, 1))


def test_shared_odors_do_not_depend_on_the_other_five():
    a = M.orn_flies(planted_long(shift5=0.0))
    b = M.orn_flies(planted_long(shift5=0.3))
    da, db = M.megamat_deviations(F=a), M.megamat_deviations(F=b)
    i12 = [a['odors'].index(o) for o in H.MEGAMAT12]
    assert np.allclose(da[:, i12], db[:, i12])
    assert a['slope'] == pytest.approx(b['slope'])
    j = a['odors'].index('aa')
    assert not np.allclose(da[:, j], 0.0)


def test_same_deviation_on_both_panels(planted):
    meg = H.abbrev_columns(M.input_fly('f4', 'megamat', F=planted) - M.lab_input('megamat'))
    a16 = H.abbrev_columns(M.input_fly('f4', 'apl16', F=planted) - M.lab_input('apl16'))
    for o in H.MEGAMAT12:
        assert np.allclose(meg[o].to_numpy(), a16[o].to_numpy())


def test_zero_deviation_input_is_the_lab_input():
    for panel in M.PANELS:
        z = M.zero_deviation_input(panel)
        assert z.equals(M.lab_input(panel).astype(float))


def test_tags_and_points():
    assert M.parse_tag(M.point_tag('A1', 'apl16', 'f8', 'M3')) == ('A1', 'apl16', 'f8', 'M3')
    with pytest.raises(ValueError):
        M.point_tag('A1', 'validation2', 'f0', 'P3')
    p1, p2 = M.planned_points(1), M.planned_points(2)
    assert len(p1) == 72 and len(p2) == 108
    assert [p.wiring for p in p1[:36]] == [w for w in M.T_WIRINGS for _ in range(9)]
    assert all(p.wiring in M.T_WIRINGS for p in p2)
    assert len({M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in p1 + p2}) == 180
    assert all(p.setting.endswith('free') for p in M.planned_points(2, free=True))
    assert M.step_of('familyfree', 'megamat') == 1 and M.step_of('family', 'apl16') == 2
