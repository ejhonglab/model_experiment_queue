"""E115 (the megamat gap against two sets of flies and the yardstick of the price in response rate ρ), on made-up
matrices: the residual correlation and the two sets' mean rank, the classes of pairs and their shares of the gap, the
carriers of a class, the size of an input change, the percentile and pooled readings, and the price and outcome
verdicts."""
import numpy as np
import pandas as pd
import pytest

from kcq import harness as H
from experiments import E115_two_yardsticks as E

ODORS = list(H.MEGAMAT12)


def pairs_series(values) -> pd.Series:
    idx = pd.MultiIndex.from_tuples(H.pair_list(ODORS), names=['odor_a', 'odor_b'])
    return pd.Series(np.asarray(values, dtype=float), index=idx)


def made_up(seed=0):
    rng = np.random.default_rng(seed)
    base = rng.random(66)
    return pairs_series(base), rng


def test_interface():
    assert E.SHARED_MIN == 0.30 and E.LEVEL_TOL == 0.03 and E.POOLED_MAX_FALL == 0.02
    assert (E.CLASS_NEAR, E.CLASS_FAR) == (10, 15)
    for f in ('residual_correlation', 'two_set_rank', 'pair_classes', 'class_shares', 'reader_row', 'resample_part_a',
              'input_size', 'fractions', 'rate_rho', 'pooled_rate_rho', 'percentile_of', 'price_verdict', 'verdicts',
              'part0', 'part_a', 'part_b', 'registry_entries'):
        assert callable(getattr(E, f))
    srcs = {c['source'] for c in E.CONTROLS}
    assert srcs == {'E93', 'E102', 'E97', 'E108', 'E109'}
    assert sum(c['source'] == 'E109' and c['kind'] != 'control search (2h, 6al)' for c in E.CONTROLS) == 30


def test_residual_correlation_and_two_set_rank():
    y, rng = made_up(1)
    m = y + rng.normal(0, 0.3, 66)
    # Remy's flies as the model: the residuals are the same, correlation 1
    assert E.residual_correlation(m, y, m) == pytest.approx(1.0)
    # Remy's flies as Yang's: Remy's residuals are all 0 (no correlation defined)
    assert np.isnan(E.residual_correlation(m, y, y))
    # the mean rank of a matrix with itself is its rank; with its reverse, constant
    t = E.two_set_rank(y, y)
    assert np.allclose(t.to_numpy(), pd.Series(y).rank().to_numpy())
    t2 = E.two_set_rank(y, -y)
    assert np.allclose(t2.to_numpy(), 33.5)
    row = E.reader_row(m, y, y)
    assert row['remy_vs_yang'] == pytest.approx(1.0) and row['rho_two_sets'] == pytest.approx(row['rho_yang'])
    assert row['n_pairs'] == 66


def test_pair_classes_and_shares():
    r = np.arange(66, dtype=float)
    y = pairs_series(r)
    remy = y.copy()
    m = y.copy()
    # pair 0: both sets at rank 1, the model at rank 66 (swap with the last pair)
    m.iloc[0], m.iloc[65] = 100.0, -1.0
    # pair 30: Yang's far from the model and from Remy's, Remy's with the model
    yy = y.copy()
    yy.iloc[30] = 200.0
    cl = E.pair_classes(m, yy, remy)
    assert cl['cls'].iloc[0] == E.CLASSES[0]
    assert cl['cls'].iloc[30] == E.CLASSES[1]
    assert set(cl['cls'].iloc[1:30]) == {E.CLASSES[2]}
    sh = E.class_shares(m, yy, cl['cls'])
    tot = sh[sh.level == 'class']['share_of_gap'].sum()
    assert tot == pytest.approx(1.0, abs=1e-9)
    # the odor rows of a class sum to twice the class's contribution (each pair counts for both odors)
    for c in E.CLASSES:
        s = sh[sh.cls == c]
        assert s[s.level == 'odor']['contribution'].sum() == pytest.approx(2 * s[s.level == 'class']['contribution'].sum())


def test_carriers():
    idx = pd.MultiIndex.from_tuples([('2h', 'eb'), ('2h', 'pa'), ('2h', 'ms'), ('eb', 'ms')], names=['odor_a', 'odor_b'])
    cl = pd.DataFrame({'cls': [E.CLASSES[1]] * 4}, index=idx)
    assert E.carriers(cl) == ['2h', 'eb', 'ms']   # 3 of 4, 2 of 4, 2 of 4
    assert E.carriers(cl.iloc[:2]) == []          # below CARRY_MIN_PAIRS


def test_input_size():
    base = pd.DataFrame(np.ones((3, 2)), index=['A', 'B', 'C'], columns=['x', 'y'])
    assert E.input_size(base * 1.1, base) == pytest.approx(0.1)
    more = pd.concat([base, pd.DataFrame([[3.0, 4.0]], index=['D'], columns=['x', 'y'])])
    assert E.input_size(more, base) == pytest.approx(5 / np.sqrt(6))


def test_percentile_and_pooled(monkeypatch):
    assert E.percentile_of(0.5, [0, 1]) == 50
    assert E.percentile_of(1, [0, 1, 2]) == 50
    data = pd.Series(np.arange(16, dtype=float), index=H.APL16_ODORS)
    monkeypatch.setattr(E.AE, 'frac_data', lambda cond, units='z': data)
    a = pd.Series(np.arange(16, dtype=float), index=H.APL16_ODORS)
    b = pd.Series(np.arange(16, dtype=float)[::-1] * 0.5, index=H.APL16_ODORS)
    assert E.rate_rho(a, 'on') == pytest.approx(1.0)
    assert E.pooled_rate_rho([a, b], 'on') == pytest.approx(E.rate_rho((a + b) / 2, 'on'))


def _controls(n=40, seed=0, size=0.1):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(dict(size=np.full(n, size), d_on_T=rng.normal(-0.03, 0.02, n)))


def test_price_verdict():
    ctl = _controls()
    lo, hi = np.percentile(ctl.d_on_T, [5, 95])
    mid = float(np.median(ctl.d_on_T))
    dirs = {'A': dict(size=0.1, d_rate_rho_on=mid, pooled_d_rate_rho_on=0.0),
            'B': dict(size=0.1, d_rate_rho_on=mid, pooled_d_rate_rho_on=-0.01)}
    v = E.price_verdict(dirs, ctl)
    assert v['inside'] and v['directions']['A']['p5'] == pytest.approx(lo) and v['directions']['B']['p95'] == pytest.approx(hi)
    dirs['B']['d_rate_rho_on'] = lo - 0.01
    assert not E.price_verdict(dirs, ctl)['inside']
    dirs['B']['d_rate_rho_on'] = mid
    dirs['B']['pooled_d_rate_rho_on'] = -0.03
    assert not E.price_verdict(dirs, ctl)['inside']
    # no comparable control: every sized control is used, and the fallback is said
    dirs['B']['pooled_d_rate_rho_on'] = 0.0
    dirs['A']['size'] = 1.0
    v = E.price_verdict(dirs, ctl)
    assert v['directions']['A']['fallback_all_sized'] and v['directions']['A']['n_used'] == 40


def test_verdicts_letters():
    def a_table(rc_on, lo_on, rc_off, lo_off, two_on, flies_on, two_off, flies_off):
        rows = []
        for c, rc, lo, two, fl in (('on', rc_on, lo_on, two_on, flies_on), ('off', rc_off, lo_off, two_off, flies_off)):
            rows.append(dict(reader='family', remy='corr', condition=c, resid_corr=rc, resid_corr_lo=lo, resid_corr_hi=rc + 0.1,
                             rho_two_sets=two, remy_vs_yang=fl))
        return pd.DataFrame(rows)
    idx = [('2h', 'eb'), ('2h', 'pa'), ('2h', 'ms')]
    cls = pd.DataFrame(dict(odor_a=[p[0] for p in idx] * 2, odor_b=[p[1] for p in idx] * 2, cls=[E.CLASSES[1]] * 6,
                            remy='corr', condition=['on'] * 3 + ['off'] * 3))
    price = dict(inside=True, directions={})
    v = E.verdicts(a_table(0.4, 0.1, 0.35, 0.05, 0.83, 0.85, 0.84, 0.86), cls, price)
    assert v['outcomes'] == ['a', 'c', 'e']
    v = E.verdicts(a_table(0.4, -0.01, 0.35, 0.05, 0.83, 0.85, 0.84, 0.86), cls.iloc[:0], dict(inside=False, directions={}))
    assert v['outcomes'] == ['b', 'd']
    v = E.verdicts(a_table(0.4, 0.1, 0.35, 0.05, 0.80, 0.85, 0.84, 0.86), cls.iloc[:0], price)
    assert v['outcomes'] == ['c'] and v['outcome_note']
