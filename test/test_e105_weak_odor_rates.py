#!/usr/bin/env python3
"""Tests for E105 (the KCs' response rates at the weak odors read like for like; experiments/E105_weak_odor_rates.py), on
made-up data: the three sets are never pooled; c* reproduces a known response rate; the pooled SD and the trial-noise
unit; the classes and W on known rows, with the 6-of-8 edge cases; the held rule and the outcomes; the blocked flies'
interval; part B's class; the lock refuses a changed definition; the reader refuses raw-only rows and any path outside
the shared directory. With the shared directory: the committed lock matches the module's definitions."""

import json

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
M = pytest.importorskip('experiments.E105_weak_odor_rates')


def _shared() -> bool:
    try:
        return (H.shared_dir() / 'data' / 'yang_weak_odor_kc').exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# the reader's refusals

def test_check_rows_refuses_raw_only_and_unlisted():
    df = pd.DataFrame(dict(recording=['a', 'b']))
    assert M.check_rows(df, ['a', 'b']) is df
    with pytest.raises(PermissionError):                     # a recording outside the processed list
        M.check_rows(df, ['a'])
    with pytest.raises(PermissionError):                     # raw only by the shared list
        M.check_rows(df, ['a', 'b'], raw_only={'b'})
    raw = df.assign(use_until_further_notice=["processed: usable", 'raw only: not used until further notice'])
    with pytest.raises(PermissionError):                     # raw only by its own mark
        M.check_rows(raw, ['a', 'b'])


def test_check_path_refuses_outside_shared(tmp_path):
    with pytest.raises(PermissionError):
        M.check_path(tmp_path / 'x.csv')
    with pytest.raises(PermissionError):
        M.read_trials('diag_mix', path=tmp_path / 'diag_mix_long.csv')
    with pytest.raises(PermissionError):
        M.raw_only_stems(tmp_path / 'recordings.csv')


# ---------------------------------------------------------------------------
# the sets are never pooled

def test_assign_set_dates():
    d = pd.Series(['2025-01-29', '2025-02-24', '2025-03-20', '2025-03-24'])
    assert list(M.assign_set(d, 'diag_mix')) == ['D1', 'D1', 'D2', 'D2']
    assert list(M.assign_set(pd.Series(['2026-01-16', '2026-02-10']), 'natmix_dilute')) == ['N', 'N']
    with pytest.raises(ValueError):                          # between the two diag_mix periods: in no set
        M.assign_set(pd.Series(['2025-03-01']), 'diag_mix')
    with pytest.raises(ValueError):
        M.assign_set(pd.Series(['2026-02-11']), 'natmix_dilute')


def _rates():
    """Made-up per-recording rates: D1 (two control flies, one blocked fly), D2 (two control flies), N (one fly)."""
    rows = []
    def add(s, rec, g, fly, rate, solvent):
        rows.append(dict(set=s, recording=rec, genotype_group=g, fly=fly, odor='ma @ -7', rate=rate, n_kcs=100))
        rows.append(dict(set=s, recording=rec, genotype_group=g, fly=fly, odor=M.SOLVENT, rate=solvent, n_kcs=100))
    add('D1', 'r1', 'control', 'f1', 0.10, 0.02)
    add('D1', 'r2', 'control', 'f1', 0.06, 0.02)             # f1's mean: 0.08, evoked 0.06
    add('D1', 'r3', 'control', 'f2', 0.04, 0.00)
    add('D1', 'r4', 'apl_output_blocked', 'f3', 0.20, 0.05)
    add('D2', 'r5', 'control', 'f4', 0.01, 0.00)
    add('D2', 'r6', 'control', 'f5', 0.03, 0.01)
    add('N', 'r7', 'control', 'f6', 0.50, 0.10)
    return M.evoked(pd.DataFrame(rows))


def test_sets_never_pooled():
    r = _rates()
    spec = dict(unit='z1', first_n=None, evoked=True, relative=False, spikes=2, drop_fly=None, pooled=False)
    c = M.set_cells(r, spec, ['ma @ -7']).set_index(['set', 'group'])
    assert sorted(set(c.index.get_level_values('set'))) == ['D1', 'D2', 'N']
    assert c.loc[('D1', 'control'), 'value'] == pytest.approx(np.mean([0.06, 0.04]))    # fly means of the evoked rates
    assert c.loc[('D2', 'control'), 'value'] == pytest.approx(np.mean([0.01, 0.02]))
    assert c.loc[('N', 'control'), 'value'] == pytest.approx(0.40)
    assert c.loc[('D1', 'control'), 'n_flies'] == 2 and c.loc[('D2', 'control'), 'n_flies'] == 2
    assert ('D2', 'apl_output_blocked') not in c.index
    # the carried reading alone pools D1 and D2, as recording means, the solvent kept
    p = M.set_cells(r, dict(spec, evoked=False, pooled=True), ['ma @ -7']).set_index(['set', 'group'])
    assert list(p.index.get_level_values('set').unique()) == [M.POOLED_SET]
    assert p.loc[(M.POOLED_SET, 'control'), 'value'] == pytest.approx(np.mean([0.10, 0.06, 0.04, 0.01, 0.03]))
    # a left-out fly leaves its set only
    lo = M.set_cells(r, dict(spec, drop_fly='f1'), ['ma @ -7']).set_index(['set', 'group'])
    assert lo.loc[('D1', 'control'), 'value'] == pytest.approx(0.04)
    assert lo.loc[('D2', 'control'), 'value'] == pytest.approx(c.loc[('D2', 'control'), 'value'])


def test_blocked_interval():
    lo, hi = M.blocked_interval(0.2, 0.05, 0.04, 0.0625)        # control interval: a factor 1.25 each way
    assert lo == pytest.approx(0.16) and hi == pytest.approx(0.25)
    assert all(np.isnan(M.blocked_interval(0.2, 0.05, 0.0, 0.06)))
    r = _rates()
    spec = dict(unit='z1', first_n=None, evoked=True, relative=False, spikes=2, drop_fly=None, pooled=False)
    c = M.set_cells(r, spec, ['ma @ -7']).set_index(['set', 'group'])
    ctrl, b = c.loc[('D1', 'control')], c.loc[('D1', 'apl_output_blocked')]
    h = 0.5 * (np.log(ctrl['hi90']) - np.log(ctrl['lo90']))
    assert b['value'] == pytest.approx(0.15)
    assert b['lo90'] == pytest.approx(0.15 * np.exp(-h)) and b['hi90'] == pytest.approx(0.15 * np.exp(h))


# ---------------------------------------------------------------------------
# the unit and c*

def _kc(values, rec='r1', odor='o1'):
    return pd.DataFrame(dict(set='APL16', recording=rec, genotype_group='control', fly='f1', roi=np.arange(len(values)),
                             odor=odor, tm=values, n=3, share_gt1=0.0, pooled_sd=1.0, tn=values))


def test_pooled_sd_known():
    t = pd.DataFrame(dict(set='D1', recording='r', genotype_group='control', fly='f', roi=0,
                          odor=['a', 'a', 'a', 'b', 'b', 'b'], repeat=[0, 1, 2, 0, 1, 2],
                          value=[1.0, 2.0, 3.0, 0.0, 0.0, 3.0]))
    # sums of squares 2 and 6 over 6 - 2 degrees of freedom
    assert float(M.pooled_sd(t).iloc[0]) == pytest.approx(np.sqrt(8.0 / 4.0))
    k = M.kc_table(t).set_index('odor')
    assert k.loc['a', 'tn'] == pytest.approx(2.0 / np.sqrt(2.0))
    assert k.loc['b', 'share_gt1'] == pytest.approx(1 / 3)
    k3 = M.kc_table(t, first_n=2)                            # the unit from the same trials
    assert float(k3['pooled_sd'].iloc[0]) == pytest.approx(np.sqrt((0.5 + 0.0) / 2.0))


def test_cstar_reproduces_a_known_rate():
    # two recordings, 100 KCs each over two odors; values 0.01 ... 1.00 in each
    v = np.linspace(0.01, 1.0, 100)
    kc = pd.concat([_kc(v, 'r1', 'o1'), _kc(v, 'r1', 'o2'), _kc(v, 'r2', 'o1'), _kc(v, 'r2', 'o2')], ignore_index=True)
    res = M.find_cstar(kc, 0.10, ['o1', 'o2'], bracket=(0.0, 2.0), tol=1e-9)
    assert res['rate_at_cstar'] == pytest.approx(0.10)       # 10 of 100 above c* in every recording and odor
    assert res['rate_below'] > 0.10
    assert 0.90 <= res['cstar'] < 0.91
    rr = M.panel_rate(M.response_rates(kc, 'trial_noise', res['cstar']), ['o1', 'o2'])
    assert rr['rate'].mean() == pytest.approx(0.10)
    with pytest.raises(RuntimeError):                        # a target the bracket cannot reach
        M.find_cstar(kc, 0.5, ['o1', 'o2'], bracket=(0.0, 0.2))


def test_response_units():
    kc = _kc(np.array([0.5, 1.5, 2.5]))
    kc['share_gt1'] = [0.0, 0.5, 1 / 3]
    assert list(M.responses(kc, 'z1')) == [False, True, True]
    assert list(M.responses(kc, 'yang')) == [False, True, False]
    assert list(M.responses(kc, 'trial_noise', 2.0)) == [False, False, True]
    with pytest.raises(ValueError):
        M.responses(kc, 'trial_noise')


# ---------------------------------------------------------------------------
# the classes, W, the held rule, the outcomes

def test_cell_class_edges():
    lo, hi = 0.04, 0.06
    assert M.cell_class([0.01] * 6 + [0.05] * 2, lo, hi) == 'short'          # 6 of 8 below
    assert M.cell_class([0.09] * 6 + [0.05] * 2, lo, hi) == 'over'
    assert M.cell_class([0.05] * 6 + [0.01] * 2, lo, hi) == 'inside'
    assert M.cell_class([lo] * 3 + [hi] * 3 + [0.0, 1.0], lo, hi) == 'inside'  # the ends are inside
    assert M.cell_class([0.01] * 5 + [0.05] * 3, lo, hi) == 'not decided'     # 5 below: median below
    assert M.cell_class([0.01] * 4 + [0.05] * 2 + [0.09] * 2, lo, hi) == 'not decided'   # median (0.03) below
    assert M.cell_class([0.01] * 3 + [0.05] * 2 + [0.09] * 3, lo, hi) == 'inside'        # neither: the median inside
    assert M.cell_class([0.01] * 5 + [0.09] * 3, lo, hi) == 'not decided'     # median between, on no side's 6
    assert M.cell_class([0.05] * 8, np.nan, hi) == 'not decided'


def test_w_class():
    assert M.w_class({'ma': 'short', '2h': 'inside', 'farn': 'inside'}) == 'short'
    assert M.w_class({'ma': 'not decided', '2h': 'short', 'farn': 'short'}) == 'short'
    assert M.w_class({'ma': 'inside', '2h': 'short', 'farn': 'short'}) == 'short'   # both short: short
    assert M.w_class({'ma': 'over', '2h': 'short', 'farn': 'inside'}) == 'not short'
    assert M.w_class({'ma': 'inside', '2h': 'short'}) == 'not short'                # set N: no farn
    assert M.w_class({'ma': 'not decided', '2h': 'short', 'farn': 'inside'}) == 'not decided'
    assert M.w_class({'ma': 'not decided', '2h': 'short'}) == 'not decided'


def test_held_w():
    h = M.held_w('short', {'unit': 'short', 'one_spike': 'short', 'fly_left_out': ['short', 'short']})
    assert h['held'] == 'short' and h['differ'] == [] and h['flips'] == []
    h = M.held_w('short', {'unit': 'short', 'fly_left_out': ['short', 'not decided']})
    assert h['held'] == 'not held' and h['differ'] == ['fly_left_out'] and h['flips'] == []
    h = M.held_w('not short', {'unit': 'short', 'first_three_trials': 'not short'})
    assert h['held'] == 'not held' and h['flips'] == ['unit']


def _held(rows):
    return pd.DataFrame([dict(estimate=e, condition=c, set=s, carrying=w, held=hd, differ=json.dumps(df),
                              flips=json.dumps(fl)) for e, c, s, w, hd, df, fl in rows])


def test_outcomes():
    rows = [(e, 'on', s, 'short', 'short', [], []) for e in M.ESTIMATES for s in M.SETS]
    rows += [(e, 'off', 'D1', 'short', 'short', [], []) for e in M.ESTIMATES]
    o = M.outcomes(_held(rows))
    assert o['a']['holds'] and not o['b']['holds'] and not o['c']['holds'] and not o['d']['holds']
    # W turns on the estimate under APL on in set D2, on the set under kiwi-control, a partner flips one cell
    rows2 = [list(r) for r in rows]
    for r in rows2:
        if r[0] == 'megamat' and r[1] == 'on' and r[2] == 'D2':
            r[3], r[4] = 'not short', 'not held'
        if r[0] == 'kiwi-control' and r[1] == 'off':
            r[4], r[5], r[6] = 'not held', ['one_spike', 'fly_left_out'], ['one_spike']
    o = M.outcomes(_held(rows2))
    assert not o['a']['holds']
    assert o['c']['where'] == ['on/D2'] and o['d']['where'] == ['on/megamat']
    assert o['e']['where'] == ['off/kiwi-control/D1: one_spike'] and not o['f']['holds']
    assert o['beside_fly_left_out_differs'] == ['off/kiwi-control/D1']
    # within D1, W differs between the conditions under validation2
    rows3 = [list(r) for r in rows]
    for r in rows3:
        if r[0] == 'validation2' and r[1] == 'on' and r[2] == 'D1':
            r[3] = 'not short'
    assert M.outcomes(_held(rows3))['f']['where'] == ['validation2']


def test_e102_class():
    assert M.e102_class(0.01, 0.02, 0.05) == 'further'      # ln 5 against ln 2.5: 1.76 times
    assert M.e102_class(0.04, 0.02, 0.05) == 'nearer'
    assert M.e102_class(0.019, 0.02, 0.05) == 'level'
    assert M.e102_class(0.0, 0.02, 0.05) == 'further'
    assert M.e102_class(0.0, 0.0, 0.05) == 'level'
    assert M.e102_class(0.02, 0.02, -0.01) == 'not read'
    # the ratio rule at its edges: 1.25 and 0.8 times the unchanged fits' log distance
    d0 = np.log(2.0)
    assert M.e102_class(0.05 * np.exp(-1.25 * d0), 0.025, 0.05) == 'further'
    assert M.e102_class(0.05 * np.exp(-0.8 * d0), 0.025, 0.05) == 'nearer'


# ---------------------------------------------------------------------------
# the lock

def _p0():
    return dict(passed=True, cstar=dict(cstar=1.93, target=0.068), carrying_unit='trial_noise', partner_unit='z1',
                unit_check=dict(rms=0.005, passed=True))


def test_lock_round_trip_and_refusals(tmp_path, monkeypatch):
    p = tmp_path / 'lock.json'
    M.write_lock(_p0(), p)
    assert M.locked(p, check_git=False)
    assert M.require_lock(p, check_git=False)['cstar'] == pytest.approx(1.93)
    assert not M.locked(p)                                   # not committed
    with pytest.raises(PermissionError):
        M.require_lock(p)
    with pytest.raises(RuntimeError):                        # other values: not overwritten
        M.write_lock(dict(_p0(), cstar=dict(cstar=2.5, target=0.068)), p)
    with pytest.raises(RuntimeError):                        # part 0 failed: no lock
        M.write_lock(dict(_p0(), passed=False), tmp_path / 'other.json')
    monkeypatch.setattr(M, 'MIN_WIRINGS', 5)                 # a changed definition
    assert not M.locked(p, check_git=False)
    with pytest.raises(PermissionError):
        M.require_lock(p, check_git=False)
    monkeypatch.undo()
    d = json.loads(p.read_text())
    d['values']['cstar'] = 1.5                               # a tampered value
    p.write_text(json.dumps(d))
    assert not M.locked(p, check_git=False)


def test_committed_lock_matches_definitions():
    if not M.LOCK.exists():
        pytest.skip('the lock is written by part 0')
    d = json.loads(M.LOCK.read_text())
    assert d['sha256'] == M.definitions_sha(M.definitions())
    assert M.locked(M.LOCK, check_git=False)


# ---------------------------------------------------------------------------
# with the shared directory

@pytest.mark.skipif(not _shared(), reason='needs the shared directory')
def test_sets_counts_from_the_data():
    t = pd.concat([M.read_trials('diag_mix'), M.read_trials('natmix_dilute')], ignore_index=True)
    c = M.set_counts(t)
    for k, v in M.SET_COUNTS.items():
        if k[0] != 'APL16':
            assert c[k] == v
    # no recording is in two sets
    assert (t.groupby('recording')['set'].nunique() == 1).all()


# ---------------------------------------------------------------------------
# part D's helpers

def test_stim_input_rules():
    k = M.model_label('2h', 'megamat')
    assert M.stim_input('2h @ -8', 'megamat', 'as_blank') is None
    assert M.stim_input('2h @ -8', 'megamat', 'as_above') == [k]
    assert M.stim_input('ma @ -8', 'megamat', 'as_blank') == [M.model_label('ma', 'megamat')]
    assert M.stim_input('2h @ -7 + farn @ -2.5', 'kiwi-control', 'as_blank') == [
        M.model_label('2h', 'kiwi-control'), M.model_label('farn', 'kiwi-control')]
    assert M.stim_input('2h @ -5 + 1o3ol @ -5', 'megamat', 'as_blank') == [M.BESIDE_STIM[0][2]]
    for s in (M.SOLVENT, 'ban', 'Bmix28', 'air', 'eb @ -7', '1o3ol @ -5'):
        assert M.stim_input(s, 'megamat', 'as_above') is None


def test_log_inverse():
    x, y = np.array([0.5, 1.0, 2.0]), np.array([0.01, 0.02, 0.04])
    assert M._log_inverse(x, y, 0.02) == (pytest.approx(1.0), '')
    f, c = M._log_inverse(x, y, 0.03)
    assert f == pytest.approx(np.sqrt(2.0)) and c == ''
    assert M._log_inverse(x, y, 0.001)[1] == 'below 0.5'
    assert M._log_inverse(x, y, 0.1)[1] == 'above 2'


def test_session_columns_weight_by_trials():
    lab = M.model_label('2h', 'kiwi-control')
    counts = pd.DataFrame({lab: [0.0, 3.0], 'other': [1.0, 1.0]})
    c, recs, g, n_blank = M._session_columns([['2h @ -7', 6], [M.SOLVENT, 3]], 'kiwi-control', 'as_blank', counts)
    assert g == 3 and n_blank == 1 and c.shape[1] == 2            # 2h's six trials: two units of three
    assert [r['column'] for r in recs] == [lab, None]
