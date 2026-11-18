#!/usr/bin/env python3
"""Tests for E108 (paa's input at the flies' concentration; experiments/E108_paa_extrapolation.py): the rules X1, X2 and
X3 and the score that chooses between them on made-up triples (X3 left one odor out), the controls' construction, the
unit conversion and the gate and class on made-up readings; with the model's input files, the inputs (only paa's column
differs from the lab's, the controls as defined); and part 0's identity on the real wiring (slow, skipped with
QUICK=1)."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E108_paa_extrapolation')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# the rules

def test_rule_x1_linear_and_capped():
    lo, hi = np.array([0.0, 10.0, 50.0, 30.0]), np.array([10.0, 30.0, 90.0, 10.0])
    np.testing.assert_allclose(X.rule_x1(lo, hi, 100.0), [20.0, 50.0, 100.0, -10.0])   # capped; no floor


def test_rule_x2_saturating_toward_the_ceiling():
    C = 100.0
    lo, hi = np.array([0.0, 50.0, 20.0]), np.array([50.0, 75.0, 20.0])
    # f = 0.5 -> 50 + 0.5 * 50 = 75; f = 0.5 -> 75 + 0.5 * 25 = 87.5; f = 0 -> unchanged
    np.testing.assert_allclose(X.rule_x2(lo, hi, C), [75.0, 87.5, 20.0])
    # never above the ceiling when rising from below it, and below X1 (saturating)
    rng = np.random.default_rng(0)
    lo = rng.uniform(0, 60, 50)
    hi = lo + rng.uniform(0, 30, 50)
    x2 = X.rule_x2(lo, hi, C)
    assert (x2 <= C + 1e-12).all() and (x2 <= X.rule_x1(lo, hi, np.inf) + 1e-12).all() and (x2 >= hi - 1e-12).all()
    # the distance to the ceiling falls geometrically: (C - x2) = (C - hi)^2 / (C - lo)
    np.testing.assert_allclose(C - x2, (C - hi) ** 2 / (C - lo))
    # a falling glomerulus keeps falling by the same rule; no room left: unchanged
    np.testing.assert_allclose(X.rule_x2(np.array([60.0]), np.array([40.0]), C), [40.0 - 20.0 * 60.0 / 40.0])
    np.testing.assert_allclose(X.rule_x2(np.array([120.0]), np.array([110.0]), C), [110.0])


def test_fit_x3_recovers_its_coefficients():
    rng = np.random.default_rng(1)
    t = pd.DataFrame(dict(low=rng.uniform(-5, 50, 200), mid=rng.uniform(-5, 80, 200)))
    t['top'] = 3.0 + 1.4 * t['mid'] - 0.2 * t['low']
    coef = X.fit_x3(t)
    np.testing.assert_allclose(coef, [3.0, 1.4, -0.2], atol=1e-9)
    np.testing.assert_allclose(X.rule_x3(np.array([10.0]), np.array([20.0]), coef), [3.0 + 28.0 - 2.0])
    with pytest.raises(ValueError):
        X.fit_x3(t.iloc[:2])


def _triples(top_fn, odors=('o1', 'o2', 'o3'), n_glom=12, seed=2):
    rng = np.random.default_rng(seed)
    rows = []
    for o in odors:
        for g in range(n_glom):
            low = rng.uniform(0, 40)
            mid = low + rng.uniform(0, 30)
            rows.append(dict(odor=o, glomerulus=f'g{g}', conc_low=-5.0, low=low, mid=mid, top=top_fn(low, mid)))
    return pd.DataFrame(rows)


def test_rule_scores_linear_ramps_favour_x1():
    t = _triples(lambda lo, mi: 2 * mi - lo)
    sc = X.rule_scores(t, cap=200.0)                             # above every top (at most 2 x 70), so X1 is not capped
    assert sc.loc['X1', 'rmse'] < 1e-9
    assert sc.loc['X3', 'rmse'] < 1e-6                       # X3 holds X1 (a=0, b=2, c=-1) and is fitted without the odor
    assert sc.loc['X2', 'rmse'] > 1e-3
    assert X.carrying_rule(sc) in ('X1', 'X3')
    assert int(sc.loc['X1', 'n_triples']) == len(t) and int(sc.loc['X1', 'n_odors']) == 3


def test_rule_scores_saturating_ramps_favour_x2():
    C = 100.0
    t = _triples(lambda lo, mi: mi + (mi - lo) / (C - lo) * (C - mi))
    sc = X.rule_scores(t, cap=C)
    assert sc.loc['X2', 'rmse'] < 1e-9
    assert X.carrying_rule(sc) == 'X2'


def test_x3_is_scored_with_its_odor_left_out():
    # one odor follows another law: X3 left that odor out cannot have learned it, so its error there is not zero
    t = pd.concat([_triples(lambda lo, mi: 1.5 * mi, odors=('o1', 'o2', 'o3')),
                   _triples(lambda lo, mi: 0.5 * mi, odors=('odd',), seed=3)], ignore_index=True)
    p = X.predictions(t, cap=1e6)
    odd = p[p['odor'] == 'odd']
    assert np.sqrt(np.mean((odd['pred_X3'] - odd['top']) ** 2)) > 1.0
    coef_without = X.fit_x3(t[t['odor'] != 'odd'])
    np.testing.assert_allclose(odd['pred_X3'], X.rule_x3(odd['low'], odd['mid'], coef_without))
    by = X.scores_by_odor(t, cap=1e6).set_index('odor')
    assert by.loc['odd', 'X3'] > by.loc['o1', 'X3']


def test_carrying_rule_is_the_lowest():
    sc = pd.DataFrame(dict(rmse=[3.0, 2.0, 2.5, 9.0]), index=['X1', 'X2', 'X3', 'none (top = mid)'])
    assert X.carrying_rule(sc) == 'X2'


# ---------------------------------------------------------------------------
# the controls

def test_derangement_moves_everything():
    for s in X.PERM_SEEDS:
        p = X.derangement(36, s)
        assert sorted(p) == list(range(36)) and not (p == np.arange(36)).any()
    assert not np.array_equal(X.derangement(36, 1), X.derangement(36, 2))


def test_paa_columns():
    g = pd.Index([f'g{i}' for i in range(10)], name='glomerulus')
    rng = np.random.default_rng(4)
    lo = pd.Series(rng.uniform(0, 40, 10), index=g)
    hi = pd.Series(lo.to_numpy() + rng.uniform(-10, 40, 10), index=g)
    C = 60.0
    t = X.paa_columns(lo, hi, C, np.array([0.0, 2.0, -1.0]))
    np.testing.assert_allclose(t['X1'], np.minimum(2 * hi - lo, C))
    np.testing.assert_allclose(t['X3'], 2 * hi - lo)            # X3 with X1's coefficients, uncapped
    np.testing.assert_allclose(t['W5'], lo)
    d = t['X1'].to_numpy() - hi.to_numpy()
    for s, name in zip(X.PERM_SEEDS, X.PERMS):
        np.testing.assert_allclose(t[name], np.minimum(hi.to_numpy() + d[X.derangement(10, s)], C))
    assert np.isclose(t['S1'].sum(), t['X1'].sum())
    np.testing.assert_allclose(t['S1'] / hi, t['S1'].iloc[0] / hi.iloc[0])   # the pattern unchanged


def test_to_spike_deltas():
    fm = pd.DataFrame(dict(date=['d1', 'd1', 'd2'], fly_num=[1, 1, 2], panel='p', odor1='o @ -3', roi=['a', 'b', 'a'],
                           value=[1.0, -0.5, 2.0]))
    mx = X.fly_maxima(fm)
    assert dict(zip(mx['date'], mx['fly_max'])) == {'d1': 1.0, 'd2': 2.0}
    s = X.to_spike_deltas(fm, mx, dict(avg_flymax=3.0, slope=100.0))
    np.testing.assert_allclose(s['rate'], [300.0, -150.0, 300.0])


# ---------------------------------------------------------------------------
# the gate and the class on made-up readings

def _made_up(d_on=0.02, d_off=0.0, d_pattern=0.0, flagged=False, w5=-0.03, perm=0.0, d_block=0.0, d_rate=0.0):
    """Readings of I0 and of 'X1' (and the controls) on the eight wirings and their population: X1 changes ρ by
    `d_on` / `d_off` everywhere and the pattern part by `d_pattern`."""
    cols = {}
    base = dict(all_pattern=0.048, all_loss=0.073, mean_effect=0.178, block=0.03, vcf_on=0.2, vcf_off=0.3,
                struct_wt_on=0.5, struct_wt_off=0.5)
    for t in X.TARGETS:
        base[f'rho_on_{t}'] = 0.69
        base[f'rho_off_{t}'] = 0.76
        base[f'rho15_on_{t}'] = 0.5
        base[f'rho15_off_{t}'] = 0.5
        base[f'{t}_pattern'] = 0.048
        base[f'{t}_loss'] = 0.073
        for k in ('amplitude', 'shape', 'sd_ratio', 'rho_effect', 'pearson_effect'):
            base[f'{t}_{k}'] = 0.1
    rows = []
    change = {'X1': (d_on, d_off, d_pattern, d_block), 'W5': (w5, w5, 0.0, 0.0), 'I0': (0, 0, 0, 0)}
    change.update({p: (perm, perm, 0.0, 0.0) for p in X.PERMS})
    for name, (a, b, c, k) in change.items():
        for unit in X.ALL_WIRINGS + ('all8',):
            for ro in X.READOUTS:
                r = dict(base, input=name, unit=unit, readout=ro)
                for t in X.TARGETS:
                    r[f'rho_on_{t}'] += a
                    r[f'rho_off_{t}'] += b
                    r[f'{t}_pattern'] += c
                r['block'] += k
                rows.append(r)
    rd = pd.DataFrame(rows)
    runs = pd.DataFrame([dict(input=n, unit=w, flagged=(flagged and n == 'X1'), rate_rho_on=0.6 + (d_rate if n == 'X1' else 0),
                              rate_rho_off=0.6, frac_paa_on=0.05, frac_paa_off=0.1)
                         for n in change for w in X.ALL_WIRINGS])
    return rd, runs


def test_gate_better():
    rd, runs = _made_up()
    g = X.gate(rd, runs, 'X1')
    assert g['gain'] and g['gain_on'] and not g['gain_off']
    assert g['n_wirings_gain_on'] == 8 and g['apl_effect_ok'] and g['realism_ok'] and g['controls_ok']
    assert X.class_of(g) == 'better' and X.outcome_of(g) == 'a'


def test_gate_gain_threshold_and_other_slack():
    rd, runs = _made_up(d_on=0.009)
    assert X.class_of(X.gate(rd, runs, 'X1')) == 'no gain'
    rd, runs = _made_up(d_on=0.02, d_off=-0.006)
    g = X.gate(rd, runs, 'X1')
    assert not g['gain_on'] and X.class_of(g) == 'no gain'
    assert X.outcome_of(g) == 'c'                                # W5 loses 0.03
    rd, runs = _made_up(d_on=0.0, w5=0.0)
    assert X.outcome_of(X.gate(rd, runs, 'X1')) == 'd'


def test_gate_price_names_the_constraint():
    rd, runs = _made_up(d_pattern=0.004)
    g = X.gate(rd, runs, 'X1')
    assert g['gain'] and not g['apl_effect_ok'] and X.class_of(g) == 'a price' and X.outcome_of(g) == 'b'
    assert any(k.startswith('pattern part') for k in g['fails'])
    rd, runs = _made_up(flagged=True)
    assert X.class_of(X.gate(rd, runs, 'X1')) == 'a price'
    rd, runs = _made_up(d_rate=-0.06)
    assert 'response rate ρ, APL on' in X.gate(rd, runs, 'X1')['fails']


def test_gate_pending_items_and_controls():
    rd, runs = _made_up(d_block=-0.011)
    g = X.gate(rd, runs, 'X1')
    assert 'the block' in g['pending_fails'] and X.class_of(g) == 'a price' and X.class_of(g, with_pending=False) == 'better'
    rd, runs = _made_up(perm=0.011)                              # the permuted inputs gain more than half
    g = X.gate(rd, runs, 'X1')
    assert g['controls_ok'] is False and X.class_of(g) == 'a price'


def test_wirings_needed():
    assert X.wirings_needed(8) == 6 and X.wirings_needed(4) == 3


# ---------------------------------------------------------------------------
# with the model's input files

def test_inputs_change_only_paa():
    lab = X.lab_input()
    p = X.paa_rates()
    assert list(p.index) == list(lab.index) and np.array_equal(p['hi'].to_numpy(), lab[X.PAA_HI_COL].to_numpy())
    coef = np.array([0.0, 2.0, -1.0])
    t = X.paa_columns(p['lo'], p['hi'], X.cap(), coef)
    X._CACHE['paa_table'] = t
    try:
        for name in X.RULES + X.CONTROLS:
            orn = X.input_of(name)
            other = [c for c in lab.columns if c != X.PAA_HI_COL]
            assert orn[other].equals(lab[other])
            np.testing.assert_allclose(orn[X.PAA_HI_COL], t[name])
        assert X.input_of('I0').equals(lab)
        assert X.cap() == float(lab.to_numpy().max())
    finally:
        X._CACHE.pop('paa_table', None)
        X._CACHE.pop('inputs', None)


@pytest.mark.skipif(QUICK or not _shared(), reason='slow (simulations); QUICK=1')
def test_identity_real_wiring(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(X, '_commits_clean', lambda: {})
    res = X.identity(units=('real',))
    assert res['units']['real']['count_diffs'] == 0 and res['units']['real']['breadth_equal']


def test_commits_clean_allows_only_unloaded_dirty_files(monkeypatch):
    commits = dict(queue_commit='abc', queue_dirty=True, pyolfsysm_dirty=False, al_analysis_dirty=False, hong2p_dirty=False)
    monkeypatch.setattr(X.H, 'code_commits', lambda: dict(commits))
    # the dirty file stands for one this process does not load, whatever other tests of the suite imported before
    unloaded = 'experiments/E97_glomerulus_claw_shares.py'
    loaded = X._loaded_queue_files
    monkeypatch.setattr(X, '_loaded_queue_files', lambda: loaded() - {unloaded})
    monkeypatch.setattr(X, '_queue_dirty_paths', lambda: [unloaded])
    cc = X._commits_clean()
    assert cc['queue_dirty_paths'] == 'experiments/E97_glomerulus_claw_shares.py'
    for bad in ('experiments/E108_paa_extrapolation.py', 'kcq/harness.py', 'experiments/E83_cross_wiring_fit.py'):
        monkeypatch.setattr(X, '_queue_dirty_paths', lambda bad=bad: [bad])
        with pytest.raises(RuntimeError):
            X._commits_clean()
    monkeypatch.setattr(X.H, 'code_commits', lambda: dict(commits, pyolfsysm_dirty=True))
    monkeypatch.setattr(X, '_queue_dirty_paths', lambda: [])
    with pytest.raises(RuntimeError):
        X._commits_clean()


def test_recorded_input_diff(tmp_path):
    """An input rebuilt on another host is compared with the recorded file: a last-bit difference is small, a changed
    value large, and other glomeruli or odors are infinitely far."""
    odors = ['pa @ -3', 'paa @ -3']
    built = pd.DataFrame([[10.0, 1.0 / 3.0], [234.41349774532097, -5.5]], index=['D', 'DA2'], columns=odors)
    lay = H.abbrev_columns(built)
    lay.index.name = 'glomerulus'
    p = tmp_path / 'X3.csv'
    p.write_text(lay.to_csv())
    assert X.recorded_input_diff(built, p) < 1e-12
    nudged = built.copy()
    nudged.iloc[1, 0] += 1e-13
    assert X.recorded_input_diff(nudged, p) <= X.INPUT_ATOL
    changed = built.copy()
    changed.iloc[0, 1] += 0.01
    assert X.recorded_input_diff(changed, p) > X.INPUT_ATOL
    assert X.recorded_input_diff(built.iloc[::-1], p) == float('inf')
