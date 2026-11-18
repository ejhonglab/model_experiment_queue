#!/usr/bin/env python3
"""Tests for E110 (KC-KC inhibition held at zero in the best family; experiments/E110_no_kckc.py) on made-up tables: the
verdict's classes ("not needed", "needed", "better without", "not determined") and the 6-of-8 rule, part C's reading,
part E's alternatives, the starts file (read back by E10_joint_fit.read_params), the counted start, the population and
paired readings, part 0's comparisons, the bands; and the fits' tags and commands against E72's. No model simulation."""

import json

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E110_no_kckc')

W8 = X.WIRINGS


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# the verdict

def _ch(**rows) -> pd.DataFrame:
    """A changes table (index the measures' keys; 'pop', 'mean_paired' and the eight wirings), zero except `rows`:
    key -> (pop, [eight per-wiring changes]) (mean_paired their mean) or (pop, per-wiring list, mean_paired)."""
    t = pd.DataFrame(0.0, index=[m.key for m in X.MEASURES], columns=['pop', 'mean_paired'] + list(W8))
    for k, v in rows.items():
        pop, per = v[0], list(v[1])
        t.loc[k, 'pop'] = pop
        t.loc[k, list(W8)] = per
        t.loc[k, 'mean_paired'] = v[2] if len(v) > 2 else float(np.mean(per))
    return t


def _bands(outside: dict = None) -> pd.DataFrame:
    """All runs inside every band except `outside`: band key -> wirings outside."""
    b = pd.DataFrame(True, index=list(W8), columns=[x.key for x in X.BANDS])
    for k, ws in (outside or {}).items():
        b.loc[list(ws), k] = False
    return b


def test_not_needed():
    v = X.verdict(_ch(pattern_cal=(0.002, [0.002] * 8), rho_on=(-0.004, [-0.004] * 8), mean_raw=(+0.02, [0.02] * 8)),
                  _bands(), c_flat=False)
    assert v['cls'] == 'not needed' and not v['not_needed_fails'] and not v['needed'] and not v['c_pending']


def test_better_without_needs_six_of_eight_and_the_population():
    six = [-0.004] * 6 + [0.0, 0.0]
    v = X.verdict(_ch(pattern_cal=(-0.004, six)), _bands(), c_flat=False)
    assert v['cls'] == 'better without' and 'pattern_cal' in v['gains'] and v['gains']['pattern_cal']['n_wirings'] == 6
    five = [-0.004] * 5 + [0.0] * 3
    assert X.verdict(_ch(pattern_cal=(-0.004, five)), _bands(), c_flat=False)['cls'] == 'not needed'
    # the population must gain too
    assert X.verdict(_ch(pattern_cal=(-0.002, six)), _bands(), c_flat=False)['cls'] == 'not needed'
    # one ρ higher by 0.010 or more
    v = X.verdict(_ch(rho_off=(0.012, [0.011] * 6 + [0.0, 0.0])), _bands(), c_flat=False)
    assert v['cls'] == 'better without' and set(v['gains']) == {'rho_off'}
    assert X.verdict(_ch(rho_off=(0.009, [0.011] * 8)), _bands(), c_flat=False)['cls'] == 'not needed'


def test_better_without_needs_nothing_failing():
    gain = [-0.004] * 8
    # a margin fails on the population only: no "better without", no "needed"
    v = X.verdict(_ch(pattern_cal=(-0.004, gain), rho_on=(-0.006, [-0.006, -0.006] + [0.0] * 6)), _bands(), c_flat=False)
    assert v['cls'] == 'not determined' and 'rho_on' in v['not_needed_fails'] and not v['needed']
    # one run outside a band
    v = X.verdict(_ch(pattern_cal=(-0.004, gain)), _bands({'band_opk_off': ['M3']}), c_flat=False)
    assert v['cls'] == 'not determined' and v['not_needed_fails']['band_opk_off']['n_runs_outside'] == 1


def test_needed_six_of_eight_rule():
    six = [0.005] * 6 + [0.0, 0.0]
    v = X.verdict(_ch(pattern_cal=(0.010, six)), _bands(), c_flat=False)
    assert v['cls'] == 'needed' and set(v['needed']) == {'pattern_cal'} and v['n_needed'] == 6
    assert any('pattern part (readout)' in s for s in v['groups']['the main targets'])
    five = [0.005] * 5 + [0.0] * 3
    v = X.verdict(_ch(pattern_cal=(0.010, five)), _bands(), c_flat=False)
    assert v['cls'] == 'not determined' and not v['needed'] and 'pattern_cal' in v['not_needed_fails']
    # six wirings worse but not the population: not "needed"
    assert X.verdict(_ch(pattern_cal=(0.002, six, 0.002)), _bands(), c_flat=False)['cls'] == 'not determined'
    # "not needed" holds on 6 of 8 wirings too: two worse wirings are fine, three are not
    two, three = [0.005] * 2 + [0.0] * 6, [0.005] * 3 + [0.0] * 5
    assert X.verdict(_ch(pattern_cal=(0.002, two)), _bands(), c_flat=False)['cls'] == 'not needed'
    v = X.verdict(_ch(pattern_cal=(0.002, three)), _bands(), c_flat=False)
    assert v['cls'] == 'not determined' and v['not_needed_fails']['pattern_cal']['n_wirings_worse'] == 3


def test_needed_margins_and_their_groups():
    v = X.verdict(_ch(rho_off=(-0.03, [-0.03] * 8), rate_rho_off=(-0.08, [-0.08] * 7 + [0.0]),
                      mean_cal=(-0.02, [-0.02] * 6 + [0.0] * 2), loss_raw=(0.006, [0.006] * 8)), _bands(), c_flat=False)
    assert v['cls'] == 'needed' and set(v['needed']) == {'rho_off', 'rate_rho_off', 'mean_cal', 'loss_raw'}
    g = v['groups']
    assert any('APL off' in s for s in g["the two ρ's"])
    assert len(g['the APL-off measures']) == 2                   # ρ with APL off and response rate ρ with APL off
    assert len(g['the main targets']) == 2                       # the mean APL effect (readout) and the raw-count loss
    # a margin's direction: a higher mean APL effect or a higher ρ never fails
    assert X.verdict(_ch(mean_raw=(0.05, [0.05] * 8), rho_on=(0.02, [0.0] * 8)), _bands(), c_flat=False)['cls'] == 'not needed'


def test_bands_six_of_eight():
    v = X.verdict(_ch(), _bands({'band_opk_off': W8[:6]}), c_flat=False)
    assert v['cls'] == 'needed' and 'band_opk_off' in v['needed']
    assert 'the APL-off measures' in v['groups'] and 'realism bands' in v['groups']
    v = X.verdict(_ch(), _bands({'band_rate_on': W8[:5]}), c_flat=False)
    assert v['cls'] == 'not determined'


def test_mean_paired_alone_blocks_not_needed():
    v = X.verdict(_ch(pattern_cal=(0.001, [0.0] * 8, 0.004)), _bands(), c_flat=False)
    assert v['cls'] == 'not determined' and v['not_needed_fails']['pattern_cal']['mean_paired'] == pytest.approx(0.004)


def test_part_c_flat_makes_not_determined_and_pending():
    assert X.verdict(_ch(), _bands(), c_flat=True)['cls'] == 'not determined'
    v = X.verdict(_ch(), _bands(), c_flat=None)
    assert v['cls'] == 'not needed' and v['c_pending']
    assert set(X.CLASSES) == {'not needed', 'needed', 'better without', 'not determined'}


def test_wirings_needed():
    assert X.wirings_needed(8) == 6 and X.wirings_needed(4) == 3


# ---------------------------------------------------------------------------
# part C and part E

def test_c_reading():
    rows = pd.DataFrame(dict(wiring=list(X.C_WIRINGS), c_kckc=[25.0, 1.0, 0.5, 30.0], c_loss=[0.0560, 0.0500, 0.0571, 0.0700],
                             fit_loss=[0.0561, 0.0520, 0.0570, 0.0600], b_loss=[0.0600, 0.0600, 0.0560, 0.0710]))
    c = X.c_reading(rows)
    p = c['per'].set_index('wiring')
    assert p.loc['real', 'returns'] and not p.loc['real', 'flat']            # back to the fit, B far
    assert p.loc['P1', 'more_evaluations'] and not p.loc['P1', 'returns']    # stays without KC-KC, lower than B
    assert p.loc['P3', 'flat'] and not p.loc['P3', 'returns']                # within 0.002 of both
    assert not p.loc['M2', 'returns']                                        # KC-KC high but the loss far from the fit
    assert c['complete'] and c['n_flat'] == 1 and c['c_flat'] is False
    rows2 = rows.assign(c_loss=[0.0561, 0.0520, 0.0571, 0.0600], b_loss=[0.0570, 0.0530, 0.0560, 0.0610], c_kckc=5.0)
    assert X.c_reading(rows2)['c_flat'] is True and X.c_reading(rows2)['n_flat'] == 4
    rows3 = rows2.copy()
    rows3.loc[3, 'c_loss'] = np.nan
    c3 = X.c_reading(rows3)
    assert not c3['complete'] and c3['c_flat'] is None


def test_e_alternatives():
    assert X.e_alternatives({1.5: 0.050, 3.0: 0.0531}, 0.056) == [1.5]
    assert X.e_alternatives({1.5: 0.053, 3.0: np.nan}, 0.056) == [1.5]
    assert X.e_alternatives({1.5: 0.0535, 3.0: 0.06}, 0.056) == []


# ---------------------------------------------------------------------------
# the starts file and the counted start

def _a_row(w='real', kc=0.0):
    p = dict(fixed_thr=220.4, wAPLKC_scale=0.0501, wKCAPL_scale=2.89e-05, wAPLPN_scale=0.376, wPNAPL_scale=0.000358,
             wKCKC_scale=kc)
    return dict(tag=f'A_{w}', wiring=w, params_on=json.dumps(p), **p), p


def test_starts_file_roundtrip(tmp_path):
    J = pytest.importorskip('experiments.E10_joint_fit')
    path = tmp_path / 'starts.json'
    row, p = _a_row()
    fit = dict(tag='rate0.068_aplboth_kckc0_k2_cap3_B', fixed_thr=180.0, wAPLKC_scale=0.06, wKCAPL_scale=3e-05,
               wAPLPN_scale=0.3, wPNAPL_scale=0.0004, kc_kc_inh=0.0)
    cur = X.merge_starts(path, {'A_real': X.a_start(row), 'C_real': X.c_start(fit, 'B_real')})
    assert set(cur) == {'A_real', 'C_real'}
    got = J.read_params(f'{path}:A_real')                     # E10's reader; the string metadata is ignored
    assert got == pytest.approx(p) and set(got) == set(X.PARAM_KEYS)
    c = J.read_params(f'{path}:C_real')
    assert c['wKCKC_scale'] == X.C_KCKC_START and c['fixed_thr'] == 180.0 and c['wAPLPN_scale'] == 0.3
    # the same entry again is fine; a changed one is refused; keys are never removed
    X.merge_starts(path, {'A_real': X.a_start(row)})
    with pytest.raises(RuntimeError):
        X.merge_starts(path, {'A_real': X.start_entry(dict(p, fixed_thr=221.0))})
    row2, _ = _a_row('P1')
    cur = X.merge_starts(path, {'A_P1': X.a_start(row2)})
    assert set(cur) == {'A_real', 'C_real', 'A_P1'}
    assert not list(tmp_path.glob('*.lock')) and not list(tmp_path.glob('.*.tmp'))
    # part A's start has KC-KC 0
    with pytest.raises(ValueError):
        X.a_start(_a_row(kc=26.4)[0])


def test_choose_counted():
    rows = [{X.CHOOSE_COL: 0.060, 'loss': 4.0}, {X.CHOOSE_COL: 0.055, 'loss': 4.5}]
    assert X.choose_counted(rows) == 1
    assert X.choose_counted(rows[:1]) == 0
    assert X.choose_counted([{X.CHOOSE_COL: 0.05}, {X.CHOOSE_COL: 0.05}]) == 0


# ---------------------------------------------------------------------------
# the population and paired readings

def _counts(seed: int, n_kc: int = 40) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    idx = pd.Index(range(n_kc), name='kc_id')
    return pd.DataFrame(rng.poisson(1.5, (n_kc, len(H.APL16_ODORS))).astype(float), index=idx, columns=list(H.APL16_ODORS))


def _run(seed):
    on, off = _counts(seed), _counts(seed + 100) + 1.0
    return {'on': ((on >= 2).astype(float), on), 'off': ((off >= 2).astype(float), off)}


def test_readings_population_and_paired(monkeypatch):
    E83 = X.E83
    monkeypatch.setattr(X, 'READOUTS', ('raw',))
    monkeypatch.setattr(X, '_dist', lambda o, ro: E83.pair_distances(o, 'raw'))
    monkeypatch.setattr(X, '_save_dist_cache', lambda: None)
    monkeypatch.setattr(X, 'scores_of', lambda d: dict(mean_on=float(d['on'].mean()), mean_off=float(d['off'].mean()),
                                                         first_on=float(E83._cs(d['on']).iloc[0])))
    outs, rows = {}, []
    for i, w in enumerate(W8):
        for kind, s in (('family', i), ('B', 50 + i)):
            outs[(kind, w)] = _run(s)
            rows.append(dict(kind=kind, wiring=w, flagged=(kind == 'B' and w == 'M3'), rate_on=0.05 + 0.001 * s,
                             rate_rho_on=0.6))
    outs[('C', 'real')] = _run(7)
    rows.append(dict(kind='C', wiring='real', flagged=False, rate_on=0.06, rate_rho_on=0.6))
    runs = pd.DataFrame(rows)
    rd = X.readings_table(outs, runs)
    # populations only where all their wirings have a run of the kind
    assert set(rd[rd.kind == 'B'].unit) == set(W8) | set(X.UNITS)
    assert set(rd[rd.kind == 'C'].unit) == {'real'}
    # the population is the mean matrix over the wirings
    d = {w: E83.pair_distances(outs[('B', w)], 'raw') for w in W8}
    pop = E83.population(list(d.values()))
    r = rd[(rd.kind == 'B') & (rd.unit == 'all8')].iloc[0]
    assert r.mean_on == pytest.approx(float(pop['on'].mean())) and r.first_on == pytest.approx(float(pop['on'].iloc[0]))
    assert r.n_flagged == 1 and r.run_rate_on == pytest.approx(np.mean([0.05 + 0.001 * (50 + i) for i in range(8)]))
    # the paired change from the family's same unit
    f = rd[(rd.kind == 'family') & (rd.unit == 'all8')].iloc[0]
    assert r.d_mean_on == pytest.approx(r.mean_on - f.mean_on)
    rw = rd[(rd.kind == 'B') & (rd.unit == 'P2')].iloc[0]
    fw = rd[(rd.kind == 'family') & (rd.unit == 'P2')].iloc[0]
    assert rw.d_first_on == pytest.approx(rw.first_on - fw.first_on)
    # changes_table: the population's change, the mean paired change and each wiring's
    ms = [X.Measure('m', 'mean distance, APL on', 'mean_on', 'raw', -1, 0.01, ('main',))]
    ch = X.changes_table(rd, 'B', measures=ms)
    per = [float(rd[(rd.kind == 'B') & (rd.unit == w)].mean_on.iloc[0] - rd[(rd.kind == 'family') & (rd.unit == w)].mean_on.iloc[0])
           for w in W8]
    assert ch.loc['m', 'pop'] == pytest.approx(r.d_mean_on)
    assert ch.loc['m', 'mean_paired'] == pytest.approx(np.mean(per))
    assert ch.loc['m', 'M2'] == pytest.approx(per[W8.index('M2')])


def test_pending_table():
    rd = pd.DataFrame([dict(kind=k, unit=u, readout='raw', block=b, vcf_on=0.3, vcf_off=0.2, vcf_diff=0.1, struct_wt_on=0.5)
                       for k, b0 in (('family', 0.04), ('B', 0.0)) for u, b in [(w, b0) for w in W8] + [('all8', b0)]])
    t = X.pending_table(rd)
    assert t.loc['block', 'pop'] == pytest.approx(-0.04) and bool(t.loc['block', 'fails'])
    assert not bool(t.loc['vcf_on', 'fails']) and t.loc['vcf_diff', 'slack'] == X.STRUCTURE_SLACK


def test_bands_of():
    ob = (0.155, 0.223)
    b = X.bands_of(0.068, 0.15, 0.6, 0.10, 0.10, ob)
    assert all(b.values())
    b = X.bands_of(0.068, 0.25, 0.30, 0.10, 0.30, ob)
    assert not b['band_ratio'] and not b['band_rate_rho_on'] and not b['band_opk_off'] and b['band_rate_on'] and b['band_opk_on']
    assert not X.bands_of(0.09, 0.18, 0.6, 0.1, 0.1, ob)['band_rate_on']


# ---------------------------------------------------------------------------
# part 0's comparisons

def test_exact_check_and_same_command():
    h0 = pd.Series(dict(eval=0, fixed_thr=102.2, loss=4.81, apl_cal_rmse=0.0560, rate_on=0.0668, time_s=9.0))
    best = pd.Series(dict(eval=120, fixed_thr=102.2, loss=4.81, apl_cal_rmse=0.0560, rate_on=0.0668, time_s=11.0))
    row = pd.Series(dict(loss=4.81, fit_apl_cal_rmse=0.0560, fit_rate_on=0.0668, k=2.0))
    e = X.exact_check(h0, best, row)
    assert e['exact'] and e['n_row_columns'] == 3
    h1 = h0.copy()
    h1['apl_cal_rmse'] = 0.05600001
    e = X.exact_check(h1, best, row)
    assert not e['exact'] and 'apl_cal_rmse' in e['history_differs'] and 'apl_cal_rmse' in e['row_differs']
    assert X.same_command(['a', '--kckc-fixed', '0'], ['a', '--kckc-fixed', '0.0'])
    assert not X.same_command(['a', '--kckc-fixed', '0'], ['a', '--kckc-fixed', '1'])
    assert not X.same_command(['a', '--tag', '_x'], ['a', '--tag', '_y'])


# ---------------------------------------------------------------------------
# tags and commands

def test_tags_match_e72():
    assert X.e10_tag('real', 2, None, '_cap3', rate=0.068) == 'rate0.068_aplboth_k2_cap3'
    assert X.e10_tag('P1', 2, None, '_cap3', rate=0.068) == 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_k2_cap3'
    assert X.e10_tag('M3', 2, None, '_cap3', rate=0.068) == 'wiring-curveball_within_type_matched-seed3_rate0.068_aplboth_k2_cap3'
    assert X.e10_tag('real', 2, 26.43035318456299, '_cap3_check0', rate=0.068) == 'rate0.068_aplboth_kckc26.4304_k2_cap3_check0'
    assert X.e10_tag('P2', 1.5, 0.0, '_cap3_E', rate=0.068) == 'wiring-claw_permutation_within_type-seed2_rate0.068_aplboth_kckc0_k1.5_cap3_E'


def test_fit_specs_and_params_commands():
    b = X.fit_specs('B')
    assert [s.name for s in b if s.part == 'B'] == [f'B_{w}' for w in W8]
    assert [s.wiring for s in b if s.part == 'B2'] == list(X.B2_WIRINGS)
    assert all(s.maxfev == 200 for s in b if s.part == 'B2') and all(s.maxfev == 150 for s in b if s.part == 'B')
    assert all(s.kckc_fixed == 0.0 for s in b)
    assert [s.name for s in X.fit_specs('B', ['P3'])] == ['B_P3']
    assert [s.wiring for s in X.fit_specs('C')] == list(X.C_WIRINGS) and all(s.kckc_fixed is None for s in X.fit_specs('C'))
    e = X.fit_specs('E')
    assert sorted({s.k for s in e}) == [1.5, 3.0] and {s.wiring for s in e} == set(X.E_WIRINGS)
    assert {s.wiring for s in X.fit_specs('E-rest', exponents=(3.0,))} == set(X.E_REST)
    cmd = X.fit_command(X.fit_specs('B', ['P1'])[0], python='py')
    s = ' '.join(cmd)
    assert s.startswith('py -u experiments/E10_joint_fit.py -k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl '
                        '--apl-readout both --tag _cap3_B --maxfev 150 --experiment E110_no_kckc/fits/B_P1')
    assert '--wiring claw_permutation_within_type --wiring-seed 1' in s and '--kckc-fixed 0.0' in s
    assert s.endswith('--init-params results/E110_no_kckc/starts.json:A_P1')
    c = ' '.join(X.fit_command(X.fit_specs('C', ['real'])[0], python='py'))
    assert '--kckc-fixed' not in c and c.endswith('--init-params results/E110_no_kckc/starts.json:C_real')
    m = X.wiring_args('M2')
    assert m[0] == '--wiring-file' and m[1].endswith('E59_coinput_shuffle/wPNKC_curveball_within_type-seed2_matched.parquet')
    assert m[2:] == ['--wiring-tag', 'curveball_within_type_matched-seed2'] and X.wiring_args('real') == []


@pytest.mark.skipif(not _shared(), reason='needs the shared directory')
def test_exact_commands_and_family_tags():
    """The second starts and part D start from each wiring's own fit; the family's tags are E10's (E72's fits.csv)."""
    m = X.E83.members()
    for w in W8:
        assert X.e10_tag(w, 2, None, '_cap3') == m[w][1]
    s = X.fit_specs('B', ['M1'])[1]
    cmd = X.fit_command(s, python='py')
    i = cmd.index('--init-from')
    assert cmd[i:i + 5] == ['--init-from', m['M1'][1], '--init-from-experiment', m['M1'][0], '--init-from-exact']
    d = X.fit_specs('D', ['real'])[0]
    assert d.kckc_fixed == pytest.approx(0.5 * X.E83.saved_fit('real')['p']['wKCKC_scale'])
    tag, cmd = X.check_specs()['exact_real']
    assert float(cmd[cmd.index('--kckc-fixed') + 1]) == X.E83.saved_fit('real')['p']['wKCKC_scale'] and '--smoke' in cmd
