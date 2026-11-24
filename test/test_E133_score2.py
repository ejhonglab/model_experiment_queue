#!/usr/bin/env python3
"""Tests for E133 step 2's score (experiments/E133_score2.py; model_experiment_queue): the class rules on made-up d's at
and around every threshold (G/2, G/4, 0.004, 0.003, 3 of the 4 wirings, every half), "no class", the provisional
reading, step 3's selection and its tie rule, G and the d table on made-up readings, the whole score on made-up runs
(provisional and final files), and that no lane loads the module. No simulations. Needs the shared directory (E10 reads
the data-target registry at import); skipped without it."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
S = pytest.importorskip('experiments.E133_score2')
E = S.E
L, HB = S.LEFT, S.HEMI
G0 = 0.016                                      # G/2 0.008, G/4 0.004: the fractions bind over the floors


def _hr(d, sign, G=G0, halves=None, d_w=None):
    """host_reading with the four wirings at `d` (or `d_w`) and every half at `halves` (default `d`)."""
    d_w = [d] * 4 if d_w is None else d_w
    hv = {t: (d if halves is None else halves.get(t, d)) for t in S.HALVES}
    return S.host_reading(d_w, hv, sign, G)


# ---------------------------------------------------------------------------
# closes / opens and a part, at and around each threshold

def test_closes_at_and_around_G_half():
    assert _hr(-0.008, -1)['full']                       # at -G/2: closes (at most -G/2)
    r = _hr(-0.0079, -1)
    assert not r['full'] and r['part']                   # just short: a part's size only
    assert _hr(+0.008, +1)['full']                       # the hemibrain's mirror opens
    assert not _hr(+0.0079, +1)['full']
    assert not _hr(-0.008, +1)['full'] and not _hr(-0.008, +1)['part']    # the wrong way on the hemibrain
    assert not _hr(+0.008, -1)['full']                   # the left raised: not closing


def test_floor_0004_binds_when_G_is_small():
    G = 0.006                                            # G/2 = 0.003 < 0.004
    assert _hr(-0.004, -1, G)['full']
    r = _hr(-0.0039, -1, G)
    assert not r['full'] and r['part']                   # over G/4 = 0.0015 and 0.003
    assert _hr(+0.004, +1, G)['full'] and not _hr(+0.0039, +1, G)['full']


def test_part_at_and_around_G_quarter():
    r = _hr(-0.004, -1)                                  # at G/4
    assert r['part'] and not r['full']
    assert not _hr(-0.0039, -1)['part']
    assert _hr(+0.004, +1)['part'] and not _hr(+0.0039, +1)['part']


def test_part_floor_0003_binds_when_G_is_small():
    G = 0.006                                            # G/4 = 0.0015 < 0.003
    assert _hr(-0.003, -1, G)['part'] and not _hr(-0.0029, -1, G)['part']
    assert _hr(+0.003, +1, G)['part'] and not _hr(+0.0029, +1, G)['part']


def test_three_of_four_wirings():
    r = _hr(None, -1, halves={t: -0.01 for t in S.HALVES}, d_w=[-0.02, -0.02, -0.02, +0.004])
    assert r['n_toward'] == 3 and r['d'] == pytest.approx(-0.0140) and r['full']
    r = _hr(None, -1, halves={t: -0.01 for t in S.HALVES}, d_w=[-0.03, -0.03, +0.004, +0.002])
    assert r['n_toward'] == 2 and r['d'] < -G0 / 2 and not r['full'] and not r['part']
    # a wiring at exactly 0 is not lower
    r = _hr(None, -1, halves={t: -0.01 for t in S.HALVES}, d_w=[-0.02, -0.02, 0.0, 0.0])
    assert r['n_toward'] == 2 and not r['full']
    # the hemibrain: higher on 3 of 4
    r = _hr(None, +1, halves={t: 0.01 for t in S.HALVES}, d_w=[0.02, 0.02, 0.02, -0.004])
    assert r['n_toward'] == 3 and r['full']


def test_every_half():
    assert _hr(-0.01, -1, halves={'A': -0.001})['full']
    for v in (0.0, +0.001):
        r = _hr(-0.01, -1, halves={'C': v})
        assert not r['every_half'] and not r['full'] and not r['part']
    assert not _hr(+0.01, +1, halves={'D': 0.0})['full']


def test_d_is_the_exact_mean_and_incomplete():
    r = _hr(-0.008, -1, d_w=[-0.008, -0.008, -0.008, -0.008])
    assert r['d'] == -0.008 and r['full']
    r = _hr(-0.01, -1, d_w=[-0.006, -0.010, -0.010, -0.006])        # mean -0.008 exactly
    assert r['d'] == -0.008 and r['full']
    r = S.host_reading([-0.01, -0.01, -0.01], {t: -0.01 for t in S.HALVES}, -1, G0)
    assert not r['complete'] and not r['full'] and not r['part'] and not r['under_none']
    r = S.host_reading([-0.01] * 4, {t: -0.01 for t in S.HALVES[:3]}, -1, G0)
    assert not r['complete']
    assert not S.host_reading([-0.01] * 4, {t: -0.01 for t in S.HALVES}, -1, float('nan'))['complete']


# ---------------------------------------------------------------------------
# an arm's class

def test_arm_class_carries_part_one_sided_none_unclear():
    full_L, full_H = _hr(-0.010, -1), _hr(+0.010, +1)
    part_L, part_H = _hr(-0.005, -1), _hr(+0.005, +1)
    assert S.arm_class(full_L, full_H) == ('carries', None)
    assert S.arm_class(part_L, part_H) == ('a part', None)
    assert S.arm_class(full_L, part_H) == ('a part', None)           # closes, and a part on the other host
    assert S.arm_class(full_L, _hr(+0.004, +1)) == ('a part', None)  # the other host at G/4 exactly: a part
    # one-sided: under G/4 on the other host, either way
    for dh in (+0.0039, -0.0039, 0.0):
        assert S.arm_class(full_L, _hr(dh, +1)) == ('one-sided', 'closes (left host)')
    for dl in (-0.0039, +0.0039):
        assert S.arm_class(_hr(dl, -1), full_H) == ('one-sided', 'opens (hemibrain host)')
    # at G/4 on the other host but short of the part rule (2 of 4 wirings): neither a part nor under G/4
    two = _hr(None, +1, halves={t: 0.004 for t in S.HALVES}, d_w=[0.012, 0.012, -0.004, -0.004])
    assert two['d'] == pytest.approx(0.004) and not two['part'] and not two['under_quarter']
    assert S.arm_class(full_L, two)[0] == 'unclear'
    # the other host moved the wrong way by G/4 or more: unclear (|d|), one-sided under the directional reading
    for dh in (-0.004, -0.010):
        assert S.arm_class(full_L, _hr(dh, +1))[0] == 'unclear'
        assert S.arm_class(full_L, _hr(dh, +1), directional=True) == ('one-sided', 'closes (left host)')
    # none: under 0.003 on both hosts, either way
    assert S.arm_class(_hr(+0.0029, -1), _hr(-0.0029, +1)) == ('none', None)
    assert S.arm_class(_hr(-0.0029, -1), _hr(+0.0029, +1)) == ('none', None)
    assert S.arm_class(_hr(+0.003, -1), _hr(0.0, +1))[0] == 'unclear'
    assert S.arm_class(_hr(-0.003, -1), _hr(0.0, +1))[0] == 'unclear'      # a part's size on one host only
    # incomplete
    inc = S.host_reading([-0.01] * 3, {t: -0.01 for t in S.HALVES}, -1, G0)
    assert S.arm_class(inc, full_H) == ('incomplete', None)


def _readings(spec):
    """{arm: (d on the left, d on the hemibrain)} -> arm_readings-like records."""
    out = {}
    for arm, (dl, dh) in spec.items():
        hr = {L: _hr(dl, -1), HB: _hr(dh, +1)}
        cls, side = S.arm_class(hr[L], hr[HB])
        out[arm] = dict(rule_class=cls, side=side, hosts=hr)
    return out


SPEC = {'C': (-0.010, +0.012), 'S': (-0.005, +0.006), 'A': (-0.015, +0.001), 'B': (+0.001, -0.001),
        'K': (-0.005, -0.005), 'ALL': (-0.020, +0.020)}
WANT = {'C': 'carries', 'S': 'a part', 'A': 'one-sided', 'B': 'none', 'K': 'unclear', 'ALL': 'carries'}


def test_classes_no_class_provisional_and_G():
    rd = _readings(SPEC)
    c = S.classes(rd, {}, G0)
    assert {a: r['cls'] for a, r in c.items()} == WANT
    c = S.classes(rd, {'C': ['S2_C_malecns-left_seed2']}, G0)
    assert c['C']['cls'] is None and c['C']['would_be'] == 'carries' and c['C']['no_class_runs'] == ['S2_C_malecns-left_seed2']
    assert c['S']['cls'] == 'a part'
    p = S.classes(rd, {}, G0, provisional=True)
    assert all(r['cls'] is None for r in p.values()) and p['ALL']['provisional_reading'] == 'carries'
    for g in (0.0, -0.01, float('nan')):
        assert all(r['cls'] is None for r in S.classes(rd, {}, g).values())


def test_step3_selection_and_ties():
    c = S.classes(_readings(SPEC), {}, G0)
    assert S.step3_selection(c) == ['ALL', 'A', 'C']               # largest |d|: 0.020, 0.015, 0.012
    assert S.step3_selection(c, include_all=False) == ['A', 'C', 'S']
    assert S.step3_selection(c, k=1) == ['ALL']
    # none and unclear never go; an arm without a class does not, its would-be list keeps it
    nc = S.classes(_readings(SPEC), {'ALL': ['x']}, G0)
    assert S.step3_selection(nc) == ['A', 'C', 'S']
    wb = {a: dict(r, would_be=r.get('cls') or r.get('would_be')) for a, r in nc.items()}
    assert S.step3_selection(wb, key='would_be') == ['ALL', 'A', 'C']
    # ties: the same larger |d| -> the other |d|; then the class; then the arm order
    tie = S.classes(_readings({'C': (-0.010, +0.005), 'S': (-0.010, +0.006), 'A': (-0.004, +0.010)}), {}, G0)
    assert [tie[a]['cls'] for a in 'CSA'] == ['a part', 'a part', 'a part']
    assert S.step3_selection(tie) == ['S', 'C', 'A']
    tie = S.classes(_readings({'K': (-0.010, +0.010), 'B': (-0.010, +0.010)}), {}, G0)
    assert S.step3_selection(tie) == ['B', 'K']
    tie = S.classes(_readings({'C': (-0.010, +0.001), 'S': (-0.010, +0.001)}), {}, G0)
    assert S.step3_selection(tie) == ['C', 'S']


# ---------------------------------------------------------------------------
# G, d, the whole score on made-up readings

def _row(v):
    r = {c: v for c in S.CAL_STATS}
    r.update({c: v for c in S.RAW_STATS})
    return r


def _step1(shift=G0):
    """readings1-like rows: each calyx's seed s at base + 0.001 s (the left `shift` above the hemibrain)."""
    rows = []
    for c, base in ((HB, 0.050), (L, 0.050 + shift), (E.RIGHT, 0.055)):
        for s in (0,) + E.SEEDS1:
            for ro in E.READOUTS:
                rows.append(dict(tag=E.tag_of(E.Point(1, c, s, None)), calyx=c, seed=s, readout=ro, **_row(base + 0.001 * s)))
    return pd.DataFrame(rows)


def _step2(rd1, spec=SPEC, skip=()):
    own = rd1.set_index(['calyx', 'seed', 'readout'])
    rows = []
    for arm, (dl, dh) in spec.items():
        for h, d in ((L, dl), (HB, dh)):
            for s in S.SEEDS:
                tag = E.tag_of(E.Point(2, h, s, arm))
                if tag in skip:
                    continue
                for ro in E.READOUTS:
                    rows.append(dict(tag=tag, arm=arm, calyx=h, seed=s, readout=ro,
                                     **_row(float(own.loc[(h, s, ro), 'all_pattern']) + d)))
    return pd.DataFrame(rows)


def test_gap_and_local_gap():
    G = S.gap(_step1())
    assert G['all']['G'] == pytest.approx(G0) and G['all']['local'] == pytest.approx(G0)
    assert G['C']['G'] == pytest.approx(G0) and G['all']['n_left'] == 8
    rd1 = _step1()
    rd1 = rd1[~((rd1['calyx'] == L) & (rd1['seed'] == 7))]
    assert np.isnan(S.gap(rd1)['all']['G']) and S.gap(rd1)['all']['local'] == pytest.approx(G0)


def test_d_table():
    rd1 = _step1()
    dt = S.d_table(_step2(rd1), rd1)
    x = dt[(dt['arm'] == 'C') & (dt['host'] == L) & (dt['readout'] == 'cal') & (dt['stat'] == 'all_pattern')].iloc[0]
    assert x['n'] == 4 and x['d'] == pytest.approx(-0.010) and x['sd_wirings'] == pytest.approx(0.0, abs=1e-12)
    assert x['sd_redraws'] == pytest.approx(0.001 * np.std(np.arange(1, 8), ddof=1))
    assert x['sd_eight'] == pytest.approx(0.001 * np.std(np.arange(0, 8), ddof=1))
    y = dt[(dt['arm'] == 'ALL') & (dt['host'] == HB) & (dt['readout'] == 'raw') & (dt['stat'] == 'all_pattern')].iloc[0]
    assert y['d'] == pytest.approx(0.020)
    assert set(dt['stat']) >= {'all_pattern_p42', 'all_pattern_p78', 'all_pattern_p3', 'A_pattern'}
    # a missing run: that wiring's d is NaN and n is 3
    dt = S.d_table(_step2(rd1, skip=('S2_C_malecns-left_seed2',)), rd1)
    x = dt[(dt['arm'] == 'C') & (dt['host'] == L) & (dt['readout'] == 'cal') & (dt['stat'] == 'all_pattern')].iloc[0]
    assert x['n'] == 3 and np.isnan(x['d_seed2']) and x['d'] == pytest.approx(-0.010)


def test_verdict2_on_made_up_readings():
    rd1 = _step1()
    dt = S.d_table(_step2(rd1), rd1)
    v = S.verdict2(dt, S.gap(rd1), {}, provisional=False)
    assert {a: r['cls'] for a, r in v['arms'].items()} == WANT
    assert v['step3']['selected'] == ['ALL', 'A', 'C'] and v['step3']['selected_without_ALL'] == ['A', 'C', 'S']
    assert v['carries'] == ['C', 'ALL']
    assert v['arms']['C']['groups']['p42'][L] == pytest.approx(-0.010) and v['arms']['C']['raw_pattern'][HB] == pytest.approx(0.012)
    assert 'cls_directional' not in v['arms']['K'] and v['arms']['A']['side'] == 'closes (left host)'
    p = S.verdict2(dt, S.gap(rd1), {}, provisional=True)
    assert all(r['cls'] is None for r in p['arms'].values()) and p['step3']['selected'] is None
    nc = S.verdict2(dt, S.gap(rd1), {'A': ['S2_A_hemibrain-argmax_seed1']}, provisional=False)
    assert nc['arms']['A']['cls'] is None and nc['arms']['A']['would_be'] == 'one-sided'
    assert nc['step3']['selected'] == ['ALL', 'C', 'S'] and nc['step3']['would_be_without_no_class'] == ['ALL', 'A', 'C']


def _runs(rd, step):
    """runs-table rows (realism, parameters, sizes) for the runs of readings `rd`."""
    rows = []
    for r in rd[rd['readout'] == 'cal'].to_dict('records'):
        rows.append(dict(tag=r['tag'], step=step, calyx=r['calyx'], seed=r['seed'], arm=r.get('arm', ''), rate_on=0.068,
                         rate_off=0.15, ratio=2.2, rate_rho_on=0.5, ratio_g=2.4, rate_on_g=0.05, rate_on_apbp=0.07,
                         rate_on_ab=0.08, rate_off_g=0.12, rate_off_apbp=0.15, rate_off_ab=0.17, max_count_on=12,
                         max_count_off=20, n_kc_ge20_on=0, n_kc_ge20_off=3, thr_factor=1.0, apl_factor=1.0,
                         fixed_thr=200.0 * (1.1 if step == 2 else 1.0), wAPLKC_scale=2.0, wKCAPL_scale=0.01, wAPLPN_scale=1.0,
                         wPNAPL_scale=0.01, wKCKC_scale=0.5, in_bands=True, corner=False, guard_max_count=False,
                         no_class=False, realism_flag='', size_C_off_share=0.19 if r.get('arm') == 'C' else np.nan,
                         **{b: True for b in E.BAND_COLS}))
    return pd.DataFrame(rows)


def test_runs2_table_beside_own():
    rd1 = _step1()
    r2 = S.runs2_table(_runs(_step2(rd1), 2), _runs(rd1, 1))
    assert np.allclose(r2['fixed_thr_vs_own'], 1.1) and np.allclose(r2['wAPLKC_scale_vs_own'], 1.0)
    assert set(r2['own_tag']) == {E.tag_of(E.Point(1, h, s, None)) for h in S.HOSTS for s in S.SEEDS}
    assert (r2['ratio_g_own'] == 2.4).all()


def test_no_class_tags_take_own_runs_too():
    rd1 = _step1()
    r2 = _runs(_step2(rd1), 2)
    r1 = _runs(rd1, 1)
    r2.loc[r2['tag'] == 'S2_B_malecns-left_seed3', 'no_class'] = True
    nc = S.no_class_tags(r2, r1)
    assert nc['B'] == ['S2_B_malecns-left_seed3'] and nc['C'] == []
    r1.loc[r1['tag'] == 'S1_hemibrain-argmax_seed2', 'no_class'] = True
    nc = S.no_class_tags(r2, r1)
    assert all('S1_hemibrain-argmax_seed2' in nc[a] for a in S.ARMS)
    # a run of the right without a class is not a host's
    r1.loc[r1['tag'] == 'S1_hemibrain-argmax_seed2', 'no_class'] = False
    r1.loc[r1['tag'] == 'S1_malecns-right_seed1', 'no_class'] = True
    assert S.no_class_tags(r2, r1)['C'] == []


def test_score2_files_partial_then_final(tmp_path, monkeypatch):
    rd1 = _step1()
    rd1.to_csv(tmp_path / 'readings1.csv', index=False)
    _runs(rd1, 1).to_csv(tmp_path / 'runs1.csv', index=False)
    full = _step2(rd1)
    state = dict(skip=('S2_K_hemibrain-argmax_seed3', 'S2_ALL_malecns-left_seed0'))

    def fake_collect(points=None):
        rd = full[~full['tag'].isin(state['skip'])]
        outs = {t: None for t in rd['tag'].unique()}
        return outs, _runs(rd, 2), [dict(tag=t, why='no row') for t in state['skip']], ['S2_PA_malecns-left_seed1']

    monkeypatch.setattr(S, 'collect2', fake_collect)
    monkeypatch.setattr(S, 'readings2', lambda outs, runs: full[full['tag'].isin(list(outs))].reset_index(drop=True))
    assert S.score2(partial=False, check_own=False, out=tmp_path) is None
    assert not list(tmp_path.glob('*2*.json'))
    v = S.score2(partial=True, check_own=False, out=tmp_path)
    assert v['provisional'] and v['n_runs_in'] == 46 and all(r['cls'] is None for r in v['arms'].values())
    assert v['arms']['C']['provisional_reading'] == 'carries' and v['arms']['K']['provisional_reading'] == 'incomplete'
    for name in S.OUTPUTS:
        assert (tmp_path / S.out_name(name, True)).exists() and not (tmp_path / name).exists()
    md = (tmp_path / 'step2_table_provisional.md').read_text()
    assert 'PROVISIONAL' in md.splitlines()[0] and 'carries |' not in md.replace('provisional reading: carries', '')
    state['skip'] = ()
    v = S.score2(partial=False, check_own=False, out=tmp_path)
    assert not v['provisional'] and {a: r['cls'] for a, r in v['arms'].items()} == WANT
    assert v['step3']['selected'] == ['ALL', 'A', 'C'] and v['other_step2_rows'] == ['S2_PA_malecns-left_seed1']
    for name in S.OUTPUTS:
        assert (tmp_path / name).exists() and not (tmp_path / S.out_name(name, True)).exists()
    vj = json.loads((tmp_path / 'verdict2.json').read_text())
    assert vj['arms']['A']['side'] == 'closes (left host)' and vj['thresholds']['G_half'] == pytest.approx(G0 / 2)
    md = (tmp_path / 'step2_table.md').read_text()
    assert md.startswith('# E133 step 2 (generated)') and '| C | carries |' in md and 'off_share 0.19' in md
    assert 'Step-2 runs without a class or with a realism flag: none.' in md
    assert len(pd.read_csv(tmp_path / 'd2.csv')) == len(S.ARMS) * len(S.HOSTS) * (len(S.CAL_STATS) + len(S.RAW_STATS))


def test_out_name_and_constants():
    assert S.out_name('d2.csv', True) == 'd2_provisional.csv' and S.out_name('verdict2.json', False) == 'verdict2.json'
    assert {pt.arm for pt in E.planned_step2()} == set(S.ARMS) and set(S.HOSTS) == set(E.HOSTS)
    assert S.SEEDS == (0, 1, 2, 3) and S.SIGN[L] < 0 < S.SIGN[HB]
    # no output name the runner or step 1 writes or globs (runs_<step>_*.csv, readings1.csv, verdict1.json)
    for name in S.OUTPUTS:
        for prov in (False, True):
            n = S.out_name(name, prov)
            assert not n.startswith('runs_') and n not in ('readings1.csv', 'runs1.csv', 'verdict1.json')


# ---------------------------------------------------------------------------
# the import graph: no lane loads this module

def test_no_lane_loads_the_score():
    """E108's `_still_clean` stops a lane when a file it loaded changes: the runner, its preloads and the builders must
    not load experiments/E133_score2.py, and no other module of experiments/ or kcq/ names it."""
    code = ('import sys; import experiments.E133_calyx_transplants as m; m._preload(); import experiments.E133_builders; '
            "print('experiments.E133_score2' in sys.modules)")
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=str(H.HERE), timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().splitlines()[-1] == 'False'
    here = Path(H.HERE)
    named = [p.relative_to(here).as_posix() for d in ('experiments', 'kcq') for p in (here / d).rglob('*.py')
             if p.name != 'E133_score2.py' and 'E133_score2' in p.read_text()]
    assert named == []
