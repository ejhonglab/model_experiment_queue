#!/usr/bin/env python3
"""Tests for E103 (the flies' APL across odors in the held-out vial series; experiments/E103_apl_across_odors_heldout.py),
on made-up stimulus files, calyx means, maps and rows: blocks of three, four and five presentations; two
paraffin-oil-only blocks and none; a block left out; the first two presentations against all; the vial labels; the
date notes' order of the three-vial odor; the split count of 5 and 8 and the vectorized split rule against E98's; the
least detectable difference; Student's t interval; the classes inside / above / below, the 6 of 8 and the held class
against the reading-dependent one; X2 / X3 over the estimates; X1's three sides; Y1 on vials with a known persistent
SD; M's share on maps built to be the odor's and the block's; Z's sign rule; the refusal of a raw-only row, of a
recording outside the lists and of any map of the series before part A's tables and lock exist; and, with the shared
directory, the eligibility table rebuilt from the series' stimulus files against the entry's lists."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
M = pytest.importorskip('experiments.E103_apl_across_odors_heldout')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        d = H.shared_dir()
        return ((d / 'data' / 'yang_apl_recordings' / 'recordings.csv').exists()
                and (d / 'results' / 'E98_apl_response_strength' / 'a1_responses.csv').exists())
    except Exception:
        return False


# ---------------------------------------------------------------------------
# made-up stimulus files (pin 37 the solvent; 38 to 40 pa's three vials at 10^-5.34, 41 ms, 42 paa, 43 pa vial 1)

PINS = ['NA'] * 36 + ['paraffin oil', 'pentyl acetate1(-5.340)', 'pentyl acetate2(-5.340)', 'pentyl acetate3(-5.340)',
                      'methyl salicylate(-2)', 'phenylacetaldehyde1(-3.431)', 'pentyl acetate1(-4.386)']


def _yaml(tmp_path, name, seq, nstim, block_out=None):
    import yaml
    d = dict(nblock=len(seq) // nstim, nstim=nstim, odor_pin=PINS, odor_stim_seq=list(seq))
    if block_out is not None:
        d['block_out'] = list(block_out)
    f = tmp_path / f'{name}.yaml'
    f.write_text(yaml.safe_dump(d))
    return f


def test_vial_labels():
    assert M.vial_label('pentyl acetate2') == ('pa_vial2', 2)
    assert M.vial_label('pentyl acetate1') == ('pa', 1)
    assert M.vial_label('phenylacetaldehyde3') == ('paa_vial3', 3)
    assert M.vial_label('ms') == ('ms', 1) and M.vial_label('paa') == ('paa', 1)
    assert M.vial_label('pin12') == ('pin12', 1)                         # only a full odor name carries a vial digit


def test_blocks_of_three_with_two_solvent_blocks(tmp_path):
    # pfo + two presentations; the all-solvent blocks are blocks 3 and 5
    seq = [37, 38, 38, 37, 39, 39, 37, 37, 37, 37, 41, 41, 37, 37, 37]
    t = M.table_from_yaml(_yaml(tmp_path, 'r3', seq, 3), 'r3')
    assert list(t.groupby('block1')['pfo_block'].first()) == [False, False, True, False, True]
    assert set(t.loc[t['odor'] != 'pfo', 'odor']) == {'pa', 'pa_vial2', 'ms'}
    assert t.loc[t['odor'] == 'pa_vial2', 'vial'].eq(2).all()
    m = np.zeros(len(seq))
    m[6:9] = [0.01, 0.02, 0.04]                                         # solvent block 3, positions 1 to 3
    m[12:15] = [0.03, 0.04, 0.06]                                       # solvent block 5
    m[1:3] = [0.30, 0.32]                                               # pa vial 1 at -5.34
    om = M.responses_from_means(m, t, 'order_matched', n_first=M.N_FIRST).set_index('odor')
    base2, base3 = np.mean([0.02, 0.04]), np.mean([0.04, 0.06])         # the mean of the two blocks at each position
    assert om.loc['pa', 'response'] == pytest.approx(np.mean([0.30 - base2, 0.32 - base3]))
    assert om.loc['pa', 'n_reps'] == 2


def test_blocks_of_four_without_a_solvent_block(tmp_path):
    seq = [37, 38, 38, 38, 37, 39, 39, 39, 37, 40, 40, 40]
    t = M.table_from_yaml(_yaml(tmp_path, 'r4', seq, 4), 'r4')
    assert not t['pfo_block'].any()
    m = np.arange(len(seq)) / 100.0
    with pytest.raises(RuntimeError):                                    # no carrying reading without a solvent block
        M.responses_from_means(m, t, 'order_matched', n_first=M.N_FIRST)
    fib = M.responses_from_means(m, t, 'first_in_block', n_first=M.N_FIRST).set_index('odor')
    base = np.mean([m[0], m[4], m[8]])                                  # every block's position-1 solvent presentation
    assert fib.loc['pa', 'response'] == pytest.approx(np.mean([m[1], m[2]]) - base)
    assert fib.loc['pa', 'n_reps'] == 2 and fib.loc['pa', 'rep3'] == pytest.approx(m[3] - base)


def test_blocks_of_five_first_two_against_all(tmp_path):
    seq = [37, 38, 38, 38, 38] + [37] * 5 + [37, 40, 40, 40, 40]
    t = M.table_from_yaml(_yaml(tmp_path, 'r5', seq, 5), 'r5')
    m = np.zeros(len(seq))
    m[5:10] = [0.0, 0.01, 0.02, 0.03, 0.04]
    m[1:5] = [0.2, 0.3, 0.6, 0.9]
    two = M.responses_from_means(m, t, 'order_matched', n_first=2).set_index('odor')
    every = M.responses_from_means(m, t, 'order_matched', n_first=None).set_index('odor')
    assert two.loc['pa', 'response'] == pytest.approx(np.mean([0.2 - 0.01, 0.3 - 0.02]))
    assert every.loc['pa', 'response'] == pytest.approx(np.mean([0.2 - 0.01, 0.3 - 0.02, 0.6 - 0.03, 0.9 - 0.04]))
    assert two.loc['pa', 'n_reps'] == 2 and every.loc['pa', 'n_reps'] == 4
    assert two.loc['pa', 'rep4'] == pytest.approx(0.9 - 0.04)              # every presentation kept in the rep columns


def test_a_block_left_out_and_a_stimulus_read_at_its_first_block(tmp_path):
    # pa vial 1 at -5.34 in blocks 1 and 4 (two passes); block 2 (ms) left out by the caller (EXTRA_OUT's kind)
    seq = [37, 38, 38, 37, 41, 41, 37, 37, 37, 37, 38, 38, 37, 39, 39]
    t = M.table_from_yaml(_yaml(tmp_path, 'rl', seq, 3), 'rl', extra_out=(2,))
    assert t.loc[t['block1'] == 2, 'out'].all()
    m = np.arange(len(seq)) / 100.0
    r = M.responses_from_means(m, t, 'order_matched', n_first=2)
    assert 'ms' not in set(r['odor'])
    assert int(r.loc[r['odor'] == 'pa', 'block1'].iloc[0]) == 1          # the first of its two blocks
    both = M.responses_from_means(m, t, 'order_matched', which='all', n_first=None)
    assert sorted(both.loc[both['odor'] == 'pa', 'block1']) == [1, 4]


def test_note_order_of_the_three_vial_odor(tmp_path):
    note = '\n'.join(['20211022 GeroAPL298', '37 PFO 0', '38 PA2 -5.340 10', '43 PA3 -4.386 15', '',
                      '1FLfun1 (Bias Voltage 36)', 'z at 15 um', '1FRfun2', 'CHANGE PA CONTROLS TO PHE!!!!',
                      '38: Phe2 -4.386', '1FLfun4', 'RUN IT AGAIN AS 1FRfun9', '1FRfun5 (Bias)', '1FRfun5 (again)',
                      'Change to PA!!!!!', '1FRfun8'])
    p = tmp_path / '20211022.txt'
    p.write_text(note)
    assert M.note_vial_order(p) == {'1FLfun1': 'pa', '1FRfun2': 'pa', '1FLfun4': 'paa', '1FRfun5': 'paa', '1FRfun8': 'pa'}
    # without a pin table: before the first change, the other odor of its target
    p2 = tmp_path / '20211029.txt'
    p2.write_text('\n'.join(['3FLfun1', '3FLfun2', 'CHANGE TO PHE!!!!', '3FLfun3', 'Change tank, vent 2s PFO', '3FLfun5']))
    assert M.note_vial_order(p2) == {'3FLfun1': 'pa', '3FLfun2': 'pa', '3FLfun3': 'paa', '3FLfun5': 'paa'}


# ---------------------------------------------------------------------------
# the split rule, the least detectable difference, the intervals

def test_split_count_of_five_and_eight_and_e98s_rule():
    rng = np.random.default_rng(3)
    a, b = rng.normal(0, 0.4, 8), rng.normal(0.5, 0.4, 5)
    share, n = M.split_share(a, b)
    assert n == 1287
    sa = pd.Series(a, index=[f'f{i}' for i in range(8)])
    sb = pd.Series(b, index=[f'g{i}' for i in range(5)])
    r = M.split_rule(sa, sb)
    assert r['n_splits'] == 1287 and r['share'] == pytest.approx(share) and r['differ'] == (share < M.SPLIT_SHARE)
    # identical sets never differ; sets far apart always do
    assert M.split_share(a, a[:5])[0] >= M.SPLIT_SHARE
    assert M.split_share(a, a[:5] + 10.0)[0] == pytest.approx(1 / 1287)


def test_least_detectable_difference():
    sd = 0.46                                                            # the entry's P1: a factor of about 1.4
    a = M.normal_scores(8) * sd
    ld = M.least_detectable(a, n_new=5)
    assert ld['n_splits'] == 1287 and ld['sd_log2'] == pytest.approx(sd)
    assert ld['factor_normal_approx'] == pytest.approx(2 ** (1.6449 * sd * np.sqrt(1 / 8 + 1 / 5)))
    assert 1.25 < ld['factor_up'] < 1.55 and ld['factor_down'] == pytest.approx(1 / ld['factor_up'], rel=1e-5)
    assert ld['factor_normal_scores'] == pytest.approx(ld['factor_up'], rel=1e-5)
    # at the returned difference the rule says differ, just below it it does not
    b0 = a.mean() + M.normal_scores(5) * sd
    d = ld['ldd_up_log2']
    assert M.split_share(a, b0 + d)[0] < M.SPLIT_SHARE <= M.split_share(a, b0 + d - 1e-4)[0]


def test_t_interval():
    from scipy.stats import t as student_t
    x = np.log([0.8, 1.1, 0.9, 1.3, 0.7])
    g, lo, hi = M.t_interval(x)
    h = student_t.ppf(0.95, 4) * x.std(ddof=1) / np.sqrt(5)
    assert g == pytest.approx(np.exp(x.mean())) and lo == pytest.approx(np.exp(x.mean() - h))
    assert hi == pytest.approx(np.exp(x.mean() + h))
    assert np.isnan(M.t_interval([0.1])[1])


# ---------------------------------------------------------------------------
# classes

def test_class_inside_above_below_and_six_of_eight():
    assert M.class_ab([0.5, 1.0, 2.0], 0.9, 1.1) == 'inside'
    assert M.class_ab([1.2, 1.5], 0.9, 1.1) == 'above'
    assert M.class_ab([0.5, 0.8], 0.9, 1.1) == 'below'
    assert M.class_ab([0.5, 1.5], 0.9, 1.1) == 'mixed'
    assert M.class_ab([1.2], 1.3, 1.6) == 'below'                       # an interval above 1: 'below', not E98's 'flatter'
    with pytest.raises(ValueError):
        M.class_ab([0.0], 0.9, 1.1)
    assert M.comparison_class(['above'] * 6 + ['inside'] * 2) == 'above'
    assert M.comparison_class(['above'] * 5 + ['inside'] * 3) == 'not decided'


def test_held_and_reading_dependent():
    cats = {'carrying': 'above', 'first_in_block': 'above', 'without_marked': 'above', 't_interval': 'above',
            'settled': 'above'}
    partners = [k for k in cats if k != 'carrying']
    assert M.held_class(cats, partners) == 'above'
    assert M.held_class(dict(cats, t_interval='inside'), partners) == 'reading-dependent'
    assert M.held_class(dict(cats, carrying='not decided', first_in_block='not decided', without_marked='not decided',
                             t_interval='not decided', settled='not decided'), partners) == 'reading-dependent'
    assert M.x_over_estimates({'a': 'above', 'b': 'reading-dependent', 'c': 'inside'}) == 'inside'
    assert M.x_over_estimates({'a': 'above', 'b': 'above'}) == 'above'
    assert M.x_over_estimates({'a': 'above', 'b': 'reading-dependent'}) == 'not decided'


def _x1_rows(spec):
    rows = []
    for rid in M.X1_RATIOS:
        r12, series, differ = spec.get(rid, (1.0, 1.0, False))
        rows.append(dict(ratio_id=rid, differ=differ, **M.x1_direction(rid, series, r12)))
    return pd.DataFrame(rows)


def test_x1_three_sides():
    assert M.x1_side(_x1_rows({'P1': (0.8, 1.5, True)})) == "the 12's side"           # one differs
    toward = {'P1': (0.8, 2.0, True), 'P2': (0.9, 2.5, True), 'Q1': (1.1, 2.4, True), 'F1': (0.6, 0.2, True)}
    assert M.x1_side(_x1_rows(toward)) == "the 2022 flies' side"
    # four differ, but one lies nearer to the 12 than to set S: between
    between = dict(toward, Q1=(1.1, 1.6, True))
    assert M.x1_side(_x1_rows(between)) == 'between'
    # four differ, one of them away from set S: between
    assert M.x1_side(_x1_rows(dict(toward, F1=(0.6, 1.2, True)))) == 'between'
    d = M.x1_direction('P2', 2.5, 0.9)                                   # inside set S's range for P2
    assert d['toward_s'] and d['nearer_s'] and d['log_dist_s'] == 0.0


# ---------------------------------------------------------------------------
# Y1, M, Z

def test_y1_variance_split_by_hand_and_on_a_known_persistent_sd():
    rows = pd.DataFrame(dict(day=['d1', 'd1', 'd2', 'd2'], recording=['a', 'b', 'c', 'd'], odor='pa', conc=-4.0,
                             pair='2/3', log2_ratio=[1.0, 1.0, -1.0, -1.0]))
    y = M.y1_variance_split(rows)
    assert y['sd_within_groups'] == 0.0 and y['persistent_sd'] == pytest.approx(np.sqrt(2.0))
    rng = np.random.default_rng(0)
    rows = []
    for g in range(600):
        eff = rng.normal(0, 0.3)
        for k in range(4):
            rows.append(dict(day=f'd{g}', recording=f'r{g}_{k}', odor='pa', conc=-4.0, pair='2/3',
                             log2_ratio=eff + rng.normal(0, 0.2)))
    y = M.y1_variance_split(pd.DataFrame(rows))
    assert y['persistent_sd'] == pytest.approx(0.3, abs=0.03)
    assert y['sd_single'] == pytest.approx(np.sqrt(0.3 ** 2 + 0.2 ** 2), abs=0.03)
    # no persistent part: floored at 0 or near it
    flat = pd.DataFrame([dict(r, log2_ratio=rng.normal(0, 0.2)) for r in rows])
    assert M.y1_variance_split(flat)['persistent_sd'] < 0.05
    assert M.vial_yardstick(0.2) == pytest.approx(0.56)


def _smooth(rng, shape, s=3.0):
    from scipy.ndimage import gaussian_filter
    x = gaussian_filter(rng.normal(size=shape), s)
    return x / x.std()


def _m_maps(odor_specific: bool, seed: int = 1):
    rng = np.random.default_rng(seed)
    shape = (96, 96)
    common = _smooth(rng, shape, 6.0)
    odor_pat = {o: _smooth(rng, shape) for o in ('pa', 'ms', 'paa')}
    pres, keys = {}, []
    for o in ('ms', 'paa', 'pa'):
        for c in range(3):
            keys.append((o, c))
    stims = keys + [('pa_vial2', c) for c in range(3)] + [('pa_vial3', c) for c in range(3)]
    for k in stims:
        base = k[0].split('_vial')[0]
        pat = odor_pat[base] if (odor_specific or '_vial' not in k[0]) else _smooth(rng, shape)
        amp = 1.0 + 0.5 * k[1]
        pres[k] = [3.0 * amp * common + pat + 0.3 * rng.normal(size=shape) for _ in range(4)]
    vial_pairs = [(('pa_vial2', c), ('pa_vial3', c)) for c in range(3)]
    odor_pairs = [(('pa_vial2', c), ('ms', c)) for c in range(3)] + [(('pa_vial3', c), ('paa', c)) for c in range(3)]
    return pres, keys, vial_pairs, odor_pairs, np.ones(shape, bool)


def test_m_share_on_maps_built_to_be_the_odors_and_the_blocks():
    pres, common, vp, op, mask = _m_maps(odor_specific=True)
    r = M.m_calyx(pres, common, vp, op, mask, blur=2.0)
    # removing one common map from three odors' maps leaves the odors' remainders anticorrelated; the share holds
    assert r['r_rep'] > 0.8 and r['r_vial'] > 0.8 and r['r_vial'] > r['r_odor'] and r['share'] > 0.8
    pres, common, vp, op, mask = _m_maps(odor_specific=False)
    r = M.m_calyx(pres, common, vp, op, mask, blur=2.0)
    assert r['r_rep'] > 0.8 and abs(r['share']) < 0.25
    assert M.m_share(0.6, 0.1, 0.6) == pytest.approx(1.0) and M.m_share(0.1, 0.1, 0.6) == pytest.approx(0.0)
    assert M.m_class(0.6, 0.9, 7) == "the odor's" and M.m_class(0.6, 0.9, 6) == 'not decided'
    assert M.m_class(-0.1, 0.2, 2) == "the block's" and M.m_class(0.2, 0.6, 3) == 'not decided'


def test_z_sign_rule():
    up = [0.2] * 7 + [-0.1]
    assert M.z_rule(up, 1.1, 1.6) == 'replicates'
    assert M.z_rule(up, 0.95, 1.6) == 'not decided'                      # the interval holds 1
    assert M.z_rule([0.2] * 4 + [-0.2] * 4, 0.8, 1.2) == 'does not'
    assert M.z_rule([0.2] * 5 + [-0.2] * 3, 1.05, 1.3) == 'not decided'


# ---------------------------------------------------------------------------
# refusals and the lock

def test_refuses_a_raw_only_row():
    with pytest.raises(PermissionError):
        M.check_row(dict(stem='20211022_1FLfun1', use_until_further_notice='raw only: not used until further notice'))
    M.check_row(dict(stem='20211022_1FLfun3', use_until_further_notice="processed: usable (Yang's processed products)"))


@pytest.mark.skipif(not _shared(), reason='shared directory not available')
def test_refuses_recordings_outside_the_lists():
    for st in ('20211022_1FLfun1', '20211109_5FRfun2', '20221031_1FLfun1'):     # raw only, raw only, another set
        with pytest.raises(PermissionError):
            M.series_rows([st])


def test_refuses_any_map_of_the_series_before_part_a(tmp_path, monkeypatch):
    res, lock = tmp_path / 'results', tmp_path / 'lock.json'
    res.mkdir()
    assert M.is_series('20211022_1FLfun3') and not M.is_series('20210714_2FRfun1')
    with pytest.raises(PermissionError):
        M.require_part_a('20211022_1FLfun3', results=res, lock=lock, check_git=False)
    M.require_part_a('20210714_2FRfun1', results=res, lock=lock, check_git=False)   # not the series: no guard
    tables = {}
    for name in M.PART_A_TABLES:
        (res / name).write_text(f'# caveat\n{name}\n')
        tables[name] = M._sha256(res / name)
    lock.write_text(M.json.dumps(dict(tables=tables)))
    assert M.part_a_locked(res, lock, check_git=False)
    M.require_part_a('20211022_1FLfun3', results=res, lock=lock, check_git=False)
    assert not M.part_a_locked(res, lock, check_git=True)                         # a lock outside git is not committed
    (res / M.PART_A_TABLES[0]).write_text('changed\n')
    with pytest.raises(PermissionError):
        M.require_part_a('20211022_1FLfun3', results=res, lock=lock, check_git=False)
    # the map reader itself refuses before touching any file
    monkeypatch.setattr(M, 'part_a_locked', lambda **kw: False)
    with pytest.raises(PermissionError):
        M.frame_means(dict(stem='20211022_1FLfun3', folder='/nonexistent',
                           use_until_further_notice="processed: usable"))
    with pytest.raises(PermissionError):
        M.part_b()


def test_tables_carry_the_caveat_lines(tmp_path):
    df = pd.DataFrame(dict(a=[1.0, 2.5], b=['x', 'y']))
    p = tmp_path / 't.csv'
    M.write_table(df, p)
    lines = p.read_text().splitlines()
    assert lines[0] == f'# {M.QUALITY_CAVEAT}' and lines[1].startswith('# unmeasured')
    pd.testing.assert_frame_equal(M.read_table(p), df)


# ---------------------------------------------------------------------------
# with the shared directory: the series' stimulus files against the entry's lists (yaml only, no map)

@pytest.mark.skipif(not _shared(), reason='shared directory not available')
def test_eligibility_from_the_stimulus_files_equals_the_lists():
    el = M.eligibility_table()
    s = el[el['stem'].isin(M.SERIES)]
    assert len(s) == 30 and s['matches'].astype(bool).all()
    assert set(el.loc[~el['stem'].isin(M.SERIES), 'raw_only'].str[:8]) == {'raw only'}
    lost = s.set_index('stem')['not_imaged_after_blocks_out']
    assert lost['20211022_1FLfun2'] == 'pa_vial2 -4.863' and lost['20211022_1FLfun5'] == 'paa_vial2 -3.908'
    vi = M.vials_identity(el)
    assert vi['agrees'].all()


# ---------------------------------------------------------------------------
# parts B to E on made-up rows and maps

def _resp_rows(rec, fly, plane, vals, baseline='order_matched', blocks=None, n_reps=2):
    rows = []
    for (o, c), v in vals.items():
        rows.append(dict(recording=rec, fly=fly, side='L', pulse_s=0.5, plane_um=float(plane), role='series', odor=o,
                         conc=float(c), vial=2 if '_vial2' in o else 3 if '_vial3' in o else 1, stimulus=f'{o} {c:g}',
                         block=(blocks or {}).get((o, c), 1) - 1, block1=(blocks or {}).get((o, c), 1), baseline=baseline,
                         blank_subtracted=False, blank=np.nan, response=v, rep1=v, rep2=v, n_reps=n_reps, source='maps'))
    return rows


def test_every_presentation_weights_blocks_by_presentations():
    rows = _resp_rows('r1', 'f1', 18, {('pa', -5.34): 0.2}, blocks={('pa', -5.34): 2}, n_reps=4)
    rows += _resp_rows('r1', 'f1', 18, {('pa', -5.34): 0.5}, blocks={('pa', -5.34): 18}, n_reps=2)
    out = M.every_presentation(pd.DataFrame(rows))
    assert len(out) == 1 and out['response'].iloc[0] == pytest.approx((0.2 * 4 + 0.5 * 2) / 6)
    assert out['n_blocks'].iloc[0] == 2 and out['block1'].iloc[0] == 2


def test_three_vial_odor_and_their_average():
    assert M.three_vial_odor(['ms', 'paa', 'pa', 'pa_vial2', 'pa_vial3']) == 'pa'
    assert M.three_vial_odor(['ms', 'paa', 'pa']) is None
    vals = {('pa', -5.34): 0.3, ('pa_vial2', -5.34): 0.6, ('pa_vial3', -5.34): 0.9, ('ms', -2.0): 0.4}
    out = M.average_three_vials(pd.DataFrame(_resp_rows('r1', 'f1', 18, vals))).set_index('odor')
    assert out.loc['pa', 'response'] == pytest.approx(0.6) and out.loc['pa_vial2', 'response'] == pytest.approx(0.6)
    assert out.loc['ms', 'response'] == pytest.approx(0.4)


def test_reading_recordings_planes_marked_flies_and_the_solvent_block():
    rows = _resp_rows('20211107_1FLfun2', '2021-11-07_fly1', 18, {('ms', -2.0): 0.3})           # marked
    rows += _resp_rows('20211107_1FLfun3', '2021-11-07_fly1', 18, {('ms', -2.0): 0.3})
    rows += _resp_rows('20211109_5FRfun1', '2021-11-09_fly5', 18, {('ms', -2.0): 0.3})
    rows += _resp_rows('20211109_5FRfun3', '2021-11-09_fly5', 15, {('ms', -2.0): 0.3})
    rows += _resp_rows('20211029_3FLfun2', '2021-10-29_fly3', 18, {('ms', -2.0): 0.3}, baseline='first_in_block')
    r = pd.DataFrame(rows)
    assert M.reading_recordings(r, 'carrying', 'order_matched') == ['20211107_1FLfun2', '20211107_1FLfun3', '20211109_5FRfun1']
    assert M.reading_recordings(r, 'without_marked', 'order_matched') == ['20211107_1FLfun3', '20211109_5FRfun1']
    assert M.reading_recordings(r, 'without 2021-11-07_fly1', 'order_matched') == ['20211109_5FRfun1']
    assert M.reading_recordings(r, 'first_in_block', 'first_in_block') == ['20211029_3FLfun2']   # no solvent block
    assert M.reading_recordings(r, 'carrying', 'order_matched', plane=None)[-1] == '20211109_5FRfun3'


def test_y1_pairs_and_the_ten22_paa_runs():
    vals = {('paa', -4.386): 0.2, ('paa_vial2', -4.386): 0.4, ('paa_vial3', -4.386): 0.01, ('ms', -2.0): 0.5}
    blocks = {('paa', -4.386): 13, ('paa_vial2', -4.386): 1, ('paa_vial3', -4.386): 4, ('ms', -2.0): 12}
    rows = _resp_rows('20211022_1FLfun4', '2021-10-22_fly1', 15, vals, blocks=blocks)                # vial 1 the other half
    rows += _resp_rows('20211107_1FRfun3', '2021-11-07_fly1', 15, {**vals, ('paa_vial3', -4.386): 0.8}, blocks=blocks)
    t = M.y1_rows(pd.DataFrame(rows), 'order_matched').set_index(['recording', 'pair'])
    a = t.loc[('20211022_1FLfun4', '1/2')]
    assert not a['carries'] and a['log2_ratio'] == pytest.approx(-1.0) and a['block_distance'] == 12
    assert t.loc[('20211022_1FLfun4', '2/3'), 'carries'] and not t.loc[('20211022_1FLfun4', '2/3'), 'entered']   # the floor
    b = t.loc[('20211107_1FRfun3', '2/3')]
    assert b['carries'] and b['log2_ratio'] == pytest.approx(-1.0) and b['block_distance'] == 3
    assert t.loc[('20211107_1FRfun3', '1/3'), 'carries']


def _bm(odor_specific=True, seed=2, far=False):
    pres, keys, vp, op, mask = _m_maps(odor_specific, seed)
    slot = {'ms': 1, 'paa': 2, 'pa_vial2': 3, 'pa_vial3': 4, 'pa': 5}      # interleaved, every pair within four blocks
    blocks = {k: slot[k[0]] + 5 * k[1] for k in pres}
    if far:                                                          # vial 3 far from everything
        for c in range(3):
            blocks[('pa_vial3', c)] = 40 + c
    conc = {0: -5.34, 1: -4.863, 2: -4.386}
    return {(k[0], conc[k[1]]): dict(block1=blocks[k], maps=pres[k]) for k in pres}, mask


def test_m_recording_pairs_filters_and_class():
    bm, mask = _bm(True)
    stim, pairs = M.m_recording(bm, mask, blur=2.0)
    assert set(pairs['kind']) == {'vial', 'odor'} and (pairs['kind'] == 'vial').sum() == 3
    assert (pairs['kind'] == 'odor').sum() == 3 * 2 * 2                    # 3 ranks x 2 vials x 2 other odors
    assert stim['r_rep'].min() > 0.8
    pairs = pairs.assign(recording='r', fly='f1', side='L')
    cal = M.m_calyces(pairs)
    assert cal['share'].iloc[0] > 0.8 and cal['vial_over_odor'].iloc[0]
    # pairs farther apart than four blocks do not enter
    bmf, _ = _bm(True, far=True)
    _, pf = M.m_recording(bmf, mask, blur=2.0)
    cf = M.m_calyces(pf.assign(recording='r', fly='f1', side='L'))
    assert cf['n_vial_pairs'].iloc[0] == 0 and np.isnan(cf['share'].iloc[0])
    # a stimulus with r_rep below 0.3 does not enter
    low = pairs.copy()
    low.loc[low['kind'] == 'vial', 'r_rep_a'] = 0.2
    assert M.m_calyces(low)['n_vial_pairs'].iloc[0] == 0
    # the class over calyces of eight flies: all shares near 1
    many = pd.concat([pairs.assign(fly=f'f{i}') for i in range(8)], ignore_index=True)
    res = M.m_result(M.m_calyces(many))
    assert res['n_calyces'] == 8 and res['cls'] == "the odor's"
    bmb, maskb = _bm(False)
    _, pb = M.m_recording(bmb, maskb, blur=2.0)
    resb = M.m_result(M.m_calyces(pd.concat([pb.assign(recording='r', fly=f'f{i}', side='L') for i in range(8)])))
    assert resb['share'] < 0.25


def test_alternating_calyces_and_drift():
    runs = pd.DataFrame(dict(recording=[f'r{i}' for i in range(9)], fly=['a'] * 4 + ['b'] * 3 + ['c'] * 2,
                             side='L', plane_um=[15, 18, 15, 18, 15, 18, 18, 18, 15], run=[1, 2, 3, 4, 1, 2, 3, 1, 3]))
    assert M.alternating_calyces(runs) == ['a_L']
    rep = pd.DataFrame(dict(block_distance=[16, 16, 8], log_later_over_first=[-0.32, -0.32, np.nan]))
    assert M.drift_per_block(rep) == pytest.approx(-0.02)


def test_widen_and_x_classes_on_made_up_wirings():
    w = M.widen({'carrying': {'P1': (1.0, 2.0)}}, 1.0)
    assert w['carrying']['P1'] == pytest.approx((0.5, 4.0))
    wir = [f'w{i}' for i in range(8)]
    on = pd.DataFrame([dict(wiring=x, stimulus=s, win=v, pre=1.0) for x in wir
                       for s, v in (('paa', 3.0), ('ms', 2.0), ('ms2', 1.5))])     # paa over ms 2, over ms2 4 (linear)
    ivs = {'carrying': {'P1': (1.5, 2.5), 'P2': (1.8, 2.2)}, 'first_in_block': {'P1': (1.5, 2.5), 'P2': (1.8, 2.2)},
           't_interval': {'P1': (0.5, 1.0), 'P2': (0.6, 1.1)}}
    readings = [('carrying', 'carrying', 'blank'), ('first_in_block', 'first_in_block', 'blank'),
                ('settled', 'carrying', 'settled')]
    cl, wcl, s = M.x_classes(ivs, readings, levels={'blank': on, 'settled': on}, num='paa',
                             estimates={'e1': 'ms', 'e2': 'ms2'})
    assert s['X3']['held_per_estimate'] == {'e1': 'inside', 'e2': 'above'} and s['X3']['over_estimates'] == 'inside'
    cl2, _, s2 = M.x_classes(ivs, readings + [('t_interval', 't_interval', 'blank')], levels={'blank': on, 'settled': on},
                             num='paa', estimates={'e1': 'ms', 'e2': 'ms2'})
    assert s2['X3']['held_per_estimate']['e1'] == 'reading-dependent'      # the t interval puts it above
    assert s2['X3']['over_estimates'] == 'not decided'


def test_x1_rows_on_made_up_sets():
    flies = [f'2021-11-0{i}_fly1' for i in range(1, 6)]
    pr = pd.DataFrame([dict(ratio_id=rid, fly=f, recording=f, entered=True, log_ratio=np.log(2.0) + 0.01 * k)
                       for rid in M.RATIO_IDS for k, f in enumerate(flies)])
    fs = pd.DataFrame([dict(ratio_id=rid, numerator='a', denominator='b', ratio=2.0, lo90=1.9, hi90=2.1, sd_log2_flies=0.01,
                            n_flies=5, n_recordings=5, num_how='imaged', den_how='imaged') for rid in M.RATIO_IDS])
    r12 = pd.DataFrame([dict(ratio_id=rid, ratio=1.0, lo90=0.9, hi90=1.1) for rid in M.RATIO_IDS]).set_index('ratio_id')
    pf12 = {rid: pd.Series(np.log(1.0) + 0.02 * np.arange(8), index=[f'g{i}' for i in range(8)]) for rid in M.RATIO_IDS}
    x = M.x1_rows(fs, pr, r12, pf12)
    assert x['differ'].all() and (x['n_splits'] == 1287).all()
    # every one differs; toward set S for all but F1 (set S lower) and S1/S2 (not in X1): between
    assert M.x1_side(x) == 'between'
