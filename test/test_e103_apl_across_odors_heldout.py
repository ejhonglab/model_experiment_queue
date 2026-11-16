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
