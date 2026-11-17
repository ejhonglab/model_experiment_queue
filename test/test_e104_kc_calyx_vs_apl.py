#!/usr/bin/env python3
"""Tests for E104 (the flies' KCs beside their APL along concentration; experiments/E104_kc_calyx_vs_apl.py), on
made-up rows and tables: the reader refuses raw-only rows, recordings outside the lists and set K without part A's
lock; the lock's sha256 round trip and its refusals (another sha, other definitions, not committed); the committed lock
matches the module's definitions; the groups of set K's rows; the entry's counts on a made-up eligibility table; the
marked-run rules; the ladder's parameters at m = 0, 1 and between; the readings of a made-up recording; and, with the
shared directory, the eligibility table's counts and lists and the reader against E100's saved set S values."""

import json

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
M = pytest.importorskip('experiments.E104_kc_calyx_vs_apl')


def _shared() -> bool:
    try:
        return (H.shared_dir() / 'data' / 'yang_apl_recordings' / 'recordings.csv').exists()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# the reader's refusals

def test_check_row_refuses_raw_only():
    M.check_row({'stem': 'x', 'use_until_further_notice': "processed: usable (Yang's processed products only)"})
    with pytest.raises(PermissionError):
        M.check_row({'stem': 'x', 'use_until_further_notice': 'raw only: not used until further notice'})
    with pytest.raises(PermissionError):
        M.check_row({'stem': 'x'})


def _recs(rows):
    return pd.DataFrame([dict(stem=s, date=int(s[:8]), fly=1, hemisphere='L', set='K', folder='/nowhere',
                              use_until_further_notice=u) for s, u in rows])


def test_recording_rows_refuses_outside_and_raw_only():
    recs = _recs([('20221019_1FLfun2', 'processed: usable'), ('20221019_1FLfun3', 'raw only: not used'),
                  ('20221031_1FLfun1', 'processed: usable')])
    allowed = {'20221019_1FLfun2': 'K6F', '20221019_1FLfun3': 'K6F'}
    out = M.recording_rows(['20221019_1FLfun2'], recs=recs, allowed=allowed)
    assert out.loc['20221019_1FLfun2', 'role'] == 'K6F' and out.loc['20221019_1FLfun2', 'fly_id'] == '2022-10-19_fly1'
    with pytest.raises(PermissionError):                  # outside the lists
        M.recording_rows(['20221031_1FLfun1'], recs=recs, allowed=allowed)
    with pytest.raises(PermissionError):                  # raw only, even when listed
        M.recording_rows(['20221019_1FLfun3'], recs=recs, allowed=allowed)


def test_allowed_stems_are_the_lists():
    a = M.allowed_stems()
    assert len(M.K6F) == 55 and len(M.K1020) == 10 and len(M.K7F) == 22 and len(M.S_STEMS) == 29
    assert set(a) == set(M.K6F) | set(M.K1020) | set(M.K7F) | set(M.S_STEMS)
    assert not set(M.K7F_LEFT_OUT_STEMS) & set(a)          # the left calyx of 2022-01-25 fly 1 is refused
    assert all(M.is_set_k(s) for s in M.K6F + M.K1020 + M.K7F + M.K7F_LEFT_OUT_STEMS)
    assert not any(M.is_set_k(s) for s in M.S_STEMS)


def test_require_part_a_refuses_set_k_without_lock(tmp_path):
    lock = tmp_path / 'lock.json'
    with pytest.raises(PermissionError):
        M.require_part_a('20221019_1FLfun2', lock=lock, check_git=False)
    with pytest.raises(PermissionError):
        M.require_part_a('20220125_1FRfun1', lock=lock, check_git=False)      # SytGCaMP7f flies are set K too
    M.require_part_a('20221031_1FLfun1', lock=lock, check_git=False)           # set S needs no lock
    M.write_lock(lock, eligibility=tmp_path / 'none.csv')
    M.require_part_a('20221019_1FLfun2', lock=lock, check_git=False)           # with a valid lock it passes


def test_read_table_k_refuses(monkeypatch, tmp_path):
    row = dict(stem='20221019_1FLfun2', use_until_further_notice='processed: usable', folder=str(tmp_path))
    monkeypatch.setattr(M, 'PART_A_LOCK', tmp_path / 'absent.json')
    with pytest.raises(PermissionError):
        M.read_table_k(row)                                                     # set K without the lock
    monkeypatch.setattr(M, 'part_a_locked', lambda **kw: True)
    with pytest.raises(PermissionError):
        M.read_table_k(dict(row, stem='20221019_1FLfun3'))                      # raw only, not in the lists
    with pytest.raises(PermissionError):
        M.read_table_k(dict(row, use_until_further_notice='raw only'))          # listed, but the row says raw only
    with pytest.raises(PermissionError):
        M.read_table_k(dict(row, stem='20220125_1FLfun1'))                      # the left-out calyx


# ---------------------------------------------------------------------------
# the lock

def test_lock_round_trip_and_refusals(tmp_path):
    lock = tmp_path / 'lock.json'
    d = M.write_lock(lock, eligibility=tmp_path / 'none.csv')
    assert d['sha256'] == M.definitions_sha(M.part_a_definitions())
    assert M.part_a_locked(lock=lock, check_git=False)
    assert not M.part_a_locked(lock=lock, check_git=True)          # not committed (tmp_path is not in a repository)
    j = json.loads(lock.read_text())
    bad = dict(j, sha256='0' * 64)
    (tmp_path / 'bad_sha.json').write_text(json.dumps(bad))
    assert not M.part_a_locked(lock=tmp_path / 'bad_sha.json', check_git=False)
    other = json.loads(json.dumps(j['definitions']))
    other['classes']['apl_steeper'] = 0.9                          # another threshold, with its own correct sha
    (tmp_path / 'other.json').write_text(json.dumps(dict(j, definitions=other, sha256=M.definitions_sha(other))))
    assert not M.part_a_locked(lock=tmp_path / 'other.json', check_git=False)
    with pytest.raises(RuntimeError):                              # a lock with other definitions is not overwritten
        M.write_lock(tmp_path / 'other.json', eligibility=tmp_path / 'none.csv')


def test_definitions_sha_stable_and_sensitive(monkeypatch):
    a = M.definitions_sha(M.part_a_definitions())
    assert a == M.definitions_sha(M.part_a_definitions()) and len(a) == 64
    monkeypatch.setattr(M, 'Q_APL_STEEPER', 0.79)
    assert M.definitions_sha(M.part_a_definitions()) != a
    monkeypatch.undo()
    monkeypatch.setattr(M, 'MARKED', dict(M.MARKED, x='dim'))
    assert M.definitions_sha(M.part_a_definitions()) != a


def test_committed_lock_matches_the_definitions():
    if not M.PART_A_LOCK.exists():
        pytest.skip('the lock is not written yet')
    j = json.loads(M.PART_A_LOCK.read_text())
    assert j['sha256'] == M.definitions_sha(j['definitions']) == M.definitions_sha(M.part_a_definitions())
    assert M.part_a_locked(check_git=False)


# ---------------------------------------------------------------------------
# the eligibility table

def test_k_group():
    p = 'processed: usable'
    assert M.k_group(dict(date=20221019, fly=1, hemisphere='L', use_until_further_notice=p)) == 'K6F'
    assert M.k_group(dict(date=20221020, fly=1, hemisphere='R', use_until_further_notice=p)) == 'K1020'
    assert M.k_group(dict(date=20220125, fly=1, hemisphere='L', use_until_further_notice=p)) == 'K7F_left_out'
    assert M.k_group(dict(date=20220125, fly=1, hemisphere='R', use_until_further_notice=p)) == 'K7F'
    assert M.k_group(dict(date=20221019, fly=1, hemisphere='L', use_until_further_notice='raw only')) == 'raw_only'
    assert M.k_group(dict(date=20230101, fly=1, hemisphere='L', use_until_further_notice=p)) == 'other'


def _made_up_eligibility():
    rows = []

    def add(group, flies_calyces, n, hold=None, series=None, planes=(15.0,)):
        for i in range(n):
            f, cx = flies_calyces[i % len(flies_calyces)]
            r = dict(stem=f'{group}_{i}', group=group, fly=f, calyx=f'{f}_{cx}', plane_um=planes[i % len(planes)],
                     plane_source='recordings.csv', has_yang_table=True, n_solvent_only_blocks=2)
            for o in M.ODORS:
                r[f'holds_{o}'] = bool(hold is None or i < hold[o])
                r[f'series_{o}'] = bool(series is None or i < series[o])
            rows.append(r)
    k6 = [(f'f{k}', s) for k in range(5) for s in ('L', 'R')][:9]          # 5 flies, 9 calyces
    add('K6F', k6, 55, hold=dict(paa=45, pa=46, ms=37), series=dict(paa=45, pa=36, ms=37),
        planes=(12.0, 15.0, 18.0, 21.0))
    add('K1020', [('g', 'L')], 10)
    add('K7F', [('a', 'L'), ('b', 'L'), ('c', 'L'), ('d', 'R')], 22)
    add('K7F_left_out', [('a', 'R')], 5)
    add('S', [('s1', 'L'), ('s2', 'R'), ('s3', 'L'), ('s3', 'R')], 29)     # 3 flies, 4 calyces
    return pd.DataFrame(rows)


def test_eligibility_counts_on_a_made_up_table():
    el = _made_up_eligibility()
    c = M.eligibility_counts(el, dict(n_recordings=20, n_flies=7))
    assert c['passed'], {g: r for g, r in c['groups'].items() if not r['equal']}
    assert c['beside']['K6F_n_series'] == dict(paa=45, pa=36, ms=37)
    assert len(c['beside']['K6F_holding_one_concentration_only']['pa']) == 10
    # one recording fewer, or one without Yang's table, or another count of pa: the identity fails
    c = M.eligibility_counts(el.drop(index=el.index[el['group'] == 'K6F'][-1]), dict(n_recordings=20, n_flies=7))
    assert not c['passed'] and not c['groups']['K6F']['equal'] and c['groups']['K1020']['equal']
    el2 = el.copy()
    el2.loc[el2.index[el2['group'] == 'K7F'][0], 'has_yang_table'] = False
    assert not M.eligibility_counts(el2, dict(n_recordings=20, n_flies=7))['groups']['K7F']['equal']
    el3 = el.copy()
    el3.loc[el3.index[(el3['group'] == 'K6F') & el3['holds_pa']][0], 'holds_pa'] = False
    assert not M.eligibility_counts(el3, dict(n_recordings=20, n_flies=7))['passed']
    assert not M.eligibility_counts(el, dict(n_recordings=19, n_flies=7))['groups']['L']['equal']
    assert not M.eligibility_counts(el)['passed']                               # set L's counts missing


def test_mark_reason():
    assert M.mark_reason('date note: ABORTED MIDWAY?') == 'aborted'
    assert M.mark_reason('date note: signigicant warpping') == 'significant warping'      # Yang's spelling
    assert M.mark_reason('date note: significant warpping, | readjust') == 'significant warping'
    assert M.mark_reason('date note: seq 3 | Getting dim, bubbles form') == 'dim'
    assert M.mark_reason('date note: warpped for about 6um') == ''
    assert M.mark_reason('dimension') == '' and M.mark_reason(np.nan) == ''


# ---------------------------------------------------------------------------
# the ladder and the readings

def test_ladder_params():
    p = dict(wAPLKC_scale=3.0, wAPLPN_scale=0.7, fixed_thr=12.5)
    assert M.ladder_params(p, 1.0) == p
    assert M.ladder_params(p, 0.0) == dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    q = M.ladder_params(p, 0.25)
    assert q['wAPLKC_scale'] == pytest.approx(0.75) and q['wAPLPN_scale'] == pytest.approx(0.175)
    assert q['fixed_thr'] == 12.5 and p['wAPLKC_scale'] == 3.0
    with pytest.raises(KeyError):
        M.ladder_params(dict(fixed_thr=1.0), 0.5)


def test_recording_readings_made_up():
    # two solvent-only blocks (first and last) around two odor blocks of a solvent and three presentations
    seq = [('solvent', np.nan)] * 4 + [('solvent', np.nan)] + [('pa', -5.0)] * 3 + [('solvent', np.nan)] \
        + [('pa', -4.0)] * 3 + [('solvent', np.nan)] * 4
    vals = [0.1, 0.2, 0.3, 0.4, 0.0, 1.0, 1.2, 1.4, 0.05, 2.0, 2.2, 2.4, 0.3, 0.4, 0.5, 0.6]
    fr = pd.DataFrame(dict(frame=range(16), block=[k // 4 for k in range(16)], position=[k % 4 + 1 for k in range(16)],
                           odor=[o for o, _ in seq], log10_conc=[c for _, c in seq], yang_peak=vals))
    fr['is_solvent'] = fr['odor'] == 'solvent'
    fr['stim'] = ['solvent' if s else f'{o} @ {c:g}' for s, o, c in zip(fr['is_solvent'], fr['odor'], fr['log10_conc'])]
    per, st = M.recording_readings(fr)
    st = st.set_index('stim')
    # order-matched: minus the mean of the two solvent-only blocks at the same position (0.3, 0.4, 0.5 at 2, 3, 4)
    assert st.loc['pa @ -5', 'yang_order_matched'] == pytest.approx(np.mean([1.0 - 0.3, 1.2 - 0.4, 1.4 - 0.5]))
    assert st.loc['pa @ -4', 'yang_minus_block'] == pytest.approx(np.mean([2.0, 2.2, 2.4]) - 0.05)
    assert st.loc['pa @ -5', 'yang_first'] == pytest.approx(1.0 - 0.3)
    assert st.loc['pa @ -4', 'n_presentations'] == 3 and st.loc['pa @ -4', 'first_block'] == 2


# ---------------------------------------------------------------------------
# with the shared directory (no response table of set K is opened)

@pytest.mark.skipif(not _shared(), reason='needs the shared directory')
def test_eligibility_on_the_copy():
    el = M.eligibility_table()
    c = M.eligibility_counts(el, M.set_l_counts())
    assert c['passed'], c
    assert M.eligibility_identity(el)['passed']


@pytest.mark.skipif(not _shared(), reason='needs the shared directory')
def test_reader_reproduces_e100_set_s():
    res, _ = M.s_identity()
    assert res['passed'] and res['max_abs_diff'] == 0.0 and res['n_recordings'] == 29


# ---------------------------------------------------------------------------
# part A on made-up sets

def test_fall_class():
    assert M.fall_class(True, 0.5, 0.4, 0.7) == 'APL steeper'
    assert M.fall_class(True, 1.5, 1.3, 1.8) == 'KCs steeper'
    assert M.fall_class(True, 0.9, 0.7, 1.1) == 'not decided'           # they differ but Q is near 1
    assert M.fall_class(False, 1.0, 0.7, 1.4) == 'one function'
    assert M.fall_class(False, 1.0, 0.6, 1.4) == 'not decided'          # Q's interval reaches outside 0.67 to 1.5
    assert M.fall_class(False, 0.5, 0.3, 0.8) == 'not decided'          # steep, but the sets do not differ
    assert M.fall_class(False, np.nan, np.nan, np.nan) == 'not decided'


def test_split_rule_on_known_sets():
    a = pd.Series(np.log([0.3, 0.32, 0.28, 0.31, 0.29, 0.3, 0.33, 0.27, 0.3, 0.31]), index=[f'a{i}' for i in range(10)])
    b = pd.Series(np.log([0.9, 0.85, 0.95, 0.88, 0.92]), index=[f'b{i}' for i in range(5)])
    r = M.split_rule(a, b)
    assert r['differ'] and r['n_splits'] == 3003 and r['share'] < 0.01
    same = M.split_rule(a.iloc[:5], pd.Series(a.iloc[5:].to_numpy(), index=[f'c{i}' for i in range(5)]))
    assert not same['differ']
    # the vectorized share is E98's `pulse_permutation`'s
    E98 = pytest.importorskip('experiments.E98_apl_response_strength')
    s = pd.concat([a.iloc[:4], b.iloc[:3]])
    assert M.split_rule(a.iloc[:4], b.iloc[:3])['share'] == pytest.approx(E98.pulse_permutation(s, list(b.index[:3]))['share'])


def test_q_of_fall_and_set_value():
    a = pd.Series(np.log([0.4, 0.5, 0.45, 0.55]), index=list('abcd'))
    k = pd.Series(np.log([0.8, 0.9, 0.85]), index=list('xyz'))
    q = M.q_of_fall(a, k, n_boot=2000)
    assert q['q'] == pytest.approx(np.exp(a.mean() - k.mean()))
    assert q['q_lo90'] <= q['q'] <= q['q_hi90'] and q == M.q_of_fall(a, k, n_boot=2000)     # seeded
    v = M.set_value(a, n_boot=2000)
    assert v['value'] == pytest.approx(np.exp(a.mean())) and v['n_flies'] == 4
    q2 = M.q_of_fall(pd.Series([0.4, 0.6]), pd.Series([0.8, 0.8]), arithmetic=True, n_boot=500)
    assert q2['q'] == pytest.approx(0.5 / 0.8)


def test_held_x1_x3():
    assert M.held_class({'carrying': 'APL steeper', 'p1': 'APL steeper'}) == 'APL steeper'
    assert M.held_class({'carrying': 'APL steeper', 'p1': 'not decided'}) == 'not decided'
    assert M.x1_verdict({'G1': 'APL steeper', 'G2': 'APL steeper', 'G3': 'not decided'}) == 'APL steeper'
    assert M.x1_verdict({'G1': 'APL steeper', 'G2': 'APL steeper', 'G3': 'KCs steeper'}) == 'not decided'
    assert M.x1_verdict({'G1': 'one function', 'G2': 'one function', 'G3': 'not decided'}) == 'one function'
    assert M.x1_verdict({'G1': 'one function', 'G2': 'one function', 'G3': 'APL steeper'}) == 'not decided'
    assert M.x1_verdict({'G1': 'APL steeper', 'G2': 'one function', 'G3': 'not decided'}) == 'not decided'
    assert M.x3_verdict({'paa': 4, 'pa': 5, 'ms': 2})['verdict'] == 'holds'
    assert M.x3_verdict({'paa': 3, 'pa': 5, 'ms': 2})['verdict'] == 'fails'


def _vals(set_, flies, falls_by_fly, odor='paa', concs=(-5.0, -4.0), extra=None):
    """Made-up per-recording values: per fly one recording with the odor's two concentrations, the higher at 1.0 and
    the lower at the fly's fall."""
    rows = []
    for f, fall in zip(flies, falls_by_fly):
        for c, v in zip(concs, (fall, 1.0)):
            rows.append(dict(set=set_, recording=f'{set_}_{f}', fly=f, calyx='', plane_um=18.0, plane_uncertain=False,
                             marked='', odor=odor, conc=c, stim=f'{odor} @ {c:g}', carrying=v, block_solvent=v,
                             first_presentation=v, raw=v + 1, **(extra or {})))
    return pd.DataFrame(rows)


def test_fall_reading_made_up():
    rng = np.random.default_rng(1)
    apl = _vals('S', [f's{i}' for i in range(10)], np.exp(np.log(0.4) + 0.05 * rng.standard_normal(10)))
    kc = _vals('K6F', [f'k{i}' for i in range(5)], np.exp(np.log(0.9) + 0.05 * rng.standard_normal(5)))
    vals = pd.concat([apl, kc], ignore_index=True)
    spec = M.reading_specs(['k0'])['carrying']
    r = M.fall_reading(vals, M.FALLS[0], dict(spec, apl_sets=('S',)))
    assert r['cls'] == 'APL steeper' and r['apl_n_flies'] == 10 and r['kc_n_flies'] == 5
    kc2 = _vals('K6F', [f'k{i}' for i in range(5)], np.exp(np.log(0.4) + 0.05 * rng.standard_normal(5)))
    r2 = M.fall_reading(pd.concat([apl, kc2], ignore_index=True), M.FALLS[0], dict(spec, apl_sets=('S',)))
    assert r2['cls'] == 'one function'
    # the floor: a recording whose lower response is below MIN_RESPONSE enters no geometric fall
    f = M.recording_fall(_vals('K6F', ['k'], [0.01]), ('paa', -5.0), ('paa', -4.0), 'carrying')
    assert np.isnan(f['log_fall']) and f['ratio'] == pytest.approx(0.01)
    g = M.recording_fall(_vals('K6F', ['k'], [0.5], concs=(-6.0, -4.0)), ('paa', -5.0), ('paa', -4.0), 'carrying')
    assert g['num_how'] == 'interpolated' and g['num'] == pytest.approx(0.75)


def test_reading_specs_and_select():
    specs = M.reading_specs(['f1', 'f2'])
    assert list(specs)[:6] == ['carrying', 'set_S_alone', 'block_solvent', 'first_presentation', 'plane_18',
                               'without_marked']
    assert 'kc_fly_left_out:f1' in specs and 'k1020_as_kc' in specs and len(specs) == 9
    v = pd.DataFrame(dict(set=['K6F'] * 4, fly=['f1', 'f1', 'f2', 'f2'], plane_um=[18.0, 15.0, 18.0, 18.0],
                          plane_uncertain=[False, False, True, False], marked=['', '', '', 'dim']))
    assert len(M.select(v, ('K6F',), specs['plane_18'], 'kc')) == 2
    assert len(M.select(v, ('K6F',), specs['without_marked'], 'kc')) == 3
    assert len(M.select(v, ('K6F',), specs['kc_fly_left_out:f1'], 'kc')) == 2
    assert len(M.select(v, ('K6F',), specs['kc_fly_left_out:f1'], 'apl')) == 4      # leaves APL sets alone


# ---------------------------------------------------------------------------
# part B on made-up curves

def test_isotonic_and_interp_monotone():
    assert np.allclose(M.isotonic([1, 3, 2, 4]), [1, 2.5, 2.5, 4])
    assert np.all(np.diff(M.isotonic([5, 1, 4, 2, 8, 3])) >= 0)
    assert np.all(np.diff(M.isotonic([1, 3, 2, 0], increasing=False)) <= 0)
    x = [0, 1, 2, 3]
    assert M.interp_monotone(x, [0, 2, 1, 3], 1.5) == pytest.approx(1.5)
    assert np.isnan(M.interp_monotone(x, [0, 1, 2, 3], 4.0))
    vals = [M.interp_monotone(x, [0, 2, 1, 3], t) for t in np.linspace(0, 3, 31)]
    assert np.all(np.diff(vals) >= -1e-12)


def test_d_star_on_a_made_up_curve():
    d = [0.0, 0.5, 1.0, 1.5, 2.0]
    ratio = [1.0, 1.5, 2.0, 2.5, 3.0]
    r = M.d_star(d, ratio, 2.0, 1.5, 2.5)
    assert r['d_star'] == pytest.approx(1.0) and r['d_star_lo'] == pytest.approx(0.5) and r['d_star_hi'] == pytest.approx(1.5)
    r = M.d_star(d, ratio, 2.75, 2.5, 3.5)                    # the interval's top beyond the ladder
    assert r['d_star'] == pytest.approx(1.75) and r['d_star_hi'] == np.inf and r['censored']
    r = M.d_star(d[::-1], ratio[::-1], 2.0)                    # order does not matter
    assert r['d_star'] == pytest.approx(1.0)
    r = M.d_star(d, [1.0, 1.6, 1.4, 2.5, 3.0], 1.5)           # a non-monotone point is pooled first
    assert 0.5 <= r['d_star'] <= 1.5


def test_y_classes():
    assert M.y_class((0.2, 0.4), (0.5, 0.8)) == 'too weak'
    assert M.y_class((0.9, 1.0), (0.5, 0.8)) == 'too strong'
    assert M.y_class((0.4, 0.6), (0.5, 0.8)) == 'inside'
    assert M.y_class((0.7, 0.7), (0.5, np.inf)) == 'inside'
    assert M.y_class((0.3, 0.3), (0.5, np.inf)) == 'too weak'
    assert M.y_over_wirings(['too weak'] * 6 + ['inside'] * 2) == 'too weak'
    assert M.y_over_wirings(['too weak'] * 5 + ['inside'] * 3) == 'not decided'
    assert M.y_over_estimates({'a': 'too weak', 'b': 'inside'}) == 'inside'
    assert M.y_over_estimates({'a': 'too weak', 'b': 'too weak'}) == 'too weak'
    assert M.y_over_estimates({'a': 'too weak', 'b': 'not decided'}) == 'not decided'
    assert M.y_verdict({'2h': 'too weak', 'ma': 'too weak', 'farn': 'inside'}) == 'too weak'
    assert M.y_verdict({'2h': 'too weak', 'ma': 'too weak', 'farn': 'too strong'}) == 'not decided'
    assert M.y_verdict({'2h': 'inside', 'ma': 'inside', 'farn': 'too weak'}) == 'inside'
    # Y's held rule on made-up readings
    assert M.held_class({'carrying': 'too weak', 'one_spike': 'too weak', 'settled': 'inside'}) == 'not decided'


def test_flies_rate_ratio_made_up():
    rows = []
    for k, (g, fly, r) in enumerate([('control', 'c1', 0.05), ('control', 'c1', 0.07), ('control', 'c2', 0.06),
                                     ('control', 'c3', 0.06), ('apl_output_blocked', 't1', 0.15),
                                     ('apl_output_blocked', 't1', 0.17), ('apl_output_blocked', 't2', 0.19)]):
        rows.append(dict(recording=f'r{k}', genotype_group=g, fly=fly, odor='2h @ -7', rate=r))
        rows.append(dict(recording=f'r{k}', genotype_group=g, fly=fly, odor='pfo @ 0', rate=0.01))
    rates = pd.DataFrame(rows)
    r = M.flies_rate_ratio(rates, '2h @ -7', n_boot=2000)
    assert r['ratio'] == pytest.approx(0.17 / 0.06) and r['lo90'] < r['ratio'] < r['hi90']
    assert r['n_ctrl_flies'] == 3 and r['n_tnt_recordings'] == 3
    s = M.flies_rate_ratio(rates, '2h @ -7', solvent=True, n_boot=2000)
    assert s['ratio'] == pytest.approx(0.16 / 0.05)


def test_add_d_made_up():
    rows = []
    for m, paa, w2h, blank, rate in ((0.0, 10.0, 9.0, 8.0, 0.2), (0.5, 12.0, 9.5, 8.0, 0.12), (1.0, 14.0, 10.0, 8.0, 0.08)):
        for s, win, rt in (('paa', paa, 0.1), ('2h [x/diag] @ -6', w2h, rate), ('blank @ 0', blank, 0.0)):
            rows.append(dict(wiring='w', m=m, source='', stimulus=s, win=win, pre=7.0, pre_settled=7.5, rate=rt,
                             rate_1spike=min(1.0, 2 * rt), max_spike_count=3.0))
    x = M.add_d(pd.DataFrame(rows))
    g = x[x['stimulus'] == '2h [x/diag] @ -6'].set_index('m')
    assert g.loc[1.0, 'd'] == pytest.approx((10.0 - 8.0) / (14.0 - 8.0))
    assert g.loc[0.5, 'd'] == pytest.approx(0.5 * 1.5 / 6.0) and g.loc[0.0, 'd'] == 0.0
    assert g.loc[1.0, 'ratio'] == pytest.approx(0.2 / 0.08) and g.loc[0.0, 'ratio'] == pytest.approx(1.0)
    assert g.loc[1.0, 'd_settled'] == pytest.approx((10.0 - 7.5) / (14.0 - 7.5))


def test_widened():
    lo, hi = M.widened(0.2, 0.3, 0.45)
    assert lo < 0.2 and hi > 0.45
