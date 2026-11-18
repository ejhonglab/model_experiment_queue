#!/usr/bin/env python3
"""Tests for E109 (the input of paa, B-myr and 2-but searched for the two ρ's on held-out wirings and recordings;
experiments/E109_input_search_rho.py): the objective (the two ρ's of a population on made-up counts), the pattern-part
penalty, the candidates' scores, the start rule (E108's class on made-up classes.csv tables), the gate and the classes on
made-up readings, the readings' plan; with the shared directory (no simulation), the objective at identity (E83's
pattern parts), the bound against the move planner's probes, E107's reference read from the copy of Remy's files and the
validation2 flies' B-myr against E108's."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E109_input_search_rho')

QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        H.data_targets_path().read_bytes()
        return True
    except (FileNotFoundError, OSError):
        return False


needs_shared = pytest.mark.skipif(not _shared(), reason='the shared directory is not available')


# ---------------------------------------------------------------------------
# the objective on made-up matrices

def _counts(rng, n_kc=60, odors=('a', 'b', 'c', 'd', 'e', 'f')):
    return pd.DataFrame(rng.poisson(1.5, size=(n_kc, len(odors))).astype(float), columns=list(odors))


def test_two_rhos_population_is_rho_of_the_mean_matrix():
    from scipy import stats
    rng = np.random.default_rng(0)
    on = {w: _counts(rng) for w in ('w1', 'w2', 'w3')}
    off = {w: _counts(rng) for w in ('w1', 'w2', 'w3')}
    pop_on = X.population_distances(on)
    pop_off = X.population_distances(off)
    # the population is the mean of the wirings' cosine distances
    d1 = X.E83._cs(H.readout_pair_distances(on['w1'], 'raw'))
    expect = pd.concat([X.E83._cs(H.readout_pair_distances(on[w], 'raw')) for w in on], axis=1).mean(axis=1)
    pd.testing.assert_series_equal(pop_on, X.E83._cs(expect), check_names=False)
    assert len(pop_on) == 15 and d1.index.equals(pop_on.index)
    data_on = pd.Series(rng.uniform(size=15), index=pop_on.index)
    data_off = pd.Series(rng.uniform(size=15), index=pop_on.index)
    r = X.two_rhos_population(on, off, data_on, data_off)
    assert r['on'] == pytest.approx(stats.spearmanr(pop_on.to_numpy(), data_on.reindex(pop_on.index).to_numpy())[0])
    assert r['off'] == pytest.approx(stats.spearmanr(pop_off.to_numpy(), data_off.reindex(pop_off.index).to_numpy())[0])
    assert r['sum'] == pytest.approx(r['on'] + r['off'])
    # (responses, counts) tuples give the same as counts
    r2 = X.two_rhos_population({w: (None, c) for w, c in on.items()}, off, data_on, data_off)
    assert r2 == r
    # the data at the model's own distances: ρ = 1
    assert X.two_rhos_population(on, off, pop_on, pop_off)['sum'] == pytest.approx(2.0)


def test_pattern_penalty():
    assert X.pattern_penalty(0.0570, 0.0569) == 0.0                       # within the margin
    assert X.pattern_penalty(0.0569 + 0.00149, 0.0569) == 0.0             # just inside it
    assert X.pattern_penalty(0.0569 + 0.0025, 0.0569) == pytest.approx(10 * 0.001 / 0.0569)
    assert X.pattern_penalty(0.03, 0.0569) == 0.0                         # a lower pattern part costs nothing
    assert X.PATTERN_MARGIN == 0.0015 and X.PENALTY == 10.0


def test_score_and_candidate_scores():
    v = dict(rho_on=0.70, rho_off=0.75, rho_sum=1.45, pattern=0.060)
    base = {'w1': (0.10, 0.20), 'w2': (0.10, 0.20)}
    s = X.score_of(v, {'w1': (0.10, 0.20), 'w2': (0.10, 0.20)}, base, rms0=0.057)
    assert s['rate_penalty'] == 0.0 and s['pattern_penalty'] == pytest.approx(10 * 0.0015 / 0.057)
    assert s['score'] == pytest.approx(1.45 - s['pattern_penalty'])
    # a response fraction 25 % above the wiring's own: 10 x (0.25 - 0.15) per condition
    s = X.score_of(dict(v, pattern=0.057), {'w1': (0.125, 0.20), 'w2': (0.10, 0.20)}, base, rms0=0.057)
    assert s['rate_penalty'] == pytest.approx(1.0) and s['score'] == pytest.approx(0.45)

    # candidates: the odor's column of each incumbent replaced; the objective sees the replaced matrices
    rng = np.random.default_rng(1)
    odors = ['2h @ -3', 'IaA @ -3', 'pa @ -3']
    inc = {}
    for w in ('w1', 'w2'):
        c_on, c_off = _counts(rng, odors=odors), _counts(rng, odors=odors)
        inc[w] = {'on': (c_on.gt(1.5).astype(float), c_on), 'off': (c_off.gt(1.5).astype(float), c_off)}
    n_kc = inc['w1']['on'][1].shape[0]
    res = {w: (np.ones((n_kc, 2)), np.full((n_kc, 2), 7.0), np.zeros((n_kc, 2)), np.full((n_kc, 2), 3.0)) for w in inc}
    seen = []

    def obj(on, off):
        seen.append((on, off))
        return dict(rho_on=float(on['w1']['IaA @ -3'].mean()), rho_off=float(off['w2']['IaA @ -3'].mean()),
                    rho_sum=float(on['w1']['IaA @ -3'].mean() + off['w2']['IaA @ -3'].mean()), pattern=0.05)
    rows = X.candidate_scores(inc, 'IaA', res, obj, {'w1': (0.5, 0.5), 'w2': (0.5, 0.5)}, rms0=0.05)
    assert len(rows) == 2
    for r in rows:
        assert r['rho_sum'] == pytest.approx(10.0)
        # fractions 1 and 0 against 0.5: 10 x (1 - 0.15) per condition and wiring
        assert r['rate_penalty'] == pytest.approx(4 * 10 * 0.85)
        assert r['fracs'] == {'w1': (1.0, 0.0), 'w2': (1.0, 0.0)}
    # the other columns are the incumbents'
    on, off = seen[0]
    pd.testing.assert_series_equal(on['w1']['2h @ -3'], inc['w1']['on'][1]['2h @ -3'])


# ---------------------------------------------------------------------------
# the start rule (E108's class)

def _classes(tmp_path, **row):
    base = dict(input='X3', carrying=True, complete=True, gain=True, fails='', pending_fails='')
    df = pd.DataFrame([dict(input='X1', carrying=False, complete=True, gain=False), dict(base, **row)])
    df.to_csv(tmp_path / 'classes.csv', index=False)
    return tmp_path


def test_e108_class_waits(tmp_path):
    st = X.e108_class(tmp_path)
    assert not st['written'] and 'classes.csv' in st['waiting_for']
    p = _classes(tmp_path, complete=False, **{'provisional, E106 note pending': True})
    st = X.e108_class(p)
    assert not st['written'] and st['carrying_rule'] == 'X3' and 'X3' in st['waiting_for']
    p = _classes(tmp_path, **{'class': 'better', 'class_without_pending': 'better', 'provisional, E106 note pending': True})
    st = X.e108_class(p)
    assert not st['written'] and st['cls'] == 'better' and st['provisional'] == 'provisional, E106 note pending'
    assert '--e108-class' in st['waiting_for']
    p = _classes(tmp_path, **{'class': 'a price', 'class_without_pending': 'better', 'provisional, E106 note pending': False})
    st = X.e108_class(p)
    assert st['written'] and st['cls'] == 'a price' and st['waiting_for'] is None


def test_start_rule(tmp_path):
    p = _classes(tmp_path, complete=False, **{'provisional, E106 note pending': True})
    # the six searches of the three odors refuse, naming what they wait for
    for n in X.MAIN:
        with pytest.raises(X.StartRefused, match="waits for E108's class"):
            X.require_e108_class(n, results=p)
    # the control searches do not depend on E108
    for n in X.CONTROL_SEARCHES:
        assert X.require_e108_class(n, results=p) == dict(applies=False, start='lab')
    # a confirmation does not help while the runs are not complete
    with pytest.raises(X.StartRefused, match='--e108-class given'):
        X.require_e108_class('A1', dict(cls='better', source='note'), results=p)
    # provisional but complete: refused without, accepted with a matching confirmation, refused with another class
    p = _classes(tmp_path, **{'class': 'better', 'class_without_pending': 'better', 'provisional, E106 note pending': True})
    with pytest.raises(X.StartRefused, match='provisional'):
        X.require_e108_class('B2', results=p)
    st = X.require_e108_class('B2', dict(cls='better', source="E108's note"), results=p)
    assert st['written'] and st['start'] == 'e108' and st['written_in'] == "E108's note"
    with pytest.raises(X.StartRefused, match='is not E108'):
        X.require_e108_class('B2', dict(cls='no gain', source='note'), results=p)
    # written: "better" starts paa from E108's input, any other class from the lab's
    p = _classes(tmp_path, **{'class': 'better', 'class_without_pending': 'better'})
    assert X.require_e108_class('C1', results=p)['start'] == 'e108'
    p = _classes(tmp_path, **{'class': 'a price', 'class_without_pending': 'better'})
    assert X.require_e108_class('C1', results=p)['start'] == 'lab'
    p = _classes(tmp_path, **{'class': 'no gain', 'class_without_pending': 'no gain'})
    assert X.require_e108_class('D1', results=p)['start'] == 'lab'


def test_start_input_lab():
    orn, info = X.start_input(None)
    assert info['start'] == 'lab' and orn.equals(H.apl16_orn_deltas())
    orn, info = X.start_input(dict(start='lab'))
    assert info['start'] == 'lab'


def test_searches_and_seeds():
    assert X.SEARCH_ODORS == ('paa', 'B-myr', '2-but') and X.CONTROL_ODORS == ('2h', '6al')
    assert [X.SEARCHES[n]['half'] for n in X.MAIN] == ['A', 'A', 'B', 'B', 'C', 'D']
    assert [X.SEARCHES[n]['half'] for n in X.CONTROL_SEARCHES] == ['A', 'B']
    assert len({c['seed'] for c in X.SEARCHES.values()}) == len(X.SEARCHES)
    assert (X.POP, X.GENS, X.SIGMA0, X.BOUNDS) == (32, 8, 0.3, (0.5, 2.0))
    assert X.GATE_SIZE == 0.6 and X.INVERT_SIZE == -0.6 and X.PERM_SEEDS == (1, 2, 3, 4)


# ---------------------------------------------------------------------------
# the gate and the classes on made-up readings

COLS_RAW = [f'rho_{c}_{t}' for c in ('on', 'off') for t in X.TARGETS] + [f'{t}_loss' for t in X.TARGETS] + \
    ['block', 'vcf_on', 'vcf_off', 'struct_wt_on', 'struct_wt_off', 'paa_bmyr_on', 'paa_bmyr_off', 'mean_effect']
COLS_CAL = [f'{t}_{k}' for t in X.TARGETS for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson_effect',
                                                     'rho_effect', 'loss')] + ['mean_effect']


def _changes():
    ch = ['I0']
    for n in X.SEARCHES:
        ch += [f'{n}_s0.6', f'{n}_s1']
    for n in X.MAIN:
        ch += [f'{n}_perm{s}' for s in X.PERM_SEEDS] + [f'{n}_inv']
    return ch


def _made_up(delta=None, wiring_delta=None, flagged=()):
    """Readings and runs with every value 0 at I0, plus `delta[change][column]` on T and every wiring (or
    `wiring_delta[change][wiring][column]` on single wirings)."""
    delta, wiring_delta = delta or {}, wiring_delta or {}
    rows, runs = [], []
    for ch in _changes():
        for unit in ('T',) + X.T_WIRINGS:
            for ro, cols in (('raw', COLS_RAW), ('cal', COLS_CAL)):
                r = dict(change=ch, unit=unit, readout=ro, **{c: 0.0 for c in cols})
                for c, v in delta.get(ch, {}).items():
                    if c in r:
                        r[c] += v
                for c, v in wiring_delta.get(ch, {}).get(unit, {}).items():
                    if c in r:
                        r[c] += v
                rows.append(r)
        for w in X.T_WIRINGS:
            runs.append(dict(change=ch, wiring=w, flagged=(ch, w) in flagged, rate_rho_on=0.5, rate_rho_off=0.4))
    return pd.DataFrame(rows), pd.DataFrame(runs)


def _all_rho(c, v):
    return {f'rho_{c}_{t}': v for t in X.TARGETS}


def _main(d, names=X.MAIN, size='0.6'):
    return {f'{n}_s{size}': dict(d) for n in names}


def test_gate_passes_and_price():
    rd, runs = _made_up(_main(_all_rho('on', 0.02)))
    g = X.gate(rd, runs)
    assert g['complete'] and g['gain_on'] and not g['gain_off'] and g['gain']
    assert g['directions']['A']['n_wirings_gain_on'] == 4 and g['mean_gain_on'] == pytest.approx(0.02)
    assert g['controls']['on']['ok'] and g['cls'] == 'passes' and g['outcome'] == 'a'
    # the pattern part 0.004 higher against the unseen half (in one direction): a price, named with its size
    d = _main(_all_rho('on', 0.02))
    for n in ('B1', 'B2'):
        d[f'{n}_s0.6'].update({'A_pattern': 0.004})
    rd, runs = _made_up(d)
    g = X.gate(rd, runs)
    assert g['cls'] == 'a price' and g['fails'] == {'pattern part (calibrated, unseen half)': pytest.approx(0.004)}
    # a run outside the bands: a price
    rd, runs = _made_up(_main(_all_rho('on', 0.02)), flagged=[('A1_s0.6', 'P3')])
    g = X.gate(rd, runs)
    assert g['cls'] == 'a price' and g['fails']['realism (runs outside the bands)'] == 1
    # a pending item only: a price, and passes without the pending items
    d = _main(dict(_all_rho('on', 0.02), block=-0.02))
    rd, runs = _made_up(d)
    g = X.gate(rd, runs)
    assert g['cls'] == 'a price' and g['class_without_pending'] == 'passes' and 'the block' in g['pending_fails']


def test_gate_gain_conditions():
    # below the threshold: carries below the gate (the mean gain 0.007 >= 0.005, the controls hold)
    rd, runs = _made_up(_main(_all_rho('on', 0.007)))
    g = X.gate(rd, runs)
    assert not g['gain'] and g['cls'] == 'carries below the gate' and g['outcome'] == 'c'
    # APL off needs 0.015
    rd, runs = _made_up(_main(_all_rho('off', 0.012)))
    assert not X.gate(rd, runs)['gain_off']
    rd, runs = _made_up(_main(_all_rho('off', 0.016)))
    assert X.gate(rd, runs)['cls'] == 'passes'
    # the other ρ lower by more than 0.005: no gain
    rd, runs = _made_up(_main(dict(_all_rho('on', 0.02), **_all_rho('off', -0.006))))
    g = X.gate(rd, runs)
    assert not g['gain_on'] and not g['gain_on_parts']['other']
    # only two of the four wirings gain: no gain
    wd = {f'{n}_s0.6': {'P3': _all_rho('on', -0.015), 'P4': _all_rho('on', -0.015)} for n in X.MAIN}
    rd, runs = _made_up(_main(_all_rho('on', 0.02)), wd)
    g = X.gate(rd, runs)
    assert g['directions']['A']['n_wirings_gain_on'] == 2 and not g['gain_on']
    assert g['cls'] == 'carries below the gate'
    # the second split the other way: no gain
    d = _main(_all_rho('on', 0.02), X.FIRST_SPLIT)
    d.update(_main(_all_rho('on', -0.01), X.SECOND_SPLIT))
    rd, runs = _made_up(d)
    g = X.gate(rd, runs)
    assert not g['gain_on'] and not g['gain_on_parts']['second_split'] and g['gain_on_parts']['directions']
    # one direction only: no gain
    d = _main(_all_rho('on', 0.02), ('A1', 'A2', 'C1', 'D1'))
    rd, runs = _made_up(d)
    assert not X.gate(rd, runs)['gain_on']
    # a direction is the mean of its two draws
    d = _main(_all_rho('on', 0.02))
    d['A1_s0.6'] = _all_rho('on', 0.035)
    d['A2_s0.6'] = _all_rho('on', 0.005)
    rd, runs = _made_up(d)
    assert X.gate(rd, runs)['directions']['A']['d_rho_on'] == pytest.approx(0.02)


def test_gate_controls():
    # the permuted factors gain half the candidate's or more: does not carry
    d = _main(_all_rho('on', 0.02))
    for n in X.MAIN:
        for s in X.PERM_SEEDS:
            d[f'{n}_perm{s}'] = _all_rho('on', 0.012)
    rd, runs = _made_up(d)
    g = X.gate(rd, runs)
    assert g['gain_on'] and not g['controls']['on']['ok'] and g['cls'] == 'does not carry' and g['outcome'] == 'd'
    assert g['controls']['on']['permuted'] == pytest.approx(0.012)
    # the inverted factors, the control search likewise
    for key in ('inv', 'ctl'):
        d = _main(_all_rho('on', 0.02))
        if key == 'inv':
            d.update({f'{n}_inv': _all_rho('on', 0.011) for n in X.MAIN})
        else:
            d.update({f'{n}_s0.6': _all_rho('on', 0.011) for n in X.CONTROL_SEARCHES})
        rd, runs = _made_up(d)
        assert X.gate(rd, runs)['cls'] == 'does not carry'
    # just under half holds
    d = _main(_all_rho('on', 0.02))
    d.update({f'{n}_s0.6': _all_rho('on', 0.0099) for n in X.CONTROL_SEARCHES})
    rd, runs = _made_up(d)
    assert X.gate(rd, runs)['cls'] == 'passes'
    # nothing gains: does not carry
    rd, runs = _made_up({})
    assert X.gate(rd, runs)['cls'] == 'does not carry'
    # the control searches read with the same gate (no second split, no controls)
    rd, runs = _made_up({f'{n}_s0.6': _all_rho('on', 0.02) for n in X.CONTROL_SEARCHES})
    g = X.gate(rd, runs, main=X.CONTROL_SEARCHES, second=(), with_controls=False)
    assert g['gain_on'] and g['controls'] is None and g['cls'] == 'passes'


def test_gate_incomplete_and_tables(tmp_path, monkeypatch):
    rd, runs = _made_up(_main(_all_rho('on', 0.02)))
    g = X.gate(rd[rd.change != 'B2_s0.6'], runs)
    assert not g['complete'] and 'missing' in g['note']
    with pytest.raises(ValueError):
        X.class_of(g)
    monkeypatch.setenv('E109_OUT_DIR', str(tmp_path))
    gt, verdict = X.gate_tables(rd, runs)
    assert verdict['cls'] == 'passes' and (tmp_path / 'gate.csv').exists() and (tmp_path / 'gate.json').exists()
    assert set(gt['which']) == {'six', 'ctl'} and gt['carrying'].sum() == 1


# ---------------------------------------------------------------------------
# the readings' plan

def test_change_meta_and_plan(monkeypatch):
    monkeypatch.setattr(X, 'search_finished', lambda n: n in ('A1', 'ctlA'))
    m = X.change_meta()
    assert set(m) == {'A1_s0.6', 'A1_s1', 'A1_perm1', 'A1_perm2', 'A1_perm3', 'A1_perm4', 'A1_inv', 'ctlA_s0.6', 'ctlA_s1'}
    assert m['A1_inv']['size'] == -0.6 and m['A1_perm3']['seed'] == 3 and m['A1_s0.6']['unseen'] == 'B'
    plan = X.read_plan()
    assert len(plan) == 4 * len(m)
    first = [p['change'] for p in plan if p['wiring'] == 'P3']
    assert first[:2] == ['A1_s0.6', 'ctlA_s0.6'] and set(first[-2:]) == {'A1_s1', 'ctlA_s1'}
    assert [p['wiring'] for p in plan][::len(m)] == list(X.T_WIRINGS)


def test_recurring_entries_on_made_up_factors(monkeypatch):
    lab = H.abbrev_columns(H.apl16_orn_deltas())
    rng = np.random.default_rng(3)
    F = {}
    for i, n in enumerate(X.MAIN):
        f = pd.DataFrame(1.0, index=lab.index, columns=lab.columns)
        for o in X.SEARCH_ODORS:
            f[o] = np.exp(rng.uniform(np.log(0.5), np.log(2.0), len(lab.index)))
        f.loc['DM1', 'paa'] = 2.0 if i < 5 else 0.6        # five of six up, mean log2 = (5 - 0.74) / 6 > 0.5
        f.loc['VA2', 'B-myr'] = 0.5                         # all six down
        F[n] = f
    fs = pd.DataFrame([dict(odor=o, glomerulus=g, flies_mean=1.0, n_flies=6, n_positive=6, flies_sd_log2=0.5,
                            flies_source='made up') for o in X.SEARCH_ODORS for g in lab.index])
    monkeypatch.setattr(X, 'flies_spread', lambda: fs)
    ent, null = X.recurring_entries(F, n_null=5, write=False)
    hit = ent[ent.meets]
    assert {('paa', 'DM1'), ('B-myr', 'VA2')} <= set(zip(hit.odor, hit.glomerulus))
    r = ent[(ent.odor == 'B-myr') & (ent.glomerulus == 'VA2')].iloc[0]
    assert r.n_agree == 6 and r.mean_log2 == pytest.approx(-1.0) and r.mean_log2_read == pytest.approx(-0.6)
    assert r.read_in_flies_sd == pytest.approx(-1.2) and r.lab_input == pytest.approx(lab.loc['VA2', 'B-myr'])
    assert set(ent.odor) == set(X.SEARCH_ODORS) and len(null) == 5
    with pytest.raises(RuntimeError, match='six'):
        X.recurring_entries({k: v for k, v in F.items() if k != 'D1'}, write=False)


# ---------------------------------------------------------------------------
# with the shared directory (no simulation)

@needs_shared
def test_objective_at_identity():
    """S's population at the saved fits: the pattern part under the calibrated readout against halves A and B is E83's
    IDENTITY_OBJECTIVE (the searches' starting point), and the two ρ's are those of the mean matrix."""
    E83 = X.E83
    saved = {w: E83.saved_fit(w)['saved'] for w in X.S_WIRINGS}
    for h in ('A', 'B'):
        v = X.RhoObjective(h)(on={w: saved[w]['on'][1] for w in X.S_WIRINGS}, off={w: saved[w]['off'][1] for w in X.S_WIRINGS})
        assert abs(v['pattern'] - E83.IDENTITY_OBJECTIVE[h]) <= 5e-5
        assert 0.6 < v['rho_on'] < 0.75 and 0.7 < v['rho_off'] < 0.8
        # the same as E106's two ρ's of the population read from E83's raw distances
        pop = E83.population([E83.pair_distances(saved[w], 'raw') for w in X.S_WIRINGS])
        dt = X.dist_targets()[h]
        assert v['rho_on'] == pytest.approx(X.H.spearman(pop['on'], dt['on']), abs=1e-12)


@needs_shared
def test_bound_reproduces_the_planners_probes():
    bt = X.bound_table()
    for (kind, unit, which, how, t), (pon, poff) in X.PROBES.items():
        if kind != 'bound':
            continue
        r = bt[(bt.unit == unit) & (bt.odors == which) & (bt.how == how) & (bt.target == t)].iloc[0]
        assert abs(r.rho_on - pon) <= 0.0015 and abs(r.rho_off - poff) <= 0.0015, (unit, which, how)
    assert set(bt[bt.odors == 'three'].n_pairs) == {42} and set(bt[bt.odors == 'controls'].n_pairs) == {29}


@needs_shared
def test_e107_reference_from_the_copy():
    from experiments import E107_validation2_panel as E107
    keep = E107.remy_dir
    ref = X.e107_reference(X.ALL_WIRINGS)
    r = X.remy_readings(ref['pop']['on'])
    assert E107.remy_dir is keep
    assert r['set'] == '136' and r['n_all'] == 136 and r['n_pb'] == 31
    assert abs(r['rho_all'] - 0.676) <= 0.0015 and abs(r['rho_pb'] - 0.542) <= 0.0015
    assert abs(r['d_paa_bmyr_on'] - 0.501) <= 0.0015 and abs(r['remy_d_paa_bmyr'] - 0.873) <= 0.0015


@needs_shared
def test_validation2_flies_bmyr_is_e108s():
    from experiments import E108_paa_extrapolation as E108
    V = X.validation2_fly_values()
    np.testing.assert_array_equal(V['B-myr'], E108.bmyr_flies()['values'])
    assert V['paa'].shape == V['B-myr'].shape and np.isfinite(V['paa']).sum() > 0
