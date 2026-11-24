"""Tests of data_prep/stage_population_data.py (E134's recorded populations and the lock): the helpers on small made-up
tables, the lock's refusal, and, where the share and the built populations are on hand, the populations against the
existing loaders' counts."""
import json

import numpy as np
import pandas as pd
import pytest

from data_prep import stage_population_data as D


def test_breadth_bins_scaled():
    assert D.breadth_bins(16) == ((0, 0), (1, 1), (2, 3), (4, 7), (8, 16))
    assert D.breadth_bins(17) == ((0, 0), (1, 1), (2, 3), (4, 8), (9, 17))
    assert D.breadth_bins(22) == ((0, 0), (1, 2), (3, 5), (6, 10), (11, 22))
    assert D.breadth_bins(8) == ((0, 0), (1, 1), (2, 2), (3, 3), (4, 8))
    assert D.breadth_bins(12) == ((0, 0), (1, 1), (2, 2), (3, 5), (6, 12))
    for n in range(4, 30):
        b = D.breadth_bins(n)
        assert b[0] == (0, 0) and b[-1][1] == n
        assert all(lo <= hi for lo, hi in b)
        assert all(b[i + 1][0] == b[i][1] + 1 for i in range(len(b) - 1))
    assert D.breadth_labels(16) == ('0', '1', '2-3', '4-7', '8-16')
    with pytest.raises(ValueError):
        D.breadth_bins(3)


def test_robust_trial_sd_planted():
    rng = np.random.default_rng(0)
    R, U, O = 3, 50, 17
    sd = rng.uniform(0.5, 3.0, U)
    signal = np.zeros((U, O))
    signal[:, :3] = 40.0                                        # a few strong odors per unit do not inflate the MAD
    X = signal[None] + rng.normal(size=(R, U, O)) * sd[None, :, None]
    est = D.robust_trial_sd(X)
    assert np.median(np.abs(est / sd - 1)) < 0.15
    X[:, 0, :] = np.nan
    X[0, 1, :] = np.nan                                         # one trial missing: still defined
    X[1:, 2, :] = np.nan                                        # one trial left of every odor: undefined
    est = D.robust_trial_sd(X)
    assert np.isnan(est[0]) and np.isfinite(est[1]) and np.isnan(est[2])


def test_long_round_trip():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(3, 4, 5))
    X[1, 2, 3] = np.nan
    odors = list('abcde')
    d = D.long_from_trials('f1', ['u0', 'u1', 'u2', 'u3'], odors, X)
    assert len(d) == X.size - 1
    Y, units = D.trials_from_long(d, odors)
    assert units == ['u0', 'u1', 'u2', 'u3']
    np.testing.assert_array_equal(np.isnan(Y), np.isnan(X))
    np.testing.assert_allclose(Y[np.isfinite(Y)], X[np.isfinite(X)])


def _planted_long(n_flies=3, n_units=20, odors=D.APL16_ODORS, seed=0, missing=None):
    rng = np.random.default_rng(seed)
    rows = []
    for f in range(n_flies):
        X = rng.gamma(1.0, 1.0, size=(3, n_units, len(odors)))
        if missing and f in missing:
            X[:, :, [odors.index(o) for o in missing[f]]] = np.nan
        rows.append(D.long_from_trials(f'fly{f}', [str(u) for u in range(n_units)], odors, X))
    return pd.concat(rows, ignore_index=True)


SPEC = dict(units=D.UNITS_BASELINE, criterion=1.9684, window={'response': (2.0, 8.0), 'pre': (-5.0, 0.0)}, notes=[],
            recording_site=D.RECORDING_SITE['kc_soma'])


def test_population_from_long_tables_and_responses():
    d = _planted_long(missing={1: ['paa', 'euc']})
    key = ('kc_soma', 'apl16', 'control')
    p = D.population_from_long(key, d, trials=True, spec=SPEC)
    assert p.flies == ('fly0', 'fly1', 'fly2')
    assert list(p.tables['fly0'].columns) == list(D.APL16_ODORS)
    assert p.tables['fly1'][['paa', 'euc']].isna().all().all()
    assert p.fly_odors('fly1') == tuple(o for o in D.APL16_ODORS if o not in ('paa', 'euc'))
    m = d[d.fly == 'fly0'].groupby(['unit', 'odor']).value.mean().unstack()
    np.testing.assert_allclose(p.tables['fly0'].loc[m.index, m.columns].to_numpy(), m.to_numpy())
    r = p.responses('fly1')
    assert r[['paa', 'euc']].isna().all().all()
    assert set(np.unique(r.to_numpy()[np.isfinite(r.to_numpy())])) <= {0.0, 1.0}
    assert ((p.tables['fly0'] > 1.9684).astype(float) == p.responses('fly0')).all().all()
    assert p.has_trials and p.trials is not None
    assert p.trials['fly0'].columns.names == ['odor', 'repeat']
    assert p.breadth_bins() == D.breadth_bins(16)
    assert p.breadth_bins('fly1') == D.breadth_bins(14)
    # trial means only
    d2 = d.assign(repeat=-1).groupby(['fly', 'unit', 'odor', 'repeat'], as_index=False).value.mean()
    p2 = D.population_from_long(key, d2, trials=True, spec=SPEC)
    assert not p2.has_trials and p2.trials is None


def test_hal_measures_read_the_population():
    from kcq import population_fit as P
    d = _planted_long(n_flies=4)
    p = D.population_from_long(('kc_soma', 'apl16', 'control'), d, spec=SPEC)
    flies = P.fly_populations(p)
    ref = P.fly_reference(flies, list(D.APL16_ODORS), D.breadth_bins(16))
    assert set(ref.parts) == {'strength/rate', 'strength/mean', 'sharing/who', 'sharing/how_strongly', 'breadth/bins'}
    assert len(ref.parts['breadth/bins'].columns) == 5


def test_mixture_names_are_the_models():
    assert D._olfdb_mix_name('kmix @ -2') == 'kmix-2'
    assert D._olfdb_mix_name('kmix @ 0') == 'kmix0'
    assert D._olfdb_mix_name('cmix @ -1.0') == 'cmix-1'
    assert D._olfdb_mix_name('IaA @ -3.7') == 'IaA'
    from kcq import harness as H
    for o in ('kmix0 @ 0', 'kmix-1 @ 0', 'kmix-2 @ 0', 'IaA @ -3.7'):
        assert H.odor_abbrev(o) in D.KMIX_ODORS
    for o in ('cmix0 @ 0', 'cmix-2 @ 0', '1o3ol @ -3'):
        assert H.odor_abbrev(o) in D.CMIX_ODORS


def test_v2_stimulus_names_and_mask():
    names, given = D._v2_stims(['PAA @ -4.0', 'pfo @ 0.0', 'geos @ -5.0', 'geos @ -4.0', '1-prop @ -3.0',
                                'euc @ -2.0', 'MethOH @ -2.0'])
    assert names == ['paa', None, None, 'geos', '1-3ol', 'euc', None]
    assert set(n for n in names if n) <= set(D.VALIDATION2_ODORS)
    np.testing.assert_array_equal(D.v2_mask([0.0, 1.0, np.nan, 3.0]), [False, True, False, True])


def test_orn_long_trial_noise_units(tmp_path, monkeypatch):
    rng = np.random.default_rng(2)
    rows = []
    for fly in (1, 2):
        for g, sd in (('DM1', 0.1), ('VA2', 0.4)):
            for o in D.MEGAMAT_ODORS:
                mu = 2.0 if o == '2h' else 0.0
                for r in range(3):
                    rows.append(dict(date='2023-04-22', fly_num=fly, roi=g, panel='megamat', odor1=f'{o} @ -3',
                                     repeat=r, value=mu + rng.normal() * sd))
    pd.DataFrame(rows).to_csv(tmp_path / 'ds1_long.csv', index=False)
    monkeypatch.setattr(D, 'olfdb_path', lambda ds: tmp_path / f'ds{ds}_long.csv')
    d = D.orn_long(1, 'megamat')
    assert sorted(d.fly.unique()) == ['2023-04-22/1', '2023-04-22/2']
    sds = d.groupby('unit').noise_sd.median()
    assert 0.05 < sds['DM1'] < 0.2 and 0.2 < sds['VA2'] < 0.8
    m = d[d.odor == '2h'].groupby('unit').value.mean()
    assert m['DM1'] > m['VA2']                                  # the same dF/F is more SDs above noise in DM1


def test_z_bridge_planted(tmp_path):
    rng = np.random.default_rng(3)
    n = 2000
    v = rng.normal(size=n)
    r = pd.DataFrame(dict(fly='f', cell=np.arange(n), trial=0, stim='2h @ -3.0', dfc=v * 2.0, sd=1.0,
                          z_respvec=v * 0.5))
    r.to_parquet(tmp_path / 'frames_reduced.parquet')
    b = D.z_bridge(tmp_path, criterion=1.0)                     # dfc / sd > 1 <=> v > 0.5 <=> z > 0.25
    assert abs(b['z'] - 0.25) < 0.02
    assert abs(b['rate_noise'] - np.mean(v > 0.5)) < 1e-12


def test_lock_refuses_a_changed_file(tmp_path, monkeypatch):
    key = ('kc_soma', 'apl16', 'control')
    d = _planted_long()
    d.to_parquet(D.pop_path(key, tmp_path), index=False)
    good = D.sha256_file(D.pop_path(key, tmp_path))
    spec = {k: (list(v) if isinstance(v, tuple) else v) for k, v in SPEC.items()}
    entry = dict(file=dict(path='x', sha256=good), sources=[], units=spec['units'], criterion=spec['criterion'],
                 window={k: list(v) for k, v in SPEC['window'].items()}, notes=[], recording_site='somata',
                 exclusions=[{'what': 'planted'}])
    lock = dict(populations={D.pop_name(key): entry})
    monkeypatch.setattr(D, 'read_lock', lambda path=None: lock)
    D._SHA_CACHE.clear()
    p = D.load_population(*key, pop_dir=tmp_path)
    assert p.exclusions == [{'what': 'planted'}] and p.criterion == 1.9684
    entry['file']['sha256'] = '0' * 64
    D._SHA_CACHE.clear()
    with pytest.raises(ValueError):
        D.load_population(*key, pop_dir=tmp_path)
    with pytest.raises(KeyError):
        D.load_population('kc_soma', 'mixtures', 'wild_type', pop_dir=tmp_path)


def test_require_lock_needs_a_committed_lock(tmp_path, monkeypatch):
    monkeypatch.setattr(D, 'read_lock', lambda path=None: None)
    with pytest.raises(RuntimeError):
        D.require_lock()


# ---------------------------------------------------------------------------
# on the share: the built populations against the existing loaders' counts

EXPECTED_FLIES = {
    ('orn', 'megamat', 'wild_type'): 9, ('pn_dendrite', 'megamat', 'wild_type'): 7,
    ('pn_bouton', 'megamat', 'wild_type'): 5, ('kc_claw', 'megamat', 'wild_type'): 4,
    ('kc_soma', 'megamat', 'wild_type'): 4, ('kc_claw', 'megamat', 'same_fly'): 1,
    ('kc_soma', 'megamat', 'same_fly'): 1, ('orn', 'validation2', 'wild_type'): 6,
    ('kc_soma', 'validation2', 'wild_type'): 20, ('orn', 'kmix', 'wild_type'): 9, ('kc_soma', 'kmix', 'wild_type'): 7,
    ('orn', 'cmix', 'wild_type'): 9, ('kc_soma', 'cmix', 'wild_type'): 4, ('kc_soma', 'apl16', 'control'): 11,
    ('kc_soma', 'apl16', 'tnt'): 8,
}


def _built() -> bool:
    return all(D.pop_path(k).exists() for k in D.POPULATIONS)


@pytest.mark.skipif(not _built(), reason='the populations are not built on this host (--build)')
def test_built_populations_flies_and_names():
    from kcq import harness as H
    assert set(EXPECTED_FLIES) == set(D.POPULATIONS)
    model_gloms = set(H.apl16_orn_deltas().index)
    for key, n in EXPECTED_FLIES.items():
        p = D.load_population(*key, check_lock=D.read_lock() is not None)
        assert len(p.flies) == n, key
        assert not any('2022-10-21__fly02' in f for f in p.flies)
        assert set(p.odors) == set(D.PANEL_ODORS[key[1]])
        if p.named_units:
            units = set().union(*(set(t.index) for t in p.tables.values()))
            assert len(units & model_gloms) >= 25, key
        assert p.criterion == D.response_criterion(*key)


@pytest.mark.skipif(not _built(), reason='the populations are not built on this host (--build)')
def test_checks_against_the_existing_loaders_passed():
    c = json.loads((D.POP_DIR / 'checks.json').read_text())
    for name in ('orn_megamat_vs_pebbled', 'per_roi_vs_e119', 'apl16_vs_trialmean', 'megamat_vs_loader',
                 'v2_rdms_cosine', 'v2_rdms_correlation'):
        assert c[name]['passed'], (name, c[name])
