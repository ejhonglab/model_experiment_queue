"""Tests of E100's rules on made-up data (experiments/E100_apl_weak_odors.py): the solvent rules, the shuffled-block reader,
W1's count rule and the rule over input estimates, W2's fly rule, W3's pass / fail and its power on planted maps, the
concentration matching to one log unit (inclusive), the input estimates per session, the reader of E98's B1 runs, and
part 0's identities that need no new recording (on the shared directory's data, skipped where it is absent)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from experiments import E100_apl_weak_odors as M


# ---------------------------------------------------------------------------
# concentration matching and the input estimates

def test_concentration_match_inclusive_one_log_unit():
    assert M.concentration_match(-3.0, [-4.0]) == -4.0            # exactly 1 log unit: in (inclusive)
    assert M.concentration_match(-3.0, [-4.01]) is None
    assert M.concentration_match(-7.0, [-6.0, -3.0]) == -6.0
    assert M.concentration_match(-5.5, [-6.0, -5.0]) == -5.0      # a tie goes to the higher concentration
    assert M.concentration_match(-2.0, [-2.0, -3.0]) == -2.0
    assert M.concentration_match(-2.0, []) is None


def _session(cols, gl=('DM4', 'VM5d', 'DC3')):
    rng = np.random.default_rng(len(cols))
    df = pd.DataFrame(rng.uniform(0, 100, size=(len(gl), len(cols))), index=pd.Index(gl, name='glomerulus'),
                      columns=pd.MultiIndex.from_tuples(cols, names=['panel', 'odor']))
    return df


def test_input_estimates_per_session_and_panel():
    inputs = {
        'kiwi-control': _session([('glomeruli_diagnostics', '2h @ -6'), ('control', '2h @ -5'),
                                  ('glomeruli_diagnostics', 'ms @ -3'), ('control', 'ms @ -3')]),
        'megamat': _session([('glomeruli_diagnostics', '2h @ -6'), ('megamat', '2h @ -3')]),
        'validation2': _session([('validation2', 'paa @ -4'), ('glomeruli_diagnostics', 'paa @ -5')]),
    }
    e = M.input_estimates(('2h', -7.0), inputs)
    # kiwi-control's nearest is -6 (1 unit; -5 is 2 away), megamat's -6; validation2 has no 2h
    assert sorted(zip(e['session'], e['column'])) == [('kiwi-control', '2h @ -6'), ('megamat', '2h @ -6')]
    assert (e['status'] == 'substitute').all()
    # the two panels of one session at the nearest concentration are separate estimates
    e = M.input_estimates(('ms', -3.0), inputs)
    assert sorted(e['panel_name']) == ['control', 'glomeruli_diagnostics'] and (e['status'] == 'exact').all()
    assert e['estimate'].is_unique
    # paa at 10^-3 takes the panel's 10^-4 (validation2's own panel), not the diagnostics' 10^-5
    e = M.input_estimates(('paa', -3.0), inputs)
    assert list(e['column']) == ['paa @ -4'] and list(e['panel_name']) == ['validation2']
    assert M.input_estimates(('2h', -8.5), inputs).empty
    # stimuli that the rule gives the same column share its label
    assert set(M.input_estimates(('2h', -7.0), inputs)['estimate']) == set(M.input_estimates(('2h', -6.0), inputs)['estimate'])


def test_stim_labels_round_trip():
    for s in [('2h', -6.0), ('ms', -5.5), ('4h-ol', -6.0), ('paa', -2.954)]:
        assert M.parse_stim(M.stim_label(s)) == s
    assert M.stim_label(('2h', -6.0)) == '2h @ -6'


# ---------------------------------------------------------------------------
# the shuffled blocks and the solvent

def _frames():
    # three shuffled blocks of two odors, the solvent last in each; no block column (inferred)
    seq = [('2h', -6), ('ma', -7), ('pfo', np.nan), ('ma', -7), ('2h', -6), ('pfo', np.nan), ('2h', -6), ('ma', -7),
           ('paraffin oil', np.nan)]
    return pd.DataFrame(dict(recording='r1', frame=range(len(seq)), odor=[o for o, _ in seq],
                             log10_conc=[c for _, c in seq]))


def test_shuffled_blocks_inferred_from_the_solvent():
    f = M.shuffled_blocks(_frames())
    assert list(f['block']) == [0, 0, 0, 1, 1, 1, 2, 2, 2]
    assert list(f['position']) == [1, 2, 3] * 3
    assert f['is_solvent'].sum() == 3 and (f.loc[f['is_solvent'], 'stim'] == M.SOLVENT).all()
    assert f['solvent_last'].all() and f['complete'].all() and not f['duplicate'].any()
    assert list(f.loc[f['stim'] == '2h @ -6', 'rep']) == [0, 1, 2]


def test_shuffled_blocks_flags_a_bad_block():
    fr = _frames()
    fr['block'] = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    fr.loc[7, ['odor', 'log10_conc']] = ['2h', -6]          # block 2: 2h twice, ma missing
    fr.loc[5, ['odor', 'log10_conc']] = ['ma', -7]          # block 1: no solvent
    f = M.shuffled_blocks(fr)
    by_block = f.groupby('block')[['solvent_last', 'complete', 'duplicate']].first()
    assert by_block.loc[0].tolist() == [True, True, False]
    assert by_block.loc[1, 'solvent_last'] == False and by_block.loc[1, 'complete'] == False  # noqa: E712
    assert by_block.loc[2, 'duplicate'] == True and by_block.loc[2, 'complete'] == False  # noqa: E712


def _resp():
    f = M.shuffled_blocks(_frames())
    # odor responses 1.0 (2h) and 0.5 (ma) plus block offsets; solvent 0.2, 0.3, 0.4
    off = {0: 0.0, 1: 0.1, 2: 0.2}
    f['response'] = [0.2 + off[b] if s else (1.0 if st == '2h @ -6' else 0.5) + off[b]
                     for s, st, b in zip(f['is_solvent'], f['stim'], f['block'])]
    return f


def test_solvent_rules():
    r = _resp()
    assert M.recording_solvent(r) == pytest.approx(0.3)                 # the mean of the three blocks' solvents
    bs = M.block_solvent(r)
    assert bs.tolist() == pytest.approx([0.2, 0.3, 0.4])
    per, stim = M.solvent_subtracted(r)
    assert stim.loc['2h @ -6', 'raw'] == pytest.approx(1.1)
    assert stim.loc['2h @ -6', 'minus_recording'] == pytest.approx(0.8)
    assert stim.loc['2h @ -6', 'minus_block'] == pytest.approx(0.8)      # each block's offset cancels
    assert per['minus_block'].groupby(per['stim']).std().max() == pytest.approx(0.0, abs=1e-12)
    assert stim.loc['ma @ -7', 'n_presentations'] == 3
    # a block without a solvent gets NaN, not another block's
    r2 = r[~((r['block'] == 1) & r['is_solvent'])]
    assert np.isnan(M.block_solvent(r2).loc[1])
    # large: the solvent (0.3) against a third of the largest mean response (1.1 / 3 = 0.367)
    assert not M.large_solvent(r)
    assert M.large_solvent(r, share=0.25)


def test_fly_ratio_geometric_or_arithmetic():
    pr = pd.DataFrame(dict(fly=['a', 'a', 'b', 'c'], num=[0.5, 0.5, 0.25, 1.0], den=[1.0, 1.0, 1.0, 1.0]))
    r = M.fly_ratio_w(pr, n_boot=200)
    assert r['mean'] == 'geometric' and r['n_flies'] == 3
    assert r['ratio'] == pytest.approx(np.exp(np.mean(np.log([0.5, 0.25, 1.0]))))
    assert r['lo'] <= r['ratio'] <= r['hi']
    # a response at or below the solvent's makes the mean arithmetic
    pr.loc[3, 'num'] = -0.1
    r = M.fly_ratio_w(pr, n_boot=200)
    assert r['mean'] == 'arithmetic' and r['ratio'] == pytest.approx(np.mean([0.5, 0.25, -0.1]))
    # a denominator at or below MIN_RESPONSE leaves the recording out
    pr.loc[2, 'den'] = M.MIN_RESPONSE
    r = M.fly_ratio_w(pr, n_boot=200)
    assert r['n_left_out'] == 1 and r['n_flies'] == 2


def test_calyx_responses_and_experiment_xml(tmp_path):
    stack = np.zeros((2, 4, 4))
    stack[0, 1:3, 1:3] = 2.0
    stack[1] = 1.0
    mask = np.zeros((4, 4), bool)
    mask[1:3, 1:3] = True
    assert M.calyx_responses(stack, mask).tolist() == [2.0, 1.0]
    x = tmp_path / 'Experiment.xml'
    x.write_text('<ThorImageExperiment><Wavelengths><Wavelength name="ChanA" exposureTimeMS="0"/></Wavelengths>'
                 '<LSM pixelX="256" pixelY="224" pixelSizeUM="0.265" channel="1"/></ThorImageExperiment>')
    d = M.experiment_xml(x)
    assert d == dict(pixel_um=0.265, channels=['ChanA'], channel_mask=1, pixel_x=256, pixel_y=224)


# ---------------------------------------------------------------------------
# the verdict rules

def test_ratio_over_estimates():
    assert M.ratio_over_estimates({'kiwi': 'flatter', 'megamat': 'inside', 'v2': 'steeper'}) == 'inside'
    assert M.ratio_over_estimates({'kiwi': 'flatter', 'megamat': 'flatter'}) == 'flatter'
    assert M.ratio_over_estimates({'kiwi': 'steeper'}) == 'steeper'
    assert M.ratio_over_estimates({'kiwi': 'flatter', 'megamat': 'steeper'}) == 'not decided'
    assert M.ratio_over_estimates({'kiwi': 'flatter', 'megamat': 'not decided'}) == 'not decided'
    assert M.ratio_over_estimates({}) == 'not decided'


def test_w1_count_rule():
    six = ['2h', 'ma', '4h-ol', 'farn', '2-but', 'paa']
    assert M.w1_verdict(dict(zip(six, ['flatter'] * 4 + ['inside'] * 2))) == 'flatter'
    assert M.w1_verdict(dict(zip(six, ['flatter'] * 3 + ['inside'] * 3))) == 'not decided'
    assert M.w1_verdict(dict(zip(six, ['inside'] * 5 + ['not decided']))) == 'inside'
    assert M.w1_verdict(dict(zip(six, ['steeper'] * 4 + ['flatter'] * 2))) == 'steeper'
    # the two rules together: per ratio over the estimates, then the count
    per = {o: {'kiwi': 'flatter', 'megamat': 'flatter', 'v2': 'flatter'} for o in six}
    per['2h'] = {'kiwi': 'flatter', 'megamat': 'inside', 'v2': 'flatter'}
    per['ma'] = {'kiwi': 'flatter', 'megamat': 'steeper', 'v2': 'flatter'}
    per['farn'] = {'kiwi': 'inside', 'megamat': 'inside', 'v2': 'inside'}
    classes = {o: M.ratio_over_estimates(v) for o, v in per.items()}
    assert classes == {'2h': 'inside', 'ma': 'not decided', '4h-ol': 'flatter', 'farn': 'inside', '2-but': 'flatter',
                       'paa': 'flatter'}
    assert M.w1_verdict(classes) == 'not decided'                     # 3 flatter, 2 inside, 1 not decided


def test_w2_four_flies():
    assert M.w2_verdict('flatter', 3) == 'described'
    assert M.w2_verdict('flatter', 4) == 'flatter'
    assert M.w2_verdict('inside', 9) == 'inside'


def _exceed(shares):
    """pair -> per-calyx shares of 10 recordings exceeding."""
    rows = []
    for pair, per_calyx in shares.items():
        for ci, sh in enumerate(per_calyx):
            k = int(round(sh * 10))
            rows += [dict(pair=pair, calyx=f'c{ci}', recording=f'r{j}', exceeds=j < k) for j in range(10)]
    return pd.DataFrame(rows)


def test_w3_pass_fail_not_decided():
    strong = [0.6] * 5 + [0.0] * 3                     # more than half in 5 of 8 calyces
    weak = [0.6] * 4 + [0.0] * 4                       # only 4 calyces
    quarter = [0.3] * 5 + [0.0] * 3                    # more than a quarter in 5, not more than half
    none = [0.2] * 8
    assert M.w3_verdict(_exceed(dict(a=strong, b=strong, c=strong)))['verdict'] == 'passes'
    assert M.w3_verdict(_exceed(dict(a=strong, b=strong, c=weak)))['verdict'] == 'not decided'
    assert M.w3_verdict(_exceed(dict(a=none, b=none, c=weak)))['verdict'] == 'fails'
    assert M.w3_verdict(_exceed(dict(a=none, b=quarter, c=none)))['verdict'] == 'not decided'
    # exactly half is not more than half
    assert M.w3_pair_calyces(_exceed(dict(a=[0.5] * 8)), M.W3_PASS_SHARE) == 0
    # only the pairs A4 kept count
    v = M.w3_verdict(_exceed(dict(a=strong, b=strong, c=strong, d=none)), pairs=['a', 'b', 'd'])
    assert v['verdict'] == 'not decided' and v['n_pairs'] == 3 and v['n_pairs_pass'] == 2
    # A4 kept no pair: W3 neither passes nor fails
    assert M.w3_verdict(_exceed(dict(a=none)), pairs=[])['verdict'] == 'no pair'


# ---------------------------------------------------------------------------
# W3's power on planted maps

def test_noise_sd_gives_the_reliability():
    rng = np.random.default_rng(0)
    sig = rng.normal(size=200000)
    sd1 = M.noise_sd(sig, 0.6, n_presentations=1)
    a = sig + rng.normal(0, sd1, sig.size)
    b = sig + rng.normal(0, sd1, sig.size)
    assert np.corrcoef(a, b)[0, 1] == pytest.approx(0.6, abs=0.01)
    assert M.noise_sd(sig, 0.6, 2) == pytest.approx(sd1 / np.sqrt(2))


def test_pair_maps_true_labels():
    rng = np.random.default_rng(1)
    G = rng.uniform(0.1, 1.0, size=(5, 144))
    Wt = rng.uniform(0, 10, size=(5, 2))
    inside = np.ones(144, bool)
    inside[:10] = False
    perms = np.vstack([np.arange(5), rng.permutation(5)])
    D = M.pair_maps(Wt, perms, G, inside)
    a, b = Wt[:, 0] @ G, Wt[:, 1] @ G
    assert np.allclose(D[0], a / a[inside].mean() - b / b[inside].mean())
    assert D.shape == (2, 144)


def _planted_state(n_glom=12, n_perm=99, seed=0, strength=1.0, shared_layout=True):
    """Eight made-up calyces sharing one layout of glomerulus fields (each glomerulus a Gaussian blob at the same place,
    with a little jitter per calyx; with `shared_layout` False each calyx has a layout of its own), one slab at 18 um and
    two blurs; the pair's two stimuli drive glomerulus 0 and glomerulus 1 alone (times `strength`) over a flat
    background."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:12, :12]
    centers = rng.uniform(1, 10, size=(n_glom, 2))
    grids = {}
    for ci, c in enumerate(M.CALYCES):
        if not shared_layout:
            centers = rng.uniform(1, 10, size=(n_glom, 2))
        cj = centers + rng.normal(0, 0.3, size=centers.shape)
        G = np.stack([np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / 4.0).ravel() for y, x in cj]) + 0.05
        perms = np.vstack([np.arange(n_glom)] + [rng.permutation(n_glom) for _ in range(n_perm)])
        side = 'L' if c.endswith('_L') else 'R'
        if side == 'L':        # a left calyx is the mirror image (posterior view)
            G = G.reshape(n_glom, 12, 12)[:, :, ::-1].reshape(n_glom, 144)
        slabs = {18.0: dict(G={2.0: G, 4.0: G}, inside=np.ones(144, bool), anat=None)}
        grids[c] = dict(names=[f'g{i}' for i in range(n_glom)], perms=perms, slabs=slabs, side=side, seed=ci)
    w = np.full((n_glom, 2), 1.0)
    w[0, 0] += 20.0 * strength
    w[1, 1] += 20.0 * strength
    rates = pd.DataFrame(w, index=[f'g{i}' for i in range(n_glom)], columns=['A [s]', 'B [s]'])
    return dict(rates=rates, grids=grids, columns={('A - B', 's'): ('A [s]', 'B [s]')})


def test_detect_planted_and_null():
    st = _planted_state()
    g = st['grids']['flywire_R']
    Wt = M._wt(st['rates'], g['names'], ['A [s]', 'B [s]'])
    P = M.pair_maps(Wt, g['perms'], g['slabs'][18.0]['G'][2.0], g['slabs'][18.0]['inside'])
    ins = np.ones(144, bool)
    best, p95, ex = M.detect({(18.0, 2.0): P}, {(18.0, 2.0): ins}, P[0].copy(), ins)
    assert best == pytest.approx(1.0) and ex
    rng = np.random.default_rng(3)
    hits = [M.detect({(18.0, 2.0): P}, {(18.0, 2.0): ins}, rng.normal(size=144), ins)[2] for _ in range(200)]
    assert np.mean(hits) < 0.15                      # about 5 % by construction


def test_power_on_planted_maps(monkeypatch):
    monkeypatch.setattr(M, 'A4_PLANE_DEPTHS_UM', (18.0,))
    st = _planted_state(strength=1.0)
    d = M.power('A - B', 's', st, reliabilities=[0.6], n_draws=3)
    # 8 truths x 3 draws x 7 scored calyces, less the two hemibrain calyces scored against each other
    assert len(d) == 3 * (8 * 7 - 2)
    assert d['exceeds'].mean() > 0.9
    # each fly's map follows its own calyx's layout, unrelated to the scored calyces': the null, about 5 %
    own = _planted_state(strength=1.0, shared_layout=False)
    d0 = M.power('A - B', 's', own, reliabilities=[0.6], n_draws=3)
    assert d0['exceeds'].mean() < 0.2
    # a weaker reliability gives less power
    d1 = M.power('A - B', 's', st, reliabilities=[0.05], n_draws=3)
    assert d1['exceeds'].mean() < d['exceeds'].mean()


# ---------------------------------------------------------------------------
# E98's B1 reader on a made-up folder

def test_read_e98_b1_made_up(tmp_path):
    for cond, win in (('on', 20.0), ('acute_off', 0.0)):
        pd.DataFrame(dict(win=[win, win + 1], pre=[8.0, 8.0]), index=pd.Index(['pa', 'paa [megamat/diag] @ -5'])
                     ).to_csv(tmp_path / f'real_{cond}_apl.csv')
        pd.DataFrame(dict(pa=[True, False, False, False], **{'paa [megamat/diag] @ -5': [True, True, False, False]})
                     ).to_parquet(tmp_path / f'real_{cond}_responses.parquet')
    b1 = M.read_e98_b1(['real', 'P1'], dirs=[tmp_path])
    assert len(b1) == 4 and b1.attrs['missing'] == [('P1', 'acute_off'), ('P1', 'on')]
    r = b1[(b1['cond'] == 'on') & (b1['stimulus'] == 'paa [megamat/diag] @ -5')].iloc[0]
    assert r['evoked'] == pytest.approx(13.0) and r['response_rate'] == pytest.approx(0.5)
    assert (r['b1_odor'], r['b1_conc'], r['b1_source']) == ('paa', -5.0, 'megamat/diag')
    with pytest.raises(NotImplementedError, match="E98's B1"):
        M.read_e98_b1(['real'], dirs=[tmp_path / 'none'])


# ---------------------------------------------------------------------------
# part 0's identities that need no new recording (shared-directory data)

def _shared():
    try:
        from kcq import harness as H
        return H.shared_dir()
    except FileNotFoundError:
        return None


@pytest.mark.skipif(_shared() is None or not (_shared() / 'data' / 'yang_weak_odor_kc').exists(),
                    reason='shared directory not mounted')
def test_part0_availability_identity():
    res = M.availability_identity()
    assert res['passed'], res['failures']


@pytest.mark.skipif(_shared() is None or not (_shared() / 'data' / 'yang_apl_calyx_maps').exists(),
                    reason='shared directory not mounted')
def test_part0_reader_identity_2021():
    res = M.reader_identity_2021()
    assert res['passed'] and res['n_stimuli'] >= 10, res


@pytest.mark.skipif(_shared() is None or not (_shared() / 'results' / M.EXPERIMENT / 'a3_pn_rates.csv').exists(),
                    reason="A3's PN-stage run not there")
def test_part0_pn_identity():
    res = M.pn_identity()
    assert res['window_vs_e90_hz'] <= 1e-9, res
    if 'legacy_evoked_vs_e79_cache_hz' in res:
        assert res['passed'], res


def test_scaled_reliability():
    assert M.scaled_reliability(0.6, 1.0) == pytest.approx(0.6)
    r = M.scaled_reliability(0.59, 0.5)
    assert r == pytest.approx(0.25 * 0.59 / (0.25 * 0.59 + 0.41))
    # the same noise on half the signal: two presentations of it agree at r
    rng = np.random.default_rng(2)
    sig = rng.normal(size=200000)
    sd = M.noise_sd(sig, 0.59, n_presentations=1)
    a, b = 0.5 * sig + rng.normal(0, sd, sig.size), 0.5 * sig + rng.normal(0, sd, sig.size)
    assert np.corrcoef(a, b)[0, 1] == pytest.approx(r, abs=0.01)


# ---------------------------------------------------------------------------
# the readers of viki's copy (made-up data shaped like it, and the 12 released 2021 recordings)

def test_odor_abbrev_yang():
    assert M.odor_abbrev_yang('phenylacetaldehyde(-6.294)') == ('paa', -6.294)
    assert M.odor_abbrev_yang('2h(-7)') == ('2h', -7.0)
    assert M.odor_abbrev_yang('4-heptanol(-6)') == ('4h-ol', -6.0)
    assert M.odor_abbrev_yang('ms(-5.5)') == ('ms', -5.5)
    s, c = M.odor_abbrev_yang('PFO')
    assert s == M.SOLVENT and np.isnan(c)
    assert M.odor_abbrev_yang('paraffin\n  oil')[0] == M.SOLVENT


W_YAML = """nblock: 3 # total number of blocks
nstim: 4 # number of odor stimuli in EACH block
odorPulseLen : 1
odor_pin: ['NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
  'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
  'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA', 'PFO',
'2h(-7)','NA','farn(-3)','paa(-3)']
odor_stim_seq: [38, 40, 41, 37,
41, 38, 40, 37,
40, 41, 38, 37]
"""


def test_yaml_stimuli_shuffled_blocks(tmp_path):
    y = tmp_path / '20250199_1FRfun1.yaml'
    y.write_text(W_YAML)
    f = M.yaml_stimuli(y)
    assert len(f) == 12 and list(f['block']) == [0] * 4 + [1] * 4 + [2] * 4
    assert list(f['odor'][:4]) == ['2h', 'farn', 'paa', M.SOLVENT] and list(f['log10_conc'][:3]) == [-7.0, -3.0, -3.0]
    # the yaml's blocks agree with the blocks inferred from the solvent, and the design flags hold
    b = M.shuffled_blocks(f.drop(columns='block'))
    assert list(b['block']) == list(f['block']) and b['solvent_last'].all() and b['complete'].all()
    y.write_text(W_YAML.replace('nblock: 3', 'nblock: 2'))
    with pytest.raises(ValueError):
        M.yaml_stimuli(y)


def test_read_recording_refuses_unreleased():
    row = dict(stem='20250115_1FRfun1', folder='/nonexistent', set='W', date='2025-01-15', fly=1, hemisphere='R',
               plane_depth_um=15, pulse_s=1.0)
    with pytest.raises(PermissionError):
        M.read_recording(row)


def test_raw_only_rows_refused():
    raw = "raw only: not used until further notice (request in viki's terminal 2026-11-15 15:28)"
    ok = "processed: usable (Yang's processed products only)"
    assert not M.usable_row(dict(use_until_further_notice=raw))
    assert M.usable_row(dict(use_until_further_notice=ok))
    assert M.usable_row(dict(stem='x'))  # a copy from before the rule
    assert M.usable_row(dict(use_until_further_notice=float('nan')))
    # refused for being raw only even if the recording were released
    row = dict(stem=M.released_2021()[0], use_until_further_notice=raw)
    with pytest.raises(PermissionError, match='raw only'):
        M.read_recording(row)


def test_recording_table_and_flies_ratios():
    t = M.recording_table(_resp())
    assert t['solvent'] == pytest.approx(0.3) and t['n_blocks'] == 3
    assert t['solvent_drift_per_block'] == pytest.approx(0.1) and t['odor_drift_per_block'] == pytest.approx(0.1)
    assert t['block_agreement_r'] == pytest.approx(1.0)
    # two flies, three recordings: 2h over ma minus the solvent
    rows = []
    for fly, rec, k in (('a', 'r1', 1.0), ('a', 'r2', 1.0), ('b', 'r3', 2.0)):
        ps = t['per_stim'].reset_index().assign(fly=fly, recording=rec)
        ps.loc[ps['stim'] == '2h @ -6', 'minus_recording'] *= k
        rows.append(ps)
    ps = pd.concat(rows)
    r = M.flies_ratios(ps, [('2h', -6.0)], ('ma', -7.0), n_boot=100)
    # per fly 0.8 / 0.3 and 1.6 / 0.3, geometric mean over the two flies
    assert r['ratio'].iloc[0] == pytest.approx(np.sqrt((0.8 / 0.3) * (1.6 / 0.3)))
    assert r['n_flies'].iloc[0] == 2 and r['mean'].iloc[0] == 'geometric'


def test_measured_pair_grid_made_up():
    rng = np.random.default_rng(0)
    fr = M.shuffled_blocks(_frames())
    fr['response'] = 0.0
    mask = np.zeros((48, 48), bool)
    mask[4:44, 4:44] = True
    yy, xx = np.mgrid[:48, :48]
    blob_a = np.exp(-((yy - 12) ** 2 + (xx - 12) ** 2) / 50.0)
    blob_b = np.exp(-((yy - 36) ** 2 + (xx - 36) ** 2) / 50.0)
    # the solvent's map is the baseline 0.5; each odor adds 0.2 everywhere and its own blob
    stack = np.stack([0.5 + (0.2 + blob_a if s == '2h @ -6' else 0.2 + blob_b if s == 'ma @ -7' else 0)
                      + rng.normal(0, 0.01, (48, 48)) for s in fr['stim']])
    g = M.measured_pair_grid(stack, fr, mask, '2h @ -6', 'ma @ -7')
    G = g['grid'].reshape(12, 12)
    assert G[2, 2] > 0 > G[9, 9]                   # 2h's blob top left, ma's bottom right (image axes)
    assert g['split_half_r'] > 0.9
    assert g['mean_a'] == pytest.approx((0.2 + blob_a)[mask].mean(), abs=0.01)
    assert M.measured_pair_grid(stack, fr, mask, '2h @ -6', 'farn @ -2') is None


def test_score_recording_planted():
    st = _planted_state()
    g = st['grids']['flywire_R']
    Wt = M._wt(st['rates'], g['names'], ['A [s]', 'B [s]'])
    sl = g['slabs'][18.0]
    truth = M.pair_maps(Wt, g['perms'][:1], sl['G'][2.0], sl['inside'])[0]
    # a fly's right calyx, which Yang's files call "L" (the entry's orientation)
    d = M.score_recording(dict(grid=truth, inside=sl['inside']), 'L', 18.0, 'A - B', 's', st)
    assert len(d) == 8 and d.set_index('calyx').loc['flywire_R', 'exceeds']
    assert d['exceeds'].mean() > 0.7                # the other calyces share the layout


def test_w3_over_estimates():
    assert M.w3_over_estimates({'a': 'passes', 'b': 'passes'}) == 'passes'
    assert M.w3_over_estimates({'a': 'passes', 'b': 'not decided'}) == 'not decided'
    assert M.w3_over_estimates({'a': 'fails', 'b': 'fails', 'c': 'fails'}) == 'fails'
    assert M.w3_over_estimates({}) == 'not decided'


@pytest.mark.skipif(_shared() is None or not (_shared() / 'data' / 'yang_apl_recordings' / 'recordings.csv').exists(),
                    reason="viki's copy not there")
def test_readers_on_released_2021():
    recs = M.read_recordings(sets=None)
    assert recs['released'].sum() == 12 or M.YANG_RELEASED
    row = recs[recs['stem'] == '20210714_2FRfun1'].iloc[0]
    rec = M.read_recording(row, maps=True)
    fr = rec['frames']
    # the yaml reader reproduces the 2021 frames.csv (built from the same yaml)
    ref = M.read_stim_table().query("recording == '20210714_2FRfun1'").reset_index(drop=True)
    assert (ref['odor'].to_numpy() == fr['odor'].to_numpy()).all()
    assert np.allclose(ref['log10_conc'].fillna(0), fr['log10_conc'].fillna(0))
    assert (ref['block'].to_numpy() == fr['block'].to_numpy()).all()
    assert rec['xml']['pixel_um'] == pytest.approx(float(row['pixel_size_um']))
    # the pair map with the recording's solvent subtracted is E79's measured pa minus ms (all paraffin oil baseline)
    CM, _ = M._e79()
    meas = CM.measured_maps('all_pfo')['20210714_2FRfun1']
    c = meas['concs']['pa-ms']
    g = M.measured_pair_grid(rec['stack'], fr, rec['mask'], M.stim_label(('pa', c['pa'])), M.stim_label(('ms', c['ms'])))
    assert np.allclose(g['grid'], meas['grid']['pa-ms'], atol=1e-9)
    assert (g['inside'] == meas['inside']).all()
