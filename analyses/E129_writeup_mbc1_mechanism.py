"""Write-up numbers for E129's step 3 (MB-C1's counted synapses in the KC-KC matrix's place with APL off), from the
model-output archive and pyolfsysm's MB-C1 loader; no simulation.

The numbers the note gives for why MB-C1's arm fails, which until 2026-11-23 were computed by hand:
- MB-C1's synapses with the model's KCs and the share of them on the wiring's claws, with each denominator
  (`mbc1_synapses`): MB-C1's calyx synapses in all ('all_calyx'), those with the 1732 KCs of the wiring ('wiring_kcs';
  the hemibrain claw x bouton wiring, `connectome_wPNKC(prat_claws=True, prat_boutons=True)`) and those of them on the
  wiring's claws ('wiring_claws', `connectome_MBC1_weights(prat_claws=True, wPNKC=...)`), with what is off the claws
  (claw -1, other claws); 'out' is MB-C1 -> KC, 'in' KC -> MB-C1, the two MB-C1_R cells summed;
- the model KCs with and without MB-C1 synapses (`n_kcs`), from the per-KC counts saved beside each step-3 run
  (<tag>_mbc1_counts.parquet, checked equal on all eight wirings and to `connectome_MBC1_weights()`);
- per wiring and condition, each arm's response rate (KC x odor entries with two or more spikes, the model side of the
  realism bands) in the KCs without output from MB-C1, in the rest, per KC type and over all KCs (`response_rate`),
  and for MB-C1's arm the share of all responses carried by the KCs without output from MB-C1
  (`share_of_responses`);
- per wiring and condition, the shares of KCs responding to 0, 1, 2-3, 4-7 and 8 or more of the 16 odors
  (`odors_per_kc`), and their total variation distance from the flies' (`odors_per_kc_tv_from_flies`; the flies:
  control flies for APL on, APL output blocked (TNT) for APL off, the registry's noise_ctrl_ / noise_tnt_odors_per_kc_*,
  the bands' noise-SD criterion);
- MB-C1's arm's calibrated pattern part minus the family's and minus KC-KC at zero's (E110 part A), against all
  recordings and each half A to D (`pattern_part_minus_family`, `pattern_part_minus_E110_A`; group the target), per
  wiring and for the populations' readings (wiring 'S4' searched, 'T4' held out, 'all8'; the loss of the population's
  mean distances, as readings.csv and the classes read them), from readings.csv.

Arms: 'family' (the family's saved runs, E67's on the real wiring and E72's seven re-draws, the (folder, tag) of
runs.csv's family rows), 'P' and 'S' (step 1a's pooled neuron and matrix at 1 s, S1_P_<w>, S1_S_<w>), 'MB' (step 3,
S3_MB_<w>), and 'flies' for the data rows. Per-wiring rows are followed by 'mean8', 'min8' and 'max8' over the eight
wirings (the KCs are the same 1732 on every wiring, so a mean of rates over wirings is the pooled rate).

Checks (raise on failure): every archived responses table equals its spike counts at two or more spikes; the arms'
KCs and odors are the same on a wiring; the per-KC MB-C1 counts are the same on all eight wirings and equal the
loader's; each arm's per-type response rates and odors-per-KC distances equal runs.csv's (rate_<cond>_<type>,
opk_<cond>_tv), which the score computed from the same runs; the pattern-part margins equal verdict.json's (MB's
d_family_per_wiring, and its class targets' arm minus E110_A per population).

Writes writeup_mbc1_mechanism.csv (columns stat, arm, wiring, cond, group, value) and
provenance_writeup_mbc1_mechanism.json (inputs with sha256, definitions, the checks) to the results folder.

Usage (from the repository root): ./venv/bin/python analyses/E129_writeup_mbc1_mechanism.py [--out DIR]
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E129_block_survivor_static'
WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
ARMS = ('family', 'P', 'S', 'MB')
ARM_FILES = {'P': ('step1a', 'S1_P_{w}'), 'S': ('step1a', 'S1_S_{w}'), 'MB': ('step3', 'S3_MB_{w}')}
CONDS = {'on': 'on', 'off': 'acute_off'}
TYPES = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
OPK_BINS = ('0', '1', '2_3', '4_7', '8plus')
RESPONSE_MIN_SPIKES = 2
TOL = 1e-9
DEFINITIONS = {
    'mbc1_synapses': "MB-C1's calyx synapses (both MB-C1_R cells summed): group <dir>_<set>, dir 'out' (MB-C1 -> KC) or "
        "'in' (KC -> MB-C1); set 'all_calyx' (all of them), 'wiring_kcs' (with the 1732 KCs of the hemibrain claw x "
        "bouton wiring), 'wiring_claws' (of those, on the wiring's claws), 'off_claws_unassigned' (claw -1) and "
        "'off_claws_other' (claws the wiring does not have)",
    'mbc1_claw_share': "the share of MB-C1's synapses on the wiring's claws: group <dir>_of_wiring_kcs (denominator its "
        "calyx synapses with the wiring's KCs) or <dir>_of_all_calyx (denominator all of its calyx synapses)",
    'n_kcs': "model KCs: 'all', 'no_mbc1_output' (no MB-C1 -> KC synapse; also by type, _g / _apbp / _ab), "
        "'mbc1_output', 'any_mbc1_synapse' (either direction), 'both_directions'",
    'response_rate': "fraction of KC x odor entries with two or more spikes, in the KCs of 'group' ('all', "
        "'no_mbc1_output', 'mbc1_output', or a KC type 'g' / 'apbp' / 'ab')",
    'share_of_responses': "the share of the arm's responses (KC x odor entries with two or more spikes) carried by the "
        "KCs of 'group'",
    'odors_per_kc': "fraction of KCs responding (two or more spikes; the flies: the bands' noise-SD criterion) to 0, 1, "
        "2-3, 4-7 or 8 or more of the 16 odors (group '0', '1', '2_3', '4_7', '8plus')",
    'odors_per_kc_tv_from_flies': "total variation distance of the arm's odors-per-KC distribution from the flies' "
        "(control flies for 'on', APL output blocked (TNT) for 'off'), as runs.csv's opk_<cond>_tv",
    'pattern_part_minus_family': "MB-C1's arm's calibrated pattern part minus the family's (positive: MB worse), "
        "against all recordings ('all') or one half ('A' to 'D'); per wiring, or the populations' readings ('S4', 'T4', "
        "'all8')",
    'pattern_part_minus_E110_A': "the same, minus KC-KC at zero's (E110 part A)",
}
PATTERN_TARGETS = ('all', 'A', 'B', 'C', 'D')
POPULATIONS = ('S4', 'T4', 'all8')


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def archive_root() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs'


def family_members(runs: pd.DataFrame) -> Dict[str, Tuple[str, str]]:
    """{wiring: (archive folder, tag)} of the family's saved runs, from runs.csv's family rows."""
    fam = runs[runs.kind == 'family'].drop_duplicates('wiring').set_index('wiring')
    missing = set(WIRINGS) - set(fam.index)
    if missing:
        raise RuntimeError(f"runs.csv has no family row for {sorted(missing)}")
    return {w: (str(fam.loc[w, 'folder']), str(fam.loc[w, 'tag'])) for w in WIRINGS}


def run_files(arm: str, w: str, family: Dict[str, Tuple[str, str]]) -> Dict[str, Dict[str, Path]]:
    """{cond: {'responses', 'spike_counts'}: archive path} of one arm's run on wiring `w`."""
    if arm == 'family':
        folder, tag = family[w]
        d = archive_root() / folder
    else:
        step, pat = ARM_FILES[arm]
        d, tag = archive_root() / EXPERIMENT / step, pat.format(w=w)
    return {c: {k: d / f'{tag}_{f}_{k}.parquet' for k in ('responses', 'spike_counts')} for c, f in CONDS.items()}


def load_run(files: Dict[str, Dict[str, Path]]) -> Dict[str, pd.DataFrame]:
    """{cond: boolean KCs x odors responses}, from the spike counts at two or more spikes, checked equal to the
    archived responses table."""
    out = {}
    for c, f in files.items():
        counts, resp = pd.read_parquet(f['spike_counts']), pd.read_parquet(f['responses'])
        r = counts.to_numpy(dtype=float) >= RESPONSE_MIN_SPIKES
        if not (counts.index.equals(resp.index) and counts.columns.equals(resp.columns)
                and np.array_equal(r, resp.to_numpy(dtype=float) > 0.5)):
            raise RuntimeError(f"{f['responses']} is not {f['spike_counts']} at {RESPONSE_MIN_SPIKES} or more spikes")
        out[c] = pd.DataFrame(r, index=counts.index, columns=counts.columns)
    return out


def opk_fractions(resp: pd.DataFrame) -> np.ndarray:
    """Fractions of KCs responding to 0, 1, 2-3, 4-7 and 8 or more odors (`OPK_BINS`)."""
    n = resp.to_numpy().sum(axis=1)
    return np.array([np.mean(n == 0), np.mean(n == 1), np.mean((n >= 2) & (n <= 3)), np.mean((n >= 4) & (n <= 7)),
                     np.mean(n >= 8)])


def flies_opk(cond: str) -> np.ndarray:
    geno = 'ctrl' if cond == 'on' else 'tnt'
    return np.array([H.data_target(f'noise_{geno}_odors_per_kc_{b}') for b in OPK_BINS])


def mbc1_synapse_rows() -> Tuple[List[dict], pd.Series, pd.Series]:
    """MB-C1's synapse totals and claw shares (`mbc1_synapses`, `mbc1_claw_share`), and its per-KC counts (out, in) for
    the KCs with any calyx synapse with it."""
    from pyolfsysm.mb_model import connectome_MBC1_weights, connectome_wPNKC, KC_ID
    wiring = connectome_wPNKC(prat_claws=True, prat_boutons=True)
    ids = pd.Index(wiring.index.get_level_values(KC_ID).unique())
    all_out, all_in, _, _ = connectome_MBC1_weights()
    dropped: dict = {}
    claw_out, claw_in, _, _ = connectome_MBC1_weights(prat_claws=True, wPNKC=wiring, dropped_counts=dropped)
    rows = []
    for d, a, c, key in (('out', all_out, claw_out, 'wMBC1KC'), ('in', all_in, claw_in, 'wKCMBC1')):
        tot_all, tot_kcs, tot_claws = float(a.sum()), float(a.reindex(ids).fillna(0).sum()), float(c.sum())
        dr = dropped[key].sum()
        if not np.isclose(tot_claws + float(dr.sum()), tot_kcs):
            raise RuntimeError(f'{d}: on claws {tot_claws} plus dropped {dict(dr)} is not {tot_kcs}')
        for g, v in (('all_calyx', tot_all), ('wiring_kcs', tot_kcs), ('wiring_claws', tot_claws),
                     ('off_claws_unassigned', float(dr['unassigned_claw'])), ('off_claws_other', float(dr['other_claw'])),
                     ('off_claws_nan', float(dr['nan_claw']))):
            rows.append(dict(stat='mbc1_synapses', arm='MB-C1', wiring='-', cond='-', group=f'{d}_{g}', value=v))
        rows.append(dict(stat='mbc1_claw_share', arm='MB-C1', wiring='-', cond='-', group=f'{d}_of_wiring_kcs',
                         value=tot_claws / tot_kcs))
        rows.append(dict(stat='mbc1_claw_share', arm='MB-C1', wiring='-', cond='-', group=f'{d}_of_all_calyx',
                         value=tot_claws / tot_all))
    rows.append(dict(stat='n_kcs', arm='MB-C1', wiring='-', cond='-', group='wiring', value=float(len(ids))))
    return rows, all_out, all_in


def mbc1_counts(w: str) -> Tuple[pd.DataFrame, Path]:
    p = archive_root() / EXPERIMENT / 'step3' / f'S3_MB_{w}_mbc1_counts.parquet'
    return pd.read_parquet(p), p


def margin_rows(readings: pd.DataFrame, verdict: dict) -> Tuple[List[dict], float]:
    """MB's calibrated pattern part minus the family's and KC-KC at zero's, per wiring and population and target
    (`pattern_part_minus_family`, `pattern_part_minus_E110_A`), from readings.csv; returns the rows and the largest
    difference from verdict.json's values."""
    rd = readings[readings.readout == 'cal'].set_index(['kind', 'unit'])
    rows, diff = [], 0.0
    c = verdict['classes']['MB']['targets']
    for ref in ('family', 'E110_A'):
        for unit in WIRINGS + POPULATIONS:
            for t in PATTERN_TARGETS:
                v = float(rd.loc[('MB', unit), f'{t}_pattern'] - rd.loc[(ref, unit), f'{t}_pattern'])
                rows.append(dict(stat=f'pattern_part_minus_{ref}', arm='MB', wiring=unit, cond='-', group=t, value=v))
                if ref == 'family' and t == 'all' and unit in WIRINGS:
                    diff = max(diff, abs(v - float(c['all']['d_family_per_wiring'][unit])))
                if unit in ('S4', 'T4'):
                    p = c[t]['populations'][unit]
                    diff = max(diff, abs(v - (float(p['arm']) - float(p[ref]))))
    return rows, diff


def compute(runs: pd.DataFrame, readings: Optional[pd.DataFrame] = None,
            verdict: Optional[dict] = None) -> Tuple[pd.DataFrame, dict]:
    """The table (rows: stat, arm, wiring, cond, group, value) and the provenance's inputs and checks; the pattern-part
    margins only when `readings` and `verdict` are passed."""
    rows, all_out, all_in = mbc1_synapse_rows()
    family = family_members(runs)
    inputs: Dict[str, str] = {}

    counts, ref_path = mbc1_counts(WIRINGS[0])
    for w in WIRINGS:
        c, p = mbc1_counts(w)
        inputs[str(p.relative_to(H.shared_dir()))] = _sha256(p)
        if not c.equals(counts):
            raise RuntimeError(f"MB-C1's per-KC counts differ between {ref_path.name} and {p.name}")
    ids = counts.index.get_level_values(H.KC_ID_LEVEL)
    for col, w_ in (('from_mbc1', all_out), ('to_mbc1', all_in)):
        if not np.array_equal(counts[col].to_numpy(dtype=float), w_.reindex(ids).fillna(0.0).to_numpy(dtype=float)):
            raise RuntimeError(f"the saved per-KC counts' {col} are not connectome_MBC1_weights()'s")
    types = np.asarray(counts.index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)
    no_out = counts['from_mbc1'].to_numpy() == 0
    has_in = counts['to_mbc1'].to_numpy() > 0
    groups = {'all': np.ones(len(counts), bool), 'no_mbc1_output': no_out, 'mbc1_output': ~no_out}
    groups.update({TYPES[t]: types == t for t in TYPES})
    kc_rows = [('all', len(counts)), ('no_mbc1_output', no_out.sum()), ('mbc1_output', (~no_out).sum()),
               ('any_mbc1_synapse', ((~no_out) | has_in).sum()), ('both_directions', ((~no_out) & has_in).sum())]
    kc_rows += [(f'no_mbc1_output_{TYPES[t]}', (no_out & (types == t)).sum()) for t in TYPES]
    rows += [dict(stat='n_kcs', arm='MB-C1', wiring='-', cond='-', group=g, value=float(n)) for g, n in kc_rows]

    checks = dict(max_abs_diff_type_rate_vs_runs_csv=0.0, max_abs_diff_opk_tv_vs_runs_csv=0.0)
    per_wiring: List[dict] = []
    for arm in ARMS:
        for w in WIRINGS:
            files = run_files(arm, w, family)
            for f in files.values():
                for p in f.values():
                    inputs[str(p.relative_to(H.shared_dir()))] = _sha256(p)
            R = load_run(files)
            run_row = runs[(runs.kind == arm) & (runs.wiring == w)]
            if len(run_row) < 1:
                raise RuntimeError(f'runs.csv has no row for {arm} on {w}')
            run_row = run_row.iloc[0]
            for cond, r in R.items():
                if not r.index.equals(counts.index):
                    raise RuntimeError(f"{arm} on {w}, {cond}: the KCs are not those of MB-C1's per-KC counts")
                a = r.to_numpy()
                for g, m in groups.items():
                    per_wiring.append(dict(stat='response_rate', arm=arm, wiring=w, cond=cond, group=g,
                                           value=float(a[m].mean())))
                for t, col in TYPES.items():
                    d = abs(float(a[groups[col]].mean()) - float(run_row[f'rate_{cond}_{col}']))
                    checks['max_abs_diff_type_rate_vs_runs_csv'] = max(checks['max_abs_diff_type_rate_vs_runs_csv'], d)
                per_wiring.append(dict(stat='share_of_responses', arm=arm, wiring=w, cond=cond, group='no_mbc1_output',
                                       value=float(a[no_out].sum() / a.sum())))
                f_ = opk_fractions(r)
                for b, v in zip(OPK_BINS, f_):
                    per_wiring.append(dict(stat='odors_per_kc', arm=arm, wiring=w, cond=cond, group=b, value=float(v)))
                tv = 0.5 * float(np.abs(f_ - flies_opk(cond)).sum())
                per_wiring.append(dict(stat='odors_per_kc_tv_from_flies', arm=arm, wiring=w, cond=cond, group='-',
                                       value=tv))
                checks['max_abs_diff_opk_tv_vs_runs_csv'] = max(checks['max_abs_diff_opk_tv_vs_runs_csv'],
                                                                abs(tv - float(run_row[f'opk_{cond}_tv'])))
    for k, v in checks.items():
        if v > TOL:
            raise RuntimeError(f'{k} = {v}: the archived runs do not reproduce runs.csv')

    pw = pd.DataFrame(per_wiring)
    agg = pw.groupby(['stat', 'arm', 'cond', 'group'], sort=False)['value'].agg(['mean', 'min', 'max']).reset_index()
    agg_rows = [dict(stat=r.stat, arm=r.arm, wiring=f'{s}8', cond=r.cond, group=r.group, value=float(getattr(r, s)))
                for r in agg.itertuples() for s in ('mean', 'min', 'max')]
    for cond in CONDS:
        for b, v in zip(OPK_BINS, flies_opk(cond)):
            rows.append(dict(stat='odors_per_kc', arm='flies', wiring='-', cond=cond, group=b, value=float(v)))
    if readings is not None and verdict is not None:
        m_rows, m_diff = margin_rows(readings, verdict)
        checks['max_abs_diff_margins_vs_verdict'] = m_diff
        if m_diff > TOL:
            raise RuntimeError(f'the pattern-part margins differ from verdict.json by up to {m_diff}')
        m = pd.DataFrame(m_rows)
        on_w = m[m.wiring.isin(WIRINGS)]
        agg_rows += [dict(stat=r.stat, arm=r.arm, wiring=f'{s}8', cond=r.cond, group=r.group, value=float(getattr(r, s)))
                     for r in on_w.groupby(['stat', 'arm', 'cond', 'group'], sort=False)['value']
                     .agg(['min', 'max']).reset_index().itertuples() for s in ('min', 'max')]
        rows += m_rows
    table = pd.concat([pd.DataFrame(rows), pw, pd.DataFrame(agg_rows)], ignore_index=True)
    return table, dict(inputs=inputs, checks=checks, family_members={w: list(v) for w, v in family.items()})


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None, help='folder for the table (default: the results folder, which '
                                                            'also holds runs.csv)')
    args = ap.parse_args(argv)
    res = H.results_dir(EXPERIMENT)
    out = args.out or res
    out.mkdir(parents=True, exist_ok=True)
    paths = {k: res / f for k, f in (('runs', 'runs.csv'), ('readings', 'readings.csv'), ('verdict', 'verdict.json'))}
    verdict = json.loads(paths['verdict'].read_text())
    table, info = compute(pd.read_csv(paths['runs']), pd.read_csv(paths['readings']), verdict)
    csv = out / 'writeup_mbc1_mechanism.csv'
    table.to_csv(csv, index=False)
    for p in paths.values():
        info['inputs'][str(p.relative_to(H.shared_dir()))] = _sha256(p)
    prov = dict(script=H.script_path(__file__), **H.code_commits(), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
                verdict_generated=verdict.get('generated'),
                response_min_spikes=RESPONSE_MIN_SPIKES,
                mbc1_source=("pyolfsysm.mb_model.connectome_MBC1_weights() and, on the hemibrain claw x bouton wiring "
                             "connectome_wPNKC(prat_claws=True, prat_boutons=True), connectome_MBC1_weights("
                             "prat_claws=True, wPNKC=wiring, dropped_counts=...)"),
                definitions=DEFINITIONS, outputs=[csv.name], **info)
    (out / 'provenance_writeup_mbc1_mechanism.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')

    show = table[table.wiring.isin(['-', 'mean8', 'min8', 'max8'])]
    with pd.option_context('display.max_rows', 400, 'display.width', 160):
        print(show.pivot_table(index=['stat', 'arm', 'cond', 'group'], columns='wiring', values='value',
                               sort=False).round(4))
        mg = table[table.stat.str.startswith('pattern_part_minus')]
        print(mg.pivot_table(index=['stat', 'wiring'], columns='group', values='value', sort=False).round(4))
    print(f'wrote {csv}')


if __name__ == '__main__':
    main()
