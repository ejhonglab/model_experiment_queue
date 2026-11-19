#!/usr/bin/env python3
"""Yang's request (2026-11-19, viki's terminal): the heatmaps of `yang_family_similarity.py` for the best performing
model that has neither KC-KC inhibition nor breadth-tuned inhibition.

Candidates are the saved variants scored in <shared>/results/apl_effect/by_variant.csv (`kcq/apl_effect.py`). A
variant counts as plain APL (`classify`) when the archive's params.csv (<shared>/data/model_outputs/<folder>/) records
KC-KC at 0 (kc_kc_inh, wKCKC_scale, kckc, wkckc; non-numeric mode labels ignored) and the breadth exponent at 0 (k,
breadth_k), or when it belongs to an early experiment on the claw or bouton bases that never sets either (`PLAIN_EARLY`,
checked by reading each script). APL weightings that target chosen KCs or claws are left out with the breadth weighting:
E11's co-input claws, E14's oracle co-responders, E15's and E17's weight rules (each keeps only its base variant). E38's
inputs other than the lab's are left out, so that every candidate has the input that the other figures have. A variant
whose parameters are not recorded either way is left out and listed in the output.

"Best performing" is the page's rule for its best-models table and best_by_experiment.csv: the lowest APL-effect loss
(spike counts, cosine) among variants without a realism flag and outside the held wirings. The candidates table also
carries the pattern part and the three ρ's, so that the ranking on each can be read.

Outputs: <shared>/results/yang_best_plain_apl/ (candidates.csv: every plain-APL variant with its scores; unclassified.csv;
similarity.png and pairs.csv for the winner; provenance.json).
Usage: ./venv/bin/python analyses/yang_best_plain_apl.py
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from kcq import harness as H

NAME = 'yang_best_plain_apl'
KCKC_COLS = ('kc_kc_inh', 'wKCKC_scale', 'kckc', 'wkckc')
BREADTH_COLS = ('k', 'breadth_k')
# experiments on the claw / bouton bases whose scripts never set KC-KC or a breadth weighting (grep of each, 2026-11-19)
PLAIN_EARLY = {'E1_apl_knockout', 'E3_drive_dependent_gain', 'E4_substrate', 'E5_kc_heterogeneity',
               'E9_centrifugal_input', 'E11_coincidence_claw_inhibition', 'E13_bouton_coincidence',
               'E14_oracle_coresponders', 'E16_claw_saturation', 'E18_bouton_depletion', 'E24_input_scaling',
               'E26_claw_power', 'E27_glomerulus_gain', 'E32_pn_stage', 'E35_divisive_apl'}
C = 'spike_counts_cosine__'
SCORES = {C + 'apl_effect_rmse': 'loss', C + 'apl_effect_rmse_centered': 'pattern', C + 'apl_effect_mean': 'mean_effect',
          C + 'apl_effect_spearman': 'rho_change', C + 'dist_rho_on_all120': 'rho_on', C + 'dist_rho_off_all120': 'rho_off'}


def _max_numeric(row: pd.Series, cols) -> float:
    vals = [pd.to_numeric(row[c], errors='coerce') for c in cols if c in row.index]
    vals = [abs(float(v)) for v in vals if pd.notna(v)]
    return max(vals) if vals else np.nan


def archive_params() -> pd.DataFrame:
    """One row per (folder, tag) of the archive: the largest KC-KC and breadth values recorded (NaN: none)."""
    rows = []
    for f in sorted((H.shared_dir() / 'data' / 'model_outputs').glob('E*/params.csv')):
        p = pd.read_csv(f, low_memory=False)
        if 'tag' not in p.columns:
            continue
        for _, r in p.iterrows():
            rows.append(dict(folder=f.parent.name, tag=r['tag'], kckc=_max_numeric(r, KCKC_COLS),
                             breadth=_max_numeric(r, BREADTH_COLS)))
    return pd.DataFrame(rows).drop_duplicates(['folder', 'tag'])


def classify(folder: str, tag: str, kckc: float, breadth: float) -> Optional[bool]:
    """True: plain APL (no KC-KC, no breadth or targeted weighting); False: not; None: not recorded either way."""
    t = str(tag)
    if folder in ('E14_oracle_coresponders', 'E15_breadth_shared_apl', 'E17_panel_general_apl'):
        return t.endswith('_base')
    if folder == 'E11_coincidence_claw_inhibition':
        return 'coinput' not in t
    if folder == 'E12_breadth_apl':
        return t == 'k0'
    if folder == 'E20_glomerulus_gain' and re.fullmatch(r'a[\d.]+_thr\d+', t):
        return True       # part A: the bouton base with the glomerulus gain and the per-claw threshold
    if folder == 'E22_soft_knee' and re.fullmatch(r'thr\d+_w\d+', t):
        return True       # the soft knee on the bouton base (the --combo runs carry _k2_kckc100)
    if folder == 'E41_claw_subgain' and t.startswith('A_'):
        return True       # part A: the bouton base with the per-claw sub-threshold gain
    if folder == 'E38_hallem_input' and t.startswith('A_') and t != 'A_lab':
        return False      # another input (part B is the reference model, with KC-KC and breadth 2)
    if folder == 'E97_glomerulus_claw_shares':
        return False      # the best family's fits on other claw shares (KC-KC, breadth 2)
    if (not np.isnan(kckc) and kckc > 0) or (not np.isnan(breadth) and breadth > 0) or re.search(r'kckc0*[1-9]', t):
        return False
    if folder in PLAIN_EARLY and 'breadth' not in t and 'kckc' not in t:
        return True
    if np.isnan(kckc) or np.isnan(breadth):
        return None
    return True


def main() -> None:
    from analyses import yang_family_similarity as YS
    from experiments import E106_rho_gap_parts as E106
    t0 = time.time()
    bv_path = H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv'
    bv = pd.read_csv(bv_path)
    m = bv.merge(archive_params(), on=['folder', 'tag'], how='left')
    m['plain'] = [classify(r.folder, r.tag, r.kckc, r.breadth) for r in m.itertuples()]
    out = H.shared_dir() / 'results' / NAME
    out.mkdir(parents=True, exist_ok=True)
    m[m.plain.isna()][['folder', 'tag', 'kckc', 'breadth']].to_csv(out / 'unclassified.csv', index=False)
    p = m[m.plain == True].rename(columns=SCORES)  # noqa: E712
    keep = ['folder', 'tag', *SCORES.values(), 'rate_on', 'rate_off', 'realism_flag']
    p = p[keep].sort_values('loss')
    p.to_csv(out / 'candidates.csv', index=False)
    ok = p[p.realism_flag.isna() & p.rho_on.notna()]
    best = ok.iloc[0]
    print(f'{len(p)} plain-APL variants, {len(ok)} without a realism flag; {int(m.plain.isna().sum())} unclassified')
    print(ok.head(6).drop(columns='realism_flag').round(4).to_string(index=False))
    model = E106.load_pairs(best.folder, best.tag)
    t = YS.pairs_table(model, E106.data_matrices())
    if len(t) != 120:
        raise RuntimeError(f'{len(t)} pairs with all six values, not 120')
    t.rename_axis(['odor_a', 'odor_b']).to_csv(out / 'pairs.csv')
    YS.draw(t, out / 'similarity.png', name='best plain APL',
            about=f'the best model without KC–KC or breadth-tuned inhibition by the APL-effect loss '
                  f'({best.folder}:{best.tag}, raw spike counts)')
    prov = H.model_provenance()
    prov.update(script=H.script_path(__file__), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), seconds=time.time() - t0,
                by_variant=str(bv_path), winner=dict(folder=best.folder, tag=best.tag),
                rule='lowest spike_counts_cosine__apl_effect_rmse among plain-APL variants without realism_flag and '
                     'with dist_rho_on_all120 (not held)')
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    print(f'winner {best.folder}:{best.tag} -> {out}')


if __name__ == '__main__':
    main()
