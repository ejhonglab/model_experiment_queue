#!/usr/bin/env python3
"""Yang's request (2026-11-18, viki's terminal): the heatmaps of `yang_family_similarity.py` for the model with
breadth-weighted APL and without KC-KC inhibition.

The model is the best family (the breadth-2 APL-effect fit, C2 of E67, on E72's eight wirings) with its KC-KC
inhibition set to 0 and both response rates held at each wiring's fit: the somatic threshold is re-tuned with APL off to
the fit's number of APL-off responses, then the APL scales (APL->KC with KC->APL and bouton->APL tied) to the fit's
number of APL-on responses (`E77_input_screen.hold_rates`, breadth weights and APL->bouton as fitted). Removing KC-KC
alone, with nothing re-tuned, raises the response rate with APL on from 0.067 to 0.185 on the real wiring, which is why
the rates are held. APL off is the acute knockout at the held threshold. Read as `yang_family_similarity` reads the
family: cosine similarity on raw spike counts, the mean over the eight wirings of the distance, beside the KC data.

Usage (from the repository root):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u analyses/yang_family_no_kckc.py run real P1 M1 P2
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u analyses/yang_family_no_kckc.py run M2 P3 M3 P4
    ./venv/bin/python analyses/yang_family_no_kckc.py plot
Each wiring's run is saved as it finishes (spike counts and responses under outputs/yang_family_no_kckc/, its held
parameters and pair distances in <shared>/results/yang_family_no_kckc/), and a wiring already saved is skipped.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List

import pandas as pd

from kcq import harness as H

NAME = 'yang_family_no_kckc'


def results_dir() -> Path:
    return H.shared_dir() / 'results' / NAME


def runs_dir() -> Path:
    return H.outputs_dir(NAME)


def run_wiring(label: str) -> dict:
    """One wiring of the family with KC-KC at 0 and both response rates held; saves its runs, returns its row."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments.E77_input_screen import hold_rates
    t0 = time.time()
    P = E83.wiring_parts(label)
    p = dict(P['p'], wKCKC_scale=0.0)
    res = hold_rates(lambda: E83.build_simulator(label), target_off=P['target_off'], target_on=P['target_on'], params=p,
                     label=f'{NAME}:{label}')
    o = {'on': res['on'], 'off': res['off']}
    d = runs_dir()
    d.mkdir(parents=True, exist_ok=True)
    for c in ('on', 'off'):
        o[c][0].to_parquet(d / f'{label}_{c}_responses.parquet')
        o[c][1].to_parquet(d / f'{label}_{c}_spike_counts.parquet')
    dist = E83.pair_distances(o, 'raw')
    out = results_dir()
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dist).rename_axis(['odor_a', 'odor_b']).to_csv(out / f'pair_distances_{label}.csv')
    row = dict(wiring=label, folder=P['folder'], tag=P['tag'], fitted_wKCKC_scale=P['p']['wKCKC_scale'],
               target_on=res['target_on'], target_off=res['target_off'], n_on=res['n_on'], n_off=res['n_off'],
               size=P['size'], rate_on=res['n_on'] / P['size'], rate_off=res['n_off'] / P['size'],
               thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], at_bound_thr=res['at_bound_thr'],
               at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], params_on=res['params'],
               params_off=res['params_off'], seconds=time.time() - t0)
    (out / f'held_{label}.json').write_text(json.dumps(row, indent=1, default=float) + '\n')
    print(f"{label}: thr x{res['thr_factor']:.3f}, APL x{res['apl_factor']:.3f}, responses {res['n_on']} / {res['n_off']} "
          f"(targets {res['target_on']} / {res['target_off']}), {res['n_sims']} sims, {row['seconds'] / 60:.1f} min",
          flush=True)
    return row


def load_population(labels: List[str]) -> Dict[str, pd.Series]:
    """{'on', 'off'}: the mean over the wirings of their saved pair distances."""
    from experiments import E106_rho_gap_parts as E106
    ds = []
    for lab in labels:
        f = results_dir() / f'pair_distances_{lab}.csv'
        if not f.exists():
            raise FileNotFoundError(f'{lab} not run yet ({f})')
        df = pd.read_csv(f, index_col=[0, 1])
        ds.append({c: df[c] for c in ('on', 'off')})
    return E106.mean_pairs(ds)


def plot() -> None:
    from scipy.stats import spearmanr
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E106_rho_gap_parts as E106
    from analyses import yang_family_similarity as YS
    labels = list(E83.ALL_WIRINGS)
    t = YS.pairs_table(load_population(labels), E106.data_matrices())
    if len(t) != 120:
        raise RuntimeError(f'{len(t)} pairs with all six values, not 120')
    out = results_dir()
    t.rename_axis(['odor_a', 'odor_b']).to_csv(out / 'pairs.csv')
    YS.draw(t, out / 'similarity.png', name='no KC–KC',
            about='the best family without KC–KC inhibition (breadth-2 APL, KC–KC 0, both response rates held at the '
                  'fits\'; mean over the 8 wirings, raw spike counts)')
    held = pd.DataFrame([json.loads((out / f'held_{lab}.json').read_text()) for lab in labels])
    held.drop(columns=['params_on', 'params_off']).to_csv(out / 'held.csv', index=False)
    prov = H.model_provenance()
    prov.update(script=H.script_path(__file__), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), wirings=labels,
                change='wKCKC_scale = 0; fixed_thr re-tuned with APL off to the fit\'s APL-off responses, then the tied '
                       'APL scales to its APL-on responses (E77_input_screen.hold_rates)')
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    rho = lambda a, b: spearmanr(t[a], t[b]).statistic
    print(held[['wiring', 'fitted_wKCKC_scale', 'rate_on', 'rate_off', 'thr_factor', 'apl_factor']].round(4).to_string())
    print(f'ρ on {rho("model_on", "flies_on"):.3f}, off {rho("model_off", "flies_off"):.3f}, change '
          f'{rho("model_change", "flies_change"):.3f} (Pearson {t.model_change.corr(t.flies_change):.3f}); mean change '
          f'{t.model_change.mean():+.3f} (flies {t.flies_change.mean():+.3f})')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('wirings', nargs='+')
    sub.add_parser('plot')
    a = ap.parse_args()
    if a.cmd == 'run':
        for lab in a.wirings:
            if (results_dir() / f'held_{lab}.json').exists():
                print(f'{lab}: already saved, skipped', flush=True)
                continue
            run_wiring(lab)
    else:
        plot()


if __name__ == '__main__':
    main()
