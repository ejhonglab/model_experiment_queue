"""Odor-odor similarity of Remy's PN bouton recordings in the calyx (megamat, 17 odors): peak responses and a movie.

Reads Remy's per-ROI bouton files read-only (/local/tensor/Remy/odor_space_collab/analysis_outputs/data_by_acq/
pn_boutons/volumetric_without_neuropil_correction/<acq>/trials/xrds_bc_trials.nc) for the five flies her analyses
use. 2022-10-21 fly02 (megamat1_001__dsub3) is left out: her dashboards drop it, it matches none of their flies
and its median cellprob is 0.03 (viki's user, 2026-11-20: "don't analyze the one fly she seemed to drop").

For each fly, each odor's response vector is the mean over its 3 trials of each ROI's F_zscore (Remy's baseline:
the median of -5 to 0 s). Similarities between the 17 odors' vectors are cosine similarity and Pearson
correlation across ROIs, computed per fly and then averaged over flies.

Writes, to <out> (default: the shared directory's results/pn_bouton_odor_rdms/):
- peak_{cosine,pearson}_flymean.csv and peak_{cosine,pearson}_fly<k>.csv: the 0.5 to 2.5 s window mean (Remy's
  "mean_peak" window);
- peak_rdms.png: fly-mean cosine and Pearson side by side;
- movie_{cosine,pearson}.npy (frames x 17 x 17, fly mean) and times.npy, from each frame's trial-mean vectors;
- movie.mp4 and movie.gif: the two fly-mean matrices over -1 to 8 s, one 0.2 s frame per movie frame.

Usage (from the repository root; a few minutes, 2 threads are plenty):
    ./venv/bin/python analyses/pn_bouton_odor_rdms.py [--out DIR]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

SRC = Path('/local/tensor/Remy/odor_space_collab/analysis_outputs/data_by_acq/pn_boutons/'
           'volumetric_without_neuropil_correction')
EXCLUDED = ['2022-10-21__fly02__megamat1_001__dsub3']
# megamat odors, grouped as Remy's rdm_explorer orders them
ODORS = ['2h', 'IaA', 'pa', '2-but', 'eb', 'ep', 'aa', 'va', 'B-cit', 'Lin', '6al', 't2h', '1-8ol', '1-5ol',
         '1-6ol', 'benz', 'ms']
PEAK = (0.5, 2.5)
MOVIE = (-1.0, 8.0)


def odor_means(path: Path, var: str = 'F_zscore') -> tuple[np.ndarray, np.ndarray]:
    """(odors x ROIs x time) trial means of one acquisition, in ODORS order, and the frame times."""
    ds = xr.open_dataset(path)
    stim = np.array([str(s).split(' @')[0].strip() for s in ds['stim'].values])
    x = ds[var].values  # trials x cells x time
    out = np.stack([x[stim == o].mean(0) for o in ODORS])
    return out, ds['time'].values


def sims(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cosine similarity and Pearson correlation between the rows of v (odors x ROIs)."""
    n = v / np.linalg.norm(v, axis=1, keepdims=True)
    return n @ n.T, np.corrcoef(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='/local/tensor/Tom/claude_kc_remapping/results/pn_bouton_odor_rdms')
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    acqs = sorted(p.name for p in SRC.iterdir() if (p / 'trials' / 'xrds_bc_trials.nc').exists() and p.name not in EXCLUDED)
    assert len(acqs) == 5, acqs
    peak = {'cosine': [], 'pearson': []}
    mov = {'cosine': [], 'pearson': []}
    times = None
    for k, acq in enumerate(acqs, 1):
        m, t = odor_means(SRC / acq / 'trials' / 'xrds_bc_trials.nc')
        w = (t >= PEAK[0] - 1e-9) & (t <= PEAK[1] + 1e-9)
        c, p = sims(m[:, :, w].mean(2))
        for name, s in (('cosine', c), ('pearson', p)):
            peak[name].append(s)
            pd.DataFrame(s, index=ODORS, columns=ODORS).to_csv(out / f'peak_{name}_fly{k}.csv')
        f = (t >= MOVIE[0] - 1e-9) & (t <= MOVIE[1] + 1e-9)
        if times is None:
            times = t[f]
        assert np.allclose(times, t[f])
        cs, ps = zip(*(sims(m[:, :, i]) for i in np.where(f)[0]))
        mov['cosine'].append(np.stack(cs)); mov['pearson'].append(np.stack(ps))
    for name in peak:
        pd.DataFrame(np.mean(peak[name], 0), index=ODORS, columns=ODORS).to_csv(out / f'peak_{name}_flymean.csv')
        np.save(out / f'movie_{name}.npy', np.mean(mov[name], 0))
    np.save(out / 'times.npy', times)
    (out / 'flies.json').write_text(json.dumps({'flies': acqs, 'excluded': EXCLUDED, 'var': 'F_zscore',
                                                'peak_window_s': PEAK, 'source': str(SRC)}, indent=1))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from kcq.harness import oklab_diverging_cmap
    cmap = oklab_diverging_cmap()

    def panel(ax, s, title):
        im = ax.imshow(s, cmap=cmap, vmin=-1, vmax=1)
        ax.set_xticks(range(17)); ax.set_xticklabels(ODORS, rotation=90, fontsize=7)
        ax.set_yticks(range(17)); ax.set_yticklabels(ODORS, fontsize=7)
        ax.set_title(title, fontsize=9)
        return im

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
    for ax, name, lab in zip(axs, ('cosine', 'pearson'), ('cosine similarity', 'Pearson correlation')):
        im = panel(ax, np.mean(peak[name], 0), f'PN boutons, {lab}\n(mean over 5 flies, {PEAK[0]}-{PEAK[1]} s)')
    fig.colorbar(im, ax=axs, shrink=0.8, label='similarity across bouton ROIs')
    fig.savefig(out / 'peak_rdms.png', dpi=150); plt.close(fig)

    M = {name: np.mean(mov[name], 0) for name in mov}
    fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    ims = [panel(ax, M[name][0], '') for ax, name in zip(axs, ('cosine', 'pearson'))]
    fig.colorbar(ims[0], ax=axs, shrink=0.8, label='similarity across bouton ROIs')

    def frame(i):
        for im, ax, name, lab in zip(ims, axs, ('cosine', 'pearson'), ('cosine similarity', 'Pearson correlation')):
            im.set_data(M[name][i])
            ax.set_title(f'PN boutons, {lab} (5 flies)\nt = {times[i]:+.1f} s from odor onset', fontsize=9)
        return ims
    anim = animation.FuncAnimation(fig, frame, frames=len(times), blit=False)
    anim.save(out / 'movie.mp4', writer=animation.FFMpegWriter(fps=5), dpi=120)
    anim.save(out / 'movie.gif', writer=animation.PillowWriter(fps=5), dpi=80)
    plt.close(fig)
    print('wrote', out, len(times), 'frames')


if __name__ == '__main__':
    main()
