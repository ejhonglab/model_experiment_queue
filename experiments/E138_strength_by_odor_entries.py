#!/usr/bin/env python3
"""E138, the odors that carry the strength error and the glomeruli that feed them: the strength view's entries odor by
odor at the KC somata of every panel, on the flies' halves, along the megamat stages, and the glomeruli behind the
odors that carry it (no simulation).

The entry is item E138 of experiment_priorities.md (proposed by viki 2026-11-24 22:54 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-25 02:08 under the idle-host rule; claimed by atlas 02:11).
"The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings (E83's), each read as one fly; the
views, one-fly units, readings and flies are E134's (kcq/population_fit.py, kcq/stage_readings.py,
data_prep/stage_population_data.py and its lock experiments/E134_lock.json). The strength view ("strength by odor") has
two parts, each odor's response rate and its mean response as multiples of the panel's mean (part 'rate' and part
'mean'); a part's error is the root mean square over odors of the model's entry minus the flies' mean entry, in the
flies' one-fly units, and the view's error the root mean square of its two parts' errors. Nothing is run, fitted or
searched: E134's, E107's and E90's saved runs are read again.

The guard: no residual (model minus flies) is computed before this experiment's classes lock exists, is committed with
no uncommitted change and equals its copy in the shared results (`require_classes_lock`); and E134's lock guards every
reading of the flies (`stage_readings.require_lock`). Every part but 'classes' calls both first.

Parts (CLI `--part`):
- classes (design item 1; no residual): the three classes of glomeruli, written to results/E138_strength_by_odor_entries/
  classes_lock.json and its committed copy experiments/E138_classes_lock.json (then commit it), with the rules of the
  classes read here ("concentrated", (d)) and the feeders' rule, fixed before any residual:
  - 'ir': the glomeruli fed by ORNs whose ligand-specific receptor is an ionotropic receptor (IR), from Task et al. 2022
    (eLife 11:e72599), the table of receptors and co-receptors by glomerulus as drosolf ships it (data/task22_table3.csv;
    the same table as al_analysis's data/sent_to_nikola/2025-02-25/task22_table3.csv), whose IR glomeruli are those of
    Silbering et al. 2011 (J Neurosci 31:13357); the co-receptors (Ir8a, Ir76b, Ir25a) beside, not a class;
  - 'not_dendrite28': the model's input glomeruli absent from the 28 glomeruli of the PN-dendrite recordings (E135's
    lock, `glomeruli`; E90's `shared_glomeruli`);
  - 'few_hemibrain': the model's input glomeruli with fewer boutons, or fewer claws, than the median input glomerulus in
    the family's hemibrain wiring (E130's `step0_wiring_only('hemibrain-first')`, the claw x bouton matrix of E83's
    'real'), each sub-list beside the union.
- entries (design items 2 and 7): per panel and condition, each odor's entry in both strength parts for the family as
  fitted and as held and for the beside models: the model's entry (each wiring, the wirings' mean views 'common'),
  the flies' mean entry and range over flies, the residual in one-fly units, each odor's share of each part's squared
  error and of the view's, and the view's error without that odor's entry; the strength split (E136's
  `strength_split`: level, spread and order) beside every error. Panels (`PANEL_SPECS`): megamat (E134's step 3,
  17 odors), the APL16 panel with APL on against the control flies and APL off against APL output blocked (TNT)
  (E90's part B, the family's own fit panel, and E90's D1 fits beside), validation2 (E107's part B for the family, E134's
  step 3 for the beside models), and the kiwi and control mixtures (E134's step 3) only once the held-out panels read
  like for like (E137) has its lock committed (`require_e137_lock`; until then the mixture rows are skipped and the
  outputs say so). The readings of the flies (`sides_of`; column 'reading'): 'e134', E134's (the model read against
  all the panel's recorded flies, the reference those of E134's `fly_side_read`), which decides the classes on the
  KC-soma panels (`MAIN_READING`); for validation2 'fly_scale' beside it (hal's `partial_reference` with rule
  'fly_scale' on all 20 flies and 22 odors, each fly on the odors it was given, 338684e); for the mixtures each of
  E137's five readings of the flies (its `reading_population`, the model read as its re-score reads it,
  `read_model_panel`: the curve on the reading's flies, their trial noise added), 'as_uploaded' (olfdb dataset 8 as
  uploaded, E134's reading) deciding and the four others beside.
- halves (design items 3 and 7): the same on the two seeded halves of each panel's flies (E136's `split_names`: hal's
  `split_halves(n, 1, 0)`, the first part the 'fit' half; the APL16 panel's control half paired with the TNT half of the
  same name, as the E134 lock's split_halves rule pairs them), each half read as the lock's split_halves rule says
  (curve, noise, reference and one-fly units from the half), on the odors of the all-flies reading (validation2's
  E134 reading splits its six flies, the 'fly_scale' reading the 20; a mixture half's trial noise is the half's,
  `mixture_noise`); the correlation over odors between the halves of the residual entries and of the shares. Remy's 14 older megamat flies with partial panels
  are not in viki's data module (stage_population_data.py leaves them out), so nothing goes beside.
- stages (design items 4 and 7): megamat along the stages: each odor's strength entries at the flies' ORNs, PN
  dendrites, PN boutons, KC claws and KC somata (mean and range over flies) beside the model's (the family as fitted and
  held, the beside models), each stage in its default reading.
- feeders (design item 5): for the odors that carry each panel's error (`carrying_odors`), the glomeruli that carry their
  ORN input (`feeders`: the model's input on its 36 glomeruli, and the panel's ORN recordings on every recorded
  glomerulus beside), each flagged for being among the model's 36, among the dendrites' 28, in each class of the lock,
  labelled by the bouton line (`BOUTON_LINE`, GH146-Gal4 by the thesis's methods; its published expression by
  glomerulus is unresolved, so the flag is left empty and never inferred from the recordings, `bouton_line_labels`),
  and its boutons and
  claws, counts and shares, in the hemibrain and in the male CNS's two calyces (Pratyush's tables through olfsysm's
  'malecns-left' / 'malecns-right' loaders, E130's `step0_wiring_only`).
- score: the classes read beforehand: "concentrated" per panel (`concentrated`) and (d) (`d_reading`); classes.json,
  tables.md, provenance_score.json.
- regression (design item 6): HELD. The regression of the family's residuals on the share of each odor's input on a
  class, read on validation2's odors against random classes, fits a quantity to a held-out panel, which waits for Tom's
  or Yang's word (CLAUDE.md, the entry's Status). The command line refuses it unless `--authorized-by` names that word
  (who, when, where); `class_regression` is the function, never run on the real residuals here.

The classes read here (fixed beforehand, in the lock's 'rules'):
- "concentrated": a panel's two largest odors (by their share of the squared strength error on all flies, the family as
  fitted, the wirings' mean views) carry more than half of the squared strength error on both halves of the flies;
- (d): the male CNS gives those glomeruli (the feeders of the panel's two largest odors) more boutons or claws: on both
  male calyces their summed share of the calyx's boutons, or of its claws, is above the hemibrain's. Read on megamat as
  the entry says, and given for every panel beside.
(a) and (b) need item 6 and are not read.

Concentrations: megamat at 10^-3 on both sides; the APL16 entries involve paa, euc and -aPine, whose model input is the
ORN data within 1 log10 unit of the flies' concentrations; validation2 is one panel at the same concentrations; the
mixtures are used as measured.

Holds (CLAUDE.md): only Pratyush's hemibrain and male CNS (olfsysm's loaders of his tables); no olfsysm change; no fit to
KC recordings; item 6 held. Nothing is published.

Outputs (<shared>/results/E138_strength_by_odor_entries/): classes_lock.json, classes_glomeruli.csv; entries.csv,
entries_split.csv, entries_summary.csv; halves.csv, halves_split.csv, halves_summary.csv, halves_agreement.csv;
stages.csv, stages_split.csv; feeders.csv; classes.json, tables.md; provenance_<part>.json.

Usage (from the repository root; atlas is shared, so two threads under nice):
    ./venv/bin/python -u experiments/E138_strength_by_odor_entries.py --part classes   # then commit the lock copy
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E138_strength_by_odor_entries.py \\
        --part entries   # likewise halves, stages, feeders, score
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import socket
import subprocess
import time
from collections import namedtuple
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E138_strength_by_odor_entries'
SCRIPT = 'experiments/E138_strength_by_odor_entries.py'
CLASSES_LOCK = 'classes_lock.json'
CLASSES_LOCK_COPY = H.HERE / 'experiments' / 'E138_classes_lock.json'
E137_LOCK_COPY = H.HERE / 'experiments' / 'E137_lock.json'

S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
STRENGTH_PARTS: Tuple[str, ...] = ('strength/rate', 'strength/mean')
PART_KEYS: Dict[str, str] = {'strength/rate': 'rate', 'strength/mean': 'mean'}
HALVES: Tuple[str, ...] = ('fit', 'held')
HALF_SEED = 0
CLASSES: Tuple[str, ...] = ('ir', 'not_dendrite28', 'few_hemibrain')
CALYCES: Tuple[str, ...] = ('hemibrain-first', 'malecns-left', 'malecns-right')
HEMIBRAIN = CALYCES[0]
MALE_CALYCES: Tuple[str, ...] = CALYCES[1:]
# the bouton recordings' driver line: Remy's PN boutons were imaged in GH146-Gal4/20xUAS-IVS-jGCaMP8m (II)
# (~/src/oconnell_thomas_2026_thesis, oconnell_thomas_2026.tex, the bouton figure's legend and methods)
BOUTON_LINE = 'GH146-Gal4'
BOUTON_LINE_SOURCE = ("the thesis's methods (oconnell_thomas_2026.tex: 'PN boutons were imaged in "
                      "GH146-Gal4/20xUAS-IVS-jGCaMP8m (II)'; '~70% of PNs (labeled by the GH146-Gal4 driver)')")
BOUTON_LINE_EXPRESSION = None     # unresolved: no published per-glomerulus table on hand (`bouton_line_labels`)
# Task et al. 2022's table of receptors and co-receptors by glomerulus, as drosolf ships it
TASK_CITATION = ('Task D, Lin C-C, Vulpe A, Afify A, Ballou S, Brbić M, Schlegel P, Raji J, Jefferis GSXE, Li H, Menuz K, '
                 'Potter CJ (2022) Chemoreceptor co-expression in Drosophila melanogaster olfactory neurons. eLife '
                 '11:e72599')
SILBERING_CITATION = ('Silbering AF, Rytz R, Grosjean Y, Abuin L, Ramdya P, Jefferis GSXE, Benton R (2011) '
                      'Complementary function and integrated wiring of the evolutionarily distinct Drosophila olfactory '
                      'subsystems. J Neurosci 31:13357-13375')
COUNTS: Tuple[str, ...] = ('n_boutons', 'n_claws')
# the class rules (fixed in the lock before any residual)
CONCENTRATED_TOP = 2
CONCENTRATED_SHARE = 0.5
# the feeders of an odor: the smallest set of glomeruli, by descending positive input, carrying this share of the odor's
# summed positive input
FEEDER_SHARE = 2.0 / 3.0
# the odors that carry a panel's error: its CONCENTRATED_TOP largest, and any other with at least CARRY_MIN_SHARE of the
# squared view error, at most CARRY_MAX per reading of the family
CARRY_MIN_SHARE = 0.10
CARRY_MAX = 5
# item 6 (held): random classes of the same size
N_RANDOM_CLASSES = 1000
REGRESSION_SEED = 0
REGRESSION_PASS_R = 0.5
REGRESSION_PASS_Q = 0.95

# the groups of saved runs (`group_runs`): the family as fitted and as held, and the beside models
BESIDE_MODELS: Tuple[str, ...] = ('published_hemibrain', 'uniform_draw', 'bouton_baseline', 'reference',
                                  'firing_ceiling')
E90_D1_MODELS: Tuple[str, ...] = ('fit1', 'P1_fc', 'P2_fc', 'M1_fc')
FAMILY_GROUPS: Tuple[str, ...] = ('family_fitted', 'family_held')
# the panels read at the KC somata: (input panel, read panel, conditions, groups)
PanelSpec = namedtuple('PanelSpec', 'input_panel panel conditions groups')
PANEL_SPECS: Tuple[PanelSpec, ...] = (
    PanelSpec('megamat', 'megamat', ('on',), FAMILY_GROUPS + BESIDE_MODELS),
    PanelSpec('apl16', 'apl16', ('on', 'off'), ('family_fitted', 'firing_ceiling_d1')),
    PanelSpec('validation2', 'validation2', ('on',), FAMILY_GROUPS + BESIDE_MODELS),
    PanelSpec('mixtures', 'kmix', ('on',), FAMILY_GROUPS + BESIDE_MODELS),
    PanelSpec('mixtures', 'cmix', ('on',), FAMILY_GROUPS + BESIDE_MODELS),
)
MIXTURE_PANELS: Tuple[str, ...] = ('kmix', 'cmix')
# the reading of the flies whose entries decide the classes, per panel (the others beside): E134's for the KC-soma
# panels; for the mixtures E137's 'as_uploaded' (olfdb dataset 8 as uploaded, E134's reading, with its trial noise)
MAIN_READING: Dict[str, str] = {'megamat': 'e134', 'apl16': 'e134', 'validation2': 'e134', 'kmix': 'as_uploaded',
                                'cmix': 'as_uploaded'}
MIXTURES_PENDING = ("the kiwi and control mixtures wait on the held-out panels read like for like (E137): its lock "
                    "(experiments/E137_lock.json) is not committed, so their rows are not read")
STAGES_MEGAMAT: Tuple[str, ...] = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')


def out_dir() -> Path:
    """<shared>/results/E138_strength_by_odor_entries/ (created)."""
    return H.results_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def _json_default(x):
    if hasattr(x, 'item'):
        return x.item()
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')
    os.replace(tmp, path)


def _write_csv(df: pd.DataFrame, name: str) -> Path:
    p = out_dir() / name
    tmp = p.with_name(f'.{p.name}.{_host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, p)
    return p


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git with no uncommitted change."""
    try:
        if subprocess.run(['git', '-C', str(H.HERE), 'ls-files', '--error-unmatch', str(path)],
                          capture_output=True).returncode != 0:
            return False
        return subprocess.run(['git', '-C', str(H.HERE), 'diff', '--quiet', 'HEAD', '--', str(path)],
                              capture_output=True).returncode == 0
    except OSError:
        return False


# ---------------------------------------------------------------------------
# design item 1: the classes lock

def model_glomeruli() -> List[str]:
    """The model's 36 input glomeruli (the lab input's rows, `H.apl16_orn_deltas`; every panel's input is on them)."""
    return [str(g) for g in H.apl16_orn_deltas().index]


def task_table_path() -> Path:
    """drosolf's copy of Task et al. 2022's table (data/task22_table3.csv)."""
    from drosolf import orns
    return Path(orns.task_csv)


def ir_glomeruli(path: Optional[Path] = None) -> dict:
    """The glomeruli fed by ORNs whose ligand-specific receptor is an IR: Task et al. 2022's table (`task_table_path`),
    every glomerulus with a receptor in its 'receptors' column starting 'Ir' (Ir64a: DC4, DP1m; Ir75a: DP1l; Ir75b:
    DL2d; Ir75c: DL2v; Ir84a: VL2a; Ir31a: VL2p; Ir41a, Ir92a, Ir76a, Ir75d, the VP glomeruli's Ir21a, Ir40a, Ir68a,
    Ir93a), which are the IR glomeruli of Silbering et al. 2011; with the file, its sha256 and both citations, and the
    co-receptor columns (orco, ir8a, ir76b, ir25a) of every glomerulus beside (not a class: Ir25a is in most ORNs,
    Task et al. 2022)."""
    path = task_table_path() if path is None else Path(path)
    t = pd.read_csv(path)
    recs = {str(r.glomerulus): [x.strip() for x in str(r.receptors).split(',')] for r in t.itertuples()}
    ir = sorted(g for g, rs in recs.items() if any(x.startswith('Ir') for x in rs))
    co = {str(r.glomerulus): {c: str(getattr(r, c)) for c in ('orco', 'ir8a', 'ir76b', 'ir25a')} for r in t.itertuples()}
    return dict(glomeruli=ir, receptors={g: recs[g] for g in ir}, coreceptors=co, file=str(path),
                sha256=_sha256(path), citation=TASK_CITATION, also=SILBERING_CITATION,
                rule="a glomerulus is IR-fed when a receptor in Task et al. 2022's 'receptors' column starts 'Ir' (its "
                     "ORNs' ligand-specific receptor is an IR); the co-receptor columns are beside, not a class")


def dendrite_glomeruli() -> dict:
    """The 28 glomeruli of the PN-dendrite recordings: E135's lock ('glomeruli', the 28 that the GH146 recordings, the
    pebbled recordings and the model's input share, E90's `shared_glomeruli`), with the lock's sha256."""
    p = H.shared_dir() / 'results' / 'E135_antennal_lobe_dendrites' / 'lock.json'
    g = [str(x) for x in json.loads(p.read_text())['glomeruli']]
    if len(g) != 28:
        raise RuntimeError(f"E135's lock has {len(g)} glomeruli, not 28")
    return dict(glomeruli=g, file=str(p), sha256=_sha256(p),
                rule="E135's lock 'glomeruli' (the 28 that GH146, the pebbled ORN recordings and the model's input "
                     "share; E90's shared_glomeruli)")


_W_CACHE: Dict[str, pd.DataFrame] = {}


def calyx_wiring(calyx: str) -> pd.DataFrame:
    """A calyx's claw x bouton wiring without a simulation (E130's `step0_wiring_only`: 'hemibrain-first' is the family's
    wiring, E83's 'real', checked equal to `J.bouton_base`'s in E130's check.json step0_wiring_same; the male CNS's two
    calyces are Pratyush's tables through olfsysm's 'malecns-left' / 'malecns-right')."""
    if calyx not in _W_CACHE:
        from experiments import E130_malecns_family as E130
        _W_CACHE[calyx] = E130.step0_wiring_only(calyx)
    return _W_CACHE[calyx]


def wiring_counts(W: pd.DataFrame) -> pd.DataFrame:
    """Per glomerulus of a claw x bouton wiring: its boutons (columns) and claws (the claws whose bouton is the
    glomerulus's), as counts and as shares of the calyx's (columns n_boutons, n_claws, share_boutons, share_claws)."""
    gl = pd.Index([str(x) for x in W.columns.get_level_values('glomerulus')])
    nb = pd.Series(1, index=gl).groupby(level=0).sum()
    nc = pd.Series(W.to_numpy(dtype=float).sum(axis=0), index=gl).groupby(level=0).sum()
    out = pd.DataFrame({'n_boutons': nb.astype(int), 'n_claws': nc.astype(int)})
    out['share_boutons'] = out.n_boutons / out.n_boutons.sum()
    out['share_claws'] = out.n_claws / out.n_claws.sum()
    out.index.name = 'glomerulus'
    return out


def few_hemibrain(counts: pd.DataFrame, universe: Sequence[str]) -> dict:
    """The glomeruli of `universe` with fewer boutons, or fewer claws, than the median glomerulus of `universe` in the
    hemibrain wiring's `counts` (a glomerulus with no bouton counts 0)."""
    c = counts.reindex(list(universe)).fillna(0)
    med = {k: float(c[k].median()) for k in COUNTS}
    fb = [g for g in universe if c.loc[g, 'n_boutons'] < med['n_boutons']]
    fc = [g for g in universe if c.loc[g, 'n_claws'] < med['n_claws']]
    return dict(glomeruli=[g for g in universe if g in set(fb) | set(fc)], few_boutons=fb, few_claws=fc,
                few_both=[g for g in fb if g in set(fc)], median_boutons=med['n_boutons'], median_claws=med['n_claws'],
                counts={g: {k: int(c.loc[g, k]) for k in COUNTS} for g in universe},
                rule='fewer boutons than the median of the 36 input glomeruli, or fewer claws than their median, in '
                     "the family's hemibrain wiring (E130 step0_wiring_only('hemibrain-first'), E83's 'real'); the union "
                     'is the class, each sub-list and their intersection (few_both) beside')


def class_rules() -> dict:
    """The rules of the classes read and of the tables' choices, fixed in the lock before any residual."""
    return dict(
        strength='the strength view of kcq/population_fit.py: parts rate and mean (multiples of the panel mean); a '
                 "part's error is the RMS over odors of (model entry - flies' mean entry) / the part's one-fly unit; "
                 "the view's error the RMS of its two parts' errors",
        odor_share="an odor's share of a part's squared error is its squared residual over the part's sum; its share of "
                   "the view's squared error is (e_rate^2 / n_rate + e_mean^2 / n_mean) / 2 over the view's error^2; the "
                   "error without an odor's entry removes that odor's term and keeps the number of odors (the planner's "
                   'lead)',
        rows="the family's rows are each wiring and the wirings' mean views ('common', whose squared error splits "
             "exactly into odor shares); the shares that decide are those of 'common', the mean over wirings of each "
             "wiring's share beside",
        reading="each stage in its default reading (stage_readings.DEFAULT_READINGS); the flies as E134's step 1 reads "
                "them (E134 fly_side_read: validation2's KC somata on the largest odor set at least 6 flies share)",
        halves="E136's seeded split (population_fit.split_halves(n, 1, 0)[0]; the first part 'fit', the rest 'held') of "
               "the flies of the all-flies reading, on its odors; the APL16 control half paired with the TNT half of "
               "the same name; each half read as the E134 lock's split_halves rule says (curve, noise, reference and "
               "one-fly units from the half)",
        concentrated=dict(rule="a panel's two largest odors carry more than half of the squared strength error on both "
                               'halves of the flies',
                          largest="the two odors with the largest share of the squared view error on all flies, the "
                                  "family as fitted, row 'common'",
                          on_halves="their summed share of the view's squared error on each half (each half's own "
                                    "reading, row 'common') above 0.5",
                          top=CONCENTRATED_TOP, share=CONCENTRATED_SHARE,
                          beside=["each half's own two largest and their share", 'the family as held, read the same way',
                                  "each part's shares"],
                          per='panel and condition (APL16 with APL on and with APL off)'),
        carrying=dict(rule=f"a panel's odors that carry its error: the family's (as fitted, and as held where read) "
                           f"{CONCENTRATED_TOP} largest on all flies and any other with at least {CARRY_MIN_SHARE} of "
                           f"the squared view error, at most {CARRY_MAX} per reading"),
        feeders=dict(rule=f"an odor's feeders: the smallest set of glomeruli, by descending positive input, carrying "
                          f'{FEEDER_SHARE:.4f} of its summed positive input',
                     sources=["the model's input on its 36 glomeruli (the panel's input; the class reading)",
                              "the panel's ORN recordings, every recorded glomerulus (beside)"]),
        d=dict(rule='the male CNS gives those glomeruli more boutons or claws',
               glomeruli="the union of the model-input feeders of the panel's two largest odors (the class "
                         "'concentrated' names)",
               reading="on both male calyces the glomeruli's summed share of the calyx's boutons, or of its claws, is "
                       "above the hemibrain's (shares, since the calyces differ in size)",
               read_on='megamat, as the entry says; every panel beside'),
        not_read=['(a) and (b): they need design item 6, held'],
        item6='held for Tom\'s or Yang\'s word (the regression fits to a held-out panel); refused by the command line')


def bouton_line_labels() -> dict:
    """The bouton line's published expression by glomerulus: unresolved. The line is GH146-Gal4 (BOUTON_LINE_SOURCE);
    no per-glomerulus table of its published expression is on hand that can be named with certainty (FlyBase: about 90
    PNs in about 40 of the about 50 glomeruli; a per-glomerulus table of GH146 and MZ699 innervation is in an eLife
    2014 figure supplement, PMC4270039, not read here). The flag is left empty, never inferred from the recordings."""
    return dict(line=BOUTON_LINE, source=BOUTON_LINE_SOURCE, labels=BOUTON_LINE_EXPRESSION,
                status='unresolved: no published per-glomerulus expression table of GH146-Gal4 on hand that can be '
                       'named with certainty; the flag is left empty (never inferred from the recordings)',
                leads=['FlyBase FBti0016783: GH146 in ~90 PNs innervating ~40 of ~50 glomeruli (no list)',
                       'eLife 2014 (PMC4270039), Figure 1 figure supplement: a per-glomerulus table of GH146-GAL4 and '
                       'MZ699-GAL4 innervation (an image; not read)'])


def classes_contents() -> dict:
    """The classes lock's contents: the universe (the model's 36 input glomeruli), the three classes with their sources
    and rules, the bouton line, and the rules of the classes read (`class_rules`). No residual is read."""
    uni = model_glomeruli()
    ir = ir_glomeruli()
    missing = [g for g in uni if g not in ir['coreceptors']]
    if missing:
        raise RuntimeError(f"Task et al. 2022's table lacks the input glomeruli {missing}")
    den = dendrite_glomeruli()
    hemi = wiring_counts(calyx_wiring(HEMIBRAIN))
    few = few_hemibrain(hemi, uni)
    classes = {
        'ir': dict(glomeruli=[g for g in uni if g in set(ir['glomeruli'])], all_ir=ir['glomeruli'],
                   receptors=ir['receptors'], source=dict(file=ir['file'], sha256=ir['sha256'],
                                                          citation=ir['citation'], also=ir['also']),
                   rule=ir['rule'], unresolved=[]),
        'not_dendrite28': dict(glomeruli=[g for g in uni if g not in set(den['glomeruli'])], dendrite28=den['glomeruli'],
                               source=dict(file=den['file'], sha256=den['sha256']), rule=den['rule']),
        'few_hemibrain': few,
    }
    return dict(experiment=EXPERIMENT, script=SCRIPT, written=_now(), host=_host(), commits=H.code_commits(),
                universe=dict(glomeruli=uni, rule="the model's 36 input glomeruli (H.apl16_orn_deltas' rows)"),
                classes=classes, coreceptors={g: ir['coreceptors'][g] for g in uni}, bouton_line=bouton_line_labels(),
                rules=class_rules(), e134_lock_sha256=_sha256(H.HERE / 'experiments' / 'E134_lock.json'),
                note='written before any residual (model minus flies) was computed; every part but classes refuses '
                     'until this file is committed (require_classes_lock)')


def classes_table(lock: Mapping, counts: Optional[Mapping[str, pd.DataFrame]] = None) -> pd.DataFrame:
    """One row per glomerulus of the universe: its membership of each class (and the sub-lists), its co-receptors, and
    (with `counts`) its boutons and claws in each calyx."""
    cl = lock['classes']
    rows = []
    for g in lock['universe']['glomeruli']:
        r = dict(glomerulus=g, in_model36=True, in_dendrite28=g in set(cl['not_dendrite28']['dendrite28']),
                 **{c: g in set(cl[c]['glomeruli']) for c in CLASSES},
                 few_boutons=g in set(cl['few_hemibrain']['few_boutons']),
                 few_claws=g in set(cl['few_hemibrain']['few_claws']),
                 receptors=';'.join(cl['ir']['receptors'].get(g, [])),
                 **{f'co_{k}': v for k, v in lock['coreceptors'][g].items()})
        for c, df in (counts or {}).items():
            for k in df.columns:
                r[f'{c}__{k}'] = df[k].get(g, 0)
        rows.append(r)
    return pd.DataFrame(rows)


def all_counts() -> Dict[str, pd.DataFrame]:
    """`wiring_counts` of every calyx (CALYCES)."""
    return {c: wiring_counts(calyx_wiring(c)) for c in CALYCES}


def part_classes() -> dict:
    """Design item 1: classes_lock.json in the shared results and its committed copy experiments/E138_classes_lock.json
    (commit it before any other part); classes_glomeruli.csv. Refuses to overwrite a committed lock."""
    if CLASSES_LOCK_COPY.exists() and _committed(CLASSES_LOCK_COPY):
        raise PermissionError(f'{CLASSES_LOCK_COPY.name} is committed already; the classes are fixed')
    lock = classes_contents()
    text = json.dumps(lock, indent=1, default=_json_default) + '\n'
    (out_dir() / CLASSES_LOCK).write_text(text)
    CLASSES_LOCK_COPY.write_text(text)
    _write_csv(classes_table(lock, all_counts()), 'classes_glomeruli.csv')
    for c in CLASSES:
        print(f"{c}: {lock['classes'][c]['glomeruli']}", flush=True)
    print(f'wrote {out_dir() / CLASSES_LOCK} and {CLASSES_LOCK_COPY} (commit the copy before any other part)',
          flush=True)
    return lock


def require_classes_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed classes lock: experiments/E138_classes_lock.json exists, is tracked with no uncommitted change,
    and has the same bytes as results/.../classes_lock.json. PermissionError otherwise (no residual is read before)."""
    copy = CLASSES_LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no classes lock yet (--part classes, then commit it); no residual is read '
                              'before it is committed')
    if not _committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = (out_dir() / CLASSES_LOCK) if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    lock['sha256'] = _sha256(copy)
    return lock


def require_e137_lock() -> dict:
    """E137's committed lock (experiments/E137_lock.json, tracked, no uncommitted change), before any mixture row is
    read; PermissionError otherwise."""
    if not E137_LOCK_COPY.exists() or not _committed(E137_LOCK_COPY):
        raise PermissionError(MIXTURES_PENDING)
    return json.loads(E137_LOCK_COPY.read_text())


def require_all() -> Tuple[dict, dict]:
    """(this experiment's classes lock, E134's lock), both required before any residual."""
    from kcq import stage_readings as R
    return require_classes_lock(), R.require_lock()


# ---------------------------------------------------------------------------
# the entries (pure functions)

def strength_reference(fp: Sequence, odors: Sequence[str]):
    """hal's `fly_reference` of the strength parts only, on `odors`, with E134's breadth bins for that many odors."""
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    return P.fly_reference(list(fp), list(odors), B.breadth_edges(len(odors)), parts=list(STRENGTH_PARTS))


def views_of(pops: Sequence, ref) -> List[Dict[str, pd.Series]]:
    """Each population's strength parts read as `ref` reads its flies (hal's `reference_views`: `population_views` on
    the reference's odors, or a partial-panel reference's strength rule)."""
    from kcq import population_fit as P
    out = []
    for p in pops:
        v = P.reference_views(p, ref)
        out.append({q: v[q] for q in STRENGTH_PARTS if q in v})
    return out


def strength_split_of(pops: Sequence, ref, names: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """E136's `strength_split` (level, spread and order of each strength part, in its one-fly unit; the same formula)
    on the parts as `ref` reads its flies (`views_of`, so that a partial-panel reference's rule applies): one row per
    population, 'mean' (their mean) and 'common' (their mean views) when there are two or more, column 'row' first."""
    from kcq import population_fit as P
    views = views_of(pops, ref)
    names = list(names) if names is not None else [p.name or str(i) for i, p in enumerate(pops)]
    sets = list(zip(names, views)) + ([('common', P.mean_views(views))] if len(views) >= 2 else [])
    rows = {}
    for name, v in sets:
        rec = {}
        for q in [q for q in STRENGTH_PARTS if q in ref.mean]:
            key = q.replace('/', '_')
            m_all = ref.mean[q].to_numpy(dtype=float)
            a_all = v[q].reindex(ref.mean[q].index).to_numpy(dtype=float)
            unit = float(ref.one_fly[q])
            ok = np.isfinite(a_all) & np.isfinite(m_all)
            a, m = a_all[ok], m_all[ok]
            if len(a) < 2:
                rec.update({f'{key}_{k}': np.nan for k in ('level', 'spread', 'order', 'sd_ratio', 'r')})
                continue
            sa, sm = float(a.std()), float(m.std())
            r = float(np.corrcoef(a, m)[0, 1]) if sa > 0 and sm > 0 else np.nan
            cross = 2.0 * sa * sm - 2.0 * float(np.mean((a - a.mean()) * (m - m.mean())))
            rec.update({f'{key}_level': abs(float(np.mean(a - m))) / unit, f'{key}_spread': abs(sa - sm) / unit,
                        f'{key}_order': float(np.sqrt(max(cross, 0.0))) / unit,
                        f'{key}_sd_ratio': sa / sm if sm > 0 else np.nan, f'{key}_r': r})
        rows[name] = rec
    df = pd.DataFrame.from_dict(rows, orient='index')
    if len(views) >= 2:
        mean = df.loc[names].mean().to_frame('mean').T
        df = pd.concat([df.loc[names], mean, df.loc[['common']]])
    return df.rename_axis('row').reset_index()


def odor_entries(pops: Sequence, ref, names: Optional[Sequence[str]] = None, common: bool = True) -> pd.DataFrame:
    """Per population (and the populations' mean views, row 'common', when there are two or more), odor: for each
    strength part (rate, mean) the model's entry, the flies' mean entry and its range over flies, the residual (model
    minus flies) and the residual in one-fly units ('<part>_e'), its share of the part's squared error and the part's
    error with and without it; and the odor's share of the view's squared error, the view's error and the view's error
    without that odor's entry (its term removed, the number of odors kept).

    Args:
        pops: hal's `Population`s (one per wiring or model).
        ref: the flies' reference (`strength_reference`).
        names: the rows' names (default the populations').
        common: add the row 'common' (the wirings' mean views).

    Returns:
        One row per (row, odor), columns 'row', 'odor', then per part p in ('rate', 'mean'): p_model, p_flies,
        p_flies_min, p_flies_max, p_resid, p_e, p_share, p_error, p_error_without; then view_term, view_share,
        view_error, view_error_without, n_flies.
    """
    from kcq import population_fit as P
    parts = [q for q in STRENGTH_PARTS if q in ref.mean]
    views = views_of(pops, ref)
    names = list(names) if names is not None else [p.name or str(i) for i, p in enumerate(pops)]
    sets = list(zip(names, views))
    if common and len(views) >= 2:
        sets.append(('common', P.mean_views(views)))
    odors = list(ref.odors)
    rows = []
    for name, v in sets:
        cols: Dict[str, np.ndarray] = {}
        terms = np.zeros(len(odors))
        for q in parts:
            k = PART_KEYS[q]
            m = ref.mean[q].reindex(odors).to_numpy(dtype=float)
            F = ref.parts[q].reindex(columns=odors).to_numpy(dtype=float)
            a = v[q].reindex(odors).to_numpy(dtype=float)
            u = float(ref.one_fly[q])
            d = a - m
            e = d / u
            ok = np.isfinite(e)
            n = max(int(ok.sum()), 1)
            e2 = np.where(ok, e ** 2, 0.0)
            err2 = e2.sum() / n
            with np.errstate(invalid='ignore', divide='ignore'):
                cols.update({f'{k}_model': a, f'{k}_flies': m, f'{k}_flies_min': np.nanmin(F, axis=0),
                             f'{k}_flies_max': np.nanmax(F, axis=0), f'{k}_resid': d, f'{k}_e': e,
                             f'{k}_share': np.where(ok, e2 / e2.sum() if e2.sum() > 0 else np.nan, np.nan),
                             f'{k}_error': np.full(len(odors), np.sqrt(err2)),
                             f'{k}_error_without': np.sqrt(np.maximum(err2 - e2 / n, 0.0))})
            terms = terms + e2 / n / len(parts)
        view2 = terms.sum()
        cols['view_term'] = terms
        cols['view_share'] = terms / view2 if view2 > 0 else np.full(len(odors), np.nan)
        cols['view_error'] = np.full(len(odors), np.sqrt(view2))
        cols['view_error_without'] = np.sqrt(np.maximum(view2 - terms, 0.0))
        df = pd.DataFrame(cols)
        df.insert(0, 'odor', odors)
        df.insert(0, 'row', name)
        df['n_flies'] = len(ref.flies)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def view_error_of(entries: pd.DataFrame, row: str = 'common') -> float:
    """The strength view's error of one row of `odor_entries`."""
    x = entries[entries.row == row]
    return float(x.view_error.iloc[0]) if len(x) else np.nan


def top_odors(entries: pd.DataFrame, row: str = 'common', k: int = CONCENTRATED_TOP, col: str = 'view_share') -> List[str]:
    """The `k` odors with the largest `col` in one row of `odor_entries` (ties in odor order)."""
    x = entries[entries.row == row].dropna(subset=[col])
    return list(x.sort_values(col, ascending=False, kind='mergesort').odor[:k])


def summed_share(entries: pd.DataFrame, odors: Sequence[str], row: str = 'common', col: str = 'view_share') -> float:
    """The summed `col` of `odors` in one row of `odor_entries`."""
    x = entries[(entries.row == row) & entries.odor.isin(list(odors))]
    return float(x[col].sum()) if len(x) else np.nan


def concentrated(all_flies: pd.DataFrame, halves: Mapping[str, pd.DataFrame], row: str = 'common',
                 top: int = CONCENTRATED_TOP, share: float = CONCENTRATED_SHARE) -> dict:
    """The class "concentrated" of one panel (the lock's rule): the `top` odors with the largest share of the squared
    view error on all flies (`all_flies`, `odor_entries` of the family as fitted) carry more than `share` of it on each
    half (`halves`: half -> `odor_entries` on that half's reading). Beside: each half's own `top` odors and their
    share, and each part's summed share of the same odors.

    Returns:
        dict with 'odors', 'share_all', 'share_halves', 'cls' ('concentrated' / 'not concentrated'; 'not readable'
        without both halves), 'halves_own' and 'part_shares'.
    """
    od = top_odors(all_flies, row, top)
    sh = {h: summed_share(e, od, row) for h, e in halves.items()}
    ok = len(sh) == len(HALVES) and all(np.isfinite(v) for v in sh.values())
    cls = ('concentrated' if all(v > share for v in sh.values()) else 'not concentrated') if ok else 'not readable'
    own = {h: dict(odors=top_odors(e, row, top), share=summed_share(e, top_odors(e, row, top), row))
           for h, e in halves.items()}
    parts = {h: {k: summed_share(e, od, row, f'{k}_share') for k in PART_KEYS.values()}
             for h, e in [('all', all_flies), *halves.items()]}
    return dict(odors=od, share_all=summed_share(all_flies, od, row), share_halves=sh, cls=cls, halves_own=own,
                part_shares=parts, rule=f'the {top} largest on all flies carry more than {share} on both halves')


def halves_agreement(a: pd.DataFrame, b: pd.DataFrame, row: str = 'common') -> dict:
    """The correlation over odors of two halves' entries (`odor_entries`): per part, Pearson r and Spearman ρ of the
    residuals (model minus the half's flies) and of the residuals in one-fly units squared (the squared error's terms),
    and of the view's shares; on the odors both define."""
    from scipy import stats
    x = a[a.row == row].set_index('odor')
    y = b[b.row == row].set_index('odor')
    od = [o for o in x.index if o in y.index]
    out = {}
    for col in [f'{k}_resid' for k in PART_KEYS.values()] + [f'{k}_share' for k in PART_KEYS.values()] + ['view_share']:
        u, v = x.loc[od, col].to_numpy(dtype=float), y.loc[od, col].to_numpy(dtype=float)
        ok = np.isfinite(u) & np.isfinite(v)
        if ok.sum() < 3 or np.std(u[ok]) == 0 or np.std(v[ok]) == 0:
            out[col] = dict(pearson=np.nan, spearman=np.nan, n=int(ok.sum()))
            continue
        out[col] = dict(pearson=float(np.corrcoef(u[ok], v[ok])[0, 1]),
                        spearman=float(stats.spearmanr(u[ok], v[ok]).statistic), n=int(ok.sum()))
    return out


def feeders(inputs: pd.Series, share: float = FEEDER_SHARE) -> pd.DataFrame:
    """The glomeruli that carry an odor's input: the smallest set, by descending positive input, whose summed input is
    at least `share` of the odor's summed positive input (ties in the given order); one row each with its input, share
    and cumulative share and rank. Empty when no input is positive."""
    x = pd.Series(inputs, dtype=float).dropna()
    pos = x[x > 0]
    if not len(pos):
        return pd.DataFrame(columns=['glomerulus', 'input', 'share', 'cum_share', 'rank'])
    s = pos.sort_values(ascending=False, kind='mergesort')
    sh = s / s.sum()
    cum = sh.cumsum()
    n = int(np.searchsorted(cum.to_numpy(), share - 1e-12) + 1)
    n = min(n, len(s))
    return pd.DataFrame(dict(glomerulus=[str(g) for g in s.index[:n]], input=s.to_numpy()[:n], share=sh.to_numpy()[:n],
                             cum_share=cum.to_numpy()[:n], rank=np.arange(1, n + 1)))


def d_reading(glomeruli: Sequence[str], counts: Mapping[str, pd.DataFrame], hemibrain: str = HEMIBRAIN,
              male: Sequence[str] = MALE_CALYCES) -> dict:
    """(d), the lock's reading: whether the male CNS gives `glomeruli` more boutons or claws: on every male calyx their
    summed share of the calyx's boutons, or of its claws, is above the hemibrain's. Per calyx the summed shares and
    counts, and per male calyx the ratios to the hemibrain's.

    Returns:
        dict with 'glomeruli', 'calyces' ({calyx: {share_boutons, share_claws, n_boutons, n_claws}}), 'more'
        ({male calyx: {boutons, claws, either}}), 'cls' ('more' / 'not more'; 'not readable' without glomeruli).
    """
    gl = [str(g) for g in glomeruli]
    if not gl:
        return dict(glomeruli=[], calyces={}, more={}, cls='not readable')
    cal = {}
    for c, df in counts.items():
        d = df.reindex(gl).fillna(0)
        cal[c] = {k: float(d[k].sum()) for k in ('share_boutons', 'share_claws', 'n_boutons', 'n_claws')}
    more = {}
    for c in male:
        b = cal[c]['share_boutons'] > cal[hemibrain]['share_boutons']
        k = cal[c]['share_claws'] > cal[hemibrain]['share_claws']
        more[c] = dict(boutons=bool(b), claws=bool(k), either=bool(b or k),
                       ratio_boutons=cal[c]['share_boutons'] / cal[hemibrain]['share_boutons']
                       if cal[hemibrain]['share_boutons'] > 0 else np.nan,
                       ratio_claws=cal[c]['share_claws'] / cal[hemibrain]['share_claws']
                       if cal[hemibrain]['share_claws'] > 0 else np.nan)
    return dict(glomeruli=gl, calyces=cal, more=more, cls='more' if all(m['either'] for m in more.values()) else 'not more')


def class_regression(residuals: pd.Series, class_share: pd.Series, read_residuals: pd.Series, read_share: pd.Series,
                     random_shares: Optional[Sequence[Tuple[np.ndarray, np.ndarray]]] = None) -> dict:
    """Design item 6 (held; never run on the real residuals without Tom's or Yang's word): the residuals (flies minus
    model, per odor) regressed by least squares on the share of each odor's input on a class (`class_share`, with an
    intercept), the fit's prediction on the read panel's odors (`read_share`) correlated with their residuals
    (`read_residuals`, Pearson r); with `random_shares` (one pair per random class: its shares on the fit odors and on
    the read odors, in the orders of `residuals` and `read_residuals` restricted to the odors with a share) the same r
    for each random class and the quantile of the class's r among them.

    Returns:
        dict with 'slope', 'intercept', 'r', and with random classes 'r_random_q95' and 'quantile' and 'passes'
        (r >= REGRESSION_PASS_R and above REGRESSION_PASS_Q of the random classes).
    """
    def fit_r(xf, yf, xr, yr) -> Tuple[float, float, float]:
        A = np.column_stack([np.ones(len(xf)), xf])
        coef, *_ = np.linalg.lstsq(A, yf, rcond=None)
        pred = coef[0] + coef[1] * xr
        r = float(np.corrcoef(pred, yr)[0, 1]) if np.std(pred) > 0 and np.std(yr) > 0 else np.nan
        return float(coef[1]), float(coef[0]), r

    od = [o for o in residuals.index if o in class_share.index]
    odr = [o for o in read_residuals.index if o in read_share.index]
    slope, icpt, r = fit_r(class_share[od].to_numpy(float), residuals[od].to_numpy(float),
                           read_share[odr].to_numpy(float), read_residuals[odr].to_numpy(float))
    out = dict(slope=slope, intercept=icpt, r=r, n_fit=len(od), n_read=len(odr))
    if random_shares is not None:
        rr = np.array([fit_r(np.asarray(f, float), residuals[od].to_numpy(float), np.asarray(g, float),
                             read_residuals[odr].to_numpy(float))[2] for f, g in random_shares])
        q95 = float(np.nanquantile(rr, REGRESSION_PASS_Q))
        out.update(r_random_q95=q95, quantile=float(np.nanmean(rr < r)),
                   passes=bool(np.isfinite(r) and r >= REGRESSION_PASS_R and r > q95))
    return out


# ---------------------------------------------------------------------------
# the saved runs and the readings

def group_runs(input_panel: str, group: str) -> List:
    """A group's saved runs on an input panel, in WIRINGS order for the family (so that each wiring's seed index is
    E134's): megamat, the family's and the beside models' step-3 runs (E134); validation2, the family's E107 runs and
    the beside models' step-3 runs; APL16, E90's part B (the family) and part D1 (the four firing-ceiling fits,
    'firing_ceiling_d1'). A run whose files are missing is left out and reported."""
    from kcq import stage_readings as R
    loaders: List[Tuple[str, Callable]] = []
    if group in FAMILY_GROUPS:
        reading = group.split('_')[1]
        for w in WIRINGS:
            if input_panel == 'apl16':
                if reading != 'fitted':
                    raise ValueError('the APL16 panel is the family\'s fit panel: no held reading')
                loaders.append((w, lambda w=w: R.load_e90_run(w, 'B')))
            elif input_panel == 'validation2':
                loaders.append((w, lambda w=w: R.load_e107_run('validation2', reading, w)))
            else:
                loaders.append((w, lambda w=w: R.load_step3_run(f'{input_panel}_{reading}_{w}')))
    elif group == 'firing_ceiling_d1':
        loaders = [(t, lambda t=t: R.load_e90_run(t, 'D1')) for t in E90_D1_MODELS]
    elif group in BESIDE_MODELS:
        loaders = [(group, lambda: R.load_step3_run(f'{input_panel}_fitted_{group}'))]
    else:
        raise ValueError(f'{group=}')
    runs = []
    for name, f in loaders:
        try:
            runs.append(f())
        except FileNotFoundError as e:
            print(f'{input_panel} {group} {name}: not found ({e})', flush=True)
    return runs


def run_names(runs: Sequence) -> List[str]:
    """Each run's row name: the wiring for the family, the model or D1 tag otherwise."""
    return [r.wiring if r.model == 'C2' else (r.tag if r.source == 'E90:D1' else r.model) for r in runs]


def soma_flies(panel: str) -> Dict[str, object]:
    """The KC somata's recorded populations a panel is read against, per model condition (viki's, through
    `stage_readings.fly_side`; None where there are none: APL off outside the APL16 panel)."""
    from kcq import stage_readings as R
    return {c: R.fly_side('kc_soma', panel, c) for c in R.CONDITIONS}


def flies_read(fl) -> Tuple[List, List[str]]:
    """The flies and odors of a recorded population as E134's step 1 reads them (`E134 fly_side_read`)."""
    from experiments import E134_population_baseline as B
    return B.fly_side_read(fl)


def half_names(names: Sequence[str], seed: int = HALF_SEED) -> Dict[str, List[str]]:
    """E136's split of a set of flies (`population_fit.split_halves(n, 1, seed)[0]`): 'fit' the first part (floor(n/2)
    flies), 'held' the rest, in the given order."""
    from kcq import population_fit as P
    a, b = P.split_halves(len(names), 1, seed)[0]
    return dict(fit=[names[i] for i in a], held=[names[i] for i in b])


def _subset(fl, names: Sequence[str]):
    from kcq import stage_readings as R
    return R.subset_flies(fl, list(names))


def read_somata(runs: Sequence, panel: str, flies: Mapping[str, object], lock: Mapping) -> Dict[str, list]:
    """Each run's KC somata on `panel` in the default reading, read against `flies` ({condition: recorded population or
    None}) alone (a FlyCache without fallback, so that a half's curve and noise are the half's): {condition: [hal's
    Population per run]}."""
    from kcq import stage_readings as R
    cache = R.FlyCache(fixed={('kc_soma', panel, c): flies.get(c) for c in R.CONDITIONS}, fallback=False)
    fam = R.read_family(list(runs), stages=['kc_soma'], lock=dict(lock), flies=cache)
    out: Dict[str, list] = {}
    for c in R.CONDITIONS:
        pops = [m.populations.get(('kc_soma', panel, c)) for m in fam.values()]
        if pops and all(p is not None for p in pops):
            out[c] = pops
    return out


# a reading of the flies of one panel on one set of flies ('all', or a half): its label ('e134' the main reading,
# 'fly_scale' validation2's partial-panel reading beside it, or one of E137's readings of the mixture flies), the set,
# the key of the model's reading (sides that read the model alike share it), how the model is read ((runs, E134's
# lock) -> {condition: populations}), the flies' reference ((condition, the model's odors) -> FlyReference), and what
# was read
Side = namedtuple('Side', 'reading flies read_key read ref info')


def _soma_read(panel: str, flmap: Mapping[str, object]) -> Callable:
    return lambda runs, lock: read_somata(runs, panel, flmap, lock)


def _strength_ref(fpmap: Mapping[str, list], odmap: Mapping[str, Sequence[str]]) -> Callable:
    return lambda c, mo: strength_reference(fpmap[c], [o for o in odmap[c] if o in set(mo)])


def _partial_ref(fps: list, odors: Sequence[str], rule: str = 'fly_scale') -> Callable:
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    return lambda c, mo: P.partial_reference(fps, rule=rule, odors=[o for o in odors if o in set(mo)],
                                             breadth_edges=D.breadth_bins, parts=list(STRENGTH_PARTS))


def soma_sides(spec: 'PanelSpec', halves: bool) -> List[Side]:
    """The readings of the flies of a KC-soma panel (megamat, APL16, validation2): 'e134' (E134's reading: the model
    read against all the panel's recorded flies, the reference those of `flies_read`; on a half, the model read
    against the half of those flies and the reference theirs), and for validation2 'fly_scale' beside (hal's
    `partial_reference` with rule 'fly_scale' on all 20 flies and 22 odors; on a half, half of the 20)."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    fl = soma_flies(spec.panel)
    for c in spec.conditions:
        if fl.get(c) is None:
            raise RuntimeError(f'no KC somata on {spec.panel} for condition {c}')
    base = {c: flies_read(fl[c]) for c in spec.conditions}
    odmap = {c: base[c][1] for c in spec.conditions}
    sides: List[Side] = []
    if not halves:
        flmap = {c: fl.get(c) for c in ('on', 'off')}
        sides.append(Side('e134', 'all', 'all', _soma_read(spec.panel, flmap),
                          _strength_ref({c: base[c][0] for c in spec.conditions}, odmap),
                          dict(flies={c: [f.name for f in base[c][0]] for c in spec.conditions}, odors=odmap,
                               model_read_on={c: list(fl[c].flies) for c in spec.conditions})))
    else:
        hn = {c: half_names([f.name for f in base[c][0]]) for c in spec.conditions}
        for h in HALVES:
            sub = {c: _subset(fl[c], hn[c][h]) for c in spec.conditions}
            sides.append(Side('e134', h, f'e134:{h}', _soma_read(spec.panel, {c: sub.get(c) for c in ('on', 'off')}),
                              _strength_ref({c: P.fly_populations(sub[c]) for c in spec.conditions}, odmap),
                              dict(flies={c: hn[c][h] for c in spec.conditions}, odors=odmap)))
    if spec.panel == 'validation2':
        od22 = list(D.panel_odors('validation2'))
        names20 = list(fl['on'].flies)
        if not halves:
            sides.append(Side('fly_scale', 'all', 'all', _soma_read(spec.panel, {'on': fl['on'], 'off': None}),
                              _partial_ref(P.fly_populations(fl['on']), od22),
                              dict(flies={'on': names20}, odors={'on': od22}, rule='fly_scale')))
        else:
            hn = half_names(names20)
            for h in HALVES:
                sub = _subset(fl['on'], hn[h])
                sides.append(Side('fly_scale', h, f'fly_scale:{h}', _soma_read(spec.panel, {'on': sub, 'off': None}),
                                  _partial_ref(P.fly_populations(sub), od22),
                                  dict(flies={'on': hn[h]}, odors={'on': od22}, rule='fly_scale')))
    return sides


def mixture_noise(panel: str, reading: str, flies: Sequence[str], lock: Mapping):
    """E137's trial noise (`TrialNoise`, as its `reading_population` sums it) of some of a mixture panel's flies under a
    reading, in the reading's primary units; None where none of them has trials."""
    from experiments import E137_heldout_readings as E
    units = E.PRIMARY_UNITS[reading]
    ss, dof, ns, with_trials = 0.0, 0, {}, []
    for f in flies:
        _, t, nz = E.reading_fly(f, reading, units, lock)
        if t is not None:
            with_trials.append(f)
        ss += nz['ss']
        dof += nz['dof']
        for k, v in nz['n'].items():
            ns.setdefault(k, []).append(v)
    if dof <= 0:
        return None
    return E.TrialNoise(sd=float(np.sqrt(ss / dof)), n={k: int(np.median(v)) for k, v in ns.items()}, dof=dof,
                        flies=with_trials)


def _mixture_read(panel: str, pop, noise) -> Callable:
    from experiments import E137_heldout_readings as E
    return lambda runs, lock: {'on': E.read_model_panel(runs, panel, pop, noise)}


def mixture_sides(panel: str, halves: bool, lock137: Mapping) -> List[Side]:
    """The readings of a mixture panel's flies ('kmix', 'cmix'): each of E137's five readings (its `READINGS`, in their
    primary units, `reading_population`), the model read as E137's re-score reads it (`read_model_panel`: the curve set
    on the reading's flies, their trial noise added); on a half, the half's flies, curve and trial noise."""
    from data_prep import stage_population_data as D
    from experiments import E137_heldout_readings as E
    from kcq import population_fit as P
    od = list(D.panel_odors(panel))
    sides: List[Side] = []
    for r in E.READINGS:
        pop, noise = E.reading_population(panel, r, lock=lock137)
        sets = {'all': (pop, noise, list(pop.flies))}
        if halves:
            hn = half_names(list(pop.flies))
            sets = {h: (_subset(pop, hn[h]), mixture_noise(panel, r, hn[h], lock137), hn[h]) for h in HALVES}
        for h, (pp, nz, names) in sets.items():
            sides.append(Side(r, h, f'{r}:{h}', _mixture_read(panel, pp, nz),
                              _strength_ref({'on': P.fly_populations(pp)}, {'on': od}),
                              dict(flies={'on': names}, odors={'on': od}, units=E.PRIMARY_UNITS[r],
                                   noise=None if nz is None else nz.as_dict())))
    return sides


def sides_of(spec: 'PanelSpec', halves: bool) -> List[Side]:
    """The readings of a panel's flies (`soma_sides`, or `mixture_sides` once E137's lock is committed)."""
    if spec.panel in MIXTURE_PANELS:
        return mixture_sides(spec.panel, halves, require_e137_lock())
    return soma_sides(spec, halves)


def _rename(pops: Sequence, names: Sequence[str]) -> list:
    return [type(p)(p.values, p.responds, name=n, named_units=p.named_units, graded=p.graded)
            for p, n in zip(pops, names)]


def _tag(df: pd.DataFrame, **cols) -> pd.DataFrame:
    for i, (k, v) in enumerate(cols.items()):
        df.insert(i, k, v)
    return df


def read_panel(spec: 'PanelSpec', lock: Mapping, halves: bool = False,
               groups: Optional[Sequence[str]] = None) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The entries of one panel's groups at the KC somata under each reading of its flies (`sides_of`), on all flies
    (`halves` False) or on each half (True): the entries (`odor_entries`), the strength splits (`strength_split_of`)
    and what was read (the flies, odors and the runs' files)."""
    sides = sides_of(spec, halves)
    ents, splits = [], []
    info = dict(panel=spec.panel, input_panel=spec.input_panel,
                sides=[dict(reading=s.reading, flies=s.flies, **s.info) for s in sides], groups={})
    for g in (groups or spec.groups):
        if g not in spec.groups:
            continue
        runs = group_runs(spec.input_panel, g)
        if not runs:
            continue
        names = run_names(runs)
        info['groups'][g] = dict(runs=[r.tag for r in runs], sources=sorted({r.source for r in runs}),
                                 files=[dict(path=f['path'], sha256=f.get('sha256')) for r in runs for f in r.files])
        reads: Dict[str, Dict[str, list]] = {}
        for side in sides:
            if side.read_key not in reads:
                reads[side.read_key] = side.read(runs, lock)
            pops = reads[side.read_key]
            for c in spec.conditions:
                if c not in pops:
                    continue
                named = _rename(pops[c], names)
                ref = side.ref(c, list(named[0].values.columns))
                e = odor_entries(named, ref, names)
                s = strength_split_of(named, ref, names)
                ents.append(_tag(e, panel=spec.panel, condition=c, group=g, reading=side.reading, flies=side.flies))
                splits.append(_tag(s, panel=spec.panel, condition=c, group=g, reading=side.reading, flies=side.flies))
        print(f'{spec.panel} {g}: read ({len(sides)} readings of the flies)', flush=True)
        del runs, reads
        gc.collect()
    return (pd.concat(ents, ignore_index=True) if ents else pd.DataFrame(),
            pd.concat(splits, ignore_index=True) if splits else pd.DataFrame(), info)


SUMMARY_KEYS: Tuple[str, ...] = ('panel', 'condition', 'group', 'reading', 'flies')


def summarize(entries: pd.DataFrame) -> pd.DataFrame:
    """Per (panel, condition, group, reading, flies, row): the view's and parts' errors, the two largest odors and
    their summed share, the error without the largest odor's entry; and a row 'mean' per group (the mean over its
    wirings of each wiring's errors, and of each wiring's shares for the two largest)."""
    keys = list(SUMMARY_KEYS)
    rows = []
    for k, x in entries.groupby(keys, sort=False):
        per = {}
        for row, y in x.groupby('row', sort=False):
            y = y.sort_values('view_share', ascending=False, kind='mergesort')
            top = list(y.odor[:CONCENTRATED_TOP])
            rows.append(dict(zip(keys, k), row=row, view_error=float(y.view_error.iloc[0]),
                             rate_error=float(y.rate_error.iloc[0]), mean_error=float(y.mean_error.iloc[0]),
                             first_odor=top[0], second_odor=top[1] if len(top) > 1 else None,
                             first_share=float(y.view_share.iloc[0]),
                             top2_share=float(y.view_share.iloc[:CONCENTRATED_TOP].sum()),
                             first_rate_share=float(y.rate_share.iloc[0]), first_mean_share=float(y.mean_share.iloc[0]),
                             view_error_without_first=float(y.view_error_without.iloc[0]),
                             n_odors=int(len(y)), n_flies=int(y.n_flies.iloc[0])))
            per[row] = y.set_index('odor')
        wir = [r for r in per if r != 'common']
        if len(wir) >= 2:
            sh = pd.concat([per[r].view_share for r in wir], axis=1).mean(axis=1).sort_values(ascending=False,
                                                                                             kind='mergesort')
            rows.append(dict(zip(keys, k), row='mean',
                             view_error=float(np.mean([per[r].view_error.iloc[0] for r in wir])),
                             rate_error=float(np.mean([per[r].rate_error.iloc[0] for r in wir])),
                             mean_error=float(np.mean([per[r].mean_error.iloc[0] for r in wir])),
                             first_odor=sh.index[0], second_odor=sh.index[1], first_share=float(sh.iloc[0]),
                             top2_share=float(sh.iloc[:CONCENTRATED_TOP].sum()), n_odors=int(len(sh)),
                             n_flies=int(per[wir[0]].n_flies.iloc[0])))
    return pd.DataFrame(rows)


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, host, time and `extra` (the locks' sha256, the inputs)."""
    _write_json(out_dir() / f'provenance_{part}.json',
                dict(part=part, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(), **extra))


def active_specs(panels: Optional[Sequence[str]] = None) -> Tuple[List['PanelSpec'], str]:
    """The panel specs to read (`PANEL_SPECS`, restricted to `panels`), the mixtures only where E137's lock is
    committed; and the mixtures' status."""
    try:
        require_e137_lock()
        status = "read on E137's lock (experiments/E137_lock.json), under each of its five readings"
    except PermissionError as e:
        status = str(e)
    specs = [s for s in PANEL_SPECS if (panels is None or s.panel in panels)
             and (s.panel not in MIXTURE_PANELS or status.startswith('read'))]
    return specs, status


# ---------------------------------------------------------------------------
# the parts

def part_entries(groups: Optional[Sequence[str]] = None, panels: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Design items 2 and 7 on all flies: entries.csv, entries_split.csv, entries_summary.csv, provenance_entries.json
    (with `panels`, entries_<panels>.csv and the rest so named, for a partial run)."""
    lock, e134 = require_all()
    t0 = time.time()
    specs, mixtures = active_specs(panels)
    print(f'mixtures: {mixtures}', flush=True)
    ents, splits, infos = [], [], []
    for spec in specs:
        e, s, info = read_panel(spec, e134, halves=False, groups=groups)
        ents.append(e)
        splits.append(s)
        infos.append(info)
    sfx = '' if panels is None else '_' + '_'.join(panels)
    ent = pd.concat(ents, ignore_index=True)
    _write_csv(ent, f'entries{sfx}.csv')
    _write_csv(pd.concat(splits, ignore_index=True), f'entries_split{sfx}.csv')
    _write_csv(summarize(ent), f'entries_summary{sfx}.csv')
    provenance(f'entries{sfx}', classes_lock_sha256=lock['sha256'], e134_lock_sha256=e134.get('sha256'),
               seconds=round(time.time() - t0, 1), panels=infos, mixtures=mixtures)
    return ent


def part_halves(groups: Optional[Sequence[str]] = None, panels: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Design items 3 and 7: the entries on each half of the flies; halves.csv, halves_split.csv, halves_summary.csv,
    halves_agreement.csv (the correlation over odors of the halves' entries), provenance_halves.json."""
    lock, e134 = require_all()
    t0 = time.time()
    specs, mixtures = active_specs(panels)
    ents, splits, infos = [], [], []
    for spec in specs:
        e, s, info = read_panel(spec, e134, halves=True, groups=groups)
        ents.append(e)
        splits.append(s)
        infos.append(info)
    sfx = '' if panels is None else '_' + '_'.join(panels)
    ent = pd.concat(ents, ignore_index=True)
    _write_csv(ent, f'halves{sfx}.csv')
    _write_csv(pd.concat(splits, ignore_index=True), f'halves_split{sfx}.csv')
    _write_csv(summarize(ent), f'halves_summary{sfx}.csv')
    _write_csv(agreement_table(ent), f'halves_agreement{sfx}.csv')
    provenance(f'halves{sfx}', classes_lock_sha256=lock['sha256'], e134_lock_sha256=e134.get('sha256'),
               seconds=round(time.time() - t0, 1), panels=infos, mixtures=mixtures,
               older_megamat_flies="Remy's 14 older megamat flies (partial panels) are not in viki's data module "
                                   '(stage_population_data.py leaves them out): nothing beside')
    return ent


def agreement_table(halves: pd.DataFrame) -> pd.DataFrame:
    """`halves_agreement` of every (panel, condition, group, reading, row) of the halves' entries."""
    rows = []
    for (panel, cond, g, rd), x in halves.groupby(['panel', 'condition', 'group', 'reading'], sort=False):
        a, b = x[x.flies == 'fit'], x[x.flies == 'held']
        for row in dict.fromkeys(x.row):
            ag = halves_agreement(a, b, row)
            rows.append(dict(panel=panel, condition=cond, group=g, reading=rd, row=row,
                             **{f'{col}__{k}': v for col, d in ag.items() for k, v in d.items()}))
    return pd.DataFrame(rows)


def part_stages(groups: Sequence[str] = FAMILY_GROUPS + BESIDE_MODELS) -> pd.DataFrame:
    """Design items 4 and 7: megamat along the stages (ORNs, PN dendrites, PN boutons, KC claws, KC somata), each in its
    default reading, the model read against all the stage's flies and the flies as E134's step 1 reads them (E134's
    score reads them so); stages.csv, stages_split.csv, provenance_stages.json."""
    from kcq import stage_readings as R
    lock, e134 = require_all()
    t0 = time.time()
    flies = R.FlyCache()
    ents, splits, info = [], [], {}
    for g in groups:
        runs = group_runs('megamat', g)
        if not runs:
            continue
        names = run_names(runs)
        fam = R.read_family(runs, lock=e134, flies=flies)
        info[g] = dict(runs=[r.tag for r in runs], stages={})
        for st in STAGES_MEGAMAT:
            key = (st, 'megamat', 'on')
            pops = [m.populations.get(key) for m in fam.values()]
            fl = flies.get(key)
            if fl is None or not pops or any(p is None for p in pops):
                continue
            fp, od = flies_read(fl)
            named = _rename(pops, names)
            ref = strength_reference(fp, [o for o in named[0].values.columns if o in set(od)])
            e = odor_entries(named, ref, names)
            s = strength_split_of(named, ref, names)
            ents.append(_tag(e, panel='megamat', condition='on', group=g, stage=st))
            splits.append(_tag(s, panel='megamat', condition='on', group=g, stage=st))
            info[g]['stages'][st] = dict(flies=[f.name for f in fp], n_odors=int(e.odor.nunique()))
        print(f'stages {g}: {time.time() - t0:.0f} s', flush=True)
        del fam, runs
        gc.collect()
    ent = pd.concat(ents, ignore_index=True)
    _write_csv(ent, 'stages.csv')
    _write_csv(pd.concat(splits, ignore_index=True), 'stages_split.csv')
    provenance('stages', classes_lock_sha256=lock['sha256'], e134_lock_sha256=e134.get('sha256'),
               seconds=round(time.time() - t0, 1), groups=info)
    return ent


def panel_input(input_panel: str) -> pd.DataFrame:
    """The model's ORN input of a panel on its 36 glomeruli, odors abbreviated (E134's `panel_input`)."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    return R._abbrev(B.panel_input(input_panel))


# the ORN recordings an odor's flies-side feeders are read from, by the panel it carries the error on
ORN_PANELS: Dict[str, Tuple[str, ...]] = {'megamat': ('megamat', 'validation2'), 'apl16': ('megamat', 'validation2'),
                                          'validation2': ('validation2', 'megamat'), 'kmix': ('kmix',),
                                          'cmix': ('cmix',)}


def orn_flies_mean(odor: str, panels: Sequence[str] = ('megamat', 'validation2')) -> Tuple[Optional[pd.Series], str]:
    """The ORN recordings' mean over flies of each recorded glomerulus's response to `odor` (trial-noise units), from
    the first panel of `panels` whose ORN flies were given it; (None, '') if none."""
    from kcq import stage_readings as R
    for p in panels:
        fl = R.fly_side('orn', p, 'on')
        if fl is None or odor not in fl.odors:
            continue
        cols = [t[odor] for t in fl.tables.values() if odor in t.columns]
        if not cols:
            continue
        m = pd.concat(cols, axis=1).mean(axis=1, skipna=True)
        m.index = [str(g) for g in m.index]
        return m, p
    return None, ''


def main_rows(entries: pd.DataFrame) -> pd.DataFrame:
    """The rows of `entries` that decide the classes: the family's groups, all flies, each panel's MAIN_READING, row
    'common'."""
    e = entries[(entries.flies == 'all') & (entries.row == 'common') & entries.group.isin(FAMILY_GROUPS)]
    return e[[MAIN_READING.get(p) == r for p, r in zip(e.panel, e.reading)]]


def carrying_odors(entries: pd.DataFrame) -> Dict[Tuple[str, str], List[str]]:
    """Per (panel, condition), the odors that carry its error (the lock's 'carrying' rule): for the family as fitted
    and as held (where read), the CONCENTRATED_TOP largest shares of the squared view error on all flies (row 'common',
    the panel's main reading) and any other with at least CARRY_MIN_SHARE, at most CARRY_MAX per reading."""
    out: Dict[Tuple[str, str], List[str]] = {}
    for (panel, cond, g), x in main_rows(entries).groupby(['panel', 'condition', 'group'], sort=False):
        y = x.sort_values('view_share', ascending=False, kind='mergesort')
        od = list(y.odor[:CONCENTRATED_TOP]) + [o for o, s in zip(y.odor[CONCENTRATED_TOP:],
                                                                  y.view_share[CONCENTRATED_TOP:])
                                                 if s >= CARRY_MIN_SHARE]
        cur = out.setdefault((panel, cond), [])
        cur += [o for o in od[:CARRY_MAX] if o not in cur]
    return out


def _frame_sha256(x: pd.DataFrame) -> str:
    return hashlib.sha256(pd.util.hash_pandas_object(x, index=True).values.tobytes()).hexdigest()


def part_feeders() -> pd.DataFrame:
    """Design item 5: for the odors that carry each panel's error (`carrying_odors`, from entries.csv), their feeders
    (`feeders`) from the model's input on its 36 glomeruli and from the panel's ORN recordings (every recorded
    glomerulus), each flagged (among the model's 36, among the dendrites' 28, each class of the lock, the bouton line's
    label: unresolved, left empty) with its boutons and claws, counts and shares, in the three calyces; feeders.csv,
    provenance_feeders.json."""
    lock, _ = require_all()
    ent = pd.read_csv(out_dir() / 'entries.csv')
    carry = carrying_odors(ent)
    counts = all_counts()
    cl = lock['classes']
    uni = set(lock['universe']['glomeruli'])
    den = set(cl['not_dendrite28']['dendrite28'])
    inputs = {p: panel_input(p) for p in ('megamat', 'apl16', 'validation2', 'mixtures')}
    rows = []
    for (panel, cond), odors in carry.items():
        inp = inputs['mixtures' if panel in MIXTURE_PANELS else panel]
        for odor in odors:
            m, op = orn_flies_mean(odor, ORN_PANELS[panel])
            for src, series, where in (('model_input', inp[odor] if odor in inp.columns else None, panel),
                                       ('orn_flies', m, op)):
                if series is None:
                    continue
                for r in feeders(series).itertuples():
                    g = r.glomerulus
                    rec = dict(panel=panel, condition=cond, odor=odor, source=src, source_panel=where, glomerulus=g,
                               input=r.input, share=r.share, cum_share=r.cum_share, rank=r.rank,
                               in_model36=g in uni, in_dendrite28=g in den,
                               **{c: g in set(cl[c]['glomeruli']) for c in CLASSES},
                               few_boutons=g in set(cl['few_hemibrain']['few_boutons']),
                               few_claws=g in set(cl['few_hemibrain']['few_claws']),
                               ir_any=g in set(cl['ir']['all_ir']), bouton_line=BOUTON_LINE,
                               bouton_line_labels=np.nan)
                    for c, df in counts.items():
                        for k in df.columns:
                            rec[f'{c}__{k}'] = df[k].get(g, 0)
                    rows.append(rec)
    df = pd.DataFrame(rows)
    _write_csv(df, 'feeders.csv')
    provenance('feeders', classes_lock_sha256=lock['sha256'], carrying={f'{p}|{c}': o for (p, c), o in carry.items()},
               bouton_line=bouton_line_labels(), feeder_share=FEEDER_SHARE,
               inputs={p: dict(sha256=_frame_sha256(x), shape=list(x.shape)) for p, x in inputs.items()})
    return df


def part_score() -> dict:
    """The classes read beforehand: "concentrated" per panel and condition (the family as fitted; as held beside; each
    panel's main reading, the others beside) and (d) per panel (read on megamat); classes.json and tables.md."""
    lock, _ = require_all()
    ent = pd.read_csv(out_dir() / 'entries.csv')
    hv = pd.read_csv(out_dir() / 'halves.csv')
    fd = pd.read_csv(out_dir() / 'feeders.csv')
    counts = all_counts()
    out = dict(concentrated={}, d={}, mixtures=active_specs()[1], not_read=lock['rules']['not_read'],
               classes_lock_sha256=lock['sha256'], main_reading=MAIN_READING)
    a_all = ent[(ent.flies == 'all') & ent.group.isin(FAMILY_GROUPS)]
    for (panel, cond, g, rd), a in a_all.groupby(['panel', 'condition', 'group', 'reading'], sort=False):
        h = {f: hv[(hv.panel == panel) & (hv.condition == cond) & (hv.group == g) & (hv.reading == rd)
                   & (hv.flies == f)] for f in HALVES}
        h = {f: v for f, v in h.items() if len(v)}
        key = f'{panel}|{cond}|{g}|{rd}'
        c = concentrated(a, h)
        c['main'] = MAIN_READING.get(panel) == rd
        out['concentrated'][key] = c
        gl = fd[(fd.panel == panel) & (fd.condition == cond) & (fd.source == 'model_input')
                & fd.odor.isin(c['odors'])].glomerulus
        out['d'][key] = dict(d_reading(list(dict.fromkeys(gl)), counts), odors=c['odors'],
                             main=c['main'], feeders_complete=bool(set(c['odors']) <= set(fd[(fd.panel == panel)
                             & (fd.condition == cond)].odor)))
    out['d_megamat'] = out['d'].get('megamat|on|family_fitted|e134', {}).get('cls')
    _write_json(out_dir() / 'classes.json', out)
    write_tables(out)
    provenance('score', classes_lock_sha256=lock['sha256'])
    return out


def write_tables(classes: Mapping) -> None:
    """tables.md: the summary per panel (the family's 'common' and 'mean' rows, the beside models), the classes."""
    d = out_dir()
    L = ['# E138: the odors that carry the strength error (KC somata, default reading)', '']
    s = pd.read_csv(d / 'entries_summary.csv')
    keep = s[(s.row.isin(['common', 'mean'])) | ((~s.group.isin(FAMILY_GROUPS)) & (s.row != 'mean'))]
    L += ['| panel | condition | group | reading | row | strength error | first odor | share | second | top-2 share | '
          'without the first |', '|---|---|---|---|---|---:|---|---:|---|---:|---:|']
    for r in keep.itertuples():
        L.append(f'| {r.panel} | {r.condition} | {r.group} | {r.reading} | {r.row} | {r.view_error:.3f} | '
                 f'{r.first_odor} | {r.first_share:.2f} | {r.second_odor} | {r.top2_share:.2f} | '
                 f"{getattr(r, 'view_error_without_first', np.nan):.3f} |")
    L += ['', '## Classes', '', '```', json.dumps(classes, indent=1, default=_json_default), '```', '']
    (d / 'tables.md').write_text('\n'.join(L) + '\n')


def part_regression(authorized_by: Optional[str]) -> dict:
    """Design item 6: refused unless `authorized_by` names Tom's or Yang's word. With it: on each half of the flies, the
    family's residuals (flies minus model, the strength view's rate part, as fitted, row 'common', the main reading)
    over the odors of megamat and the APL16 panel (APL on) regressed on the share of each odor's positive model input on
    each class of the lock (`class_regression`), read on validation2's odors, against N_RANDOM_CLASSES random classes
    of the same size among the 36 input glomeruli; regression.json. Never run here."""
    if not authorized_by:
        raise SystemExit("E138 design item 6 (the regression on a class, read on validation2) is held for Tom's or "
                         "Yang's word; pass --authorized-by \"<who, when, where>\" only once it is given")
    lock, _ = require_all()
    hv = pd.read_csv(out_dir() / 'halves.csv')
    inputs = {p: panel_input(p) for p in ('megamat', 'apl16', 'validation2')}
    uni = lock['universe']['glomeruli']
    rng = np.random.default_rng(REGRESSION_SEED)

    def share_on(cls: Sequence[str], inp: pd.DataFrame) -> pd.Series:
        pos = inp.clip(lower=0)
        return pos.loc[[g for g in cls if g in pos.index]].sum() / pos.sum()

    def stacked(cls: Sequence[str]) -> pd.Series:
        return pd.concat([share_on(cls, inputs[p]).rename(lambda o, p=p: f'{p}|{o}') for p in ('megamat', 'apl16')])

    out = dict(authorized_by=authorized_by, time=_now(), host=_host(), classes={})
    for c in CLASSES:
        cls = lock['classes'][c]['glomeruli']
        rand = [list(rng.choice(uni, size=len(cls), replace=False)) for _ in range(N_RANDOM_CLASSES)]
        res = {}
        for h in HALVES:
            x = hv[(hv.flies == h) & (hv.row == 'common') & (hv.group == 'family_fitted') & (hv.condition == 'on')
                   & (hv.reading == 'e134')]
            fit = pd.concat([-x[x.panel == p].set_index('odor').rate_resid.rename(lambda o, p=p: f'{p}|{o}')
                             for p in ('megamat', 'apl16')])
            rd = -x[x.panel == 'validation2'].set_index('odor').rate_resid
            rshares = [(stacked(rc).reindex(fit.index).to_numpy(float),
                        share_on(rc, inputs['validation2']).reindex(rd.index).to_numpy(float)) for rc in rand]
            res[h] = class_regression(fit, stacked(cls).reindex(fit.index), rd,
                                      share_on(cls, inputs['validation2']).reindex(rd.index), rshares)
        out['classes'][c] = dict(glomeruli=cls, halves=res, passes=all(r.get('passes') for r in res.values()))
    _write_json(out_dir() / 'regression.json', out)
    return out


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('classes', 'entries', 'halves', 'stages', 'feeders', 'score', 'regression'))
    ap.add_argument('--groups', nargs='+', default=None,
                    help='entries / halves / stages: the groups read (default all: the family as fitted and held, the '
                         'beside models, the firing-ceiling fits on the APL16 panel)')
    ap.add_argument('--panels', nargs='+', default=None,
                    help='entries / halves: the panels read (default all; a partial run writes files with a suffix)')
    ap.add_argument('--authorized-by', default=None,
                    help="regression only: Tom's or Yang's word that lifts the hold on design item 6 (who, when, where)")
    a = ap.parse_args(argv)
    if a.part == 'regression':
        part_regression(a.authorized_by)
    elif a.part in ('entries', 'halves'):
        dict(entries=part_entries, halves=part_halves)[a.part](a.groups, a.panels)
    elif a.part == 'stages':
        part_stages(tuple(a.groups) if a.groups else FAMILY_GROUPS + BESIDE_MODELS)
    else:
        dict(classes=part_classes, feeders=part_feeders, score=part_score)[a.part]()


if __name__ == '__main__':
    main()
