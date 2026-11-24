#!/usr/bin/env python3
"""E133's transplant builders: one feature of a donor calyx's tables given to a host calyx (step 2, and step 3's
options), and step 2′'s two re-draws of a calyx's own wiring.

Steps 0 and 1 of experiments/E133_calyx_transplants.py never import this module (its test checks it), so it can be
committed while step-1 lanes run from the same checkout without tripping E108's `_still_clean` guard (which stops a
lane when a file it has loaded changes).

The tables. A calyx's tables are the dict of E133_calyx_transplants.calyx_tables: 'calyx', 'seed', 'orn', 'base', the
claw x bouton wiring 'W' (PN->claw weights; rows the claws, index levels kc_id, claw_id, kc_type, ...; columns
glomerulus, pn_id, bouton_id), the normalized APL weights 'wnorm' ({'wAPLKC', 'wKCAPL'} per claw, indexed as W's rows;
{'wAPLPN', 'wPNAPL'} per bouton, indexed as W's columns; as E10_joint_fit.bouton_base returns them), 'breadth', the KC x
KC matrix 'kckc' (post x pre, KC body IDs, mean row sum 1; the male calyces' E130's matrices, None for the hemibrain,
whose matrix olfsysm loads; a builder that needs it reads olfsysm's `connectome_KCKC_weights(normalize=True)`) and
'meta'. A builder returns a new dict (it never changes its inputs; what it leaves unchanged is the host's own object),
with 'breadth' None when 'W' changed (the runner recomputes it, E133_calyx_transplants.breadth_run) and
meta['transplant'] the transplant's size (the entry's Realism; `size_row` flattens it into runs-table columns).

The arms (the entry's Definitions), with the readings fixed here where the entry leaves them open:
- C, claws per KC by type: within each KC type the KCs ranked by claw count take the donor's count at the same rank; a
  surplus claw, drawn at random, is switched off (its PN->claw weight 1e-6 on its own bouton, its APL weights spread
  over the KC's other claws); a missing claw is a copy of one of the KC's own (that claw's PN->claw weight raised by 1);
  each KC keeps its summed APL->KC and KC->APL weights. Readings: "the same rank" between types of different sizes is
  host rank i of n_h taking the donor's sorted value at floor((i + 1/2) n_d / n_h) (`rank_match`), ties broken at random
  by the seed; a KC's claw count is its rows with input, a copied row counting its weight; the switched-off claw's APL
  weights go to the KC's remaining claws in proportion to theirs (evenly when those are all 0); a copy changes the
  PN->claw weight only (as the entry writes it), so the copied claw's APL weights stay as they are; copies are drawn
  without replacement, a KC needing more copies than it has claws copying each claw floor(k / n) times first. Each KC
  type has its own random stream, so C on one type alone (step 3's split) is the full arm's change on that type.
- S, the glomeruli's shares of the claws: E97's `move_claws` (`experiments/E97_glomerulus_claw_shares.py`) at E97's
  ratios, read as E97 reads them: one r_g per glomerulus, the donor's share of its claws on the host's glomeruli over
  the host's (`E97.ratios`, r_g = 1 where either calyx has fewer than E97's 30 claws of g), put onto the host within
  each KC type (each type's total kept, largest-remainder rounding); not per-type shares. A moved claw keeps its own APL
  weights and takes its new bouton's.
- A, APL's synapses at the KCs: within each KC type the KCs ranked by their summed weight take the donor's amount at
  the same rank, in units of the donor's mean per KC; spread over the KC's claws in its own proportions (evenly over its
  claws with input when it had none, E99's convention for a KC without a synapse on a kept claw); KC->APL the same way,
  ranked on its own. Reading: the new totals are scaled by one factor so that the calyx's total is the host's (with
  "the donor's mean per KC" taken literally the total is off by the ratio of the rank-matched amounts' mean to the
  donor's mean, recorded as 'literal_total_ratio'; both it and the per-KC factors are the same on counted synapses and on
  normalized weights).
- B, APL's synapses at the boutons: glomerulus g's factor f_g = (m_d(g) / M_d) / (m_h(g) / M_h) (m a glomerulus's mean
  weight per bouton, M the calyx's mean), times one common factor so that the calyx's total over the glomeruli that
  change is kept; a glomerulus the donor lacks, or without APL synapses in the host, keeps its own weights exactly. So
  each glomerulus's mean as a multiple of the calyx's is the donor's up to that common factor (recorded as 'rescale').
  Bouton->APL the same way.
- K, the matrix among KCs: block (receiving type a, sending type b) scaled by s_d(a, b) / s_h(a, b) (its share of the
  summed weights), then the matrix scaled by one factor so that its sum, hence its mean row sum (1), is the host's; a
  block empty in the host stays empty (listed). The hemibrain's matrix is olfsysm's, passed as `_wKCKC` (E122's route).
- ALL: S, C, A, B and K, applied in that order, each with the same donor and seed.
The step-3 options are tokens after the arm, joined by '-' (no '_', which separates a tag's parts): 'three' (S, A, B on
the three odors' glomeruli: S with r_g = 1 outside them, as E97's rule for small glomeruli; A and B the full arm's
change on those glomeruli's claws or boutons only), 'g' / 'apbp' / 'ab' (C on that KC type alone), 'aplkc' / 'kcapl'
(A) and 'aplpn' / 'pnapl' (B), one direction alone, and 'dhb' / 'dml' / 'dmr' (the donor, for the runner: `donor_of`).
Step 2′'s re-draws (on a calyx's own tables, seed 0; no donor):
- PA: the within-type claw permutation of `seed` (pyolfsysm's `shuffle_wPNKC`, as step 1's re-draws), each claw's
  APL->KC and KC->APL weights moving with its input (the plain permutation leaves them at the claw).
- CI: E59's co-input-preserving shuffle (`experiments/E59_coinput_shuffle.py`'s `_match_one`): within-type Curveball
  of `seed`, then `match_coinput` to the calyx's own shared-KC counts (seed 100 + seed, within type, type layer held,
  E59's annealing). APL weights stay at their claws. `Refused` if it does not take the wiring.
The identity the entry asks for: with donor and host the same tables every builder returns them unchanged, bit for bit
(`identity_holds`). Its factors are built so that they are exactly 1 then (products in the same order, exact sums).

Holds: Pratyush's hemibrain and male CNS only (the tables E133_calyx_transplants builds); no olfsysm change; bouton
identities are the tables' own, never inferred.

Usage (no simulation; from the repository root):
    OMP_NUM_THREADS=1 nice -n 19 ./venv/bin/python -u -m experiments.E133_builders --sizes
    (each arm's size on each host's four wirings from the model's tables, APL's counted synapses in place of the
    normalized weights, which differ from them by one factor per vector: results/E133_calyx_transplants/
    builder_sizes.csv)
"""
from __future__ import annotations

import argparse
import math
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

ARMS = ('C', 'S', 'A', 'B', 'K', 'ALL')
ALL_ORDER = ('S', 'C', 'A', 'B', 'K')
PRIME_ARMS = ('PA', 'CI')
SWITCHED_OFF_WEIGHT = 1e-6
ACTIVE_MIN = 0.5                     # a row with input at least this is a claw (a switched-off claw is not)
SMALL_CLAWS = 30                     # E97's SMALL_CLAWS (test checks)
# E59's MATCH_KWS (test checks) and its seed of the matching (100 + the draw's seed)
COINPUT_MATCH_KWS = dict(n_steps=3_000_000, t_start=5.0, t_end=0.05, type_layer_weight=4.0)
COINPUT_SEED_OFFSET = 100
DIRECTIONS = {'A': ('wAPLKC', 'wKCAPL'), 'B': ('wAPLPN', 'wPNAPL')}
DIRECTION_TOKENS = {'aplkc': 'wAPLKC', 'kcapl': 'wKCAPL', 'aplpn': 'wAPLPN', 'pnapl': 'wPNAPL'}
TYPE_TOKENS = {'g': 'g', 'apbp': "a'b'", 'ab': 'ab'}
TYPE_COL = {v: k for k, v in TYPE_TOKENS.items()}
HEMI, LEFT, RIGHT = 'hemibrain-argmax', 'malecns-left', 'malecns-right'
DONOR_TOKENS = {'dhb': HEMI, 'dml': LEFT, 'dmr': RIGHT}
DEFAULT_DONOR = {HEMI: LEFT, LEFT: HEMI, RIGHT: HEMI}
# random streams: np.random.default_rng([seed, stream, ...])
RNG_S, RNG_C, RNG_A = 1, 2, 3
_CACHE: dict = {}


class Refused(RuntimeError):
    """A builder that cannot take the tables it was given (e.g. E59's shuffle on a wiring it does not take)."""


# ---------------------------------------------------------------------------
# arms, donors, small helpers

def parse_arm(arm: str) -> dict:
    """An arm string's parts: 'base' (one of ARMS or PRIME_ARMS), 'subset' ('three' or None), 'types' (tuple of KC
    types or None), 'directions' (tuple of wnorm keys or None), 'donor' (a calyx or None). ValueError on a token the base
    does not take."""
    parts = str(arm).split('-')
    base, toks = parts[0], parts[1:]
    if base not in ARMS + PRIME_ARMS:
        raise ValueError(f'arm {arm!r}: {base!r} not in {ARMS + PRIME_ARMS}')
    if '_' in str(arm):
        raise ValueError(f"arm {arm!r}: no '_' (it separates a tag's parts)")
    out = dict(base=base, subset=None, types=None, directions=None, donor=None)
    types, dirs = [], []
    for t in toks:
        if t == 'three' and base in ('S', 'A', 'B'):
            out['subset'] = 'three'
        elif t in TYPE_TOKENS and base == 'C':
            types.append(TYPE_TOKENS[t])
        elif t in DIRECTION_TOKENS and DIRECTION_TOKENS[t] in DIRECTIONS.get(base, ()):
            dirs.append(DIRECTION_TOKENS[t])
        elif t in DONOR_TOKENS and base not in PRIME_ARMS:
            out['donor'] = DONOR_TOKENS[t]
        else:
            raise ValueError(f'arm {arm!r}: token {t!r} not taken by {base}')
    out['types'] = tuple(types) or None
    out['directions'] = tuple(dirs) or None
    return out


def donor_of(host_calyx: str, arm: str) -> Optional[str]:
    """The donor calyx of `arm` on `host_calyx`: None for step 2′'s re-draws; the arm's donor token if any; else the
    entry's (the hemibrain for a male calyx, the male CNS's left for the hemibrain)."""
    p = parse_arm(arm)
    if p['base'] in PRIME_ARMS:
        return None
    d = p['donor'] or DEFAULT_DONOR.get(host_calyx)
    if d is None:
        raise ValueError(f'{host_calyx}: no donor for arm {arm!r}')
    if d == host_calyx:
        raise ValueError(f'arm {arm!r}: the donor is the host ({host_calyx})')
    return d


def _fsum(x) -> float:
    """The exactly rounded sum (independent of the order of the terms)."""
    return math.fsum(np.asarray(x, dtype=float).ravel().tolist())


def _same(a, b) -> bool:
    """Bit for bit: both None, or the same index (and columns) and the same bytes and dtypes."""
    if a is None or b is None:
        return a is None and b is None
    if not a.index.equals(b.index):
        return False
    if isinstance(a, pd.DataFrame) != isinstance(b, pd.DataFrame):
        return False
    if isinstance(a, pd.DataFrame):
        if not a.columns.equals(b.columns) or not a.dtypes.equals(b.dtypes):
            return False
    elif a.dtype != b.dtype:
        return False
    return np.ascontiguousarray(a.to_numpy()).tobytes() == np.ascontiguousarray(b.to_numpy()).tobytes()


def _rows(W: pd.DataFrame) -> dict:
    """Per row of a claw x bouton wiring (at most one input per row): 'kc', 'type', 'v' (values), 'wt' (the row's input
    weight), 'col' (its bouton's column, -1 without), 'code' (KC code in first-appearance order), 'ids' (KC IDs by code),
    'kc_type' (type by code), 'active' (rows with input at least ACTIVE_MIN), 'glom' (the row's glomerulus or '')."""
    v = W.to_numpy(dtype=float)
    if not np.isfinite(v).all() or (v < 0).any():
        raise ValueError('wiring values must be finite and non-negative')
    if ((v > 0).sum(axis=1) > 1).any():
        raise ValueError('a claw x bouton wiring has at most one input per claw')
    kc = np.asarray(W.index.get_level_values(H.KC_ID_LEVEL))
    typ = np.asarray(W.index.get_level_values(H.KC_TYPE_LEVEL)).astype(str)
    wt = v.max(axis=1)
    col = np.where(wt > 0, v.argmax(axis=1), -1)
    code, ids = pd.factorize(kc)
    first = np.unique(code, return_index=True)[1]
    kc_type = typ[first]
    if not np.array_equal(kc_type[code], typ):
        raise ValueError('a KC has rows of more than one type')
    gl = np.asarray([str(x) for x in W.columns.get_level_values('glomerulus')], dtype=object)
    glom = np.where(col >= 0, gl[np.maximum(col, 0)], '')
    return dict(kc=kc, type=typ, v=v, wt=wt, col=col, code=code, ids=np.asarray(ids), kc_type=kc_type,
                active=wt >= ACTIVE_MIN, glom=glom, n_kcs=len(ids))


def claw_counts(W: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(claws per KC, the KCs' types, the KCs' IDs), KCs in first-appearance order; a claw is a row with input of at
    least ACTIVE_MIN, counted by its weight (a copied claw counts twice)."""
    r = _rows(W)
    n = np.bincount(r['code'], weights=np.where(r['active'], np.rint(r['wt']), 0.0), minlength=r['n_kcs'])
    return n.astype(np.int64), r['kc_type'], r['ids']


def rank_match(host_vals, donor_vals, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """Each host value replaced by the donor's value at the same rank: host rank i (ascending, of n_h; ties in random
    order with `rng`, else in position order) takes the donor's sorted value at floor((i + 1/2) n_d / n_h). With the
    same values on both sides every value is returned unchanged (whatever the tie order)."""
    h = np.asarray(host_vals)
    d = np.sort(np.asarray(donor_vals), kind='stable')
    n_h, n_d = len(h), len(d)
    if n_h == 0:
        return h.copy()
    if n_d == 0:
        raise ValueError('rank_match: no donor values')
    if rng is None:
        order = np.argsort(h, kind='stable')
    else:
        p = rng.permutation(n_h)
        order = p[np.argsort(h[p], kind='stable')]
    j = ((2 * np.arange(n_h, dtype=np.int64) + 1) * n_d) // (2 * n_h)
    out = np.empty(n_h, dtype=d.dtype)
    out[order] = d[j]
    return out


def _sha(df: pd.DataFrame) -> str:
    from experiments.E133_calyx_transplants import sha256_frame
    return sha256_frame(df)


def _derive(tables: Mapping, part: str, donor: Optional[Mapping], seed: int, size: dict, **changed) -> dict:
    """A new dict from `tables` with the changed entries (only those that differ from the host's), 'breadth' None when
    'W' changed, and meta['transplant'] with `size` under its parts[`part`]."""
    T = dict(tables)
    wn = dict(tables['wnorm'])
    for k, v in changed.items():
        if k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL'):
            if not _same(v, wn[k]):
                wn[k] = v
        elif not _same(v, tables.get(k)):
            T[k] = v
    if any(wn[k] is not tables['wnorm'][k] for k in wn):
        T['wnorm'] = wn
    meta = dict(tables.get('meta') or {})
    if T['W'] is not tables['W']:
        T['breadth'] = None
        meta.update(W_sha256=_sha(T['W']), breadth_sha256=None)
    if T.get('kckc') is not tables.get('kckc'):
        from experiments.E77_input_screen import frame_sha256
        meta.update(kckc_sha256=frame_sha256(T['kckc']), kckc_source=f'transplant {part} on the calyx\'s matrix')
    tr = dict(meta.get('transplant') or {})
    parts = dict(tr.get('parts') or {})
    parts[part] = dict(size)
    tr.update(host=tables.get('calyx'), host_seed=tables.get('seed'), donor=None if donor is None else donor.get('calyx'),
              seed=int(seed), parts=parts)
    meta['transplant'] = tr
    T['meta'] = meta
    return T


# ---------------------------------------------------------------------------
# C, claws per KC by type

def claws_per_kc(host: Mapping, donor: Mapping, seed: int = 0, types: Optional[Sequence[str]] = None) -> dict:
    """Arm C (the module docstring); `types` restricts it to those KC types (step 3's γ alone)."""
    W = host['W']
    r = _rows(W)
    n_h, kc_type_h, _ = claw_counts(W)
    n_d, kc_type_d, _ = claw_counts(donor['W'])
    v = r['v'].copy()
    w = {k: host['wnorm'][k].to_numpy(dtype=float).copy() for k in DIRECTIONS['A']}
    for k in w:
        if not host['wnorm'][k].index.equals(W.index):
            raise ValueError(f"{k}: not indexed as the wiring's rows")
    order = np.argsort(r['code'], kind='stable')
    rows_of = np.split(order, np.cumsum(np.bincount(r['code'], minlength=r['n_kcs']))[:-1])
    target = n_h.copy()
    rec = dict(kcs=int(r['n_kcs']), claws_host=int(n_h.sum()), claws_off=0, claws_copied=0, kcs_losing=0, kcs_gaining=0,
               spread_factor_min=np.nan, spread_factor_max=np.nan, kcs_spread_evenly=0, types_skipped='')
    spread = []
    skipped = []
    all_types = sorted(set(kc_type_h))
    for ti, t in enumerate(all_types):
        if types is not None and t not in types:
            continue
        m, md = kc_type_h == t, kc_type_d == t
        col_t = TYPE_COL.get(t, t)
        rec[f'claws_per_kc_{col_t}_host'] = float(n_h[m].mean())
        rec[f'claws_per_kc_{col_t}_donor'] = float(n_d[md].mean()) if md.any() else np.nan
        if not md.any():
            skipped.append(t)
            rec[f'claws_per_kc_{col_t}_after'] = float(n_h[m].mean())
            continue
        rng = np.random.default_rng([int(seed), RNG_C, ti])
        target[m] = rank_match(n_h[m], n_d[md], rng)
        off_t = cop_t = 0
        for k in np.flatnonzero(m):
            dlt = int(target[k] - n_h[k])
            if dlt == 0:
                continue
            rows = rows_of[k]
            act = rows[r['active'][rows]]
            if dlt < 0:
                ones = act[r['wt'][act] == 1.0]
                if len(ones) < -dlt:
                    raise RuntimeError(f'KC {r["ids"][k]}: {-dlt} claws to switch off, {len(ones)} single claws')
                off = np.sort(rng.choice(ones, size=-dlt, replace=False))
                v[off, r['col'][off]] = SWITCHED_OFF_WEIGHT
                rest = np.setdiff1d(act, off)
                for kk in w:
                    s_off, s_rest = _fsum(w[kk][off]), _fsum(w[kk][rest])
                    if s_off == 0:
                        continue
                    if s_rest > 0:
                        f = (s_rest + s_off) / s_rest
                        w[kk][rest] = w[kk][rest] * f
                        spread.append(f)
                    else:
                        w[kk][rest] = s_off / len(rest)
                        rec['kcs_spread_evenly'] += 1
                    w[kk][off] = 0.0
                off_t += -dlt
                rec['kcs_losing'] += 1
            else:
                reps, extra = divmod(dlt, len(act))
                add = np.full(len(act), reps, dtype=float)
                if extra:
                    add[rng.choice(len(act), size=extra, replace=False)] += 1
                v[act, r['col'][act]] += add
                cop_t += dlt
                rec['kcs_gaining'] += 1
        rec[f'claws_per_kc_{col_t}_after'] = float(target[m].mean())
        rec[f'off_share_{col_t}'] = off_t / float(n_h[m].sum())
        rec[f'copied_share_{col_t}'] = cop_t / float(n_h[m].sum())
        rec['claws_off'] += off_t
        rec['claws_copied'] += cop_t
    rec['off_share'] = rec['claws_off'] / float(n_h.sum())
    rec['copied_share'] = rec['claws_copied'] / float(n_h.sum())
    rec['max_claw_weight'] = float(v.max())
    if spread:
        rec.update(spread_factor_min=float(min(spread)), spread_factor_max=float(max(spread)))
    rec['types_skipped'] = ' '.join(skipped)
    if types is not None:
        rec['types'] = ' '.join(types)
    changed = (rec['claws_off'] + rec['claws_copied']) > 0
    W1 = pd.DataFrame(v, index=W.index, columns=W.columns) if changed else W
    new = {k: (pd.Series(w[k], index=host['wnorm'][k].index, name=host['wnorm'][k].name) if changed else host['wnorm'][k])
           for k in w}
    return _derive(host, 'C', donor, seed, rec, W=W1, **new)


# ---------------------------------------------------------------------------
# S, the glomeruli's shares of the claws

def shares_ratios(host_W: pd.DataFrame, donor_W: pd.DataFrame, small: int = SMALL_CLAWS
                  ) -> Tuple[pd.Series, pd.Series, pd.DataFrame]:
    """E97's ratios of the donor to the host on the host's glomeruli: (r_g, whether g is small (r_g = 1 by rule), the
    counts and shares per glomerulus)."""
    from experiments import E97_glomerulus_claw_shares as E97
    gl = list(dict.fromkeys(str(x) for x in host_W.columns.get_level_values('glomerulus')))
    hc = E97.type_glomerulus_counts(host_W).sum(axis=0)
    dc = E97.type_glomerulus_counts(donor_W).sum(axis=0)
    counts = pd.DataFrame({'host': hc.reindex(gl).fillna(0), 'donor': dc.reindex(gl).fillna(0)}).astype(float)
    r, sm = E97.ratios(counts, gl, ref='host', small=small)
    tab = counts.assign(share_host=counts['host'] / counts['host'].sum(), share_donor=counts['donor'] / counts['donor'].sum())
    return r['donor'].astype(float), sm['donor'] | sm['host'], tab


def glomerulus_shares(host: Mapping, donor: Mapping, seed: int = 0, glomeruli: Optional[Sequence[str]] = None,
                      small: int = SMALL_CLAWS) -> dict:
    """Arm S (the module docstring); `glomeruli` keeps r_g = 1 outside them (step 3's three odors' glomeruli)."""
    from experiments import E97_glomerulus_claw_shares as E97
    W = host['W']
    r, sm, tab = shares_ratios(W, donor['W'], small)
    if glomeruli is not None:
        r = r.where(r.index.isin(list(glomeruli)), 1.0)
    W1, moved, _ = E97.move_claws(W, r.to_dict(), seed=[int(seed), RNG_S])
    n = float(tab['host'].sum())
    ch = r[r != 1.0]
    rec = dict(claws_host=int(n), claws_moved=int(len(moved)), moved_share=len(moved) / n, n_glomeruli=int(len(r)),
               n_glomeruli_changed=int(len(ch)), n_glomeruli_small=int(sm.sum()),
               r_min=float(ch.min()) if len(ch) else 1.0, r_max=float(ch.max()) if len(ch) else 1.0)
    if len(moved):
        after = E97.type_glomerulus_counts(W1).sum(axis=0).reindex(tab.index).fillna(0).astype(float)
        sh = after / after.sum()
        big = ~sm.reindex(tab.index).to_numpy()
        rec['share_change_max'] = float((sh - tab['share_host']).abs().max())
        rec['share_gap_to_donor_max'] = float((sh - tab['share_donor']).abs()[big].max()) if big.any() else 0.0
        rec['share_gap_to_donor_before_max'] = float((tab['share_host'] - tab['share_donor']).abs()[big].max()) if big.any() else 0.0
    if glomeruli is not None:
        rec['glomeruli'] = ' '.join(sorted(glomeruli))
    return _derive(host, 'S', donor, seed, rec, W=W1 if len(moved) else W)


# ---------------------------------------------------------------------------
# A, APL's synapses at the KCs

def apl_at_kcs(host: Mapping, donor: Mapping, seed: int = 0, directions: Optional[Sequence[str]] = None,
               glomeruli: Optional[Sequence[str]] = None) -> dict:
    """Arm A (the module docstring); `directions` restricts it to those of ('wAPLKC', 'wKCAPL'), `glomeruli` to the
    claws fed by those glomeruli (the full arm's change there, the other claws unchanged)."""
    r, rd = _rows(host['W']), _rows(donor['W'])
    rec: dict = dict(kcs=int(r['n_kcs']))
    new = {}
    sub = None if glomeruli is None else np.isin(r['glom'], list(glomeruli)) & r['active']
    all_types = sorted(set(r['kc_type']))
    for di, k in enumerate(DIRECTIONS['A']):
        if directions is not None and k not in directions:
            continue
        s = host['wnorm'][k]
        if not s.index.equals(host['W'].index) or not donor['wnorm'][k].index.equals(donor['W'].index):
            raise ValueError(f"{k}: not indexed as the wiring's rows")
        w = s.to_numpy(dtype=float)
        h = np.bincount(r['code'], weights=w, minlength=r['n_kcs'])
        d = np.bincount(rd['code'], weights=donor['wnorm'][k].to_numpy(dtype=float), minlength=rd['n_kcs'])
        matched = h.copy()
        skipped = []
        for ti, t in enumerate(all_types):
            m, md = r['kc_type'] == t, rd['kc_type'] == t
            if not md.any():
                skipped.append(t)
                continue
            matched[m] = rank_match(h[m], d[md], np.random.default_rng([int(seed), RNG_A, di, ti]))
        c = _fsum(h) / _fsum(matched)
        tot = matched * c
        pos = h > 0
        fac = np.ones_like(h)
        fac[pos] = tot[pos] / h[pos]
        full = w * fac[r['code']]
        z = (~pos) & (tot > 0)
        if z.any():
            n_act = np.bincount(r['code'], weights=r['active'].astype(float), minlength=r['n_kcs'])
            zr = z[r['code']]
            full = np.where(zr, np.where(r['active'], tot[r['code']] / np.maximum(n_act[r['code']], 1.0), 0.0), full)
        out = full if sub is None else np.where(sub, full, w)
        key = 'aplkc' if k == 'wAPLKC' else 'kcapl'
        f = fac[pos]
        rec.update({f'{key}_factor_min': float(f.min()) if f.size else 1.0,
                    f'{key}_factor_median': float(np.median(f)) if f.size else 1.0,
                    f'{key}_factor_max': float(f.max()) if f.size else 1.0,
                    f'{key}_literal_total_ratio': float((_fsum(matched) / len(matched)) / (_fsum(d) / len(d))),
                    f'{key}_kcs_from_zero': int(z.sum()), f'{key}_kcs_to_zero': int((pos & (tot == 0)).sum()),
                    f'{key}_kcs_zero_host': int((~pos).sum()), f'{key}_types_skipped': ' '.join(skipped)})
        if sub is not None:
            rec[f'{key}_claws_in_subset'] = int(sub.sum())
            rec[f'{key}_subset_total_factor'] = float(_fsum(out[sub]) / _fsum(w[sub])) if _fsum(w[sub]) > 0 else np.nan
        new[k] = pd.Series(out, index=s.index, name=s.name)
    if glomeruli is not None:
        rec['glomeruli'] = ' '.join(sorted(glomeruli))
    if directions is not None:
        rec['directions'] = ' '.join(directions)
    return _derive(host, 'A', donor, seed, rec, **new)


# ---------------------------------------------------------------------------
# B, APL's synapses at the boutons

def apl_at_boutons(host: Mapping, donor: Mapping, seed: int = 0, directions: Optional[Sequence[str]] = None,
                   glomeruli: Optional[Sequence[str]] = None) -> dict:
    """Arm B (the module docstring); `directions` restricts it to those of ('wAPLPN', 'wPNAPL'), `glomeruli` to those
    glomeruli's boutons (the full arm's change there, the others unchanged)."""
    gh = np.asarray([str(x) for x in host['W'].columns.get_level_values('glomerulus')], dtype=object)
    gd = np.asarray([str(x) for x in donor['W'].columns.get_level_values('glomerulus')], dtype=object)
    gl = list(dict.fromkeys(gh))
    rec: dict = dict(n_glomeruli=len(gl))
    new = {}
    for k in DIRECTIONS['B']:
        if directions is not None and k not in directions:
            continue
        s, sd = host['wnorm'][k], donor['wnorm'][k]
        if not s.index.equals(host['W'].columns) or not sd.index.equals(donor['W'].columns):
            raise ValueError(f"{k}: not indexed as the wiring's boutons")
        w, wd = s.to_numpy(dtype=float), sd.to_numpy(dtype=float)
        Mh, Md = _fsum(w) / len(w), _fsum(wd) / len(wd)
        f = np.ones(len(w))
        elig = np.zeros(len(w), dtype=bool)
        missing, zero = [], []
        for g in gl:
            mb, mbd = gh == g, gd == g
            if not mbd.any():
                missing.append(g)
                continue
            mh, md = _fsum(w[mb]) / mb.sum(), _fsum(wd[mbd]) / mbd.sum()
            if not mh > 0:
                zero.append(g)
                continue
            f[mb] = (md * Mh) / (mh * Md)
            elig |= mb
        den = _fsum(w[elig] * f[elig])
        if not den > 0:
            raise RuntimeError(f'{k}: the donor has no APL synapses on the shared glomeruli')
        c = _fsum(w[elig]) / den
        fac = np.where(elig, f * c, 1.0)
        full = w * fac
        out = full if glomeruli is None else np.where(np.isin(gh, list(glomeruli)), full, w)
        key = 'aplpn' if k == 'wAPLPN' else 'pnapl'
        fe = fac[elig]
        rec.update({f'{key}_factor_min': float(fe.min()) if fe.size else 1.0,
                    f'{key}_factor_max': float(fe.max()) if fe.size else 1.0, f'{key}_rescale': float(c),
                    f'{key}_kept_missing': ' '.join(missing), f'{key}_kept_zero': ' '.join(zero)})
        if glomeruli is not None:
            m = np.isin(gh, list(glomeruli))
            rec[f'{key}_subset_total_factor'] = float(_fsum(out[m]) / _fsum(w[m])) if _fsum(w[m]) > 0 else np.nan
        new[k] = pd.Series(out, index=s.index, name=s.name)
    if glomeruli is not None:
        rec['glomeruli'] = ' '.join(sorted(glomeruli))
    if directions is not None:
        rec['directions'] = ' '.join(directions)
    return _derive(host, 'B', donor, seed, rec, **new)


# ---------------------------------------------------------------------------
# K, the matrix among KCs

def hemibrain_kckc(W: pd.DataFrame) -> pd.DataFrame:
    """olfsysm's normalized hemibrain KC x KC matrix among the wiring's KCs (post x pre, mean row sum 1; E122's
    `hemibrain_kckc`), cached in the process."""
    from pyolfsysm.mb_model import connectome_KCKC_weights
    ids = pd.Index(pd.unique(W.index.get_level_values(H.KC_ID_LEVEL)), name=H.KC_ID_LEVEL)
    key = ('hemibrain_kckc', len(ids), hash(tuple(ids.tolist())))
    if key not in _CACHE:
        _CACHE[key] = connectome_KCKC_weights(ids)
    return _CACHE[key]


def kckc_of(tables: Mapping) -> pd.DataFrame:
    """The tables' KC x KC matrix: their own, or the hemibrain's from olfsysm when they have none."""
    if tables.get('kckc') is not None:
        return tables['kckc']
    if not str(tables.get('calyx', '')).startswith('hemibrain'):
        raise ValueError(f"{tables.get('calyx')}: no KC-KC matrix in the tables")
    return hemibrain_kckc(tables['W'])


def _kc_types(W: pd.DataFrame) -> pd.Series:
    r = _rows(W)
    return pd.Series(r['kc_type'], index=pd.Index(r['ids'], name=H.KC_ID_LEVEL))


def block_sums(M: pd.DataFrame, types: pd.Series) -> pd.DataFrame:
    """Each (receiving type, sending type) block's summed weight (exact sums)."""
    tp, tq = types.reindex(M.index).to_numpy(), types.reindex(M.columns).to_numpy()
    if pd.isna(tp).any() or pd.isna(tq).any():
        raise ValueError("the matrix has KCs that are not the wiring's")
    A = M.to_numpy(dtype=float)
    ts = sorted(set(tp) | set(tq))
    return pd.DataFrame([[_fsum(A[np.ix_(tp == a, tq == b)]) for b in ts] for a in ts], index=ts, columns=ts)


def kckc_blocks(host: Mapping, donor: Mapping, seed: int = 0) -> dict:
    """Arm K (the module docstring)."""
    M = kckc_of(host)
    Md = kckc_of(donor)
    th, td = _kc_types(host['W']), _kc_types(donor['W'])
    Bh, Bd = block_sums(M, th), block_sums(Md, td)
    ts = list(Bh.index)
    Bd = Bd.reindex(index=ts, columns=ts).fillna(0.0)
    Th, Td = _fsum(Bh.to_numpy()), _fsum(Bd.to_numpy())
    F = pd.DataFrame(1.0, index=ts, columns=ts)
    empty = []
    for a in ts:
        for b in ts:
            if Bh.loc[a, b] > 0:
                F.loc[a, b] = (Bd.loc[a, b] * Th) / (Bh.loc[a, b] * Td)
            elif Bd.loc[a, b] > 0:
                empty.append(f'{TYPE_COL.get(a, a)}<{TYPE_COL.get(b, b)}')
    tp, tq = th.reindex(M.index).to_numpy(), th.reindex(M.columns).to_numpy()
    ip, iq = pd.Index(ts).get_indexer(tp), pd.Index(ts).get_indexer(tq)
    Fm = F.to_numpy()[ip[:, None], iq[None, :]]
    A = M.to_numpy(dtype=float)
    g = _fsum(A) / _fsum(A * Fm)
    G = F * g
    rec = dict(rescale=float(g), blocks_empty_in_host=' '.join(empty))
    for a in ts:
        for b in ts:
            rec[f'f_{TYPE_COL.get(a, a)}_from_{TYPE_COL.get(b, b)}'] = float(G.loc[a, b])
    same = np.equal.outer(np.arange(len(ts)), np.arange(len(ts)))
    rec.update(within_type_share_host=float(Bh.to_numpy()[same].sum() / Th),
               within_type_share_donor=float(Bd.to_numpy()[same].sum() / Td))
    if (G.to_numpy() == 1.0).all():
        rec['within_type_share_after'] = rec['within_type_share_host']
        return _derive(host, 'K', donor, seed, rec)
    out = pd.DataFrame(A * G.to_numpy()[ip[:, None], iq[None, :]], index=M.index, columns=M.columns)
    Ba = block_sums(out, th)
    rec['within_type_share_after'] = float(Ba.to_numpy()[same].sum() / _fsum(Ba.to_numpy()))
    rec['mean_row_sum_after'] = float(_fsum(out.to_numpy()) / len(out))
    return _derive(host, 'K', donor, seed, rec, kckc=out)


# ---------------------------------------------------------------------------
# step 2′: the re-draws that move APL's weights with the input, and E59's co-input-preserving shuffle

def within_type_permutation(W: pd.DataFrame, seed: int) -> np.ndarray:
    """The source row of each row under pyolfsysm's `shuffle_wPNKC(W, 'claw_permutation', within_kc_type=True,
    seed=seed)` (its draws replicated: row j takes row src[j]'s input; rows without input their own)."""
    from pyolfsysm.wpnkc_shuffles import _kc_codes, _kc_groups
    values = W.to_numpy()
    row_kc, n_kcs = _kc_codes(W)
    claw_rows = np.flatnonzero(values.any(axis=1))
    claw_kc = row_kc[claw_rows]
    kc_group = _kc_groups(W, row_kc, n_kcs)
    claw_group = kc_group[claw_kc]
    rng = np.random.default_rng(int(seed))
    src = np.arange(len(values))
    for g in np.unique(kc_group):
        idx = np.flatnonzero(claw_group == g)
        p = rng.permutation(len(idx))
        src[claw_rows[idx]] = claw_rows[idx[p]]
    return src


def _own_wiring(host: Mapping, arm: str) -> None:
    if int(host.get('seed', 0) or 0) != 0:
        raise ValueError(f"arm {arm} re-draws the calyx's own wiring: pass its seed-0 tables (got seed {host.get('seed')})")


def permutation_with_apl(host: Mapping, seed: int) -> dict:
    """Step 2′'s PA: the within-type claw permutation of `seed` with each claw's APL->KC and KC->APL weights moving with
    its input; RuntimeError unless the wiring is `shuffle_wPNKC`'s exactly."""
    from pyolfsysm.wpnkc_shuffles import shuffle_wPNKC
    _own_wiring(host, 'PA')
    W0 = host['W']
    src = within_type_permutation(W0, seed)
    W1 = pd.DataFrame(W0.to_numpy()[src], index=W0.index.copy(), columns=W0.columns.copy())
    if W0.dtypes.nunique() > 1:
        W1 = W1.astype(W0.dtypes.to_dict())
    ref = shuffle_wPNKC(W0, 'claw_permutation', within_kc_type=True, seed=int(seed))
    if not _same(W1, ref):
        raise RuntimeError("PA: the replicated permutation does not give shuffle_wPNKC's wiring")
    r = _rows(W0)
    moved = src != np.arange(len(src))
    rec = dict(claws_reassigned_share=float(moved.mean()), claws_bouton_changed_share=float((r['col'][src] != r['col']).mean()))
    new = {}
    for k in DIRECTIONS['A']:
        s = host['wnorm'][k]
        if not s.index.equals(W0.index):
            raise ValueError(f"{k}: not indexed as the wiring's rows")
        w = s.to_numpy(dtype=float)
        nw = w[src]
        h = np.bincount(r['code'], weights=w, minlength=r['n_kcs'])
        hn = np.bincount(r['code'], weights=nw, minlength=r['n_kcs'])
        pos = h > 0
        f = hn[pos] / h[pos]
        key = 'aplkc' if k == 'wAPLKC' else 'kcapl'
        rec.update({f'{key}_kc_factor_min': float(f.min()), f'{key}_kc_factor_median': float(np.median(f)),
                    f'{key}_kc_factor_max': float(f.max()), f'{key}_kcs_from_zero': int(((~pos) & (hn > 0)).sum())})
        new[k] = pd.Series(nw, index=s.index, name=s.name)
    T = _derive(host, 'PA', None, seed, rec, W=W1, **new)
    T['meta'] = dict(T['meta'], wiring='claw_permutation_within_type_apl_moves', wiring_seed=int(seed))
    return T


def coinput_shuffle(host: Mapping, seed: int, match_kws: Optional[Mapping] = None) -> dict:
    """Step 2′'s CI: E59's co-input-preserving shuffle of the calyx's own wiring at `seed` (within-type Curveball, then
    `match_coinput` to the wiring's own shared-KC counts, seed COINPUT_SEED_OFFSET + seed, within type, type layer held,
    E59's annealing unless `match_kws`); `Refused` if pyolfsysm's shuffle or matching does not take the wiring."""
    from pyolfsysm.wpnkc_shuffles import match_coinput, shuffle_wPNKC, wPNKC_shuffle_summary
    _own_wiring(host, 'CI')
    W0 = host['W']
    kws = dict(COINPUT_MATCH_KWS if match_kws is None else match_kws)
    try:
        S = shuffle_wPNKC(W0, 'curveball', within_kc_type=True, seed=int(seed))
        L, stats = match_coinput(S, W0, seed=COINPUT_SEED_OFFSET + int(seed), within_kc_type=True, match_type_layer=True,
                                 **kws)
        summ = wPNKC_shuffle_summary(W0, L)
    except (ValueError, AssertionError, KeyError, IndexError) as e:
        raise Refused(f"E59's co-input-preserving shuffle does not take the wiring of {host.get('calyx')}: {e!r}") from e
    rec = {k: (v.item() if hasattr(v, 'item') else v) for k, v in stats.items()
           if k in ('objective_start', 'objective_end', 'pair_abs_diff_max', 'pair_abs_diff_mean', 'coinput_corr',
                    'type_layer_abs_diff_max', 'frac_claws_changed_vs_start', 'frac_claws_changed_vs_target', 'n_steps_made',
                    'seconds')}
    rec.update(kc_n_claws_preserved=bool(summ['kc_n_claws_preserved']),
               max_abs_glom_n_claws_change=int(summ['max_abs_glom_n_claws_change']), match_seed=COINPUT_SEED_OFFSET + int(seed),
               **{f'match_{k}': v for k, v in kws.items()})
    T = _derive(host, 'CI', None, seed, rec, W=L)
    T['meta'] = dict(T['meta'], wiring='curveball_within_type_coinput_matched', wiring_seed=int(seed))
    return T


# ---------------------------------------------------------------------------
# the entry point

BUILDERS: Dict[str, Callable] = {'C': claws_per_kc, 'S': glomerulus_shares, 'A': apl_at_kcs, 'B': apl_at_boutons,
                                 'K': kckc_blocks}


def three_glomeruli(tables: Mapping) -> List[str]:
    """The three odors' glomeruli (E133_calyx_transplants.odor_glomeruli on the tables' input: those carrying at least a
    tenth of 2-but's, paa's or B-myr's summed positive ORN input), sorted."""
    from experiments.E133_calyx_transplants import odor_glomeruli
    return sorted(set().union(*[set(map(str, g.index)) for g in odor_glomeruli(tables['orn']).values()]))


def build(arm: str, host: Mapping, donor: Optional[Mapping], seed: int = 0) -> dict:
    """The host's tables with arm `arm` of the donor's (`parse_arm`; ARMS with step 3's options) or step 2′'s re-draw
    (PA, CI: `host` the calyx's seed-0 tables, `donor` ignored); `seed` draws C's and S's claws, PA's permutation and
    CI's shuffle.

    Args:
        arm: e.g. 'C', 'ALL', 'S-three', 'C-g', 'A-aplkc', 'B-pnapl-dmr', 'PA'.
        host: the host calyx's tables at one of its four wirings (E133_calyx_transplants.calyx_tables).
        donor: the donor calyx's tables at its own wiring.
        seed: for the random draws.

    Returns:
        new tables (a dict), meta['transplant'] holding the arm, host, donor, seed and each part's size.
    """
    p = parse_arm(arm)
    b = p['base']
    if b == 'PA':
        T = permutation_with_apl(host, seed)
    elif b == 'CI':
        T = coinput_shuffle(host, seed)
    else:
        if donor is None:
            raise ValueError(f'arm {arm}: no donor tables')
        sub = three_glomeruli(host) if p['subset'] == 'three' else None
        if b == 'ALL':
            T = dict(host)
            for a in ALL_ORDER:
                T = BUILDERS[a](T, donor, seed=seed)
        elif b == 'C':
            T = claws_per_kc(host, donor, seed=seed, types=p['types'])
        elif b == 'S':
            T = glomerulus_shares(host, donor, seed=seed, glomeruli=sub)
        elif b in ('A', 'B'):
            T = BUILDERS[b](host, donor, seed=seed, directions=p['directions'], glomeruli=sub)
        else:
            T = kckc_blocks(host, donor, seed=seed)
    T['meta'] = dict(T['meta'], transplant=dict(T['meta']['transplant'], arm=str(arm)))
    return T


def identity_holds(arm: str, tables: Mapping, seed: int = 0) -> dict:
    """Whether `build(arm, tables, tables, seed)` returns `tables` unchanged (W, every wnorm vector and kckc bit for
    bit; and the breadth kept): a dict of what differs (empty when it holds)."""
    T = build(arm, tables, tables, seed)
    diff = {}
    if not _same(T['W'], tables['W']):
        diff['W'] = True
    for k in tables['wnorm']:
        if not _same(T['wnorm'][k], tables['wnorm'][k]):
            diff[k] = True
    if not _same(T.get('kckc'), tables.get('kckc')):
        diff['kckc'] = True
    if T.get('breadth') is not tables.get('breadth'):
        diff['breadth'] = True
    return diff


def size_row(meta: Mapping, prefix: str = 'size_') -> dict:
    """A transplant's size (meta['transplant'] of a builder's tables, or the tables themselves) as flat columns:
    'size_arm', 'size_donor', and '<prefix><part>_<name>' for each part's entries."""
    if 'meta' in meta:
        meta = meta['meta']
    tr = meta.get('transplant') or {}
    row = {f'{prefix}arm': tr.get('arm', ''), f'{prefix}donor': tr.get('donor') or ''}
    for part, rec in (tr.get('parts') or {}).items():
        for k, v in rec.items():
            row[f'{prefix}{part}_{k}'] = v
    return row


# ---------------------------------------------------------------------------
# sizes from the model's tables (no simulation)

def light_tables(c: str, seed: int = 0) -> dict:
    """A calyx's tables without a simulation, for the sizes: E130's wiring (`step0_wiring_only`; at seed n its
    within-type claw permutation, as `calyx_tables`), APL's counted synapses (E130's `apl_counts`) as the weights (the
    normalized ones differ by one factor per vector, which no size depends on), the male calyx's KC-KC matrix (E130's
    `malecns_kckc`; None for the hemibrain), the model's input."""
    from experiments import E130_malecns_family as E130
    from pyolfsysm.wpnkc_shuffles import shuffle_wPNKC
    W0 = E130.step0_wiring_only(c)
    a = E130.apl_counts(c, W0)
    W = W0 if not seed else shuffle_wPNKC(W0, 'claw_permutation', within_kc_type=True, seed=int(seed))
    wn = {'wAPLKC': a['apl_kc'].astype(float), 'wKCAPL': a['kc_apl'].astype(float), 'wAPLPN': a['apl_pn'].astype(float),
          'wPNAPL': a['pn_apl'].astype(float)}
    key = ('light_kckc', c)
    if c.startswith('malecns') and key not in _CACHE:
        _CACHE[key] = E130.malecns_kckc(c, E130.model_kc_ids(W0))
    return dict(calyx=c, seed=int(seed), orn=H.apl16_orn_deltas(), base={}, wnorm=wn, W=W, breadth=None,
                kckc=_CACHE.get(key), meta=dict(wiring='own' if not seed else 'claw_permutation_within_type', wiring_seed=int(seed)))


def sizes(hosts: Sequence[Tuple[str, str]] = ((HEMI, LEFT), (LEFT, HEMI), (RIGHT, HEMI), (LEFT, RIGHT)),
          seeds: Sequence[int] = (0, 1, 2, 3), arms: Sequence[str] = ARMS) -> pd.DataFrame:
    """Each arm's size on each (host, donor) and host wiring (`light_tables`), with PA on the left's seeds 1 to 3 (CI's
    is written when it runs): one row each, `size_row`'s columns."""
    rows = []
    for hc, dc in hosts:
        donor = light_tables(dc, 0)
        for s in seeds:
            host = light_tables(hc, s)
            for arm in arms:
                T = build(arm, host, donor, seed=s)
                rows.append(dict(host=hc, donor=dc, host_seed=s, arm=arm, **size_row(T)))
            if hc == LEFT and dc == HEMI and s:
                T = build('PA', light_tables(LEFT, 0), None, seed=s)
                rows.append(dict(host=hc, donor='', host_seed=s, arm='PA', **size_row(T)))
            print(f'E133 sizes: {hc} seed {s} (donor {dc})', flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sizes', action='store_true', help='write results/E133_calyx_transplants/builder_sizes.csv')
    a = ap.parse_args()
    if a.sizes:
        from experiments.E133_calyx_transplants import _host, _now, _write_json, out_dir
        df = sizes()
        p = out_dir() / 'builder_sizes.csv'
        df.to_csv(p, index=False)
        _write_json(out_dir() / 'builder_sizes.json', dict(
            generated=_now(), host=_host(), script='experiments/E133_builders.py --sizes', commits=H.code_commits(),
            n_rows=int(len(df)), note=("Each arm's size (E133_builders.size_row) on each (host, donor) and host wiring "
                                       "(seeds 0 to 3), from the model's tables without a simulation: E130's wiring and "
                                       "APL's counted synapses (light_tables), which differ from the runner's normalized "
                                       "weights by one factor per vector; no recorded size depends on that factor. PA on "
                                       "the left's seeds 1 to 3; CI's size is recorded when it runs.")))
        print(f'wrote {p} ({len(df)} rows)', flush=True)


if __name__ == '__main__':
    main()
