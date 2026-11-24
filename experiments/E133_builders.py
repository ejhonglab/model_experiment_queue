#!/usr/bin/env python3
"""E133's transplant builders (step 2): one feature of a donor calyx's tables given to a host calyx.

STUB (2026-11-24): the interface only. The bodies are written by a later code subagent; steps 0 and 1 of
experiments/E133_calyx_transplants.py never import this module, so it can be committed while step-1 lanes run from the
same checkout without tripping E108's `_still_clean` guard (which stops a lane when a file it has loaded changes).

The tables. A calyx's tables are E133_calyx_transplants.Tables: the claw x bouton wiring 'W' (PN->claw weights; rows
the claws, index levels kc_id, claw_id, kc_type, ...; columns glomerulus, pn_id, bouton_id), the normalized APL weights
'wnorm' ({'wAPLKC', 'wKCAPL'} per claw, {'wAPLPN', 'wPNAPL'} per bouton, as E10_joint_fit.bouton_base returns them), the
KC x KC matrix 'kckc' (post x pre, KC body IDs, mean row sum 1; the male calyces' from E130's matrices, the hemibrain's
from olfsysm's `connectome_KCKC_weights(normalize=True)` when a builder needs it) and the counted APL synapses 'counts'
(E130's `apl_counts`: 'apl_kc', 'kc_apl' per claw, 'apl_pn', 'pn_apl' per bouton). A builder returns new Tables (it
never changes its inputs) with 'meta' describing the transplant's size (the entry's Realism: the share of claws switched
off and copied, the range of the factors on APL's weights, the block factors); the runner recomputes the breadth when
'W' changed (E133_calyx_transplants.breadth_run) and holds both response rates as step 1 does.

The arms (the entry's Definitions, fixed before any run):
- C, claws per KC by type: within each KC type the KCs ranked by claw count take the donor's count at the same rank; a
  surplus claw, drawn at random, is switched off (PN->claw weight 1e-6, its APL weights spread over the KC's other
  claws); a missing claw is a copy of one of the KC's own (that claw's PN->claw weight raised by 1); each KC keeps its
  summed APL->KC and KC->APL weights.
- S, the glomeruli's shares of the claws: E97's `move_claws` to the donor's proportions, within KC type, with E97's rule
  for small glomeruli.
- A, APL's synapses at the KCs: within each KC type the KCs ranked by their summed APL->KC count take the donor's amount
  at the same rank, in units of the donor's mean per KC, spread over the KC's claws in its own proportions; KC->APL the
  same way.
- B, APL's synapses at the boutons: each glomerulus's mean APL->bouton weight, as a multiple of the calyx's mean, set to
  the donor's by one factor on its boutons; bouton->APL the same way; a glomerulus the donor lacks keeps its own.
- K, the matrix among KCs: each block (receiving type by sending type) scaled so that its share of the summed weights
  is the donor's, then normalized to a mean row sum of 1.
- ALL: S, C, A, B and K together, in that order.
The identity the entry asks for: with donor and host the same tables, every builder returns the tables unchanged.
"""
from __future__ import annotations

from typing import Callable, Dict, Optional

ARMS = ('C', 'S', 'A', 'B', 'K', 'ALL')
ALL_ORDER = ('S', 'C', 'A', 'B', 'K')
SWITCHED_OFF_WEIGHT = 1e-6


def build(arm: str, host: 'object', donor: 'object', seed: int = 0) -> 'object':
    """The host's Tables with arm `arm` (one of ARMS) of the donor's; `seed` draws the claws C switches off or copies.

    Args:
        arm: 'C', 'S', 'A', 'B', 'K' or 'ALL'.
        host: the host calyx's Tables at one of its four wirings (its own or a within-type claw permutation).
        donor: the donor calyx's Tables at its own wiring.
        seed: for the random draws (C's surplus and missing claws).

    Returns:
        new Tables, with 'meta' holding the transplant's size.
    """
    raise NotImplementedError('E133 step 2: the builders are not written yet')


def claws_per_kc(host: 'object', donor: 'object', seed: int = 0) -> 'object':
    """Arm C (the module docstring)."""
    raise NotImplementedError


def glomerulus_shares(host: 'object', donor: 'object', seed: int = 0) -> 'object':
    """Arm S (E97's `move_claws` to the donor's proportions within KC type)."""
    raise NotImplementedError


def apl_at_kcs(host: 'object', donor: 'object', seed: int = 0) -> 'object':
    """Arm A (APL->KC and KC->APL per KC, rank-matched within type, in units of the donor's mean per KC)."""
    raise NotImplementedError


def apl_at_boutons(host: 'object', donor: 'object', seed: int = 0) -> 'object':
    """Arm B (each glomerulus's mean APL->bouton and bouton->APL weight relative to the calyx's, set to the donor's)."""
    raise NotImplementedError


def kckc_blocks(host: 'object', donor: 'object', seed: int = 0) -> 'object':
    """Arm K (the KC-KC matrix's type blocks to the donor's shares, then a mean row sum of 1)."""
    raise NotImplementedError


BUILDERS: Dict[str, Callable] = {'C': claws_per_kc, 'S': glomerulus_shares, 'A': apl_at_kcs, 'B': apl_at_boutons,
                                 'K': kckc_blocks}


def identity_holds(arm: str, tables: 'object', seed: int = 0) -> Optional[dict]:
    """Whether `build(arm, tables, tables, seed)` returns `tables` unchanged (W, wnorm, kckc bit for bit): a dict of
    what differs (empty when it holds)."""
    raise NotImplementedError
