#!/usr/bin/env python3
"""E124: the inhibition among KCs as a part of APL, fitted with a KC-KC scale for each condition (STUB: interface only).

The entry is item E124 of experiment_priorities.md (proposed by viki 2026-11-21 08:32, moved 08:42 under the idle-host
rule; planned by Fable subagents, claude-fable-5-2). It builds on E123 (experiments/E123_kckc_form.py: the runner at fixed
scales, the classes, the family, S and T).

Parts (CLI `--part`): check, a1, b, score, a2 (fits, per wiring, --wirings), a3 (gated on a2), variants, chain.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

EXPERIMENT = 'E124_kckc_apl_borne'
DOSES = (0.25, 0.5, 0.75)
A1_ARMS = ('HBoff25', 'HBoff50', 'HBoff75')
B_ARMS = ('TS16s0', 'TS16s1', 'K4s1', 'K4s2', 'K16s1', 'K16s2')


def type_structured_partners(ref, types, k: int = 16, seed: int = 0):
    """k partners per KC with equal weights per postsynaptic type, drawn within the hemibrain's totals between KC types."""
    raise NotImplementedError


def arm_array(arm: str, ref, types):
    """The arm's KC x KC matrix (post x pre)."""
    raise NotImplementedError


def off_kckc_fraction(arm: str) -> Optional[float]:
    """The APL-off KC-KC scale as a fraction of the fit's (A1's arms), None for the others."""
    raise NotImplementedError


def part_check() -> dict:
    """Reproduces E123's held point S1_HBab_real exactly through this runner; the K4 / K16 seed-0 matrices are E123's."""
    raise NotImplementedError


def part_runs(arms: Sequence[str], wirings: Sequence[str]) -> None:
    """Parts a1 and b: the fixed-scale held points (resumable, runs_1_<host>_<lane>.csv)."""
    raise NotImplementedError


def score() -> dict:
    """readings.csv, runs.csv, verdict.json (E123's classes; part A1's dose reading; part B's reading; part A2's class)."""
    raise NotImplementedError


def fit_command(name: str, w: str) -> List[str]:
    """The E10_joint_fit.py command of one fit of part A2 ('A2_<w>') or A3 ('A3_<w>')."""
    raise NotImplementedError


def part_fits(part: str, wirings: Optional[Sequence[str]] = None, dry_run: bool = False) -> List[str]:
    """Parts a2 and a3: one E10_joint_fit.py process per wiring (claims by mkdir; skipping those done)."""
    raise NotImplementedError


def write_variants():
    """variants.csv and model_pair_distances_long.csv for the archive."""
    raise NotImplementedError


def chain(wirings: Sequence[str]) -> None:
    """check, a1, b on `wirings`, then score and variants."""
    raise NotImplementedError


if __name__ == '__main__':
    raise SystemExit('E124: stub only')
