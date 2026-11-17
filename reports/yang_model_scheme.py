#!/usr/bin/env python3
"""The scheme of the baseline model (the bouton model) for Yang's version of the summary page: her request for "a
scheme showing the original version of the model, the connectivity, the transformation at each step" (2026-11-17).

One landscape PNG: the chain ORNs -> PNs -> PN boutons -> KC claws -> KC somata -> read-out as six numbered stage
boxes, the connectivity on the arrows between them, the LNs at the antennal lobe and APL at the calyx; one panel per
stage with its equation, one for APL, how the model is tuned, a legend and a footer. The wording was set with the
request; every number on the figure has a source:

- counts, read at run time (`load_counts`): the glomeruli, PNs, boutons and claws of the model's wiring from E29's
  glomerulus_substrate.csv (rows with model_n_claws > 0; ORN input where `driven`), the boutons, claws and KCs from
  E94's a1_path.csv (row wiring 'real'; its boutons and claws must equal E29's), and from the saved baseline run
  (E13, tag bouton_base) its odors (model_pair_distances_long.csv) and its response rates with APL on and off
  (stats_variants.csv);
- the simulator's parameters and literals and two settings of the experiments, in `CONSTANTS` (value, file, what to
  look for), checked against the olfsysm checkout's olfsysm.cpp (where pip installed olfsysm from, else the checkout
  beside this repository) and against experiments/E3_drive_dependent_gain.py and E4_substrate.py
  (`check_constants`): a value that differs from the table fails; without the olfsysm source the figure is drawn with
  the table's values, and the provenance and stdout say so.

Nothing here runs or imports the simulator, not even through kcq.harness (which imports it): the shared directory is
found as kcq.harness.shared_dir finds it. Writes <shared>/results/summary_yang/model_scheme.png and
model_scheme_provenance.json beside it (the inputs with sha256, the counts, the constants and whether each was
checked, the values printed on the figure, the commits, host and time).
Usage: ./venv/bin/python reports/yang_model_scheme.py [--out-dir DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import unquote, urlparse

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT_SUBDIR = Path('results') / 'summary_yang'
PNG_NAME = 'model_scheme.png'
PROV_NAME = 'model_scheme_provenance.json'
# the results tables the counts come from, relative to the shared directory
TABLES: Dict[str, str] = {
    'glomerulus_substrate': 'results/E29_oracle_substrate/glomerulus_substrate.csv',
    'a1_path': 'results/E94_bouton_share/a1_path.csv',
    'stats_variants': 'results/E13_bouton_coincidence/stats_variants.csv',
    'pair_distances': 'results/E13_bouton_coincidence/model_pair_distances_long.csv',
}
BASE_TAG = 'bouton_base'   # the saved baseline run in E13's tables
OLFSYSM_CPP = 'libolfsysm/src/olfsysm.cpp'   # relative to the olfsysm checkout


# ---------------------------------------------------------------- the constants and where they come from

_NUM = r'([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)'


@dataclass(frozen=True)
class Constant:
    """A number of the simulator or of the experiments that the figure shows (or computes a shown number from).

    Attributes:
        name: the name in the source ('p.pn.taum', 'TARGET'), or what an unnamed literal is.
        value: the value the figure is drawn with.
        root: 'olfsysm' (the olfsysm checkout) or 'queue' (this repository).
        file: the source file, relative to the root.
        pattern: a regular expression (multiline) matching the value exactly once in the file, the value as group 1.
        look_for: what to look for, for people.
        line: the line it was on when the table was written (2026-11-17); the check records where it is now.
    """
    name: str
    value: float
    root: str
    file: str
    pattern: str
    look_for: str
    line: int


def _default(name: str, value: float, line: int) -> Constant:
    """A field of olfsysm's DEFAULT_PARAMS, `<name> = <value>;` (the bouton base sets none of these)."""
    return Constant(name, value, 'olfsysm', OLFSYSM_CPP, rf'^\s*{re.escape(name)}\s*=\s*{_NUM}\s*;',
        f'`{name} = {value:g};` in DEFAULT_PARAMS', line)


CONSTANTS: Tuple[Constant, ...] = (
    _default('p.time.start', -0.5, 197),
    _default('p.time.end', 0.75, 198),
    _default('p.time.stim.start', 0.0, 199),
    _default('p.time.stim.end', 0.5, 200),
    _default('p.orn.taum', 0.01, 203),
    _default('p.ln.tauGA', 0.1, 210),
    _default('p.ln.tauGB', 0.4, 211),
    _default('p.pn.taum', 0.01, 216),
    _default('p.pn.offset', 2.9410, 217),
    _default('p.pn.tanhsc', 5.3395, 218),
    _default('p.pn.inhsc', 368.6631, 219),
    _default('p.pn.inhadd', 31.4088, 220),
    _default('p.kc.sp_factor_pre_apl', 2.0, 301),
    _default('p.kc.taum', 0.01, 327),
    _default('p.kc.apl_taum', 0.05, 345),
    Constant('PN tanh amplitude (literal)', 200.0, 'olfsysm', OLFSYSM_CPP,
        rf'dPNdt\s*\+=\s*{_NUM}\s*\*\s*\(\(orn_delta', 'sim_PN_layer: `dPNdt += 200.0*((orn_delta...`', 2601),
    Constant('PN tanh divisor (literal)', 200.0, 'olfsysm', OLFSYSM_CPP,
        rf'p\.pn\.tanhsc\s*/\s*{_NUM}\s*\*\s*inh_PN', 'sim_PN_layer: `p.pn.tanhsc/200.0*inh_PN`', 2601),
    Constant('PN gain: weight of inhA (literal)', 0.25, 'olfsysm', OLFSYSM_CPP,
        rf'p\.pn\.inhadd\s*\+\s*{_NUM}\s*\*\s*inhA\(t\)',
        'sim_PN_layer: `inh_PN = p.pn.inhsc/(p.pn.inhadd+0.25*inhA(t)+0.75*inhB(t))`', 2606),
    Constant('PN gain: weight of inhB (literal)', 0.75, 'olfsysm', OLFSYSM_CPP,
        rf'\*\s*inhA\(t\)\s*\+\s*{_NUM}\s*\*\s*inhB\(t\)',
        'sim_PN_layer: `inh_PN = p.pn.inhsc/(p.pn.inhadd+0.25*inhA(t)+0.75*inhB(t))`', 2606),
    Constant('KC threshold: factor on the spontaneous input (literal)', 2.0, 'olfsysm', OLFSYSM_CPP,
        rf'return\s+thr_const\s*\+\s*spont_in\.array\(\)\s*\*\s*{_NUM}\s*;',
        'choose_KC_thresh_uniform: `return thr_const + spont_in.array()*2.0;`', 1489),
    Constant('TARGET', 0.05, 'queue', 'experiments/E3_drive_dependent_gain.py', rf'^TARGET\s*=\s*{_NUM}',
        '`TARGET = 0.05` (the response rate tuned for, APL on)', 48),
    Constant('n_spikes_for_response', 2, 'queue', 'experiments/E4_substrate.py',
        rf'^CLAWS\s*=\s*dict\([^)]*?\bn_spikes_for_response\s*=\s*{_NUM}',
        '`n_spikes_for_response=2` in the CLAWS dict, which the bouton base (E3\'s BOUTON) extends', 45),
)


def constant_values(constants: Sequence[Constant] = CONSTANTS) -> Dict[str, float]:
    """{name: value} of the table."""
    return {c.name: c.value for c in constants}


def olfsysm_source_dir() -> Optional[Path]:
    """The olfsysm checkout holding the simulator's source: where pip installed the venv's olfsysm from (its
    direct_url.json, as pyolfsysm.mb_model finds it, without importing olfsysm), else the checkout beside this
    repository; None if neither has libolfsysm/src/olfsysm.cpp."""
    cands = []
    try:
        info = json.loads(distribution('olfsysm').read_text('direct_url.json') or '{}')
        url = urlparse(info.get('url', ''))
        if url.scheme == 'file':
            cands.append(Path(unquote(url.path)))
    except (PackageNotFoundError, ValueError):
        pass
    cands.append(REPO.parent / 'olfsysm')
    return next((c for c in cands if (c / OLFSYSM_CPP).is_file()), None)


def check_constants(constants: Sequence[Constant], roots: Mapping[str, Optional[Path]]) -> List[dict]:
    """Checks each constant against its source file.

    Args:
        constants: the table (`CONSTANTS`).
        roots: {root name: directory, or None where it was not found}.

    Returns:
        One dict per constant: name, value, source ('<root>:<file>'), look_for, line_listed, path, and status
        'checked' (the value is in the file, at `line`) or 'source not found' (no root, or no such file; line None).

    Raises:
        ValueError: the file exists but the pattern does not match it exactly once, or the value there differs from
            the table's.
    """
    texts: Dict[Path, str] = {}
    out = []
    for c in constants:
        root = roots.get(c.root)
        path = None if root is None else Path(root) / c.file
        row = dict(name=c.name, value=c.value, source=f'{c.root}:{c.file}', look_for=c.look_for, line_listed=c.line,
            path=None if path is None else str(path))
        if path is None or not path.is_file():
            out.append(dict(row, status='source not found', line=None))
            continue
        if path not in texts:
            texts[path] = path.read_text()
        text = texts[path]
        found = list(re.finditer(c.pattern, text, re.M))
        if len(found) != 1:
            raise ValueError(f'{c.name}: {len(found)} matches in {path}, expected one (look for {c.look_for})')
        m = found[0]
        line = text.count('\n', 0, m.start(1)) + 1
        if float(m.group(1)) != c.value:
            raise ValueError(f'{c.name} is {m.group(1)} in {path}:{line}, not {c.value:g} as CONSTANTS says: the '
                'figure would be wrong; update CONSTANTS (and the wording, if the change needs it)')
        out.append(dict(row, status='checked', line=line))
    return out


# ---------------------------------------------------------------- the counts

def _whole(x, what: str) -> int:
    x = float(x)
    if not x.is_integer():
        raise ValueError(f'{what} is {x}, not a whole number')
    return int(x)


def _one_row(df: pd.DataFrame, col: str, value: str, table: str) -> pd.Series:
    rows = df[df[col] == value]
    if len(rows) != 1:
        raise ValueError(f'{table}: {len(rows)} rows with {col} == {value!r}, expected one')
    return rows.iloc[0]


def load_counts(shared: Path) -> Dict[str, object]:
    """The counts on the figure, from the results tables `TABLES` in the shared directory.

    Returns:
        n_glomeruli (of the model's wiring: model_n_claws > 0), n_driven / n_undriven (of those, with / without ORN
        input), n_pn, n_boutons, n_claws, n_kcs, n_odors (of the saved baseline run), claws_per_kc (one decimal),
        claws_per_bouton (the nearest integer), and the saved baseline run's rate_on, rate_off and target_rate.

    Raises:
        ValueError: a table lacks the row it should have, a count is not a whole number, `driven` is not boolean, the
            boutons or claws of the two wiring tables disagree, or the saved run has another number of KCs.
    """
    shared = Path(shared)
    sub = pd.read_csv(shared / TABLES['glomerulus_substrate'])
    if sub['driven'].dtype != bool:
        raise ValueError(f"glomerulus_substrate.csv: the column driven is {sub['driven'].dtype}, not boolean")
    model = sub[sub['model_n_claws'] > 0]
    n_glomeruli, n_driven = len(model), int(model['driven'].sum())
    n_pn, sub_boutons, sub_claws = (_whole(model[c].sum(), f'glomerulus_substrate.csv: the sum of {c}')
        for c in ('model_n_pn', 'model_n_boutons', 'model_n_claws'))
    path = _one_row(pd.read_csv(shared / TABLES['a1_path']), 'wiring', 'real', 'a1_path.csv')
    n_boutons, n_claws, n_kcs = (_whole(path[c], f'a1_path.csv: {c}') for c in ('n_boutons', 'n_claws', 'n_kcs'))
    if (n_boutons, n_claws) != (sub_boutons, sub_claws):
        raise ValueError(f'the wiring tables disagree: a1_path.csv (wiring real) has {n_boutons} boutons and '
            f'{n_claws} claws, glomerulus_substrate.csv (model_n_claws > 0) {sub_boutons} and {sub_claws}')
    stats = _one_row(pd.read_csv(shared / TABLES['stats_variants']), 'tag', BASE_TAG, 'stats_variants.csv')
    if _whole(stats['n_kcs'], 'stats_variants.csv: n_kcs') != n_kcs:
        raise ValueError(f"the saved baseline run has {stats['n_kcs']} KCs, a1_path.csv {n_kcs}")
    pairs = pd.read_csv(shared / TABLES['pair_distances'])
    pairs = pairs[pairs['tag'] == BASE_TAG]
    if pairs.empty:
        raise ValueError(f'model_pair_distances_long.csv: no rows with tag {BASE_TAG!r}')
    return dict(n_glomeruli=n_glomeruli, n_driven=n_driven, n_undriven=n_glomeruli - n_driven, n_pn=n_pn,
        n_boutons=n_boutons, n_claws=n_claws, n_kcs=n_kcs, n_odors=len(set(pairs['odor_a']) | set(pairs['odor_b'])),
        claws_per_kc=round(n_claws / n_kcs, 1), claws_per_bouton=int(round(n_claws / n_boutons)),
        rate_on=float(stats['rate_on']), rate_off=float(stats['rate_off']), target_rate=float(stats['target_rate']))


# ---------------------------------------------------------------- the values printed on the figure

def _g(x: float) -> str:
    """x as short as it is exact (0.5, 10, 2), with a real minus sign."""
    return f'{round(float(x), 9):g}'.replace('-', '−')


def _frac(x: float) -> str:
    f = Fraction(x).limit_denominator(100)
    if float(f) != x:
        raise ValueError(f'{x} is not a simple fraction')
    return rf'\dfrac{{{f.numerator}}}{{{f.denominator}}}'


def figure_values(counts: Mapping[str, object], consts: Mapping[str, float]) -> Dict[str, str]:
    """Every number the figure prints, as printed ({name: text}), from the counts and the constants."""
    c = consts
    return dict(
        n_glomeruli=f"{counts['n_glomeruli']:,}", n_driven=f"{counts['n_driven']:,}",
        n_undriven=f"{counts['n_undriven']:,}", n_pn=f"{counts['n_pn']:,}", n_boutons=f"{counts['n_boutons']:,}",
        n_claws=f"{counts['n_claws']:,}", n_kcs=f"{counts['n_kcs']:,}", n_odors=f"{counts['n_odors']:,}",
        claws_per_bouton=f"{counts['claws_per_bouton']:,}", claws_per_kc=f"{counts['claws_per_kc']:.1f}",
        rate_on=f"{counts['rate_on']:.3f}", rate_off=f"{counts['rate_off']:.3f}",
        tau_orn=_g(1000 * c['p.orn.taum']), tau_pn=_g(1000 * c['p.pn.taum']), tau_kc=_g(1000 * c['p.kc.taum']),
        tau_apl=_g(1000 * c['p.kc.apl_taum']), tau_ga=_g(c['p.ln.tauGA']), tau_gb=_g(c['p.ln.tauGB']),
        odor_s=_g(c['p.time.stim.end'] - c['p.time.stim.start']),
        window_from=_g(c['p.time.start'] - c['p.time.stim.start']),
        window_to=_g(c['p.time.end'] - c['p.time.stim.start']),
        tanh_amp=_g(c['PN tanh amplitude (literal)']), tanh_div=_g(c['PN tanh divisor (literal)']),
        tanhsc=f"{c['p.pn.tanhsc']:.2f}", offset=f"{c['p.pn.offset']:.2f}", inhsc=f"{c['p.pn.inhsc']:.1f}",
        inhadd=f"{c['p.pn.inhadd']:.1f}", w_a=_frac(c['PN gain: weight of inhA (literal)']),
        w_b=_frac(c['PN gain: weight of inhB (literal)']),
        thr_factor=_g(c['KC threshold: factor on the spontaneous input (literal)']),
        n_spikes=_g(c['n_spikes_for_response']),
        pct_pre_apl=_g(100 * c['TARGET'] * c['p.kc.sp_factor_pre_apl']), pct_target=_g(100 * c['TARGET']))


# ---------------------------------------------------------------- the drawing

FIG_W, FIG_H, DPI = 18.4, 10.15, 180  # inches (near 16:9, so it fills a landscape screen); dots per inch
RC = {'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans', 'axes.unicode_minus': True}
PT = dict(title=15.5, region=11.5, name=12, count=10.5, conn=10, head=10.5, text=10.5, formula=12.5, small=10,
    foot=9.5)
INK, MUTED, CONN = '#1f1f1f', '#555555', '#707070'   # text, secondary text, connectivity labels (grey italic)
EXC = '#3b3b3b'   # excitation arrows, outlines, the number circles
INH, INH_FILL = '#D55E00', '#fbe8dc'   # inhibition (APL and the LNs): the one hue
BOX_FILL, PANEL_EDGE, PANEL_FILL = '#efefec', '#c6c6c1', '#fbfbfa'
LEFT = 0.3   # left and right margin
COLS = (2.25, 3.95, 2.85, 2.85, 2.6, 2.55)   # widths of the six columns of the transformation panels
COL_GAP = 0.15
BOX_W = 1.66   # stage boxes
BOX_TOP = FIG_H - 1.15   # the top of the stage boxes
PAD = 0.16   # inside panels
LINESPACING = 1.3


def _col_x(k: int) -> float:
    return LEFT + sum(COLS[:k]) + k * COL_GAP


def _col_c(k: int) -> float:
    return _col_x(k) + COLS[k] / 2


class _Canvas:
    """Draws on a figure in inches from its bottom-left corner, measures what it drew and keeps the texts, boxes
    and arrows for the layout check."""

    def __init__(self, fig):
        self.fig, self.tr = fig, fig.dpi_scale_trans
        self.r = fig.canvas.get_renderer()
        self.texts, self.boxes, self.arrows = [], [], []
        self.problems: List[str] = []   # what did not fit where it should go

    def text(self, x, y, s, size, **kw):
        kw.setdefault('color', INK)
        t = self.fig.text(x, y, s, fontsize=size, transform=self.tr, zorder=4, **kw)
        self.texts.append(t)
        return t

    def bbox(self, artist):
        return artist.get_window_extent(self.r).transformed(self.tr.inverted())

    def width(self, s, size, **font) -> float:
        from matplotlib.font_manager import FontProperties
        w, _, _ = self.r.get_text_width_height_descent(s, FontProperties(size=size, **font), ismath=s.count('$') >= 2)
        return w / self.fig.dpi

    def line_height(self, size) -> float:
        """The height (inches) of one line of text at `size`, as its bounding box has it."""
        from matplotlib.font_manager import FontProperties
        _, h, _ = self.r.get_text_width_height_descent('lp', FontProperties(size=size), ismath=False)
        return h / self.fig.dpi

    def wrap(self, s, width, size, **font) -> List[str]:
        """`s` broken greedily into lines no wider than `width` (inches); $...$ stays whole; '\\n' breaks a line."""
        lines = []
        for para in s.split('\n'):
            line = ''
            for tok in _tokens(para):
                trial = f'{line} {tok}' if line else tok
                if line and self.width(trial, size, **font) > width:
                    lines.append(line)
                    line = tok
                else:
                    line = trial
            lines.append(line)
        for line in lines:
            if self.width(line, size, **font) > width + 1e-6:
                self.problems.append(f'too wide for {width:.2f} in at {size} pt: {line!r}')
        return lines

    def box(self, x0, y0, x1, y1, *, fill, edge, lw=1.0, rounding=0.06, margin=0.05, zorder=1):
        """A rounded box; a text inside it must keep `margin` (inches) from its edges (`_layout_problems`)."""
        from matplotlib.patches import FancyBboxPatch
        p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f'round,pad=0,rounding_size={rounding}',
            transform=self.tr, facecolor=fill, edgecolor=edge, linewidth=lw, zorder=zorder)
        self.fig.add_artist(p)
        self.boxes.append((x0, y0, x1, y1, margin))
        return p

    def circle(self, x, y, r, label, *, fill=EXC, size=10):
        from matplotlib.patches import Circle
        self.fig.add_artist(Circle((x, y), r, transform=self.tr, facecolor=fill, edgecolor='none', zorder=3))
        self.boxes.append((x - r, y - r, x + r, y + r, 0.0))
        self.text(x, y - 0.004, label, size, ha='center', va='center', color='white', weight='bold')

    def arrow(self, xy0, xy1, *, inhibit=False, lw=1.6, color=None):
        from matplotlib.patches import ArrowStyle, FancyArrowPatch
        style = (ArrowStyle('-[', widthB=0.6, lengthB=0) if inhibit else
            ArrowStyle('-|>', head_length=0.7, head_width=0.35))
        p = FancyArrowPatch(xy0, xy1, transform=self.tr, arrowstyle=style, mutation_scale=14,
            color=color or (INH if inhibit else EXC), linewidth=lw, shrinkA=0, shrinkB=0, zorder=3)
        self.fig.add_artist(p)
        self.arrows.append(p)
        return p


def _tokens(s: str) -> List[str]:
    """`s` split at spaces outside $...$."""
    out, cur, math = [], '', False
    for ch in s:
        if ch == '$':
            math = not math
        if ch == ' ' and not math:
            if cur:
                out.append(cur)
            cur = ''
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


class _Flow:
    """Text placed top-down in a column (`x`, `width`) from `top`; `y` is where the next item starts."""

    def __init__(self, cv: _Canvas, x: float, top: float, width: float):
        self.cv, self.x, self.y, self.width = cv, x, top, width

    def par(self, s, *, size=PT['text'], gap=0.08, indent=0.0, **kw):
        font = {k: kw[k] for k in ('style', 'weight') if k in kw}
        lines = self.cv.wrap(s, self.width - indent, size, **font)
        t = self.cv.text(self.x + indent, self.y, '\n'.join(lines), size, va='top', linespacing=LINESPACING, **kw)
        self.y = self.cv.bbox(t).y0 - gap
        return t

    def item(self, mark, s, *, gap=0.08):
        """A list item with a hanging indent (`mark`, e.g. '1.')."""
        self.cv.text(self.x, self.y, mark, PT['text'], va='top')
        self.par(s, gap=gap, indent=0.24)

    def formula(self, tex, *, note=None, gap=0.1, indent=0.0):
        """A formula (mathtext, at the formula size); `note` beside it if it fits, else under it."""
        t = self.cv.text(self.x + indent, self.y, tex, PT['formula'], va='top')
        b = self.cv.bbox(t)
        if b.x1 > self.x + self.width + 1e-6:
            self.cv.problems.append(f'formula {b.width:.2f} in wide, {self.width - indent:.2f} in available: {tex}')
        self.y = b.y0 - gap
        if note:
            if b.x1 + 0.2 + self.cv.width(note, PT['text']) <= self.x + self.width:
                self.cv.text(b.x1 + 0.2, (b.y0 + b.y1) / 2, note, PT['text'], va='center', color=MUTED)
            else:
                self.y += gap - 0.05
                self.par(note, color=MUTED, gap=gap)
        return t


def _short(p0: Tuple[float, float], p1: Tuple[float, float], d: float) -> Tuple[float, float]:
    """`p1` moved `d` inches back towards `p0`."""
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    n = (dx * dx + dy * dy) ** 0.5
    return p1[0] - d * dx / n, p1[1] - d * dy / n


def _panel(cv: _Canvas, x: float, top: float, width: float, mark: str, heading: str) -> _Flow:
    """A transformation panel's mark (a stage number in a circle, or the APL chip) and bold heading; returns the flow
    for its content. The panel's box is drawn later (`_panel_box`), once the column's height is known."""
    r = 0.14
    if mark == 'APL':
        w = cv.width('APL', PT['small'], weight='bold') + 0.16
        cv.box(x + PAD, top - PAD - 2 * r, x + PAD + w, top - PAD, fill=INH, edge=INH, rounding=0.05, zorder=3)
        cv.text(x + PAD + w / 2, top - PAD - r, 'APL', PT['small'], ha='center', va='center', color='white',
            weight='bold')
        indent = w + 0.14
    else:
        cv.circle(x + PAD + r, top - PAD - r, r, mark)
        indent = 2 * r + 0.14
    f = _Flow(cv, x + PAD + indent, top - PAD - r + 0.1, width - 2 * PAD - indent)
    f.par(heading, size=PT['head'], weight='bold', gap=0.14)
    f.x, f.width = x + PAD, width - 2 * PAD
    f.y = min(f.y, top - PAD - 2 * r - 0.14)
    return f


def _panel_box(cv: _Canvas, x: float, top: float, width: float, bottom: float, *, edge=PANEL_EDGE):
    cv.box(x, bottom, x + width, top, fill=PANEL_FILL, edge=edge, lw=0.9, zorder=0)


def _chain(cv: _Canvas, v: Mapping[str, str]) -> float:
    """The stage boxes, the connectivity on the arrows between them, the LNs and APL. Returns the bottom of what
    was drawn."""
    from matplotlib.lines import Line2D
    stages = [('ORNs', [f"{v['n_glomeruli']} glomeruli", f"odor responses measured for {v['n_driven']}"]),
        ('PNs', [f"{v['n_pn']} PNs", 'one rate per glomerulus']),
        ('PN boutons', [f"{v['n_boutons']} boutons"]),
        ('KC claws', [f"{v['n_claws']} claws"]),
        ('KC somata', [f"{v['n_kcs']} KCs"]),
        ('Read-out', ['spike counts', f"{v['n_odors']}\u00a0odors × {v['n_kcs']}\u00a0KCs"])]   # (no-break spaces)
    r = 0.14
    bottoms = []
    for k, (name, lines) in enumerate(stages):
        x0 = _col_c(k) - BOX_W / 2
        cv.circle(x0 + 0.1 + r, BOX_TOP - 0.12 - r, r, str(k + 1))
        cv.text(x0 + 0.1 + 2 * r + 0.1, BOX_TOP - 0.12 - r, name, PT['name'], va='center', weight='bold')
        f = _Flow(cv, x0 + 0.14, BOX_TOP - 0.12 - 2 * r - 0.12, BOX_W - 0.28)
        for line in lines:
            f.par(line, size=PT['count'], gap=0.04)
        bottoms.append(f.y)
    box_bottom = min(bottoms) - 0.1
    for k in range(6):
        x0 = _col_c(k) - BOX_W / 2
        cv.box(x0, box_bottom, x0 + BOX_W, BOX_TOP, fill=BOX_FILL, edge=EXC, lw=1.2, rounding=0.08)
    # the regions, as brackets over the boxes
    for k0, k1, label in ((0, 1, 'antennal lobe'), (2, 4, 'mushroom body')):
        xa, xb, y = _col_c(k0) - BOX_W / 2, _col_c(k1) + BOX_W / 2, BOX_TOP + 0.2
        cv.fig.add_artist(Line2D([xa, xa, xb, xb], [y - 0.09, y, y, y - 0.09], transform=cv.tr, color=MUTED,
            linewidth=1.1))
        cv.text((xa + xb) / 2, y + 0.07, label, PT['region'], ha='center', va='bottom', color=MUTED)
    # the connectivity between the stages, on the arrows
    conn = [("within each glomerulus", None), ('a bouton carries the rate of its own glomerulus', None),
        ('each claw is fed by one bouton, weight 1', f"(about {v['claws_per_bouton']} claws per bouton)"),
        ('each KC sums its own claws', f"({v['claws_per_kc']} per KC on average)"), ('spikes', None)]
    y_arrow = (BOX_TOP + box_bottom) / 2 - 0.05
    for k, (above, below) in enumerate(conn):
        xa, xb = _col_c(k) + BOX_W / 2, _col_c(k + 1) - BOX_W / 2
        cv.arrow((xa + 0.06, y_arrow), (xb - 0.06, y_arrow))
        wrap = dict(size=PT['conn'], style='italic')
        lines = cv.wrap(above, xb - xa - 0.16, **wrap)
        t = cv.text((xa + xb) / 2, y_arrow + 0.08, '\n'.join(lines), PT['conn'], ha='center', va='bottom',
            style='italic', color=CONN, linespacing=1.2)
        if cv.bbox(t).y1 > BOX_TOP + 0.05:
            cv.problems.append(f'the connectivity label {above!r} rises above the boxes')
        if below:
            lines = cv.wrap(below, xb - xa - 0.16, **wrap)
            t = cv.text((xa + xb) / 2, y_arrow - 0.08, '\n'.join(lines), PT['conn'], ha='center', va='top',
                style='italic', color=CONN, linespacing=1.2)
            if cv.bbox(t).y0 < box_bottom - 0.05:
                cv.problems.append(f'the connectivity label {below!r} drops below the boxes')
    node_h, gap = 0.62, 0.62
    node_top = box_bottom - gap
    node_y = node_top - node_h / 2
    # the LNs, between the ORNs and the PNs
    b1, b2 = _col_c(0) + BOX_W / 2, _col_c(1) - BOX_W / 2
    lx, lw_ = (b1 + b2) / 2, 0.8
    cv.box(lx - lw_ / 2, node_y - 0.22, lx + lw_ / 2, node_y + 0.22, fill=INH_FILL, edge=INH, lw=1.4, rounding=0.1)
    cv.text(lx, node_y, 'LNs', PT['name'], ha='center', va='center', weight='bold')
    cv.arrow((b1 - 0.3, box_bottom), (lx - 0.18, node_y + 0.22))
    cv.arrow((lx + 0.18, node_y + 0.22), _short((lx + 0.18, node_y + 0.22), (b2 + 0.3, box_bottom), 0.06),
        inhibit=True)
    lines = cv.wrap('LNs: driven by the mean ORN rate over all glomeruli; they divide the gain of every PN', 3.2,
        PT['conn'], style='italic')
    t = cv.text(lx, node_y - 0.22 - 0.1, '\n'.join(lines), PT['conn'], ha='center', va='top', style='italic',
        color=CONN, linespacing=1.2)
    low = cv.bbox(t).y0
    # APL, under the boutons and the claws: drive down, inhibition up, beside each other under each box
    c3, c4 = _col_c(2), _col_c(3)
    ax0, ax1 = c3 - 0.45, c4 + 0.45
    cv.box(ax0, node_top - node_h, ax1, node_top, fill=INH_FILL, edge=INH, lw=1.4, rounding=0.1)
    cv.text((ax0 + ax1) / 2, node_y + 0.11, 'APL', PT['name'], ha='center', va='center', weight='bold')
    cv.text((ax0 + ax1) / 2, node_y - 0.14, 'one neuron', PT['count'], ha='center', va='center')
    for c, out in ((c3, -1), (c4, 1)):
        cv.arrow((c + 0.14 * out, box_bottom), (c + 0.14 * out, node_top))   # drive: stage -> APL (outer)
        cv.arrow((c - 0.14 * out, node_top), (c - 0.14 * out, box_bottom - 0.06), inhibit=True)   # inhibition (inner)
    labels = [(False, 'drive: bouton→APL and claw→APL synapse counts'),
        (True, 'inhibition: APL→bouton and APL→claw synapse counts')]
    # beside APL: the two labels, each after a short arrow of its kind, and the note under them
    kx = ax1 + 0.25
    tx = kx + 0.58
    tw = max(cv.width(s, PT['conn'], style='italic') for _, s in labels)
    f = _Flow(cv, tx, node_top, tw)
    for inhibit, s in labels:
        t = f.par(s, size=PT['conn'], style='italic', color=CONN, gap=0.1)
        y = cv.bbox(t).y1 - cv.line_height(PT['conn']) / 2
        cv.arrow((kx, y), (kx + 0.44, y), inhibit=inhibit, lw=1.4)
    note = ('All four sets of weights are synapse counts of the hemibrain (right calyx), each scaled to a mean of 1, '
        'times one gain $\\gamma$.')
    f.x, f.width = kx, tx + tw - kx
    f.y -= 0.04
    f.par(note, size=PT['small'], color=INK)
    return min(low, node_top - node_h, f.y)


def _panels(cv: _Canvas, v: Mapping[str, str], top: float, foot_top: float):
    """The six transformation panels, APL's, the tuning panel and the legend, in columns from `top`; the columns
    end together, above `foot_top`."""
    x = [_col_x(k) for k in range(6)]
    tau = r'$\tau$'
    ends = {}
    f = _panel(cv, x[0], top, COLS[0], '1', 'ORN rate of glomerulus $g$')
    f.formula(r'$\tau\,\dot{r}_g = -r_g + s_g + \Delta_g$', note=f"{tau} = {v['tau_orn']} ms")
    f.par(rf"$\Delta_g$: the odor's measured change of rate, applied as a step for the {v['odor_s']} s of the odor "
        f"(smoothed), and 0 for the {v['n_undriven']} glomeruli without data.")
    f.par('$s_g$: spontaneous rate.')
    ends[0] = f.y
    f = _panel(cv, x[1], top, COLS[1], '2', 'PN rate of glomerulus $g$: saturation, with a gain set by lateral '
        'inhibition')
    f.formula(r'$\tau\,\dot{P}_g = -P_g + P^0_g$', gap=0.04)
    f.formula(rf"$+\ {v['tanh_amp']}\,\tanh\left[\dfrac{{{v['tanhsc']}}}{{{v['tanh_div']}}}\,G_{{\mathrm{{LN}}}}(t)\,"
        rf"(r_g - s_g + {v['offset']})\right]$", indent=0.2, gap=0.06)
    f.par(f"floored at 0; {tau} = {v['tau_pn']} ms", color=MUTED, gap=0.14)
    f.formula(rf"$G_{{\mathrm{{LN}}}}(t) = \dfrac{{{v['inhsc']}}}{{{v['inhadd']} + {v['w_a']}A(t) + {v['w_b']}B(t)}}$",
        gap=0.12)
    f.par(f"$A$, $B$: the LNs' response, filtered over {v['tau_ga']} s and {v['tau_gb']} s. The LNs are driven by "
        'the cube of the mean ORN rate over all glomeruli: the more of the antennal lobe an odor drives, the lower the '
        'gain of every PN.', gap=0.14)
    f.formula(rf"$P^0_g = \dfrac{{{v['inhsc']}\,s_g}}{{{v['inhadd']} + \sum_{{g\prime}}\, s_{{g\prime}}}}$",
        note="(the PN's spontaneous rate)")
    ends[1] = f.y
    f = _panel(cv, x[2], top, COLS[2], '3', 'Rate of bouton $j$, of glomerulus $g$')
    f.formula(r'$b_j = \max\left(P_g - \gamma\,w_j^{\mathrm{APL}\to b}\,a,\ 0\right)$')
    f.par("a subtraction, floored at 0; $a$: APL's activity")
    ends[2] = f.y
    f = _panel(cv, x[3], top, COLS[3], '4', 'Drive of claw $i$, on bouton $j$')
    f.formula(r'$c_i = \max\left(b_j - \gamma\,w_i^{\mathrm{APL}\to c}\,a,\ 0\right)$')
    f.par('a subtraction, floored at 0, at each claw separately')
    ends[3] = f.y
    f = _panel(cv, x[4], top, COLS[4], '5', 'KC $k$: sum, leaky integration, threshold')
    f.formula(r'$\tau\,\dot{V}_k = -V_k + \sum_{i\in k}\, c_i$', note=f"{tau} = {v['tau_kc']} ms")
    f.par(r'a spike when $V_k > \theta_k$, then $V_k$ is set to 0 (no refractory period, no adaptation)', gap=0.12)
    f.par(rf"$\theta_k = \theta + {v['thr_factor']}\times\,$(the KC's summed claw input before the odor)")
    ends[4] = f.y
    f = _panel(cv, x[5], top, COLS[5], '6', 'Read-out')
    f.par(f"$n_k$: the KC's spikes in the simulated window ({v['window_from']} to {v['window_to']} s around the "
        "odor's\u00a0onset)", gap=0.12)
    f.par(rf"a response: $n_k \geq {v['n_spikes']}$", gap=0.12)
    f.par("odor-odor correlation: the cosine similarity of the two odors' spike counts across the\u00a0KCs")
    ends[5] = f.y
    row = {0: ends[0] - PAD, 1: ends[1] - PAD, 2: min(ends[2], ends[3]) - PAD, 3: min(ends[2], ends[3]) - PAD,
        4: min(ends[4], ends[5]) - PAD, 5: min(ends[4], ends[5]) - PAD}
    for k in (0, 2, 3, 4, 5):
        _panel_box(cv, x[k], top, COLS[k], row[k])
    # under panels 3 and 4: APL; under 5 and 6: the tuning; under 1: the legend
    xa, wa = x[2], COLS[2] + COL_GAP + COLS[3]
    ta = row[2] - COL_GAP
    f = _panel(cv, xa, ta, wa, 'APL', "APL's activity $a$")
    f.formula(r'$\tau_A\,\dot{a} = -a + \dfrac{\gamma}{N_{\mathrm{KC}}}\sum_i\, w_i^{c\to\mathrm{APL}}\,c_i'
        r' + \dfrac{\gamma}{N_b}\sum_j\, w_j^{b\to\mathrm{APL}}\,b_j$', gap=0.08)
    f.par(rf"$\tau_A$ = {v['tau_apl']} ms", color=MUTED, gap=0.1)
    f.par("APL sums what its own inhibition leaves of the claws' and the boutons' drive (a feedback loop), and "
        'inhibits every claw and every bouton with the one value $a$.', gap=0.1)
    f.par(rf"$N_{{\mathrm{{KC}}}}$ = {v['n_kcs']}, $N_b$ = {v['n_boutons']}", size=PT['small'], color=MUTED)
    end_apl = f.y
    xt, wt = x[4], COLS[4] + COL_GAP + COLS[5]
    tt = row[4] - COL_GAP
    f = _Flow(cv, xt + PAD, tt - PAD, wt - 2 * PAD)
    f.par(f"Tuning: two numbers, on all {v['n_odors']} odors", size=PT['head'], weight='bold', gap=0.12)
    f.item('1.', r"$\theta$, with APL's weights at zero:" + f"\n{v['pct_pre_apl']} % of the KC × odor entries respond.")
    f.item('2.', rf"$\gamma$, the one gain on APL's four sets of weights: {v['pct_target']} % respond.", gap=0.12)
    f.par("APL off (the acute knock-out): APL's output weights are set to zero and $\\theta$ is left as it is, so "
        f"{v['pct_pre_apl']} % respond again.", gap=0.12)
    f.par(f"In the saved run: {v['rate_on']} of the entries with APL on,\n{v['rate_off']} with APL off.")
    end_tuning = f.y
    tl = row[0] - COL_GAP
    end_legend = _legend(cv, x[0] + PAD, tl - PAD, COLS[0] - 2 * PAD)
    bottom = min(ends[1], end_apl, end_tuning, end_legend) - PAD
    if bottom < foot_top:
        cv.problems.append(f'the panels end at {bottom:.2f} in, below the footer ({foot_top:.2f} in): make the figure '
            'taller')
    _panel_box(cv, x[1], top, COLS[1], bottom)
    _panel_box(cv, xa, ta, wa, bottom, edge=INH)
    _panel_box(cv, xt, tt, wt, bottom)
    _panel_box(cv, x[0], tl, COLS[0], bottom)
    return bottom


def _legend(cv: _Canvas, x: float, top: float, width: float) -> float:
    """The legend's three entries from `top`; returns where they end."""
    key, gap = 0.55, 0.16
    f = _Flow(cv, x + key + gap, top, width - key - gap)
    for kind, s in (('exc', 'excitation'), ('inh', 'inhibition (APL: subtraction; LNs: division of the gain)'),
            ('conn', 'connectivity')):
        style = {'style': 'italic', 'color': CONN} if kind == 'conn' else {}
        t = f.par(s, gap=0.14, **style)
        first = cv.bbox(t).y1 - cv.line_height(PT['text']) / 2   # the middle of its first line
        if kind == 'conn':
            cv.text(x + key / 2, first, 'Aa', PT['text'], ha='center', va='center', **style)
        else:
            cv.arrow((x, first), (x + key, first), inhibit=kind == 'inh')
    return f.y


def _layout_problems(cv: _Canvas) -> List[str]:
    """Texts that overlap each other, cross a box's edge, touch an arrow or leave the figure."""
    from matplotlib.path import Path as MPath
    # each arrow as polylines in display units (its curves flattened, CLOSEPOLY resolved)
    arrows = [[MPath(q) for q in a.get_transform().transform_path(a.get_path()).to_polygons(closed_only=False)]
        for a in cv.arrows]
    items = [(t.get_text().split('\n')[0][:40], cv.bbox(t)) for t in cv.texts if t.get_text().strip()]
    out = []
    eps = 0.012
    for i, (s, b) in enumerate(items):
        if b.x0 < 0.05 or b.y0 < 0.05 or b.x1 > FIG_W - 0.05 or b.y1 > FIG_H - 0.05:
            out.append(f'outside the figure: {s!r}')
        for s2, b2 in items[i + 1:]:
            if min(b.x1, b2.x1) - max(b.x0, b2.x0) > eps and min(b.y1, b2.y1) - max(b.y0, b2.y0) > eps:
                out.append(f'texts overlap: {s!r} and {s2!r}')
        for (x0, y0, x1, y1, m) in cv.boxes:
            overlaps = min(b.x1, x1) - max(b.x0, x0) > eps and min(b.y1, y1) - max(b.y0, y0) > eps
            inside = b.x0 >= x0 + m - eps and b.x1 <= x1 - m + eps and b.y0 >= y0 + m - eps and b.y1 <= y1 - m + eps
            if overlaps and not inside:
                out.append(f'text crosses a box edge or is within {m} in of it: {s!r}')
        disp = b.transformed(cv.tr)
        if any(q.intersects_bbox(disp, filled=False) for polys in arrows for q in polys):
            out.append(f'text touches an arrow: {s!r}')
    return out


def draw(v: Mapping[str, str], out_path: Path, *, dpi: int = DPI, check: bool = True) -> Tuple[int, int]:
    """Draws the scheme with the printed values `v` (`figure_values`) to `out_path` (a PNG).

    Args:
        v: the values, as printed.
        out_path: the PNG to write (its folder is made).
        dpi: dots per inch.
        check: raise on a layout problem (below); if False, print the problems and write the PNG anyway (for a look
            at a draft).

    Returns:
        (width, height) in pixels.

    Raises:
        RuntimeError: something does not fit where it should go; a text overlaps another, crosses a box edge, touches
            an arrow or leaves the figure; or a value of `v` is not on the figure.
    """
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    with matplotlib.rc_context(RC):
        fig = Figure(figsize=(FIG_W, FIG_H), dpi=dpi, facecolor='white')
        FigureCanvasAgg(fig)
        cv = _Canvas(fig)
        cv.text(LEFT, FIG_H - 0.25, 'The baseline model (the bouton model): its connectivity and the transformation '
            'at each step', PT['title'], va='top', weight='bold')
        low = _chain(cv, v)
        foot_y = 0.22
        _panels(cv, v, low - 0.2, foot_y + 0.3)
        cv.text(LEFT, foot_y, 'Equations as in the simulator (olfsysm) with the settings of the bouton baseline; '
            "counts from the model's own wiring tables.", PT['foot'], va='bottom', color=MUTED)
        shown = '\n'.join(t.get_text() for t in cv.texts)
        problems = cv.problems + _layout_problems(cv) + [f'not on the figure: {k} ({s})' for k, s in v.items()
            if s not in shown]
        if problems and check:
            raise RuntimeError('layout:\n' + '\n'.join(problems))
        for msg in problems:
            print('layout:', msg)
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=dpi, facecolor='white')
    return int(round(FIG_W * dpi)), int(round(FIG_H * dpi))


# ---------------------------------------------------------------- running it

def shared_dir() -> Path:
    """The shared directory, found as kcq.harness.shared_dir finds it (KC_REMAPPING_SHARED, else
    Tom/claude_kc_remapping under /local, /mnt or /media tensor), without importing the model's packages."""
    for c in [os.environ.get('KC_REMAPPING_SHARED')] + [
            f'/{root}/tensor/Tom/claude_kc_remapping' for root in ('local', 'mnt', 'media')]:
        if c and Path(c).is_dir():
            return Path(c)
    raise FileNotFoundError('shared directory not found; set KC_REMAPPING_SHARED')


def _sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(shared: Path, out_dir: Path, *, roots: Optional[Mapping[str, Optional[Path]]] = None) -> dict:
    """Reads the counts, checks the constants, draws the scheme to `out_dir`/PNG_NAME and writes PROV_NAME beside it.

    Args:
        shared: the shared directory (the tables are read from it).
        out_dir: where the PNG and its provenance go.
        roots: {'olfsysm': checkout or None, 'queue': this repository} for the constants' sources; default the
            olfsysm checkout `olfsysm_source_dir()` finds and this repository.

    Returns:
        The provenance written.
    """
    from hong2p.util import git_commits
    shared, out_dir = Path(shared), Path(out_dir)
    counts = load_counts(shared)
    if roots is None:
        roots = {'olfsysm': olfsysm_source_dir(), 'queue': REPO}
    checks = check_constants(CONSTANTS, roots)
    consts = constant_values(CONSTANTS)
    if abs(counts['target_rate'] - consts['TARGET']) > 1e-12:
        raise ValueError(f"the saved baseline run was tuned to {counts['target_rate']}, not TARGET "
            f"({consts['TARGET']})")
    unchecked = [c['name'] for c in checks if c['status'] != 'checked']
    if unchecked:
        print('simulator source not found: drawn with the values of CONSTANTS, not checked: ' + ', '.join(unchecked))
    else:
        print(f'all {len(checks)} constants checked against their sources')
    values = figure_values(counts, consts)
    png = out_dir / PNG_NAME
    width, height = draw(values, png)
    olf = roots.get('olfsysm')
    sources = {}
    for c in checks:
        if c['path'] and c['status'] == 'checked':
            sources[c['source']] = dict(path=c['path'], sha256=_sha256(Path(c['path'])))
    prov = dict(script='reports/yang_model_scheme.py',
        **git_commits({'queue': REPO, **({'olfsysm': olf} if olf is not None else {})}),
        inputs={k: dict(path=str(shared / p), sha256=_sha256(shared / p)) for k, p in TABLES.items()},
        counts=counts, constants_checked=not unchecked, constants_unchecked=unchecked,
        olfsysm_source=None if olf is None else str(olf), constant_sources=sources, constants=checks,
        values_on_figure=values,
        figure=dict(path=str(png), width_px=width, height_px=height, dpi=DPI, size_in=[FIG_W, FIG_H],
            sha256=_sha256(png)),
        generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), host=socket.gethostname())
    (out_dir / PROV_NAME).write_text(json.dumps(prov, indent=1, default=str) + '\n')
    print('values on the figure: ' + ', '.join(f'{k}={s}' for k, s in values.items()))
    print(f'written: {png} ({width} x {height} px, {png.stat().st_size / 1024:.0f} KB) and {out_dir / PROV_NAME}')
    return prov


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--out-dir', type=Path, default=None,
        help='write the PNG and its provenance here instead of <shared>/results/summary_yang/')
    args = ap.parse_args(argv)
    shared = shared_dir()
    run(shared, args.out_dir if args.out_dir is not None else shared / OUT_SUBDIR)


if __name__ == '__main__':
    main()
