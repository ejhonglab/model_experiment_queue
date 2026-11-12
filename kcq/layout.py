"""The repository's layout, and the flat one it had before 2026-11-12.

Until then every script sat at the repository's root (and, until 2026-11-11, in al_analysis's
scripts/claude_pair_remapping_experiments) and imported the others as top-level modules
(`import kc_remapping_harness as H`, `from E10_joint_fit import ...`). Since then they are modules of
the packages below, installed editable (pyproject.toml), and import each other absolutely
(`from kcq import harness as H`, `from experiments.E10_joint_fit import ...`). `package_imports`
rewrites the imports of a file from before the move, so that an old version of a script (a
baseline run of E78_claw_knee.py's check-e10) runs against the current packages.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

# the repository's root
ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ('kcq', 'experiments', 'data_prep', 'analyses', 'reports', 'archive_tools', 'coordination')
# the modules whose names changed with the move (the others kept their names)
RENAMED = {'kc_remapping_harness': 'kcq.harness', 'apl_effect_scores': 'kcq.apl_effect',
    'noise_floor_trial_distances': 'kcq.noise_floor'}


def new_module(old: str) -> Optional[str]:
    """The dotted module name today of a module that was flat (`old`, e.g. 'E10_joint_fit' ->
    'experiments.E10_joint_fit'), or None if `old` is no module of this repository."""
    if old in RENAMED:
        return RENAMED[old]
    for pkg in PACKAGES:
        if (ROOT / pkg / f'{old}.py').is_file():
            return f'{pkg}.{old}'
    return None


_IMPORT = re.compile(r'^(?P<ind>[ \t]*)import (?P<old>\w+)(?: as (?P<alias>\w+))?(?P<rest>[ \t]*(?:#.*)?)$', re.M)
_FROM = re.compile(r'^(?P<ind>[ \t]*)from (?P<old>\w+) import ', re.M)


def package_imports(src: str) -> str:
    """`src` (a script's source from before the move) with its imports of this repository's flat
    modules made absolute: `import X as Y` -> `from <pkg> import <X today> as Y`, `import X` ->
    `from <pkg> import <X today>` (bound to the old name if it changed), `from X import ...` ->
    `from <pkg>.<X today> import ...`. Other imports (one module per line, as the scripts write
    them) are left alone."""
    def one(m):
        new = new_module(m['old'])
        if new is None:
            return m[0]
        pkg, mod = new.rsplit('.', 1)
        alias = m['alias'] or (m['old'] if mod != m['old'] else None)
        return f"{m['ind']}from {pkg} import {mod}{f' as {alias}' if alias else ''}{m['rest']}"

    def frm(m):
        new = new_module(m['old'])
        return m[0] if new is None else f"{m['ind']}from {new} import "

    return _FROM.sub(frm, _IMPORT.sub(one, src))
