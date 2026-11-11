"""pytest setup for the queue's tests: the experiment scripts and the harness are flat modules at
the repository root (no package), so the root goes on sys.path for `import kc_remapping_harness`
and the like. Run from the repository root with QUICK=1 (CLAUDE.md, "Code changes")."""
import sys
from pathlib import Path

# the repository root, where the scripts are
EXP_DIR = Path(__file__).resolve().parents[1]
if str(EXP_DIR) not in sys.path:
    sys.path.insert(0, str(EXP_DIR))
