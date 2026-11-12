"""pytest setup for the queue's tests: they import the repository's packages (kcq, experiments, ...) as installed
in the venv (`./venv/bin/pip install -e .` from the repository root; pyproject.toml). Run from the repository root
with QUICK=1 (CLAUDE.md, "Code changes")."""
import importlib.metadata

try:
    importlib.metadata.distribution('model_experiment_queue')
except importlib.metadata.PackageNotFoundError:
    # without it most tests would be skipped (importorskip) and the scripts they run by path would not import
    raise ImportError("the repository's packages are not installed in this Python: run `./venv/bin/pip install -e .` "
        'from the repository root') from None
