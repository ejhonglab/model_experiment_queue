#!/usr/bin/env python3
"""Tests for archive_tools/migrate_model_output_params.py's command line: a named archive folder without a params.csv is
skipped with a line, not an error, and a name with no folder at all (a typo) exits with a message. Needs the shared
directory (the module reads it at import); skipped without it."""

import sys

import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.shared_dir()
    M = pytest.importorskip('archive_tools.migrate_model_output_params')
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory is not available', allow_module_level=True)


def test_folder_without_params_is_skipped(tmp_path, monkeypatch, capsys):
    (tmp_path / 'E0_spike_times_only').mkdir()
    monkeypatch.setattr(M, 'ARCHIVE', tmp_path)
    monkeypatch.setattr(sys, 'argv', ['migrate_model_output_params.py', '--apply', 'E0_spike_times_only'])
    M.main()
    out = capsys.readouterr().out
    assert 'E0_spike_times_only: no params.csv, nothing to fill' in out
    assert not (tmp_path / 'params_migration_log.csv').exists()


def test_missing_folder_exits(tmp_path, monkeypatch):
    monkeypatch.setattr(M, 'ARCHIVE', tmp_path)
    monkeypatch.setattr(sys, 'argv', ['migrate_model_output_params.py', 'E0_typo'])
    with pytest.raises(SystemExit) as e:
        M.main()
    assert 'no archive folder for E0_typo' in str(e.value.code)
