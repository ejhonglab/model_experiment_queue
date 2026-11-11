#!/bin/bash
# Copies the queue's tests from al_analysis's test/ folder into this repository's test/, as they
# are at a ref of al_analysis (default origin/main; fetch it first), with their path to the
# scripts rewritten: they located the scripts as <al_analysis>/scripts/claude_pair_remapping_experiments
# and find them at this repository's root instead (test/conftest.py also puts it on sys.path);
# docstrings naming that folder name the scripts alone.
# A test is the queue's if it names that folder. Used once for the move (2026-11-11) and again at
# the cutover, for tests added to al_analysis in the meantime; prints what it copied, and what it
# would overwrite with a different version (which it does, so check `git diff test/` afterwards).
#
# Usage: migration/copy_al_analysis_tests.sh [<al_analysis checkout> [<ref>]]
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
al=${1:-$(cd "$here/.." && pwd)/al_analysis}
ref=${2:-origin/main}
git -C "$al" rev-parse --verify -q "$ref^{commit}" >/dev/null || { echo "no $ref in $al" >&2; exit 1; }
mkdir -p "$here/test"
n=0
for f in $(git -C "$al" ls-tree --name-only "$ref" test/ | grep -E '^test/test_.*\.py$'); do
  src=$(git -C "$al" show "$ref:$f")
  grep -q claude_pair_remapping_experiments <<<"$src" || continue
  dst=$here/$f
  new=$(sed -e "s| / 'scripts' / 'claude_pair_remapping_experiments'||g" \
    -e 's|scripts/claude_pair_remapping_experiments/||g' \
    -e 's|scripts/claude_pair_remapping_experiments|model_experiment_queue|g' <<<"$src")
  if grep -q claude_pair_remapping_experiments <<<"$new"; then
    echo "$f: names the old folder in a way this does not rewrite; copy it by hand" >&2
    continue
  fi
  if [ -f "$dst" ] && [ "$(cat "$dst")" != "$new" ]; then
    echo "$f: replaced (differs from the copy here)"
  elif [ ! -f "$dst" ]; then
    echo "$f: added"
  fi
  printf '%s\n' "$new" > "$dst"
  n=$((n + 1))
done
echo "$n test files from $al at $ref ($(git -C "$al" rev-parse --short "$ref"))"
