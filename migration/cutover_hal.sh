#!/bin/bash
# The cutover of the KC-remapping queue from al_analysis's folder scripts/claude_pair_remapping_experiments
# to this repository, on hal, at the pause (MIGRATION.md has the whole checklist and the order).
# Each phase checks its preconditions and stops with a message rather than guess; run them in order,
# reading each one's output before the next:
#
#   merge [<al_analysis ref>]  (a) splits the folder's history again from al_analysis's main (default
#                              origin/main, fetched first) in a scratch clone, and merges what changed
#                              in the folder since this repository's history was split (87c871f) into
#                              main here; if main was rewritten meanwhile (from_remy/ and
#                              apl_pair_analysis/ removed), the new commits are first replayed onto the
#                              rewritten counterpart of 87c871f, without those two folders, and that
#                              is merged. Then copies tests added to or changed in al_analysis's test/
#                              since. Commits locally; pushing is left to the operator.
#   host                       (c) hal's own setup: moves outputs/, from_yang/ and .fable_consent.local
#                              from the old folder, makes ./venv with al_analysis, hong2p and olfsysm
#                              editable (olfsysm only when no process has it loaded), relinks the
#                              queue monitor's systemd unit to this repository and restarts it, and
#                              points the 12 hook commands in ~/src/dotfiles/claude_settings.json here
#                              (left uncommitted, for Tom).
#   olfsysm                    the olfsysm install alone, for when `host` skipped it.
#   al-analysis                (b) prints the commands for al_analysis's removal commit (run by hal's
#                              session, after `host`).
#   check                      what is in place and what is not.
#
# Usage: migration/cutover_hal.sh <phase> [args]
# CUTOVER_IGNORE_PIDS='<pid> ...' leaves those processes out of the check for processes in the old
# folder (a person's idle shell whose working directory is there, or the queue monitor, which `host`
# restarts from here itself); everything else there still stops the phase.
set -uo pipefail

here=$(cd "$(dirname "$0")/.." && pwd -P)
src=$(cd "$here/.." && pwd -P)
al=$src/al_analysis
old=$al/scripts/claude_pair_remapping_experiments
prefix=scripts/claude_pair_remapping_experiments
split_base=87c871f  # this repository's main when it was made (split of al_analysis 16f552e)
venv=$here/venv
unit=claude-queue-monitor
dotfiles_settings=$HOME/src/dotfiles/claude_settings.json

die() { echo "cutover_hal.sh: $*" >&2; exit 1; }
say() { printf '\n== %s\n' "$*"; }

ancestors() {  # pids of this shell and its ancestors
  local p=$$
  while [ -n "$p" ] && [ "$p" != 0 ]; do
    echo "$p"
    p=$(awk '/^PPid:/{print $2}' "/proc/$p/status" 2>/dev/null)
  done
}

procs_in() {  # $1 folder: pids (not ours) whose cwd is in it or whose command line names it
  local mine p pid c
  mine=" $(ancestors | paste -sd' ') ${CUTOVER_IGNORE_PIDS//,/ } "
  for p in /proc/[0-9]*; do
    pid=${p#/proc/}
    [[ $mine == *" $pid "* ]] && continue
    c=$(readlink "$p/cwd" 2>/dev/null) || continue
    case $c/ in "$1"/*) echo "$pid"; continue ;; esac
    tr '\0' ' ' < "$p/cmdline" 2>/dev/null | grep -qF "$1/" && echo "$pid"
  done | sort -un
}

show_procs() {  # pids on stdin: one line each
  local pid
  while read -r pid; do
    printf '    %s  %s\n' "$pid" "$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null | cut -c1-150)"
  done
}

olfsysm_users() {  # pids that have an olfsysm extension built in the olfsysm checkout loaded
  grep -l "$src/olfsysm/olfsysm" /proc/[0-9]*/maps 2>/dev/null | sed 's|/proc/\([0-9]*\)/maps|\1|' | sort -un
}

clean_here() {
  [ -z "$(git -C "$here" status --porcelain --untracked-files=no)" ] ||
    die "$here has uncommitted changes to tracked files; commit or stash them first"
}

phase_merge() {
  local ref=${1:-origin/main} tmp split new n r
  clean_here
  [ "$(git -C "$here" symbolic-ref -q --short HEAD)" = main ] || die "$here is not on main"
  git -C "$here" fetch -q origin || die "fetching this repository's origin failed"
  [ "$(git -C "$here" rev-parse HEAD)" = "$(git -C "$here" rev-parse origin/main)" ] ||
    die "main here is not origin/main; pull (or push) first so the merge starts from what is on GitHub"
  if [ "$ref" = origin/main ]; then
    git -C "$al" fetch -q https://github.com/ejhonglab/al_analysis.git main:refs/remotes/origin/main ||
      die "fetching al_analysis failed"
  fi
  r=$(git -C "$al" rev-parse --verify -q "$ref^{commit}") || die "no $ref in $al"
  [ -z "$(git -C "$al" status --porcelain --untracked-files=no -- "$prefix")" ] ||
    echo "warning: $old has uncommitted changes; they are not in $ref and will not move"
  [ -z "$(git -C "$al" log --oneline "$r..HEAD" -- "$prefix" 2>/dev/null)" ] ||
    echo "warning: al_analysis's HEAD has commits to the folder that are not in $ref (unpushed?)"

  say "splitting $prefix at al_analysis $ref ($(git -C "$al" rev-parse --short "$r"))"
  tmp=$(mktemp -d /tmp/mq_cutover.XXXXXX) || die "mktemp failed"
  git clone -q --shared --no-checkout "$al" "$tmp/al" || die "cloning $al failed"
  git -C "$tmp/al" sparse-checkout set "$prefix" &&
    git -C "$tmp/al" checkout -q --detach "$r" || die "checking out $r in the scratch clone failed"
  split=$(cd "$tmp/al" && git subtree split --prefix="$prefix" "$r" 2>/dev/null | tail -1)
  [[ $split =~ ^[0-9a-f]{40}$ ]] || die "git subtree split failed in $tmp/al"
  git -C "$tmp/al" update-ref refs/migration/split "$split"
  git -C "$here" fetch -q "$tmp/al" "+refs/migration/split:refs/migration/split" ||
    die "fetching the split failed"
  git -C "$here" merge-base --is-ancestor "$split_base" "$split" || die "the split ($split) does not" \
    "contain $split_base: it was not made the same way as this repository's history; stop and look"
  n=$(git -C "$here" rev-list --count "$split_base..$split")
  echo "split $split: $n commits to the folder since $split_base"
  git -C "$here" log --oneline --reverse "$split_base..$split" | sed 's/^/    /'

  if [ "$n" = 0 ]; then
    echo "nothing to merge"
  elif git -C "$here" merge-base --is-ancestor "$split_base" HEAD; then
    say "merging (main still contains $split_base)"
    git -C "$here" merge --no-ff --no-edit -m "Merge the folder's last commits in al_analysis (to $(git -C "$al" rev-parse --short "$r"))" \
      refs/migration/split || die "the merge stopped (conflicts above): resolve them, commit, and go on"
  else
    # main was rewritten (from_remy/ and apl_pair_analysis/ removed): find 87c871f's counterpart by
    # its author, author date and subject, replay the new commits onto it without those folders,
    # and merge that, so the old history (and the two folders) do not come back
    local key base wt c
    key=$(git -C "$here" log -1 --format='%ae %at %s' "$split_base")
    base=$(git -C "$here" log --format='%H %ae %at %s' HEAD | awk -v k="$key" 'substr($0, 42) == k {print $1}')
    [ "$(wc -w <<<"$base")" = 1 ] || die "main does not contain $split_base and its rewritten" \
      "counterpart is not unique ($(wc -w <<<"$base") matches of '$key'); replay by hand (MIGRATION.md)"
    say "main was rewritten: replaying the $n commits onto $(git -C "$here" rev-parse --short "$base") ($split_base rewritten)"
    wt=$tmp/replay
    git -C "$here" worktree add -q --detach "$wt" "$base" || die "git worktree add failed"
    for c in $(git -C "$here" rev-list --reverse "$split_base..$split"); do
      if ! git -C "$wt" cherry-pick --no-commit "$c" >/dev/null 2>&1; then
        # the only expected conflicts: changes to the two removed folders
        git -C "$wt" rm -r -q --cached --ignore-unmatch from_remy apl_pair_analysis >/dev/null
        if [ -n "$(git -C "$wt" diff --name-only --diff-filter=U)" ]; then
          die "replaying $(git -C "$here" log -1 --format='%h %s' "$c") conflicts outside from_remy/ and" \
            "apl_pair_analysis/ in $wt; resolve there, git -C $wt commit -C $c, and continue by hand (MIGRATION.md)"
        fi
      fi
      git -C "$wt" rm -r -q --cached --ignore-unmatch from_remy apl_pair_analysis >/dev/null
      git -C "$wt" clean -q -fd -- from_remy apl_pair_analysis
      if git -C "$wt" diff --cached --quiet; then
        echo "    $(git -C "$here" log -1 --format='%h %s' "$c" | cut -c1-90): nothing left, skipped"
        git -C "$wt" reset -q --hard
        continue
      fi
      git -C "$wt" commit -q -C "$c" || die "committing the replay of $c failed in $wt"
    done
    new=$(git -C "$wt" rev-parse HEAD)
    git -C "$here" update-ref refs/migration/replayed "$new"
    git -C "$here" worktree remove --force "$wt"
    say "merging the replayed commits"
    git -C "$here" merge --no-ff --no-edit -m "Merge the folder's last commits in al_analysis (to $(git -C "$al" rev-parse --short "$r"), replayed onto the rewritten history)" \
      refs/migration/replayed || die "the merge stopped (conflicts above): resolve them, commit, and go on"
    [ -z "$(git -C "$here" ls-tree --name-only HEAD -- from_remy apl_pair_analysis)" ] ||
      die "from_remy/ or apl_pair_analysis/ is back in the tree; remove it before pushing"
  fi

  say "tests added to or changed in al_analysis since"
  "$here/migration/copy_al_analysis_tests.sh" "$al" "$r" || die "copying the tests failed"
  if [ -n "$(git -C "$here" status --porcelain -- test)" ]; then
    git -C "$here" add test &&
      git -C "$here" commit -q -m "test/: the queue's tests as in al_analysis $(git -C "$al" rev-parse --short "$r")

migration/copy_al_analysis_tests.sh run again at the cutover: tests added to or changed in
al_analysis/test since the first copy" || die "committing the tests failed"
    git -C "$here" log -1 --stat --format='%h %s' | sed 's/^/    /'
  else
    echo "no change"
  fi
  rm -rf "${tmp:?}"
  cat <<EOF

Done; nothing is pushed. Next: compile and test from $here (QUICK=1 ./venv/bin/python -m pytest, once
the venv exists, or with al_analysis's venv), then push: git push origin HEAD:main
EOF
}

move_host_files() {
  local f e
  [ "$(stat -c %d "$old")" = "$(stat -c %d "$here")" ] || die "$old and $here are on different" \
    "filesystems; move the files by hand (rsync -a, then delete the originals)"
  for f in outputs from_yang .fable_consent.local; do
    [ -e "$old/$f" ] || continue
    if [ ! -e "$here/$f" ]; then
      mv "$old/$f" "$here/$f" || die "could not move $old/$f"
      echo "moved $f"
    elif [ -d "$old/$f" ] && [ -d "$here/$f" ]; then
      for e in "$old/$f"/* "$old/$f"/.[!.]*; do
        [ -e "$e" ] || continue
        if [ -e "$here/$f/${e##*/}" ]; then
          echo "kept $e: $here/$f/${e##*/} exists already (merge them by hand)"
        else
          mv "$e" "$here/$f/" || die "could not move $e"
        fi
      done
      rmdir "${old:?}/${f:?}" 2>/dev/null && echo "moved $f"
    else
      echo "kept $old/$f: $here/$f exists already (merge them by hand)"
    fi
  done
}

install_olfsysm() {
  local users
  users=$(olfsysm_users)
  if [ -n "$users" ]; then
    echo "olfsysm is loaded by these processes, so it is not rebuilt now (pip install -e rebuilds"
    echo "the extension in $src/olfsysm in place); run 'migration/cutover_hal.sh olfsysm' when they end:"
    show_procs <<<"$users"
    return 1
  fi
  "$venv/bin/python" -m pip install -q -e "$src/olfsysm" || die "installing olfsysm failed"
  echo "olfsysm installed editable from $src/olfsysm"
}

check_venv() {
  (cd / && "$venv/bin/python" -I - "$src" <<'EOF')
import os, sys
src = os.path.realpath(sys.argv[1])
import al_analysis, hong2p, olfsysm
from al_analysis import mb_model
for name, mod in (('al_analysis', al_analysis), ('hong2p', hong2p), ('olfsysm', olfsysm)):
    loc = getattr(mod, '__file__', None) or next(p for p in mod.__path__ if os.path.isdir(p))
    print(f'{name}: {os.path.realpath(loc)}')
p = mb_model.olfsysm_clean_problems()
print('olfsysm: ' + ('clean build of its checkout' if not p else 'NOT clean: ' + '; '.join(p)))
EOF
}

phase_host() {
  local busy
  [ -d "$old" ] || echo "($old is gone already)"
  say "processes in or from the old folder"
  busy=$(procs_in "$old")
  if [ -n "$busy" ]; then
    echo "These still run in or from $old (model runs, chat watchers, sessions):"
    show_procs <<<"$busy"
    die "let the runs finish and stop the watchers (TaskStop) and sessions started there, then re-run"
  fi
  echo "none"

  if [ -d "$old" ]; then
    say "moving outputs/, from_yang/ and .fable_consent.local"
    move_host_files
  fi

  say "venv ($venv)"
  if [ ! -e "$venv" ]; then
    command -v python3.14 >/dev/null || die "no python3.14 on PATH"
    python3.14 -m venv "$venv" || die "making the venv failed"
    "$venv/bin/python" -m pip install -q -U pip || die "upgrading pip failed"
    # al_analysis pins hong2p and olfsysm to their GitHub versions; the editable installs replace them
    "$venv/bin/python" -m pip install -q -e "$al[dev]" || die "installing al_analysis failed"
    "$venv/bin/python" -m pip install -q -e "$src/hong2p[all]" || die "installing hong2p failed"
  else
    echo "exists; not changed"
  fi
  if "$venv/bin/pip" show olfsysm 2>/dev/null | grep -q "^Editable project location: $src/olfsysm$"; then
    echo "olfsysm already editable from $src/olfsysm"
  else
    install_olfsysm || true
  fi
  check_venv || echo "the venv does not import everything yet (above)"

  say "queue monitor ($unit)"
  [ -f "$here/monitor-server/$unit.service" ] || die "no $here/monitor-server/$unit.service"
  systemctl --user stop "$unit" 2>/dev/null
  systemctl --user disable "$unit" 2>/dev/null  # also removes the link to the old folder's unit file
  systemctl --user link "$here/monitor-server/$unit.service" || die "systemctl --user link failed"
  systemctl --user daemon-reload
  systemctl --user enable --now "$unit" || die "systemctl --user enable --now $unit failed"
  sleep 3
  systemctl --user is-active "$unit" && systemctl --user show -p FragmentPath "$unit"

  say "hook commands in $dotfiles_settings"
  if [ -f "$dotfiles_settings" ]; then
    local n
    n=$(grep -c "\$HOME/src/al_analysis/$prefix/" "$dotfiles_settings")
    if [ "$n" -gt 0 ]; then
      # one exact substitution, the same in all 12 commands; the file is a symlink target
      # (~/.claude/settings.json), so it is rewritten in place and running sessions reload it
      python3 -I - "$dotfiles_settings" "\$HOME/src/al_analysis/$prefix/" '$HOME/src/model_experiment_queue/' <<'EOF' || die "rewriting the settings failed"
import json, sys
from pathlib import Path
p, old, new = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
s = p.read_text().replace(old, new)
json.loads(s)
p.write_text(s)
EOF
      echo "$n commands now run the copies in ~/src/model_experiment_queue (uncommitted in ~/src/dotfiles;"
      echo "dotfiles is Tom's repository, not one of the queue's four: ask him to commit and push it):"
      git -C "$(dirname "$dotfiles_settings")" diff --stat
    else
      echo "nothing names the old folder"
    fi
  else
    echo "no $dotfiles_settings"
  fi
  cat <<EOF

Done. Next: start the queue sessions in $here (accept the workspace-trust prompt), restart the chat
watcher from there, then al_analysis's removal commit (migration/cutover_hal.sh al-analysis).
EOF
}

phase_al_analysis() {
  local tests
  tests=$(git -C "$al" grep -l "$prefix" HEAD -- 'test/test_*.py' | sed 's/^HEAD://')
  cat <<EOF
In $al, as hal's session, once 'host' is done and no session runs in the old folder:

    cd $al
    git pull --rebase
    git rm -r -q $prefix
    git rm -q \$(git grep -l $prefix -- 'test/test_*.py')

(at al_analysis's HEAD now, $(git -C "$al" rev-parse --short HEAD), the last line removes
$(echo $tests | wc -w) tests: $(echo $tests | sed 's|test/||g'); after the pull it takes the ones there
then, which should be the ones in this repository's test/)

then edit, by hand:
- README.md: "The experiment harness (\`$prefix\`)" becomes "The KC-remapping experiment
  harness (ejhonglab/model_experiment_queue since 2026-11-11)";
- al_analysis/mb_model.py, olfsysm_checkout_path's docstring: "(see
  $prefix/join_queue.sh)" becomes "(see join_queue.sh in
  ejhonglab/model_experiment_queue)";
- .gitignore: the comment "(the experiments folder's .claude/settings.json is shared)" can become
  "(a .claude/settings.json in a checkout is shared)"; the two patterns stay;
then run QUICK=1 ./venv/bin/python -m pytest test/test_model_provenance.py, commit (one line per
file, and one for the tests moved to model_experiment_queue's test/) and push. The old folder keeps
only ignored files afterwards (caches); delete it when nothing runs there.
EOF
}

phase_check() {
  say "this repository"
  git -C "$here" log -1 --format='%h %s' | cut -c1-100
  git -C "$here" status --short | head
  say "venv"
  [ -x "$venv/bin/python" ] && check_venv || echo "no venv at $venv"
  say "old folder"
  if [ -d "$old" ]; then
    ls -A "$old" | head
    procs_in "$old" | show_procs
  else
    echo "gone"
  fi
  say "monitor"
  systemctl --user is-active "$unit"; systemctl --user show -p FragmentPath "$unit"
  say "hooks"
  grep -c "src/model_experiment_queue/" "$dotfiles_settings" 2>/dev/null | sed 's/$/ dotfiles commands name this repository/'
  grep -c "$prefix" "$dotfiles_settings" 2>/dev/null | sed 's/$/ name the old folder/'
}

case ${1:-} in
  merge) shift; phase_merge "$@" ;;
  host) phase_host ;;
  olfsysm) [ -x "$venv/bin/python" ] || die "no venv yet (run host)"; install_olfsysm && check_venv ;;
  al-analysis) phase_al_analysis ;;
  check) phase_check ;;
  *) sed -n '2,/^set -uo/{/^set -uo/d;s/^# \{0,1\}//;p}' "$0"; exit 1 ;;
esac
