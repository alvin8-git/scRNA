#!/usr/bin/env bash
# PreToolUse hook (matcher: Bash).
#
# Policy: auto-allow bash commands, but force a permission prompt ("ask") when
# the command may DELETE files. Reads the hook JSON on stdin and emits a
# permissionDecision on stdout.
#
# Fails SAFE: on any error it emits nothing and exits 0, so Claude Code falls
# back to its normal permission prompt rather than silently auto-allowing.
set -uo pipefail

input=$(cat)
cmd=$(printf '%s' "$input" | jq -r '.tool_input.command // ""' 2>/dev/null) || exit 0

emit() { # $1 = decision (allow|ask), $2 = reason
  jq -nc --arg d "$1" --arg r "$2" \
    '{hookSpecificOutput:{hookEventName:"PreToolUse",permissionDecision:$d,permissionDecisionReason:$r}}'
  exit 0
}

# 1) rm / rmdir / unlink / shred in COMMAND position only: at the start, after a
#    shell separator ( ; & | ( && || ), or after a command prefix word
#    (sudo/xargs/time/exec/nohup/env/then/do/else/git ...), with optional /path/.
#    This deliberately does NOT fire on rm as a mere argument (e.g. `grep -r rm`).
if printf '%s' "$cmd" | grep -Eq '(^|[;&|(]|&&|\|\||\b(sudo|xargs|time|command|exec|nohup|env|then|do|else|git)[[:space:]]+)[[:space:]]*(/[^[:space:]]*/)?(rm|rmdir|unlink|shred)([[:space:];&|)]|$)'; then
  emit ask "Command may delete files (rm/rmdir/unlink/shred) — confirm before running."
fi

# 2) find ... -delete   or   find ... -exec rm
if printf '%s' "$cmd" | grep -Eq '(^|[[:space:]])find[[:space:]].*(-delete|-exec[[:space:]]+(/[^[:space:]]*/)?rm)'; then
  emit ask "Command deletes via find (-delete / -exec rm) — confirm before running."
fi

# 3) git clean with a force flag (deletes untracked files)
if printf '%s' "$cmd" | grep -Eq '(^|[[:space:]])git[[:space:]]+clean[[:space:]][^;&|]*-[A-Za-z]*f'; then
  emit ask "git clean -f can delete untracked files — confirm before running."
fi

# 4) xargs ... rm (flags may sit between xargs and rm, e.g. `xargs -0 rm`)
if printf '%s' "$cmd" | grep -Eq '\bxargs\b[^;&|]*[[:space:]](/[^[:space:]]*/)?rm([[:space:];&|)]|$)'; then
  emit ask "Command pipes into rm via xargs — confirm before running."
fi

emit allow "Auto-approved: no file-deletion command detected."
