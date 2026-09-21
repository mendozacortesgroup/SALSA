#!/bin/bash
# Configure the shell environment SALSA needs.
#
#   ./setup_env.sh                      # set SALSA_DIR only
#   ./setup_env.sh /path/to/vasp/pot    # also set VASP_PP_PATH (stage 3 only)
#   ./setup_env.sh --print              # show the block, change nothing
#
# Appends a marked block to ~/.bashrc. Re-running replaces that block rather
# than adding a duplicate, so it is safe to run as often as you like. Set
# RC_FILE to write somewhere other than ~/.bashrc.

set -euo pipefail

SALSA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RC="${RC_FILE:-$HOME/.bashrc}"
VASP_POT=""
PRINT_ONLY=0

for arg in "$@"; do
    case "$arg" in
        --print) PRINT_ONLY=1 ;;
        -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
        -*)
            echo "ERROR: unknown option '$arg'." >&2
            echo "  usage: $0 [--print] [/path/to/vasp/potentials]" >&2
            exit 1
            ;;
        *)
            if [[ -n "$VASP_POT" ]]; then
                echo "ERROR: more than one pseudopotential directory given" >&2
                echo "       ('$VASP_POT' and '$arg')." >&2
                exit 1
            fi
            VASP_POT="$arg"
            ;;
    esac
done

if [[ -n "$VASP_POT" ]]; then
    if [[ ! -d "$VASP_POT" ]]; then
        echo "ERROR: '$VASP_POT' is not a directory." >&2
        exit 1
    fi
    # Store an absolute path: the value is read from ~/.bashrc in shells whose
    # working directory has nothing to do with where this script was run.
    VASP_POT="$(cd "$VASP_POT" && pwd)"
    # A VASP pseudopotential tree holds per-element subdirectories each with a POTCAR.
    if ! compgen -G "$VASP_POT/*/POTCAR" > /dev/null; then
        echo "WARNING: no */POTCAR found under '$VASP_POT'." >&2
        echo "         Expected a layout like <dir>/Ag/POTCAR, <dir>/Br/POTCAR, ..." >&2
    fi
fi

BLOCK="# >>> SALSA environment >>>
export SALSA_DIR=\"$SALSA_ROOT\""
if [[ -n "$VASP_POT" ]]; then
    BLOCK="$BLOCK
export VASP_PP_PATH=\"$VASP_POT\""
fi
BLOCK="$BLOCK
# <<< SALSA environment <<<"

if [[ $PRINT_ONLY -eq 1 ]]; then
    echo "$BLOCK"
    exit 0
fi

# Remove any previous SALSA block, then append the current one.
if grep -q '^# >>> SALSA environment >>>' "$RC" 2>/dev/null; then
    sed -i '/^# >>> SALSA environment >>>$/,/^# <<< SALSA environment <<<$/d' "$RC"
    echo "Replaced the existing SALSA block in $RC"
else
    echo "Adding a SALSA block to $RC"
fi

# Deleting the block leaves behind the blank line that separated it, so without
# this every run would add one more. Command substitution strips all trailing
# newlines; writing through cat keeps the file's inode and permissions.
if [[ -s "$RC" ]]; then
    printf '%s\n' "$(cat "$RC")" > "$RC.salsa.tmp"
    cat "$RC.salsa.tmp" > "$RC"
    rm -f "$RC.salsa.tmp"
fi
printf '\n%s\n' "$BLOCK" >> "$RC"

echo
echo "  SALSA_DIR     = $SALSA_ROOT"
[[ -n "$VASP_POT" ]] && echo "  VASP_PP_PATH  = $VASP_POT" \
                     || echo "  VASP_PP_PATH  = (unset - required only for stage 3)"
echo
echo "Run 'source $RC' or open a new shell to apply."
