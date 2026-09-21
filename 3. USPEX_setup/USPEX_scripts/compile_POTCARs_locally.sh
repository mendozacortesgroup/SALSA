#!/bin/bash

# Source of VASP PAW pseudopotentials. NOT distributed with SALSA - these are
# proprietary and separately licensed. Supply your own, via (in priority order):
#   1. the first argument to this script
#   2. $VASP_PP_PATH   (the variable ASE and pymatgen already use)
# Set it persistently with ../../setup_env.sh /path/to/potentials
global_POTCAR_Dir="${1:-${VASP_PP_PATH:-}}"

if [[ -z "$global_POTCAR_Dir" ]]; then
    echo "ERROR: no VASP pseudopotential directory given." >&2
    echo "  usage: $0 /path/to/vasp/potentials/PBE" >&2
    echo "  or:    export VASP_PP_PATH=/path/to/vasp/potentials/PBE" >&2
    echo "  or:    ../../setup_env.sh /path/to/vasp/potentials/PBE" >&2
    exit 1
fi
if [[ ! -d "$global_POTCAR_Dir" ]]; then
    echo "ERROR: '$global_POTCAR_Dir' is not a directory." >&2
    exit 1
fi

# Destination. Defaults to the in-repo POTCARs directory (which ships empty,
# with only a README) so the script works whether or not SALSA_DIR is exported.
_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
local_POTCAR_Dir="${SALSA_POTCAR_DIR:-$_script_dir/../POTCARs}"

elements=( H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og )

if [[ ! -d "$local_POTCAR_Dir" ]]; then
    mkdir -p "$local_POTCAR_Dir" || {
        echo "ERROR: could not create '$local_POTCAR_Dir'." >&2
        exit 1
    }
    echo "Made directory $local_POTCAR_Dir"
fi
# Collapse the "USPEX_scripts/../POTCARs" form so the paths we print back are
# the ones a user would actually type.
local_POTCAR_Dir="$(cd "$local_POTCAR_Dir" && pwd)"

copied=0
skipped=0
missing=0

for element in "${elements[@]}"; do
    potential_POTCAR="${global_POTCAR_Dir}/${element}/POTCAR"
    potential_new_POTCAR="${local_POTCAR_Dir}/POTCAR_${element}"
    if [[ -f "$potential_new_POTCAR" ]]; then
        echo "POTCAR_${element} already exists locally. Skipping."
        skipped=$(( skipped + 1 ))
    elif [[ -f "$potential_POTCAR" ]]; then
        # Report the copy only if it actually happened. An earlier version
        # printed "Copied" unconditionally, so a failing cp looked like success
        # for all 118 elements and the script still exited 0.
        if cp "$potential_POTCAR" "$potential_new_POTCAR"; then
            echo "Copied POTCAR_${element} to local directory"
            copied=$(( copied + 1 ))
        else
            echo "ERROR: failed to copy ${potential_POTCAR}" >&2
            exit 1
        fi
    else
        echo "Did not find ${potential_POTCAR}. Skipping."
        missing=$(( missing + 1 ))
    fi
done

echo
echo "Done: ${copied} copied, ${skipped} already present, ${missing} not found in"
echo "  ${global_POTCAR_Dir}"
echo "Destination: ${local_POTCAR_Dir}"
if (( copied == 0 && skipped == 0 )); then
    echo
    echo "WARNING: no pseudopotentials were found. This script expects one" >&2
    echo "  subdirectory per element, each containing a file named POTCAR:" >&2
    echo "    ${global_POTCAR_Dir}/Ag/POTCAR" >&2
    echo "    ${global_POTCAR_Dir}/Si/POTCAR" >&2
    echo "  which is the layout of the VASP potpaw_PBE distribution." >&2
    exit 1
fi

