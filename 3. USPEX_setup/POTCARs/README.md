# POTCARs — user-supplied

VASP PAW pseudopotentials are proprietary and are **not** distributed with
SALSA. You must supply your own licensed copies.

Populate this directory with `compile_POTCARs_locally.sh`, which copies from
your VASP installation:

```bash
export VASP_PP_PATH=/path/to/vasp/potentials/PBE   # or pass it as an argument
../USPEX_scripts/compile_POTCARs_locally.sh
```

Expected result: one file per element, named `POTCAR_<Symbol>` (`POTCAR_Ag`,
`POTCAR_Br`, ...). The source tree is expected in VASP's own layout, i.e.
`$VASP_PP_PATH/<Symbol>/POTCAR`.

**These are needed only for stage 3 (USPEX structure prediction).** Stages 1, 2
and 4, and the convex-hull analysis in v2.0.0, do not require VASP at all.
