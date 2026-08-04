# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A small MATLAB workflow to compute MOV (Meridional Overturning Volume/transport) from the SAMBA array's inverted echo sounder (IES/PIES) mooring data at the western (Samba-W) and eastern (Samba-E) boundaries. Not a package — no build system, test suite, or dependency manifest. Just two scripts run in order inside MATLAB.

## Running the workflow

Run in MATLAB, in this order, from this directory:

1. `concat_IES.m` — loads `samba_w` and `samba_e` mat-files (not present in this repo; must exist on the MATLAB path or in the working directory), restricts both to the common overlap period 2013-09-11 to 2017-07-17, and extracts per-mooring temperature/salinity arrays plus pressure/depth/longitude vectors for each site. Ends by clearing the raw loaded variables, leaving only the site-level `tem_*`, `sal_*`, `prew`/`pree`, `depw`/`depe`, `lonw`/`lone`, and `dt` in the workspace. The `save concat_IESsamba` line at the end is commented out — if step 2 is run in a fresh session rather than immediately after step 1, uncomment that save (and the corresponding `load` at the top of `mov_samba.m`) to persist/restore the workspace between steps.
2. `mov_samba.m` — takes the workspace from step 1 and computes geostrophic dynamic height anomaly (`gpan_*`) per mooring via `gsw_geo_strf_dyn_height`, referenced to 1000 dbar. Saves the result to `gpan_samba.mat`. Requires the **GSW (TEOS-10) Gibbs SeaWater MATLAB toolbox** on the MATLAB path for `gsw_geo_strf_dyn_height`. The remainder of the script (geostrophic velocity/transport calculation and plotting, using `gsw_geostrophic_velocity` and `gsw_distance`) is commented out as a sketch of the next step, not working code.

## Site/mooring naming

- Samba-W (western boundary) moorings: `A, AA, B, BB, C, CC, D`, with corresponding depths `depw = [1350 2902 3510 4173 4558 4730 4756]` m.
- Samba-E (eastern boundary) moorings: `P1, P2, P4, P5, P6, P8` (P3 and P7 are missing from the array), with corresponding depths `depe = [1266 2129 4482 4969 5185 4608]` m.
- Variable naming convention throughout: `<field>_<site>`, e.g. `tem_AA`, `sal_P4`.

## Notes for future changes

- There is no version control in this directory yet — treat any edits as the first tracked state unless the user says otherwise.
- The two scripts communicate purely through the shared MATLAB workspace (or via the commented-out `save`/`load` of `concat_IESsamba`); there's no function interface. If adding a step 3, follow the same pattern (load `gpan_samba.mat`, operate on `gpan_*`/`lonw`/`lone`/`dt`) unless refactoring the whole workflow into functions.
