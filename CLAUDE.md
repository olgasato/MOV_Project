# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A MATLAB workflow computing MOV (Meridional Overturning component of freshwater transport) and its full decomposition (mean + Mov + gyre) at 34.5°S, from the SAMBA array's inverted echo sounder (IES/PIES) mooring data at the western (Samba-W) and eastern (Samba-E) boundaries. Not a package — no build system, test suite, or dependency manifest; a sequence of MATLAB scripts run in order. Version-controlled with git, pushed to `https://github.com/olgasato/MOV_Project` (git identity: `olgasato <olga.sato@usp.br>`).

## Versioning convention

Scripts are iterated as `<name>_vN.m` rather than edited in place — e.g. `concat_IES.m` → `concat_IES_v4.m` → `concat_IES_v5.m` (current). Older versions are kept for history/comparison rather than deleted. When fixing or extending a script, save the result under the next version number rather than overwriting the file being changed.

## The two pipelines

Both start from the same Step 1 (`concat_IES_v5.m`), then diverge:

### Step 1: `concat_IES_v5.m` (current version)

Loads `samba_w.mat`/`samba_e.mat` (not in this repo — must be on the MATLAB path or in the working directory), restricts both to the common overlap period 2013-09-11 to 2017-07-17, and extracts per-mooring temperature/salinity as flat variables (`tem_A`, `tem_AA`, ..., `sal_P1`, `sal_P8`, ...) plus `prew`/`pree` (pressure), `depw`/`depe` (depth), `lonw`/`lone` (longitude, sorted), and `dt`. Saves everything to `concat_IESsamba.mat` (gitignored — regenerate by running the script, requires the raw `samba_w`/`samba_e` inputs).

Earlier versions (`concat_IES.m`, `concat_IES_v2.m`, `concat_IES_v4.m`) are kept for comparison — `v4` aggregated all sites into single `tem_w`/`sal_w` matrices, which was a bug (mixed different moorings along the time axis, incompatible with downstream per-site calls); `v5` reverted to per-site variables while keeping `v4`'s usability improvements (progress messages, `-v7.3` save).

### Step 2a: `mov_samba.m` — dynamic height only, not maintained further

Computes geostrophic dynamic height anomaly (`gpan_*`) per mooring via `gsw_geo_strf_dyn_height`, referenced to 1000 dbar. Requires the GSW (TEOS-10) toolbox. The velocity/transport part is commented out as a sketch, superseded by the `mov_samba_marion*.m` line below.

### Step 2b: `mov_samba_marion_v6.m` (current version) — full Mov decomposition

Combines the Step 1 output with Marion's independently-computed absolute velocities (external file, see Dependencies below) to compute freshwater transport at 9 sites (A, C, D, P1, P2, P4, P5, P6, P8 — a subset of the full mooring list), decomposed into three components:

- **`mean`** — net/barotropic transport
- **`mov`** — overturning component (the classic "Mov")
- **`gyre`** — local (per-mooring-pair) deviation from the zonal-mean profiles

A `total_direct` (independently computed, non-decomposed) is compared against `mean+mov+gyre` as a self-consistency check, printed via `fprintf` — should be ~0 (currently closes to floating-point noise, ~1e-16 Sv).

**Citations** (in the script header): Marion's 2021 paper for the velocity fields (https://doi.org/10.1029/2020JC016947); two papers defining the mean/Mov/gyre decomposition being implemented (https://doi.org/10.1016/S0074-6142(01)80134-0, https://doi.org/10.1029/2023JC020558). Citations only — article PDFs are deliberately not committed to this (public) repo to avoid redistributing copyrighted journal content.

**Version history** (`mov_samba_marion.m` → `_v6.m`): each version's header comment documents what changed and why — worth reading before touching this script, since several non-obvious bugs were found and fixed in sequence:
- `v2`: added `mean`/`gyre` alongside the pre-existing `mov`, plus the `total_direct` residual check.
- `v3`–`v5`: chased down why the residual didn't close to zero. Root causes, in order found: (1) `V0` (velocity zonal-mean) and `wsal` (salinity zonal-mean) used inconsistent weighting — masking `V0` the same way as `wsal` cut the residual ~85%; (2) confirmed `dx1`'s validity mask is time-invariant, ruling out the `dx2`-vs-`dx1` time-slice mismatch; (3) `V00`/`S0` were unweighted depth-means, inconsistent with `dx2`'s depth-varying weighting — redefining both as `dx2`-weighted depth-means closed the residual to ~1e-16 Sv. **Note: this changed `mov`'s numeric value** relative to the original `mov_samba_marion.m` (used for the OSM26 presentation) — not a pure refactor.
- `v6`: robustness fixes anticipating more station data — `assert(isequal(prew,pree))` (the depth axis `pre=prew` is silently assumed identical for both West and East sites), the `dx1` time-invariance check kept as a permanent `warning()` rather than a one-off, monthly-means array sizing derived from `dt` instead of hardcoded (`47` months, `init=9`/`iend=7`). **Not addressed**: the site list itself (9 sites, 8 gaps) is still hardcoded throughout — `sal00`'s dimension, the `for ii=1:8` gap loop, `cpiesW([1 5 6])`/`cpiesE([1 2 4 5 6 8])` index selection, `p_cpiesW`/`p_cpiesE` cutoff arrays. Adding/removing a mooring station requires editing all of those by hand.

## External dependencies (not in this repo)

- `~/research/sambar/renellys_sent/Marions_code/functions/positions_pies/position_West.m` and `position_East.m` — define the `cpiesW`/`cpiesE` structs (site name/lat/lon) used for inter-site distances. Added via `addpath(...)` at the top of `mov_samba_marion_v6.m`; not on the default MATLAB path (`~/matlab/startup.m` only adds `~/matlab/*`).
- `~/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT.mat` — Marion's absolute velocity/transport-per-unit-depth (`Total_TPUD1..8`) data, loaded via an explicit absolute path (`prefix` variable).
- Custom utility functions assumed on the path but not in this repo: `lowpass_filter`, `sinfitb_tot`, `unico`.
- GSW (TEOS-10) Gibbs SeaWater MATLAB toolbox (`gsw_geo_strf_dyn_height`, `gsw_distance`) — on the path via `~/matlab/startup.m`.

## Site/mooring naming

- Samba-W (western boundary) moorings: `A, AA, B, BB, C, CC, D`, with corresponding depths `depw = [1350 2902 3510 4173 4558 4730 4756]` m.
- Samba-E (eastern boundary) moorings: `P1, P2, P4, P5, P6, P8` (P3 and P7 are missing from the array), with corresponding depths `depe = [1266 2129 4482 4969 5185 4608]` m.
- Variable naming convention throughout: `<field>_<site>`, e.g. `tem_AA`, `sal_P4`.
- `mov_samba_marion*.m` uses a 9-site subset (A, C, D, P1, P2, P4, P5, P6, P8) — the ones Marion's velocity dataset covers.

## Outputs

- `concat_IESsamba.mat`, `gpan_samba.mat` — gitignored (`.mat` files excluded generally; regenerate by running the corresponding script).
- `mov_IES.png` — plot output from `mov_samba_marion*.m`, committed (not gitignored).
