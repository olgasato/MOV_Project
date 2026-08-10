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

### Step 2b: `mov_samba_marion_v7.m` (current version) — full Mov decomposition

Combines the Step 1 output with Marion's independently-computed absolute velocities (external file, see Dependencies below) to compute freshwater transport at 9 sites (A, C, D, P1, P2, P4, P5, P6, P8 — a subset of the full mooring list), decomposed into three components:

- **`mean`** — net/barotropic transport
- **`mov`** — overturning component (the classic "Mov")
- **`gyre`** — local (per-mooring-pair) deviation from the zonal-mean profiles

A `total_direct` (independently computed, non-decomposed) is compared against `mean+mov+gyre` as a self-consistency check, printed via `fprintf` — should be ~0 (currently closes to floating-point noise, ~1e-16 Sv).

**Citations** (in the script header): Marion's 2021 paper for the velocity fields (https://doi.org/10.1029/2020JC016947); two papers defining the mean/Mov/gyre decomposition being implemented (https://doi.org/10.1016/S0074-6142(01)80134-0, https://doi.org/10.1029/2023JC020558). Citations only — article PDFs are deliberately not committed to this (public) repo to avoid redistributing copyrighted journal content.

**Version history** (`mov_samba_marion.m` → `_v7.m`): each version's header comment documents what changed and why — worth reading before touching this script, since several non-obvious bugs were found and fixed in sequence:
- `v2`: added `mean`/`gyre` alongside the pre-existing `mov`, plus the `total_direct` residual check.
- `v3`–`v5`: chased down why the residual didn't close to zero. Root causes, in order found: (1) `V0` (velocity zonal-mean) and `wsal` (salinity zonal-mean) used inconsistent weighting — masking `V0` the same way as `wsal` cut the residual ~85%; (2) confirmed `dx1`'s validity mask is time-invariant, ruling out the `dx2`-vs-`dx1` time-slice mismatch; (3) `V00`/`S0` were unweighted depth-means, inconsistent with `dx2`'s depth-varying weighting — redefining both as `dx2`-weighted depth-means closed the residual to ~1e-16 Sv. **Note: this changed `mov`'s numeric value** relative to the original `mov_samba_marion.m` (used for the OSM26 presentation) — not a pure refactor.
- `v6`: robustness fixes anticipating more station data — `assert(isequal(prew,pree))` (the depth axis `pre=prew` is silently assumed identical for both West and East sites), the `dx1` time-invariance check kept as a permanent `warning()` rather than a one-off, monthly-means array sizing derived from `dt` instead of hardcoded (`47` months, `init=9`/`iend=7`). **Not addressed**: the site list itself (9 sites, 8 gaps) is still hardcoded throughout — `sal00`'s dimension, the `for ii=1:8` gap loop, `cpiesW([1 5 6])`/`cpiesE([1 2 4 5 6 8])` index selection, `p_cpiesW`/`p_cpiesE` cutoff arrays. Adding/removing a mooring station requires editing all of those by hand.
- `v7`: `v6`'s new assert fired for real — `prew` (501 points, 0-5000 dbar) and `pree` (531 points, 0-5300 dbar) had different *lengths*, not just values. Root cause: `ind=pree<=5000` was already used to trim `sal_P*`/`tem_P*` down to 501 points, but `pree` itself was never reassigned with that same `ind`. Fixed with one line (`pree=pree(ind);`). Not a real grid mismatch — both grids share spacing (10 dbar) and origin (0).

## Upstream data pipeline (outside this repo)

`mov_samba_marion*.m`'s inputs (`Total_TPUD1..8`, and transitively `samba_w.mat`/`samba_e.mat`) come from a much larger MATLAB pipeline at `~/research/sambar/renellys_sent/Marions_code/`, documented in `README_MOC` and `README_MHT` there (not in this repo). Traced through steps A–E of `README_MOC` on 2026-08-07; summary below so this doesn't need to be re-derived from scratch.

**Chain:** raw PIES travel-time calibration (`README_CALIBRATE_IES_EAST`) → GEM fields (`README_GEM`) → step A `IES_Make_Profiles_*.m` (combines travel time + GEM → T/S profiles, this is what `samba_w.mat`/`samba_e.mat` ultimately derive from) → step B `Ekman_transports.m` (CCMP wind → Ekman transport) → step C `ECCO_TPUD_ShelfBits.m` (static, time-mean shelf-transport correction from ECCO, ~-4.4 Sv west / ~0 Sv east) → step D `correct_topo_all_MSM60.m` (static, per-gap continuous topography/open-area fraction `AREA_TOPO_1..8`, from real MSM60 cruise bathymetry) → step E `Full_Depth_OverturningEstimate_Cst_for_MHT.m` (combines all of the above into `Total_TPUD1..8`, saved to `Full_Depth_OverturningEstimate_Cst_for_MHT.mat`).

**Confirmed: `Total_TPUD1..8` already includes the topography correction, absolute reference, and Ekman.** Step E's script builds `Absolute_TPUD1..8` (baroclinic shear from dynamic height + absolute reference from real PIES bottom-pressure sensors + an ECCO mean-velocity offset at 1500 dbar), multiplies each by its per-gap `AREA_TOPO_i` profile, then adds per-gap Ekman (`Ekman_TPUD1..8`) — that sum *is* `Total_TPUD1..8`. So `mov_samba_marion*.m` does not need to (and currently does not) reapply any topography correction on the velocity side — it's already baked in.

**Open question flagged, not yet resolved:** `mov_samba_marion*.m`'s own salinity masking (`p_cpiesW`/`p_cpiesE`, a binary NaN cutoff at a fixed depth per *site*) is inconsistent with the continuous, per-*gap* `AREA_TOPO_i` weighting already applied to the velocity it's multiplied against in `mov`/`mean`/`gyre`. E.g. the D–P8 gap's `AREA_TOPO_3` starts dropping below 1.0 at 940 dbar (a mid-gap ridge), far shallower than either site's own individual cutoff (D: 4850, P8: 1280) — so the velocity smoothly attenuates there while the salinity is either fully valid or fully NaN depending on the site-level cutoff. Worth investigating whether this creates a real inconsistency in the decomposition, or whether it's negligible in practice.

**Date range provenance:** the `2013-09-11`–`2017-07-17` window used throughout this repo isn't an independent choice — it's inherited directly from a hardcoded trim in `Full_Depth_OverturningEstimate.m` (`bad=find(dt_SAM<datenum(2013,9,11) | dt_SAM>datenum(2017,7,17))`), matching the East PIES array's (Samba-E / "GH") availability window at the time that pipeline was last run.

**Data currency (relevant to "add more station data") — UPDATE 2026-08-10, see below for what changed.** Originally (2026-08-07): both the East PIES calibration (`Merge_SiteP*_records.m`, hardcoded to only merge the `2013_2015`/`2015_2017` legs) and the CCMP wind data (`Reformat_WindData.m`, hardcoded to only merge `LoadWinds_2009.mat`...`LoadWinds_2017.mat`) stopped at 2017. West PIES calibration is handled externally ("Chris" provides pre-calibrated `IES_SAM.mat`) and extends to 2018-04-30; `samba_w.mat` is a byte-identical copy of `.../IES profile data/IES_Make_Profiles_plusBrazil_plusDynHgt_Olga.mat`.

**Resolved (2026-08-10), in the sibling repo `~/research/sambar/renellys_sent/Marions_code/` (now version-controlled, `https://github.com/olgasato/SAMBA_E_IES`, own `CLAUDE.md` there with full detail):**
- **East array**: found `~/research/sambar/samba_e/Daily_Tau-*.zip`, an officially published, already-calibrated dataset (SAMOC-SA/DFFE) covering **2013-09 to 2023-09** for P1/P2/P3/P3a/P4/P4a/P5/P6/P8 — cross-validated against a from-scratch recalibration of 2017-2019 (agreed to 0.0024s/0.16 dbar) and adopted as the primary calibration source. `README_MOC` step A (T/S profiles) has been rebuilt on top of it and rerun through 2023-09-24 for the 6 sites `mov_samba_marion*.m` uses (P1/P2/P4/P5/P6/P8) — Ekman/shelf/topo (steps B-D) and the final MOC estimate (step E) have not been rerun yet, so `Total_TPUD1..8` here still reflects the old 2013-2017 window.
- **West array's `samba_w.mat`/`IES_Make_Profiles_plusBrazil_plusDynHgt_Olga.mat` mystery is resolved**: it's almost certainly `DATA/West_PIES_Chris/IES_SAM.mat` (sites A/B/C/D) merged with `IES_brazil_SAM.mat` (sites AA/BB, a separate file the original `IES_USA_SAM.m` never loaded) plus a dynamic-height calculation — the naming lines up exactly. **Confirmed data-capped**: neither source file has anything past 2018-04-30 (A/B/C/D) / 2016-10-30 (AA/BB) despite `samba_w.mat`'s recent file-modification date — that regeneration was a reprocessing, not a data extension. Extending West coverage needs fresh data from Chris; nothing newer exists locally. Site CC has no usable data for a long stretch and is being treated as permanently excluded (per user).
- Also flagged there: `position_East.m`/`position_West.m`'s stale per-site positions (see the D-P8 gap example already noted above for PIES4) may equally affect `dx` in `mov_samba_marion*.m` here, via the same `cpiesE`/`cpiesW` structs — not yet checked in this repo specifically.

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
