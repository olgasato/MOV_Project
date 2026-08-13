# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A MATLAB workflow computing MOV (Meridional Overturning component of freshwater transport) and its full decomposition (mean + Mov + gyre) at 34.5°S, from the SAMBA array's inverted echo sounder (IES/PIES) mooring data at the western (Samba-W) and eastern (Samba-E) boundaries. Not a package — no build system, test suite, or dependency manifest; a sequence of MATLAB scripts run in order. Version-controlled with git, pushed to `https://github.com/olgasato/MOV_Project` (git identity: `olgasato <olga.sato@usp.br>`).

## Versioning convention

Scripts are iterated as `<name>_vN.m` rather than edited in place — e.g. `concat_IES.m` → `concat_IES_v4.m` → `concat_IES_v5.m` (current). Older versions are kept for history/comparison rather than deleted. When fixing or extending a script, save the result under the next version number rather than overwriting the file being changed.

## The two pipelines

Both start from the same Step 1 (`concat_IES_v6.m`), then diverge:

### Step 1: `concat_IES_v6.m` (current version)

Loads `samba_w.mat`/`samba_e.mat` (not in this repo — must be on the MATLAB path or in the working directory), restricts both to the common overlap period **2013-09-06 to 2022-12-11** (opened 2026-08-10; see below), and extracts per-mooring temperature/salinity as flat variables (`tem_A`, `tem_AA`, ..., `sal_P1`, `sal_P8`, ...) plus `prew`/`pree` (pressure), `depw`/`depe` (depth), `lonw`/`lone` (longitude, sorted), and `dt`. Saves everything to `concat_IESsamba.mat` (gitignored — regenerate by running the script, requires the raw `samba_w`/`samba_e` inputs).

**`v6` (2026-08-10) opened the date window.** `v5` hardcoded 2013-09-11–2017-07-17, matching the East array's availability at the time. Since then: `samba_w.mat` already covered 2009-03-18 to 2022-12-11 the whole time (unused past 2017 only because of this window); `samba_e.mat` was rebuilt from `Marions_code`'s Daily_Tau-based calibration (see "Upstream data pipeline" below) and now covers 2013-09-06 to 2023-09-24. The real overlap is bounded by `samba_w.mat`'s end date — **West is now the limiting side**, the reverse of the `v5`-era situation. Verified: `concat_IESsamba.mat` now spans 3383 days (~2.4x `v5`'s 1415), including West mooring **CC** for the first time (0% coverage in the old window; CC's instrument wasn't redeployed until mid-2019).

**Extended (2026-08-12)**: `README_MOC` steps B-E have now been rebuilt in `Marions_code` (see "Upstream data pipeline" below) — `mov_samba_marion_v11.m` consumes the new `Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022.mat` and the full Mov/mean/gyre decomposition now runs on the complete 2013-09-06 to 2022-12-10 range (3381 days).

Earlier versions (`concat_IES.m`, `concat_IES_v2.m`, `concat_IES_v4.m`, `concat_IES_v5.m`) are kept for comparison — `v4` aggregated all sites into single `tem_w`/`sal_w` matrices, which was a bug (mixed different moorings along the time axis, incompatible with downstream per-site calls); `v5` reverted to per-site variables while keeping `v4`'s usability improvements (progress messages, `-v7.3` save).

### Step 2a: `mov_samba.m` — dynamic height only, not maintained further

Computes geostrophic dynamic height anomaly (`gpan_*`) per mooring via `gsw_geo_strf_dyn_height`, referenced to 1000 dbar. Requires the GSW (TEOS-10) toolbox. The velocity/transport part is commented out as a sketch, superseded by the `mov_samba_marion*.m` line below.

### Step 2b: `mov_samba_marion_v13.m` (current version) — full Mov decomposition

Combines the Step 1 output with Marion's independently-computed absolute velocities (external file, see Dependencies below) to compute freshwater transport at 9 sites (A, C, D, P1, P2, P4, P5, P6, P8 — a subset of the full mooring list), decomposed into three components:

- **`mean`** — net/barotropic transport
- **`mov`** — overturning component (the classic "Mov")
- **`gyre`** — local (per-mooring-pair) deviation from the zonal-mean profiles

A `total_direct` (independently computed, non-decomposed) is compared against `mean+mov+gyre` as a self-consistency check, printed via `fprintf` — should be ~0 (currently closes to floating-point noise: mean=-1.06e-14 Sv, max=3.13e-13 Sv on the full 2013-2022 dataset).

**Citations** (in the script header): Marion's 2021 paper for the velocity fields (https://doi.org/10.1029/2020JC016947); two papers defining the mean/Mov/gyre decomposition being implemented (https://doi.org/10.1016/S0074-6142(01)80134-0, https://doi.org/10.1029/2023JC020558). Citations only — article PDFs are deliberately not committed to this (public) repo to avoid redistributing copyrighted journal content.

**Version history** (`mov_samba_marion.m` → `_v13.m`): each version's header comment documents what changed and why — worth reading before touching this script, since several non-obvious bugs were found and fixed in sequence:
- `v2`: added `mean`/`gyre` alongside the pre-existing `mov`, plus the `total_direct` residual check.
- `v3`–`v5`: chased down why the residual didn't close to zero. Root causes, in order found: (1) `V0` (velocity zonal-mean) and `wsal` (salinity zonal-mean) used inconsistent weighting — masking `V0` the same way as `wsal` cut the residual ~85%; (2) confirmed `dx1`'s validity mask is time-invariant, ruling out the `dx2`-vs-`dx1` time-slice mismatch; (3) `V00`/`S0` were unweighted depth-means, inconsistent with `dx2`'s depth-varying weighting — redefining both as `dx2`-weighted depth-means closed the residual to ~1e-16 Sv. **Note: this changed `mov`'s numeric value** relative to the original `mov_samba_marion.m` (used for the OSM26 presentation) — not a pure refactor.
- `v6`: robustness fixes anticipating more station data — `assert(isequal(prew,pree))` (the depth axis `pre=prew` is silently assumed identical for both West and East sites), the `dx1` time-invariance check kept as a permanent `warning()` rather than a one-off, monthly-means array sizing derived from `dt` instead of hardcoded (`47` months, `init=9`/`iend=7`). **Not addressed**: the site list itself (9 sites, 8 gaps) is still hardcoded throughout — `sal00`'s dimension, the `for ii=1:8` gap loop, `cpiesW([1 5 6])`/`cpiesE([1 2 4 5 6 8])` index selection, `p_cpiesW`/`p_cpiesE` cutoff arrays. Adding/removing a mooring station requires editing all of those by hand.
- `v7`: `v6`'s new assert fired for real — `prew` (501 points, 0-5000 dbar) and `pree` (531 points, 0-5300 dbar) had different *lengths*, not just values. Root cause: `ind=pree<=5000` was already used to trim `sal_P*`/`tem_P*` down to 501 points, but `pree` itself was never reassigned with that same `ind`. Fixed with one line (`pree=pree(ind);`). Not a real grid mismatch — both grids share spacing (10 dbar) and origin (0).
- `v8`: after `concat_IES_v6.m` opened `concat_IESsamba.mat`'s window to 2013-09-2022-12, this script would have silently broken — both `concat_IESsamba.mat` and `Full_Depth_OverturningEstimate_Cst_for_MHT.mat` (`Total_TPUD1..8`) have their own `dt`, and the second `load` was silently overwriting the first's. Harmless before only because both files happened to cover identical dates; not safe now that they don't (`Total_TPUD` still reflects the old 2013-2017 window — `README_MOC` steps B-E haven't been rebuilt). Fixed by keeping both `dt`'s separately and explicitly aligning via `[dt,i_TS,i_TPUD]=intersect(dt_TS,dt_TPUD)`, then indexing every time-dependent array by the matching `i_TS`/`i_TPUD` rather than assuming positional alignment. Verified: aligns to 1405 common days (2013-09-11 to 2017-07-16, `Total_TPUD` being the narrower/limiting side), residual still ~1e-16/1e-14 Sv. Output range is unchanged until `Total_TPUD` itself is extended, but the script no longer silently mismatches or crashes now that its two inputs differ in length.
- `v9`: `README_MOC` steps B-E were rebuilt in `Marions_code` (2026-08-12, see "Upstream data pipeline" below) — `Total_TPUD1..8` now live in `Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022.mat` (2013-09-06 to 2022-12-11). Updated the load path accordingly. Also fixed a latent bug: `pre_0`/`difpre` hardcoded the day count to `1405` (the `v8`-era dataset size) instead of deriving it from the actual data — harmless while `difpre` itself was unused downstream, but a landmine waiting to happen.
- `v10`–`v11`: `v9` ran (3381 aligned days, West now the limiting side, matching `concat_IES_v6`'s prediction) but the residual check — closed to ~1e-16 Sv since `v5` — blew up to mean=-2.24 Sv, max=219.9 Sv. Three independent, compounding bugs, all latent since `v5` but only triggered once real (not just hypothetical) time-varying, non-identical coverage gaps existed between salinity and velocity data:
  1. `dx2` (the width profile weighting `mov`/`mean_term`) was built from day 1's coverage mask only (`dx2=nansum(squeeze(dx1(:,1,:)),2)`) and silently reused as if constant across all days — true when coverage genuinely was constant (`v5`-era 2013-2017), false once real per-site gaps (CC pre-2019, P5/P6 ~58% coverage, etc.) made 22-38% of points vary day to day. Fixed (`v10`) by making `dx2` itself time-resolved (`squeeze(nansum(dx1,3))`, i.e. summing the real per-day mask instead of day 1's).
  2. `dx1` (and therefore `dx2`) was masked by salinity's NaN pattern only. `Total_TPUD`'s own coverage gaps are real and independent of salinity's — diagnosed 1.27M of 13.55M (depth,time,gap) points with valid salinity but NaN velocity. Those points kept full width in `dx2`/`wsal`'s denominator while contributing nothing to `V0`'s numerator, biasing `V0` relative to `wsal`. Fixed (`v11`) by building `geo`/`vel00` *before* `dx1`, then masking `dx1` by the union of both NaN patterns.
  3. The actual dominant cause: `S0=nansum(wsal.*dx2)/nansum(dx2)` and `V00=nansum(V0.*dx2)/nansum(dx2)` used `/` (MATLAB matrix right division), not `./` (elementwise). While `dx2` was a plain `[depth x 1]` column, `nansum(dx2)` was always a true scalar, where `/` and `./` agree — harmless. Once `dx2` became `[depth x time]` (fix #1 above), `nansum(dx2)` became `[1 x time]` too, and `A/B` between two same-size non-scalar arrays is a least-squares matrix division, not elementwise — this silently collapsed `S0` (and `V00`) to a **single best-fit scalar for the entire record** (confirmed: `size(S0)` was `[1 1]`, not `[1 3381]`) instead of one value per day, biasing every downstream quantity. Fixed (`v11`) by changing both to `./`. Also changed `mean_term`'s `sum(dx2)` to `nansum(dx2)` while in there — this repo's custom `nansum.m` (`Marions_code/functions/nan/`) returns NaN (not 0) when an entire reduced dimension is NaN, unlike plain `sum`, so a depth level with zero valid gaps on some day would have propagated NaN across that whole day's `mean_term`; didn't fire on this dataset but was inconsistent with every other reduction in the script.
  Confirmed (`v11`): residual mean=-1.06e-14 Sv, max=3.13e-13 Sv — back to the `v5`-era floating-point-noise baseline, now on the full 3381-day 2013-2022 dataset.
- `v12`: fixes the `AREA_TOPO`/salinity-cutoff inconsistency (see the "Upstream data pipeline" section below for the full quantified writeup — this bullet covers the `mov_samba_marion*.m`-side implementation only). Switched `geo`/`vel00` from `Total_TPUD1..8` to `Absolute_TPUD1..8_preTopo` (Step E `v2`'s new pre-`AREA_TOPO` output), and built `dxA = dx1 .* AREA_TOPO_gap(z)` — the real per-gap open width — used everywhere `dx1` previously weighted an average (`dx2`, `wsal`, `V0`, `aux_gyre`, `aux_total`). Two bugs surfaced while validating against a new sanity check (`vel00.*dxA` should reconstruct `Total_TPUD` exactly): (1) `V0` was `nansum(geo,3)./nansum(dxA,3)` — `geo` alone has no `AREA_TOPO` factor, mismatching a non-topo-weighted numerator against a topo-weighted denominator (residual blew back up to mean=-14.49 Sv); fixed to `nansum(vel00.*dxA,3)./nansum(dxA,3)`. (2) `Absolute_TPUD_preTopo` is captured *before* Step E's Ekman addition too (not just before `AREA_TOPO`), so it silently excluded Ekman from the entire calculation; fixed by having Step E `v2` also save `Ekman_TPUD1..8` and adding it into `geo`/`vel00` at the 7 surface rows (confirmed `AREA_TOPO_i==1` there for all 8 gaps, so no extra scaling needed). Confirmed: residual mean=-8.83e-16 Sv, max=1.71e-13 Sv; sanity check ~0 at every row. **Impact on the actual series**: `mean_term` essentially unchanged (diff ~1e-15 Sv); `mov` changes modestly (mean|diff|=0.007 Sv, max|diff|=0.12 Sv on 2018-09-27, corr=0.997 vs. `v11`); `gyre` barely moves (mean|diff|=0.001 Sv, corr=0.9998) — a real, necessary correction, but not a dramatic one numerically.
- `v13`: adds the coverage/confidence diagnostics for the 2017-08/2019-10 low-variance window (see "Known limitation" below) — two new per-day outputs, `n_gaps_valid` (how many of the 8 gaps have any valid depth that day, 0-8) and `coverage_totalwidth` (`nansum(dx2)`, the same total open width `mean_term` is weighted by). Also **enabled the save line** (`save('mov_samba_marion_v13.mat', ...)`), which had been commented out since at least `v8` — this script previously computed everything but never persisted it beyond the `mov_IES.png` plot.

## Known limitation: the 2017-08/2019-10 low-variance period in `mov` is a coverage artifact, not real quiescence

**Observation (user, 2026-08-12):** looking at `mov_IES.png`, `mov`'s variance is visibly high from the start of the record to ~2017-09, much lower from ~2017-09 to ~2019-09, then high again from ~2019-09 to the end — but `mean_term` and `gyre` don't show this pattern. Investigated and confirmed quantitatively; root-caused to a real, ~2-year loss of spatial coverage in the deep/eastern part of the array, not a genuine oceanographic quiet period.

**Confirmed with `std()` over three periods** (split points chosen from the actual instrument-gap boundaries found below, not from eyeballing the plot):

| | 2013-09–2017-08 | 2017-08–2019-10 | 2019-10–2022-12 |
|---|---|---|---|
| `std(mov)` | 0.144 | **0.090** | 0.188 |
| `std(mean_term)` | 39.26 | 42.23 | 41.30 |
| `std(gyre)` | 0.080 | 0.101 | 0.118 |

`mov` drops ~40% in the middle period and recovers after; `mean_term` stays flat (even ticks up slightly) and `gyre` rises monotonically. Neither shows the dip `mov` does.

**Root cause — real instrument outages, found via per-site monthly coverage fractions:**
- **P5 and P6: essentially zero coverage from 2017-08 to 2021-09/10** (~4 years). This was already noted in passing in `Marions_code/CLAUDE.md` ("P5/P6 have no instrument 2017-2021") but its effect on `mov`'s variance specifically hadn't been traced through until now.
- **P8**: zero coverage 2018-04 to 2019-09/10.
- **P4**: zero coverage 2019-02 to 2019-09/10.
- **D**: a shorter gap, 2017-10 to 2018-04/05.
- A, C, P1, P2 have essentially continuous coverage throughout.

**Per-gap valid-day fraction in each period** (fraction of days with at least some valid depth for that gap):

| Gap | 2013-09–2017-08 | 2017-08–2019-10 | 2019-10–2022-12 |
|---|---|---|---|
| A–C | 1.00 | 1.00 | 1.00 |
| C–D | 1.00 | 0.74 | 0.99 |
| D–P8 | 1.00 | 0.08 | 0.99 |
| P8–P6 | 0.99 | **0.00** | 0.38 |
| P6–P5 | 0.99 | **0.00** | 0.38 |
| P5–P4 | 0.99 | **0.00** | 0.38 |
| P4–P2 | 1.00 | 0.68 | 0.99 |
| P2–P1 | 1.00 | 1.00 | 1.00 |

Three of the eight gaps (`P8–P6`, `P6–P5`, `P5–P4` — everything touching P5/P6) have **zero** coverage for essentially the entire 2017-08/2019-10 window. Total section width (`sum(dx2)` per day) drops from ~3.8e9 (period 1) to ~9.5e8 (period 2, a **~75% drop**) and back to ~3.6e9 (period 3) once P4/P8 return in late 2019 (P5/P6 themselves don't fully return until 2021, but that alone is apparently enough to restore most of `mov`'s variance).

**Why `mov` specifically, and not `mean_term`/`gyre`:** `mov` is the depth-integrated correlation between the zonal-mean velocity and salinity *anomaly profiles* (`vel_prime(z)`, `sal_prime(z)`) — it's fundamentally about baroclinic structure spanning the section's full depth. The three gaps lost are precisely the deep/eastern ones, where the deep return branch of the overturning circulation would be expected to show up. Losing them doesn't just remove noise — it removes the array's ability to see genuine deep baroclinic variability, so what's left (computed from a reduced, shallower/more-western subset) is quieter *because the observing system's sensitivity dropped*, not because the ocean did. `mean_term` (a single barotropic number per day, not sensitive to which particular gaps are present as long as depth-weighting is consistent) and `gyre` (per-gap local deviations, which don't need every gap simultaneously present) are structurally less exposed to this.

**Practical implication:** treat `mov` values in 2017-08 to 2019-10 as **less reliable / likely biased toward zero variance** than the rest of the record — this is a data-coverage artifact, not a real trend or event, and should be flagged as such in any analysis or plot that uses this window.

**Addressed in `v13` (2026-08-13):** `mov_samba_marion_v13.m` now saves two per-day coverage diagnostics alongside `mov`/`mean_term`/`gyre`/`total_direct` in `mov_samba_marion_v13.mat` — `n_gaps_valid` (how many of the 8 gaps have any valid depth that day, 0-8) and `coverage_totalwidth` (`nansum(dx2)`, the same total open width `mean_term` is weighted by). Confirmed these reproduce the table above automatically: `mean(n_gaps_valid)` = 7.96 / 3.51 / 6.10 across the three periods, `mean(coverage_totalwidth)` drops ~80% in the middle period. Downstream users can now screen low-coverage days programmatically (e.g. `n_gaps_valid < 6`) instead of re-deriving this by hand. The underlying data gap itself is of course not fixable — no amount of code can recover an instrument that wasn't there.

**Distinct from the `AREA_TOPO`/salinity-cutoff issue below** — that one is a *static*, depth-based masking inconsistency present throughout the *entire* record; this one is a *time-varying*, real instrument-outage issue affecting one specific ~2-year window. Both are legitimate, independent caveats.

## Upstream data pipeline (outside this repo)

`mov_samba_marion*.m`'s inputs (`Total_TPUD1..8`, and transitively `samba_w.mat`/`samba_e.mat`) come from a much larger MATLAB pipeline at `~/research/sambar/renellys_sent/Marions_code/`, documented in `README_MOC` and `README_MHT` there (not in this repo). Traced through steps A–E of `README_MOC` on 2026-08-07; summary below so this doesn't need to be re-derived from scratch.

**Chain:** raw PIES travel-time calibration (`README_CALIBRATE_IES_EAST`) → GEM fields (`README_GEM`) → step A `IES_Make_Profiles_*.m` (combines travel time + GEM → T/S profiles, this is what `samba_w.mat`/`samba_e.mat` ultimately derive from) → step B `Ekman_transports.m` (CCMP wind → Ekman transport) → step C `ECCO_TPUD_ShelfBits.m` (static, time-mean shelf-transport correction from ECCO, ~-4.4 Sv west / ~0 Sv east) → step D `correct_topo_all_MSM60.m` (static, per-gap continuous topography/open-area fraction `AREA_TOPO_1..8`, from real MSM60 cruise bathymetry) → step E `Full_Depth_OverturningEstimate_Cst_for_MHT.m` (combines all of the above into `Total_TPUD1..8`, saved to `Full_Depth_OverturningEstimate_Cst_for_MHT.mat`).

**Steps B-E rebuilt for 2013-2022 (2026-08-12).** Step B: CCMP wind data extended to 2018-2022 (fetched from RSS's public archive, the original `/data/kersale` mount being unreachable), `Ekman_transports_9pies_2009_2022.m` now covers 2009-2022. Steps C/D confirmed static/climatological, unchanged. Step E: `Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022.m` switches its West pressure/T/S source to `samba_w.mat` (from the old, 2018-04-30-capped source), pads `samba_w.mat`'s pressure grid (only to 5000 dbar) back to the shared 5150 dbar/516-level grid `topo_corr_msm60.mat`'s `AREA_TOPO_*` are fixed to, and replaces three hardcoded 2013-2017 date trims with a 3-way `intersect()` alignment. Output (`Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022.mat`) now spans 2013-09-06 to 2022-12-11 (3382 days, was 1405). A real bug was found and fixed while validating this: `Daily_Tau`'s bottom-pressure column uses `-9999` as a documented "no data" flag that `Marions_code/ies/read_daily_tau.m` never converted to NaN (unlike `tau1000`, which escaped by coincidence via a downstream range filter) — caused ~1e9-scale transport blow-ups in 4 of 8 segments until fixed at the source. Full detail in `Marions_code/CLAUDE.md`.

**Confirmed: `Total_TPUD1..8` already includes the topography correction, absolute reference, and Ekman.** Step E's script builds `Absolute_TPUD1..8` (baroclinic shear from dynamic height + absolute reference from real PIES bottom-pressure sensors + an ECCO mean-velocity offset at 1500 dbar), multiplies each by its per-gap `AREA_TOPO_i` profile, then adds per-gap Ekman (`Ekman_TPUD1..8`) — that sum *is* `Total_TPUD1..8`. So `mov_samba_marion*.m` does not need to (and currently does not) reapply any topography correction on the velocity side — it's already baked in.

**Resolved 2026-08-13 (was: "open question, partially investigated, not yet fixed"):** `mov_samba_marion*.m`'s own salinity masking (`p_cpiesW`/`p_cpiesE`, a binary NaN cutoff at a fixed depth per *site*) is inconsistent with the continuous, per-*gap* `AREA_TOPO_i` weighting already applied to the velocity it's multiplied against in `mov`/`mean`/`gyre`. E.g. the D–P8 gap's `AREA_TOPO_3` starts dropping below 1.0 at 940 dbar (a mid-gap ridge), far shallower than either site's own individual cutoff (D: 4850, P8: 1280) — so the velocity smoothly attenuates there while the salinity is either fully valid or fully NaN depending on the site-level cutoff.

Quantified per gap: within each gap's salinity-valid depth range, the fraction of nominal cross-sectional width that's actually topographically blocked (`1-AREA_TOPO`, "phantom width" — counted at full weight for `wsal`/`dx1`/`dx2` even though the real open area there is smaller or zero):

| Gap | Salinity cutoff (dbar) | Mean `AREA_TOPO` in that range | Phantom width |
|---|---|---|---|
| A–C | 4620 | 0.75 | **24.7%** |
| C–D | 4850 | 0.99 | 1.3% (negligible) |
| D–P8 | 4850 | 0.84 | **15.6%** |
| P8–P6 | 2150 | 1.00 | 0.0% |
| P6–P5 | 4560 | 1.00 | 0.0% |
| P5–P4 | 5000 (uncapped) | 0.96 | 4.4% |
| P4–P2 | 5000 (uncapped) | 0.78 | **21.9%** |
| P2–P1 | 5000 (uncapped) | 0.35 | **64.8%** |

`P2–P1` is the worst case: `AREA_TOPO` is essentially 0 by the salinity cutoff, meaning nearly two-thirds of the "width" counted toward `wsal`/`V0`'s zonal averaging in the deep part of that gap is topographically nonexistent. This does **not** break the `mean+Mov+gyre == total_direct` algebraic identity (that holds regardless of what `dx1` contains, as long as it's used consistently — see the `v10`/`v11` fixes above), but it likely biases the *absolute values* of `V0`/`wsal`/`S0`/`mov`/`mean_term` toward whatever the "phantom" depths' salinity/velocity happen to be.

**How it was fixed (2026-08-13):** Step E was rebuilt as `Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.m` (`Marions_code/Full_Depth_MOC/`), which additionally saves `Absolute_TPUD1..8_preTopo` (a copy taken right before the `AREA_TOPO_i` multiplication overwrites `Absolute_TPUD1..8` in place) and `Ekman_TPUD1..8` (needed because `Absolute_TPUD_preTopo` is captured before Ekman is added too — see `mov_samba_marion`'s `v12` note above), plus re-saves `AREA_TOPO_1..8` themselves. `mov_samba_marion_v12.m` then uses `Absolute_TPUD*_preTopo` for `vel00`/`geo` and weights `V0`/`wsal`/`dx2`/`aux_total`/`aux_gyre` by `dxA = dx1 .* AREA_TOPO_gap(z)` (the real per-gap open width) instead of the flat `dx1`. A sanity check (`vel00.*dxA` should reconstruct `Total_TPUD` exactly) confirms the plumbing is correct, and the `mean+Mov+gyre==total_direct` identity still closes to floating-point noise (re-derived from scratch with the new definitions, not just assumed). See `mov_samba_marion`'s own version history above for the two bugs found while getting there, and the resulting `mov`/`mean_term`/`gyre` impact (modest: `mov` corr=0.997 vs. the old flat-`dx1` version).

**Date range provenance:** the `2013-09-11`–`2017-07-17` window used throughout this repo isn't an independent choice — it's inherited directly from a hardcoded trim in `Full_Depth_OverturningEstimate.m` (`bad=find(dt_SAM<datenum(2013,9,11) | dt_SAM>datenum(2017,7,17))`), matching the East PIES array's (Samba-E / "GH") availability window at the time that pipeline was last run.

**Data currency (relevant to "add more station data") — UPDATE 2026-08-10, see below for what changed.** Originally (2026-08-07): both the East PIES calibration (`Merge_SiteP*_records.m`, hardcoded to only merge the `2013_2015`/`2015_2017` legs) and the CCMP wind data (`Reformat_WindData.m`, hardcoded to only merge `LoadWinds_2009.mat`...`LoadWinds_2017.mat`) stopped at 2017. West PIES calibration is handled externally ("Chris" provides pre-calibrated `IES_SAM.mat`) and extends to 2018-04-30; `samba_w.mat` is a byte-identical copy of `.../IES profile data/IES_Make_Profiles_plusBrazil_plusDynHgt_Olga.mat`.

**Resolved (2026-08-10), in the sibling repo `~/research/sambar/renellys_sent/Marions_code/` (now version-controlled, `https://github.com/olgasato/SAMBA_E_IES`, own `CLAUDE.md` there with full detail):**
- **East array**: found `~/research/sambar/samba_e/Daily_Tau-*.zip`, an officially published, already-calibrated dataset (SAMOC-SA/DFFE) covering **2013-09 to 2023-09** for P1/P2/P3/P3a/P4/P4a/P5/P6/P8 — cross-validated against a from-scratch recalibration of 2017-2019 (agreed to 0.0024s/0.16 dbar) and adopted as the primary calibration source. `README_MOC` step A (T/S profiles) has been rebuilt on top of it and rerun through 2023-09-24 for the 6 sites `mov_samba_marion*.m` uses (P1/P2/P4/P5/P6/P8). Steps B-E have since also been rebuilt (2026-08-12, see above) — `Total_TPUD1..8` now reflects 2013-2022, not 2013-2017.
- **West array's `samba_w.mat`/`IES_Make_Profiles_plusBrazil_plusDynHgt_Olga.mat` mystery is resolved**: it's almost certainly `DATA/West_PIES_Chris/IES_SAM.mat` (sites A/B/C/D) merged with `IES_brazil_SAM.mat` (sites AA/BB, a separate file the original `IES_USA_SAM.m` never loaded) plus a dynamic-height calculation — the naming lines up exactly.
- **CORRECTED 2026-08-11** (an earlier version of this note got the date range wrong — checked a different, older file and mistakenly attributed its range here without directly verifying `samba_w.mat` itself): `samba_w.mat` actually covers **2009-03-18 to 2022-12-11**, not just to 2018. Per-site real coverage: A/C 100%, D 95.7%, B 69.2%, BB 61.8%, AA 54.8%, **CC 25.1% (from 2019-06-27 on — CC is not permanently absent, just had a long equipment gap before that)**. `samba_w.mat` is already the adopted primary West source (see `Marions_code/CLAUDE.md` for full detail) — **no further West-side data work is needed to reach the same ~2022-2023 ceiling the East array's Daily_Tau rebuild reaches**. What's actually discarding this already-available range is `concat_IES*.m`'s own hardcoded `2013-09-11`/`2017-07-17` trim — opening that window (on both `concat_IES*.m` and `mov_samba_marion*.m`) is the real next step for extending this repo's output, not more data acquisition. (Raw West data reportedly exists through June 2025 per the user, not yet processed into tau1000 — not needed yet since the East array caps the combined calculation at 2023 regardless.)
- Also flagged there: `position_East.m`/`position_West.m`'s stale per-site positions (see the D-P8 gap example already noted above for PIES4) may equally affect `dx` in `mov_samba_marion*.m` here, via the same `cpiesE`/`cpiesW` structs — not yet checked in this repo specifically.

## External dependencies (not in this repo)

- `~/research/sambar/renellys_sent/Marions_code/functions/positions_pies/position_West.m` and `position_East.m` — define the `cpiesW`/`cpiesE` structs (site name/lat/lon) used for inter-site distances. Added via `addpath(...)` at the top of `mov_samba_marion_v6.m`; not on the default MATLAB path (`~/matlab/startup.m` only adds `~/matlab/*`).
- `~/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat` (was `..._Cst_for_MHT.mat` through `v8`, `..._2013_2022.mat` through `v11`) — Marion's absolute velocity/transport-per-unit-depth data, loaded via an explicit absolute path (`prefix` variable). `v12` onward reads `Absolute_TPUD1..8_preTopo`, `AREA_TOPO_1..8`, and `Ekman_TPUD1..8` from this file (not just `Total_TPUD1..8`) — see the `AREA_TOPO`/salinity-cutoff fix above.
- `~/research/sambar/renellys_sent/Marions_code/functions/nan/nansum.m`/`nanmean.m` — custom NaN-aware reductions (Jan Gläscher's NaN-suite), NOT the Statistics Toolbox's built-ins. Notably, `nansum(X,dim)` returns NaN (not 0) when the entire reduced dimension is NaN — differs from naive expectations and was the proximate trigger for a real bug in `v11` (see version history above: `A/nansum(...)` silently doing matrix division instead of elementwise once the divisor stopped being a scalar).
- Custom utility functions assumed on the path but not in this repo: `lowpass_filter`, `sinfitb_tot`, `unico`.
- GSW (TEOS-10) Gibbs SeaWater MATLAB toolbox (`gsw_geo_strf_dyn_height`, `gsw_distance`) — on the path via `~/matlab/startup.m`.

## Site/mooring naming

- Samba-W (western boundary) moorings: `A, AA, B, BB, C, CC, D`, with corresponding depths `depw = [1350 2902 3510 4173 4558 4730 4756]` m.
- Samba-E (eastern boundary) moorings: `P1, P2, P4, P5, P6, P8` (P3 and P7 are missing from the array), with corresponding depths `depe = [1266 2129 4482 4969 5185 4608]` m.
- Variable naming convention throughout: `<field>_<site>`, e.g. `tem_AA`, `sal_P4`.
- `mov_samba_marion*.m` uses a 9-site subset (A, C, D, P1, P2, P4, P5, P6, P8) — the ones Marion's velocity dataset covers.

## Outputs

- `concat_IESsamba.mat`, `gpan_samba.mat`, `mov_samba_marion_v13.mat` — gitignored (`.mat` files excluded generally; regenerate by running the corresponding script). The latter (`dt`, `mov`, `mean_term`, `gyre`, `total_direct`, `n_gaps_valid`, `coverage_totalwidth`) is new as of `v13` — the save line existed but was commented out in every version before that.
- `mov_IES.png` — plot output from `mov_samba_marion*.m`, committed (not gitignored).
