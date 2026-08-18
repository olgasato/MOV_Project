# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A MATLAB workflow computing MOV (Meridional Overturning component of freshwater transport) and its full decomposition (mean + Mov + gyre) at 34.5°S, from the SAMBA array's inverted echo sounder (IES/PIES) mooring data at the western (Samba-W) and eastern (Samba-E) boundaries. As of 2026-08-14, also has sibling MHT (heat transport, `mht_samba_marion*.m`) and MOC (volume-transport streamfunction, `moc_streamfunction_v1.m`) lines — see their own sections below. Not a package — no build system, test suite, or dependency manifest; a sequence of MATLAB scripts run in order. Version-controlled with git, pushed to `https://github.com/olgasato/MOV_Project` (git identity: `olgasato <olga.sato@usp.br>`).

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

### Step 2b: `mov_samba_marion_v15.m` (current version) — full Mov decomposition

Combines the Step 1 output with Marion's independently-computed absolute velocities (external file, see Dependencies below) to compute freshwater transport at 9 sites (A, C, D, P1, P2, P4, P5, P6, P8 — a subset of the full mooring list), decomposed into three components:

- **`mean`** — net/barotropic transport
- **`mov`** — overturning component (the classic "Mov")
- **`gyre`** — local (per-mooring-pair) deviation from the zonal-mean profiles

A `total_direct` (independently computed, non-decomposed) is compared against `mean+mov+gyre` as a self-consistency check, printed via `fprintf` — should be ~0 (currently closes to floating-point noise: mean=-1.06e-14 Sv, max=3.13e-13 Sv on the full 2013-2022 dataset).

**Citations** (in the script header): Marion's 2021 paper for the velocity fields (https://doi.org/10.1029/2020JC016947); two papers defining the mean/Mov/gyre decomposition being implemented (https://doi.org/10.1016/S0074-6142(01)80134-0, https://doi.org/10.1029/2023JC020558). Citations only — article PDFs are deliberately not committed to this (public) repo to avoid redistributing copyrighted journal content.

**Version history** (`mov_samba_marion.m` → `_v15.m`): each version's header comment documents what changed and why — worth reading before touching this script, since several non-obvious bugs were found and fixed in sequence:
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
- `v14`: visualizes the low-coverage window directly on `mov_IES.png` instead of either leaving it unmarked (`v13` — looks like a real quiet period) or masking `mov` to `NaN` there. A one-off diagnostic script (`mov_samba_marion_v13_9pies.m`, not kept) tried the latter: requiring all 9 of Marion's original sites simultaneously present leaves a blank ~4-year hole (2017-08/2021-09) and discards the ~54% of days that DO have partial, still-informative coverage — user's reaction: "ficou um 'buraco' no meio da série" (a "hole" in the middle of the series), not an improvement. `v14` instead shades the plot background wherever a 30-day rolling mean of `coverage_totalwidth` drops below 50% of the record's median (via `fill(...)`, drawn before the `mov` line so it sits behind it) — a single stable stretch, **06-Oct-2017 to 05-Oct-2019** (729 days, 21.6%), not very sensitive to the exact threshold (30-60% of median all identify essentially the same stretch). Narrower than the full P5/P6 outage (2017-08/2021-09) because `coverage_totalwidth` recovers enough once P4/P8 return in Nov 2019, even though P5/P6 themselves don't return until 2021. `mov`/`mean_term`/`gyre`/`total_direct` themselves are **not modified or masked** — every value is the same real, fully-computed estimate as `v13`; only the plot changes (adds `low_coverage`, the shading flag, to the saved `.mat` too).
- `v15`: replaces `v14`'s single binary low-coverage flag with graded shading bands by `n_sites_valid` (how many of the 9 original sites are present that day, 0-9) — follow-up question: is there a minimum site count more informative than just "low"/"not low"? A one-off test (`test_minsites.m`, not kept) swept thresholds `>=6/7/8/9` sites: `>=7` almost exactly reproduces `v14`'s `coverage_totalwidth`-based stretch (78.1% of days kept, one gap 2017-10-05 to 2019-10-11 — a useful cross-check that the two independent criteria, width-based and site-count-based, agree); `>=8` and `>=9` are nearly identical to each other (54.8% vs. 54.4% kept — P5 and P6 are almost always either both missing or both present); `>=6` keeps 90.7% of days with only two short gaps (59 and 255 days). A second test (`test_colorband.m`, not kept) turned this into 4 shading bands instead of a single threshold: **9/9 sites** (54.4% of days, unshaded), **7-8/9** (23.7%, light gray), **6/9** (12.5%, medium gray), **<6/9** (9.3%, dark gray) — revealing a previously unflagged patch of reduced (7-8/9) coverage in 2020-2022 that `v14`'s single stretch missed entirely. `mov`/`mean_term`/`gyre`/`total_direct` are still **not modified or masked**, identical to `v13`/`v14` — only the shading changed (adds `n_sites_valid` and `category`, the 4 band assignments, to the saved `.mat`).

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

**Addressed in `v13`-`v15` (2026-08-13):** `mov_samba_marion_v13.m` saves per-day coverage diagnostics alongside `mov`/`mean_term`/`gyre`/`total_direct` — `n_gaps_valid` (how many of the 8 gaps have any valid depth that day, 0-8) and `coverage_totalwidth` (`nansum(dx2)`, the same total open width `mean_term` is weighted by). Confirmed these reproduce the table above automatically: `mean(n_gaps_valid)` = 7.96 / 3.51 / 6.10 across the three periods, `mean(coverage_totalwidth)` drops ~80% in the middle period. `v14` shaded a single low-coverage stretch (06-Oct-2017 to 05-Oct-2019) directly on `mov_IES.png`; `v15` replaced that with 4 graded shading bands by `n_sites_valid` (9/7-8/6/<6 of the 9 original sites present), which also surfaced a second, previously-unflagged reduced-coverage patch in 2020-2022 — see `v15`'s entry in the version history above for the full comparison against requiring a fixed minimum site count. Downstream users can screen low-coverage days programmatically (e.g. `n_sites_valid < 7`) via the saved `.mat`, or just look at the plot. The underlying data gap itself is of course not fixable — no amount of code can recover an instrument that wasn't there.

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

## MHT (heat transport): `mht_samba_marion_v1.m`, a new sibling line to `mov_samba_marion*.m` (2026-08-13)

User's request: reuse the MOV pipeline's fixes/lessons to also compute MOC (deferred — see below) and MHT (meridional heat transport). Unlike freshwater transport, MHT has no natural reference-temperature analog to `S0` (the classic "reference temperature problem" in the MHT literature), so this does **not** reuse `mov_samba_marion*.m`'s mean+overturning+gyre decomposition math directly — it's adapted instead from Marion's own `Marions_code/MHT/MHT_Estimate_constituents.m` (the "final step" of `README_MHT`), which computes **total** heat transport (`rho*cp*integral(v*T dz)`, no anomaly/reference subtraction) plus a sensitivity/attribution breakdown (Ekman/relative/reference velocity held constant in turn), not (yet) a spatial mean/overturning/gyre split.

**Two-phase plan, `v1` covers phase 1-2 only:**
1. **Faithful reproduction** of the original 2013-2017 method, unmodified — confirmed exact match against the existing reference (`Marions_code/MHT/Wrk/MHT_estimates_Constituents.mat`): `mean(Heat_total)=0.5466` PW to 4 decimal places. Done *before* any extension, to validate the method was understood correctly. See `Marions_code/CLAUDE.md`'s README_MHT section for detail.
2. **Extended to 2013-2022** (`mht_samba_marion_v1.m`, this repo): updates the East and West seasonally-corrected temperature sources (`Marions_code/CLAUDE.md`) and switches to Step E `v2`'s `Total_TPUD1..8` (+ `_RefConst`/`_RelConst`/`_EKMANconst` variants), replacing the three hardcoded 2013-2017 date trims with a 3-way `intersect()` (same pattern as everywhere else in this project). Result: `Heat_total` etc. now span **2013-09-06 to 2022-12-11** (3382 days). Reuses `mov_samba_marion_v15`'s coverage-shading approach (same `n_sites_valid`/graded-band logic, since it's the same array and outages) on the output plot (`mht_IES.png`).
   - **Deliberate simplification, not yet resolved**: gap 3 (D-P8) originally used a satellite-altimetry-derived "interior" temperature, not a simple site average — that sub-pipeline is itself capped at 2017 and extending it is a materially bigger undertaking (new AVISO SLA data + its own GEM/seasonal-correction chain). Decided with the user to use `(Temperature_D+Temperature_P8)/2` instead for now (same convention as every other gap). **Quantified impact**: over the *same* 2013-2017 window as the reference, this shifts `mean(Heat_total)` from 0.5466 to 0.6251 PW (corr=0.99 against the reference — the *shape* of variability is well preserved, but the *absolute level* is measurably affected by ~0.08 PW). Revisit the altimetry pipeline before treating `Heat_total`'s absolute magnitude as final.

**Phase 3 (`mht_samba_marion_v2.m`, 2026-08-14): mean+overturning+gyre decomposition for heat**, derived from scratch (not ported from `mov_samba_marion` directly, given the reference-temperature difference). Also switched the velocity side to `Absolute_TPUD*_preTopo` + explicit `dxA=dx(gap)*AREA_TOPO(z)` weighting (the `mov_samba_marion_v12` fix) — same phantom-width bias applies here, same array/velocity field.

Derivation (`T0` = `dx2`-weighted section-mean temperature, same convention as `S0`):
- `Q_overturning = rho*Cp*sum_z[(V0-V00)*(wtemp-T0)*dx2]*dz` — confirmed algebraically (and empirically, by construction) **invariant to the `T0` reference choice**: shifting `T0` by a constant adds a term proportional to `sum_z[(V0-V00)*dx2]`, which is zero since `V00` is *defined* as that weighted mean.
- `Q_mean = rho*Cp*V00*T0*sum_z(dx2)*dz` — the **only** term that depends on the `T0` choice (unlike freshwater, where the analogous `1/S0` factor cancels algebraically out of `mean_term` — heat has no such cancellation, so `T0` must appear explicitly). This is exactly "the reference temperature problem" from the MHT literature: the overturning (baroclinic) heat transport is a robust, reference-independent quantity, but the net/mean (barotropic) term genuinely depends on an arbitrary reference choice.
- `Q_gyre` — computed and closes the identity, unlike MOC's (see below) since there IS a second field (temperature) to correlate against here.

**Real bug caught by cross-checking against `v1`**: first draft copied `mov_samba_marion`'s leading `-` sign into all four `Q_*` formulas by reflex. That sign is specific to freshwater's salt-flux-to-freshwater-flux conversion (`F=-(1/S0)*salt_flux`) — not a general convention. Heat transport is `Q=rho*Cp*integral(v*T dA)` directly, with no leading minus (matching `v1`/the original `MHT_Estimate_constituents.m`'s `Heat_total`). The sign bug gave `Q_total_direct=-0.57` PW — wrong sign *and*, once fixed, `Q_total_full` (`Q_total_direct` + static shelf/deep) came out `0.352` PW, closely matching `v1`'s `Heat_total` (`0.332` PW) over the same period — a good independent cross-check that both scripts are internally consistent with each other, not just self-consistent.

**Results** (2013-2022 mean): `Q_mean=-0.14` PW, `Q_overturning=0.76` PW (std `0.59` PW), `Q_gyre=-0.05` PW, `Q_total_direct=0.57` PW, `Q_total_full=0.35` PW. Residual (`Q_total_direct` vs. `Q_mean+Q_overturning+Q_gyre`) closes to ~3 W max (vs. values of order 1e14-1e15 W — floating-point noise). Static shelf/deep contributions (`Q_deep`, `Q_shelfE`, `Q_shelfW`) are kept separate from the mean/overturning/gyre split (external additions to what the CPIES array itself measures, not something this decomposition can attribute spatially) — only folded into `Q_total_full` for comparison against `v1`. Same graded coverage shading as `mov_samba_marion_v15`/`mht_samba_marion_v1` on the output plot (`mht_overturning_IES.png`, showing `Q_overturning` specifically, since that's the reference-independent, physically robust quantity).

## MOC: attempted 2026-08-14, the "reuse mov_samba_marion, drop the tracer" approach doesn't work -- mathematically, not a data problem

Original plan: define MOC as "the overturning-circulation volume-transport component" (Sv, no salinity/temperature weighting), reusing `mov_samba_marion`'s `V0`/`V00`/`dx2`/`dxA` machinery directly, over the classical streamfunction-maximum definition. Built `moc_samba_marion_v1.m` this way (simpler than `mov_samba_marion` structurally too — no T/S needed at all, so no `concat_IESsamba.mat`, no cross-source `dt` alignment, no `p_cpies` salinity cutoffs; could use Step E's full 516-level 0-5150dbar grid instead of trimming to West's 5000dbar salinity limit).

**Result: both `moc_overturning` and `moc_gyre` came out identically zero** (~1e-16 Sv, floating-point noise) — `moc_mean` alone equals `moc_total_direct` exactly (both 6.285 Sv). Re-deriving the algebra without a tracer explains why, and it's a **general mathematical fact, not specific to this array/dataset**:
- `moc_gyre`'s integrand is `nansum_gap(vel_pp.*dxA)` where `vel_pp=vel00-V0` and `V0` is *by definition* the `dxA`-weighted mean of `vel00` across gaps — so this is exactly the weighted sum of deviations from a weighted mean, which is always zero.
- `moc_overturning`'s integrand is `vel_prime.*dx2` where `vel_prime=V0-V00` and `V00` is *by definition* the `dx2`-weighted depth-mean of `V0` — same identity, same reason, zero for the same structural cause.
- Physically: `mean+overturning+gyre` only produces nonzero `overturning`/`gyre` terms for freshwater/heat because those are *correlations between two independently-varying fields* (velocity anomaly × tracer anomaly). Decomposing a single field against its own weighted mean, with no second field to correlate against, is mathematically guaranteed to zero out everything past the mean — true for any array, any dataset, any time period. This isn't a bug and isn't fixable by better data.

**Conclusion**: the classical overturning-*streamfunction* definition (cumulative integral of `vel00(x,z)` from the surface down to each depth, then take the max over depth — the RAPID-26°N-style "MOC strength") is the only one of the originally-offered options that's physically meaningful for pure volume transport. Deferred at the user's choice (2026-08-14) in favor of MHT phase 3 first; `moc_samba_marion_v1.m` and its outputs were deleted (not committed) since the deliverable it computed isn't useful — this section preserves the reasoning so the same dead end isn't rediscovered from scratch.

## MOC streamfunction: `moc_streamfunction_v1.m` (2026-08-14) -- implemented after MHT phase 3, several real fixes found along the way

`Psi(z,t) = cumsum_{z'=0}^{z} [ vel_prime(z',t) * dx2(z',t) ] dz'`, `MOC(t) = max_z Psi(z,t)` (Sv) -- the classical definition concluded above. No salinity/temperature needed at all (same structural simplification noted for the abandoned `moc_samba_marion_v1.m`); loads Step E `v2` directly (`Absolute_TPUD*_preTopo` + `AREA_TOPO` + `Ekman_TPUD`, same `v12`-style weighting as `mov_samba_marion`/`mht_samba_marion_v2`).

**Three real issues found and fixed while validating, each caught by comparing against expectations rather than trusting the first result:**

1. **Used `V0` (full zonal-mean velocity, including the barotropic/net component) instead of `vel_prime=V0-V00` for the streamfunction integrand.** First attempt gave `mean(Psi at max depth)=6.285` Sv exactly matching `mean_term`/`moc_mean` — looked like a good sanity check, but meant every individual day's profile sat on top of that day's own barotropic offset (`mean_term`'s std is ~40 Sv vs. a ~6 Sv mean, from the 2017-2019 coverage-artifact investigation), swamping the baroclinic structure. A 45-day lowpass barely helped (barotropic noise is too energetic to filter out after the fact). Fixed by integrating `vel_prime` instead of `V0` -- since `V00` is *by definition* the `dx2`-weighted depth-mean of `V0`, `sum_z(vel_prime.*dx2)=0` exactly, so `Psi` now closes to ~0 at the seafloor *by construction*, cleanly isolating the overturning signal. (Same cross-term-cancellation identity that made the abandoned `moc_samba_marion_v1.m`'s overturning term trivially zero when applied to the *whole* profile -- applied here correctly, to just the *bottom boundary condition* of a cumulative sum, which is exactly the right place for it.)
2. **Depth axis was inverted on the Hovmöller plot** (`YDir` set to `'normal'` instead of `'reverse'`) -- since `Pressure` runs `0->5150` ascending, `'normal'` put the surface at the *bottom* of the figure and the seafloor at the *top*, backwards from both the standard oceanographic convention and how a Kersalé et al. (2021)-style figure would be drawn. Caught while building the requested mean-profile-vs-depth plot, at the user's prompt.
3. **Sign convention**: this project's established `-V00*...`-style leading minus (used throughout for `mean_term`/`mov`/`Q_mean`) gives the *wrong* sign here -- with it, the time-mean profile is negative through most of the water column (peaking around -20 Sv near 1000-1500dbar); *without* it, the profile matches the standard/expected MOC shape (~0 at the surface, growing to +20 Sv around 1000-1500dbar -- the upper cell -- tapering through the NADW layer, small negative dip near the bottom hinting at an abyssal cell). Confirmed by literally plotting both signs side by side and checking against the user's physical expectation ("a camada superior deveria dar positiva") rather than guessing. Adopted the no-leading-minus convention.

**Result**: `MOC` (45-day lowpass) mean=**20.9 Sv**, std=11.3 Sv, typical depth-of-maximum ~1330-1540dbar -- comparable in magnitude to RAPID's ~17 Sv at 26°N, a reassuring cross-check given the very different array/method. Same graded coverage shading available (`n_gaps_valid` computed, not yet plotted on these particular figures). Outputs: `moc_streamfunction_hovmoller.png` (Ψ(z,t), depth-time), `moc_streamfunction_mean_profile.png` (time-mean Ψ(z) ± 1 std, the Kersalé-Figure-2-style plot), `moc_index_IES.png` (MOC(t) raw vs. lowpass), `moc_streamfunction_v1.mat`.

## MOC "Pilot" method: `moc_pilot_v1.m` (2026-08-15) -- two-endpoint (A-P1) dynamic-height-difference estimate, for comparison against the full 8-gap streamfunction

Requested to reproduce the red "Pilot" curve in Kersalé et al. (2021) Figure 2a, which estimates the MOC from the dynamic-height difference between just the two basin-edge sites instead of integrating all 8 gaps. User's explicit choice on a 3-way clarifying question: use BPR (bottom pressure) absolute referencing at **both** endpoints (site A, West edge; site P1, East edge) -- not baroclinic-only shear with a level-of-no-motion assumption.

**P1 confirmed as the true East edge of the array**: `position_East.m` gives `cpiesE(1).lon=17.5576`°E for P1, the largest (easternmost) longitude of all 8 East sites (P2=17.30°, P4=15.00°, P5=11.20°, P6=7.45°, P8≈0°) -- and `cpiesW(1).lon=-51.5`°W for A is the most-negative (westernmost) of the West sites. So A-to-P1 does span the full basin in one hop, as intended.

**Not built from Step E v2's saved output directly** -- that `.mat` only has the final per-gap `Total_TPUD`/`Absolute_TPUD_preTopo` (already resolved to 8 gaps), not the intermediate `Gpan`/`pres`/`rhob` needed to construct a *new*, single A-to-P1 "gap". Reloads the same raw sources Step E v2 does (`samba_w.mat` for A; `ies/Wrk/IES_FrSA_SAMBA_6PIES.mat` + `ies_profiles/Wrk/IES_Make_Profiles_FrSA_SAMBA_6PIES.mat` for P1; `Ekman_transports_9pies.mat`'s `ekmanN_AtoZ`, the whole-basin Ekman total already computed alongside the per-gap `ekmanN_i`) and reruns just the A/P1 slice of Step E v2's machinery.

**Formula: mirrors gap 1 (A-C) exactly, substituting P1 for C** -- rather than deriving from scratch or guessing which of the 8 gaps' (non-uniform, direction-dependent) sign/ordering conventions to copy (checked: e.g. gap1 and gap3's `RefTPUD` formulas have opposite term order/sign despite looking structurally similar, and it's not a simple "shallow vs. deep" or "west vs. east" rule -- looks hand-derived per gap in the original code, not from one universal pattern). This is a safe analogy because **P1 is actually the *shallowest* East site** (`p_cpiesE(1)=1280dbar`, shallower even than A's 1370dbar) despite being the far-East edge -- so A and C's shallow/deep relationship in gap 1 (A=1370 shallow, C=4620 deep) maps directly onto A and P1 here, with P1 now playing the "shallow reference-index" role and A playing the "deep, own-profile-for-correction" role:
- `Gpan_AtoP1_corr = Gpan_A(p_ind_A,:) - Gpan_A(p_ind_P1,:)` (A's own profile at both bottom depths, demeaned)
- `pres_AtoP1 = pres_A - Gpan_AtoP1_corr.*rhob_A` (project A's bottom pressure up to P1's shallower depth)
- `RefTPUD_pilot = (pres_AtoP1 - pres_P1) / (f*rhob_P1)`, `RelTPUD_pilot = -(1/f)*(Gpan_P1-Gpan_A)` (same east-minus-west + sign-flip convention as all 8 gaps)
- `Offset_pilot = RefTPUD_pilot - RelTPUD_pilot(p_ind_P1,:)` (applied at P1's index, the shallower site)

**No `AREA_TOPO` correction** -- `topo_corr_msm60.mat` only has per-adjacent-gap topographic width, nothing for a single A-to-P1 span. Uses the plain `gsw_distance` between A's and P1's lon/lat (6197 km). A genuine simplification vs. the full method, consistent with "Pilot" being the simpler of the two methods being compared, not an oversight.

**Streamfunction**: same construction as `moc_streamfunction_v1.m` (remove the depth-mean/barotropic component so `Psi` closes to ~0 at depth, cumsum from the surface down, 45-day lowpass, `MOC(t)=max_z Psi(z,t)`) -- since there's only one "gap" here, the depth-mean removal is a plain `nanmean` over depth of the A-P1 velocity (no `AREA_TOPO`-weighted `dx2`, width is depth-constant). Verified `mean(Psi_pilot at max depth)=-2e-16` Sv (closes as expected). Same no-leading-minus sign convention as the full streamfunction -- confirmed correct by checking the resulting mean profile independently: rises smoothly from ~0 at the surface to a peak ~1400-1600dbar, back to ~0 at the bottom, same shape as the full array's already-validated profile.

**Result**: Pilot MOC (45-day lowpass) mean=**29.8 Sv**, std=7.0 Sv, peak of the time-mean profile at 1480dbar (vs. full array's 20.1 Sv at 1320dbar) -- same overall shape and depth-of-maximum as the full array (see `moc_profile_comparison_IES.png`, a Kersalé-Figure-2-style overlay of both time-mean profiles ±1 std), but a materially higher/smoother index, and a modest correlation (0.24) between the two daily indices over the full record.

**The intensity gap is mostly explained by the 2017-2019 coverage artifact, not a real methodological difference**: comparing means with/without the already-documented low-coverage window (2017-08-01 to 2020-01-01, when P5/P6/P4/P8 were down -- see the Mov coverage-artifact section above) --

| Period | Full array | Pilot | Diff |
|---|---|---|---|
| Full record (2013-09/2022-12) | 20.9 Sv | 29.8 Sv | 8.8 Sv |
| Excluding 2017-08/2020-01 | 24.7 Sv | 28.6 Sv | **4.0 Sv** |
| Only 2017-08/2020-01 | 10.3 Sv | 32.9 Sv | 22.6 Sv |

Since the Pilot only depends on A and P1 (neither of which had an outage in that window), it stays near its normal level throughout, while the full array's index is pulled down hard by losing 4 of 8 gaps -- independently confirming that dip is instrumental, not oceanographic (the Pilot shows no corresponding dip on `moc_pilot_vs_full_IES.png`). Excluding that window closes more than half the gap (8.8→4.0 Sv); the remaining ~4 Sv likely reflects the Pilot's real simplifications (no `AREA_TOPO`, only 2 points vs. 8 spatially-resolved gaps).

Outputs: `moc_pilot_vs_full_IES.png` (MOC(t) overlay, both methods), `moc_profile_comparison_IES.png` (time-mean Ψ(z) ±1 std overlay, from the separate `moc_profile_comparison_v1.m` which just reloads both `.mat` outputs -- no recomputation), `moc_pilot_v1.mat`.

## MOC Pilot, level-of-no-motion (LNM) variant: `moc_pilot_lnm_v1.m` (2026-08-15) -- confirms the streamfunction MOC index is mathematically insensitive to the BPR-vs-LNM reference choice

Revisited the two other options from the original 3-way clarifying question before `moc_pilot_v1.m` was built ("baroclinic-only shear with a level-of-no-motion assumption" vs. the chosen BPR-absolute option) -- user asked whether it was worth investigating the LNM alternative for comparison.

`moc_pilot_lnm_v1.m` shares `moc_pilot_v1.m`'s `RelTPUD_pilot` derivation verbatim (same `Gpan_A`/`Gpan_P1`, same dt-alignment, same Ekman/distance) but needs no `pres_A`/`pres_P1`/`rhob_A`/`rhob_P1` at all -- instead of BPR-referencing, it assumes zero velocity at the deepest pressure level where neither `Gpan_A` nor `Gpan_P1` is entirely NaN across the record (found programmatically via `~all(isnan(...),2)`, lands at 5000dbar in practice -- `samba_w.mat`'s real A profile is capped there, see Step E v2's header note).

**Result: the two reference choices give numerically IDENTICAL streamfunction/MOC output** -- `corr(BPR,LNM)=1.000`, means equal to 4 decimal places (29.7565 Sv both), curves fully overlapping on `moc_pilot_lnm_vs_bpr_IES.png` (the BPR curve is plotted first and completely hidden under the LNM curve). **This is not a data coincidence, it's a structural consequence of the streamfunction definition itself**: both the BPR offset (`RefTPUD_pilot - RelTPUD_pilot(p_ind_P1,:)`) and the LNM offset (`-RelTPUD_pilot(ref_idx,:)`) are pure per-day CONSTANTS -- a single scalar added uniformly across every depth level, shifting the whole profile up/down without changing its shape. But the streamfunction construction (shared with `moc_streamfunction_v1.m`/`moc_pilot_v1.m`) explicitly removes the depth-mean/barotropic component (`vel_prime=vel-V00`) *before* integrating -- and that step cancels out any uniform vertical shift exactly, regardless of its value. So `Psi`/`MOC(t)` depend only on the *shape* of the baroclinic shear, never on which absolute-reference method anchored it. This matches known theory behind why RAPID/MOVE-style streamfunction MOC indices are relatively insensitive to reference-level uncertainty -- now empirically confirmed for this array/pipeline specifically, not just asserted from the literature.

**Practical implication**: choosing BPR over LNM for the Pilot method's *MOC-strength index* specifically turned out not to matter numerically -- but BPR remains the physically correct choice if the *absolute* velocity/transport at the reference level itself is ever needed (e.g. reporting the barotropic/net transport, not just the overturning cell shape), since LNM's "zero velocity assumed" is only ever an assumption, not a measurement, for that quantity.

Outputs: `moc_pilot_lnm_vs_bpr_IES.png` (3-way MOC(t) overlay: full array, Pilot-BPR, Pilot-LNM), `moc_pilot_lnm_v1.mat`.

## Root cause of the full-array-vs-Pilot MOC intensity gap (investigated 2026-08-15, no script changes -- diagnostics only, not committed)

User's premise: for pure volume transport (MOC only -- explicitly does NOT apply to Mov/MHT, which involve tracer correlations, not just geopotential differences), the relative/baroclinic transport across the whole section should depend *only* on the geopotential-anomaly difference between the two edge sites (A, P1), because the intermediate gaps' `Gpan` differences telescope: `(Gpan_C-Gpan_A)+(Gpan_D-Gpan_C)+...+(Gpan_P1-Gpan_P2) = Gpan_P1-Gpan_A` exactly, by algebra, regardless of what happens in between -- this is the literal reason the SAMBA array design and Kersalé et al. (2021) Figure 2 expect the "Pilot" and full-array curves to track closely. Asked to find why our two implementations diverge by as much as they do (~8.8 Sv mean, `moc_profile_comparison_IES.png`'s peak: 20.1 vs 29.7 Sv).

**Investigated in three independent, isolated diagnostics** (recomputed everything from raw sources in scratch scripts, not from the two committed scripts' saved output, to rule out a shared bug):

1. **Site A/P1 coverage**: confirmed near-total (A valid 3379/3382 common days = 99.9%, P1 valid 3382/3382 = 100%). Not the cause.
2. **The telescoping identity holds exactly, as expected**: recomputed `Gpan` for all 9 sites from scratch and summed the 8 gaps' `RelTPUD_i` -- matches `RelTPUD_pilot=-(1/f)*(Gpan_P1-Gpan_A)` to `max|diff|=2.9e-11` (floating-point noise) at *every* depth and every day. Confirms the user's "beauty of the method" premise is exactly right for the baroclinic part -- this is *not* where the discrepancy comes from.
3. **The BPR reference/offset part diverges, but as a pure depth-CONSTANT** -- `sum_i(Offset_i) - Offset_pilot = 1505.3` m²/s, identical at every pressure level checked (0 to 4800dbar). Each gap's independently-calibrated BPR sensor doesn't have to sum to the same absolute reference a direct A-P1 BPR comparison gives (no telescoping guarantee for measured, not algebraically-derived, quantities) -- but because it's depth-uniform, it's exactly the kind of quantity the streamfunction's depth-mean removal cancels (same mechanism as the BPR≡LNM finding above). **Also not the cause of the MOC/Ψ discrepancy**, even though it's a real, nonzero difference in the raw (non-streamfunction) transport.
4. **`AREA_TOPO` is the actual cause.** Each of the 8 gaps' `Total_TPUD_i` is multiplied by that gap's own depth-dependent topographic-width fraction (`topo_corr_msm60.mat`, real per-gap bathymetric blocking near the seafloor) -- the Pilot method has no such correction (no basin-spanning topographic profile exists for a single A-P1 span; documented simplification in the Pilot section above). With `AREA_TOPO` included, `sum_i(Total_TPUD_i) - Total_TPUD_pilot` stops being a depth-constant and grows sharply with depth: ~1500 m²/s near the surface, ~1430-1670 through 1200-2000dbar, up to **10,900 m²/s by ~4400dbar**. Because the streamfunction removes only *one* overall depth-mean (not a depth-local correction), this deep-ocean distortion shifts the *entire* profile baseline -- including the shallow (~1300-1500dbar) peak that dominates the MOC index -- which is what actually produces the intensity gap.

**Conclusion: no bug found.** The ~8.8 Sv full-record gap has two distinct, both-real contributors, and neither is a computational error: (a) the already-documented 2017-2019 coverage-outage window (explains about half, per the earlier full-coverage-only-days test: gap shrinks from 8.8 to ~4.0 Sv when that window or any incomplete-coverage day is excluded), and (b) `AREA_TOPO`'s absence from the Pilot method, which explains the remaining ~4 Sv persisting even on full-coverage days. The baroclinic (relative) part -- the actual "beauty of the method" the user described -- checks out exactly as expected; the residual is entirely attributable to two known, already-documented methodological asymmetries between the two implementations, not to an undiscovered bug.

## MOC Pilot v2: `moc_pilot_v2.m` (2026-08-15) -- adds an aggregate AREA_TOPO correction, closes most of the structural residual

Follow-up to the root-cause investigation above: user asked "nao tem um jeito de considerarmos essa AREA_TOPO para o Pilot?" -- since no independent bathymetry exists for a single A-to-P1 span, built a data-driven proxy instead of acquiring new data:

```
AREA_TOPO_pilot(z) = sum_i[ dx_i * AREA_TOPO_i(z) ] / sum_i[ dx_i ]
```

a width-weighted average of the 8 real, **static** (depth-only, no time dimension) per-gap `AREA_TOPO_i(z)` profiles from Step E v2, weighted by each gap's own along-section distance `dx_i`. Deliberately built from the static profiles/weights only, *not* from the full array's own day-varying `dx2(z,t)` -- using `dx2` directly would have silently imported the full array's own coverage outages into the Pilot, defeating the purpose of an independent cross-check.

**Side finding while building this**: `dx_total=sum(dx_i)` (the 8-gap along-section path, 7846km) is 26.6% *longer* than `v1`'s `dx_pilot` (direct A-P1 great-circle, 6197km) -- not a bug, a real geometric fact: two points on the same latitude (except the equator) are joined by a great circle that bows toward the pole, shorter than the constant-latitude arc the array physically follows. `v1`'s `dx_pilot` never affected `v1`'s `Psi`/MOC output (it's used to convert TPUD→velocity and immediately multiplied back by the identical value, canceling exactly -- same "any depth-constant factor cancels via the depth-mean removal" mechanism as the BPR≡LNM finding). `v2` switches to `dx_total=sum(dx_i)` throughout since `AREA_TOPO_pilot`'s weighted average needs the physically-correct along-section width.

**Result**: closes most of the "structural" (coverage-independent) part of the gap. Restricting to full-coverage (8/8 gaps) days specifically, where the earlier investigation found a persistent ~4.0 Sv residual unexplained by coverage: with `AREA_TOPO_pilot` applied, that residual drops to **0.82 Sv** (~80% closed) -- full array 24.96 Sv vs. Pilot v2 25.78 Sv on those days, vs. 24.96 vs. 28.98 for `v1`. Over the full record (all coverage levels), the gap shrinks from 8.82 Sv (`v1`) to 5.60 Sv (`v2`): full array 20.93 Sv, Pilot v1 29.76 Sv, Pilot v2 26.54 Sv. Correlation with the full array's daily index is essentially unchanged (0.234 vs. 0.240) -- expected, since `AREA_TOPO_pilot` is a static depth correction, it shifts the overall level/shape but not day-to-day variability.

**Coverage-threshold test (more precise than a fixed date window)**, using `n_gaps_valid` directly rather than the earlier "exclude 2017-08/2020-01" date range (which misses the second, smaller 2020-2022 coverage patch already flagged in the Mov section above):

| Filter | Days | Full array | Pilot v1 | Pilot v2 |
|---|---|---|---|---|
| All days | 3382 | 20.9 Sv | 29.8 Sv (diff 8.8) | 26.5 Sv (diff 5.6) |
| n_gaps≥5 | 2642 | 24.4 Sv | 28.6 Sv (diff 4.2) | 25.5 Sv (diff 1.1) |
| n_gaps=8 (perfect) | 1836 | 25.0 Sv | 29.0 Sv (diff 4.0) | 25.8 Sv (diff **0.8**) |

Confirms both real effects cleanly separated: `AREA_TOPO` alone (not coverage) explains ~80% of the gap that remains even at perfect coverage; coverage (not `AREA_TOPO`) explains the sharp additional widening once `n_gaps_valid` drops below ~5 (the remaining 5.6 Sv full-record gap is now overwhelmingly a coverage-outage effect, since it's only 0.8 Sv on good-coverage days).

Outputs: `moc_pilot_v2_vs_full_IES.png` (3-way MOC(t): full array, Pilot v1 no-topo, Pilot v2 with-topo), `moc_pilot_v2_profile_comparison_IES.png` (Kersalé-style time-mean profile overlay, full array vs. Pilot v2), `moc_pilot_v2.mat`.

## MOC redefinition: `moc_streamfunction_v2.m` + `moc_pilot_v3.m` (2026-08-15) -- v1/v2 used the WRONG MOC definition; corrected against the actual paper text

User asked to verify against the paper directly: "No artigo da Kersale, a MOC nao e a integral total na bacia, sem decomposicao em media e baroclinica?" Read `kersale_jgr2021_samba.pdf` (found at `~/Documents/references/moc/`) directly rather than relying on general RAPID-26N-style conventions, via `pdftotext -layout` (needed `brew install poppler` first, not previously on this machine).

**Confirmed the user was right -- `moc_streamfunction_v1.m`/`moc_pilot_v1.m`/`v2.m` implemented a materially different, incorrect definition.** Direct quotes from Section 2.1.1/2.1.4:
- *"Because the present study does not use a residual method, and obtains the reference (barotropic) velocity variability directly from data (bottom pressure gradients), no 'residual' zero net volume transport assumption is made here. In other words, the velocities are not uniformly adjusted to ensure zero net volume transport..."* -- directly contradicts `v1`'s `vel_prime=V0-V00` step (a RAPID-26°N-style residual/zero-net construction), which the paper explicitly disavows.
- *"the strength of the MOCup flow is defined as the basin-wide transport integrated from the surface down to the time-mean pressure interface where the zonally integrated meridional flow changes from northward to southward (1,315 dbar here at 34.5°S for both [pilot and full array] configurations)"* -- a FIXED depth found ONCE from the long-term time-mean profile, not `v1`'s per-day `max_z(Psi(z,t))` search.

**Corrected construction** (`moc_streamfunction_v2.m` for the full array, `moc_pilot_v3.m` for the Pilot -- both keep everything else from `v1`/`v2` unchanged: same `vel00`/`geo`/`dxA`/`dx2`/`V0` zonal-averaging, same BPR reference derivation, same aggregate `AREA_TOPO_pilot(z)` from `v2`):
1. `Psi_raw(z,t) = cumsum_z[ V0(z,t)*dx2(z,t) ]` -- RAW absolute-velocity cumulative transport, NO depth-mean removal. `V0`/`Absolute_TPUD*_preTopo` already contain the real BPR-derived barotropic component AND the real ECCO-based time-mean reference (Step E's "Add the time-mean velocity from OFES" step, lines ~549-693) -- exactly what the paper describes combining, no new data needed.
2. `h_star = argmax_z[ mean_t(Psi_raw(z,t)) ]` -- a single FIXED depth from the time-mean profile (found independently for the full array and for the Pilot, not forced to match).
3. `MOCup(t) = Psi_raw(h_star, t)` for every day, 45-day lowpass applied to this one time series (not per-depth-level anymore, since only one depth is needed).

**Validates cleanly against the paper's own numbers**: full array `h_star=1220dbar` (paper: 1315dbar, close given the different/extended period), mean **MOCup=18.02 Sv** (paper: 17.3 Sv full array, 17.7 Sv pilot -- both essentially matched). Time-mean profile shape (`moc_streamfunction_v2_mean_profile.png`) reproduces the paper's Figure 2a well: ~0 at surface, peak ~18 Sv at `h_star`, decreasing through NADW to a local minimum around -7 Sv near 3800-4000dbar (matches the paper's separately-reported "6.6±2.7 Sv southward, integrated to 4700dbar" reasonably), hint of AABW reversal near the very bottom.

**Bonus retroactive validation**: the paper's 6.6±2.7 Sv deep/abyssal-transport number is essentially identical to the ABANDONED `moc_samba_marion_v1.m`'s `moc_mean=6.285` Sv (deleted 2026-08-14 because its `moc_overturning`/`moc_gyre` terms were trivially zero -- see the MOC dead-end section above). That `moc_mean` term itself was a correct, useful, independently-validated quantity all along; only the "decompose into mean+overturning+gyre" framing was inapplicable to pure volume transport, not the underlying calculation.

**Pilot v3 result**: `h_star_pilot=1270dbar` (close to the full array's 1220dbar, mirroring the paper's finding that both configurations land near the same transition depth), mean **MOCup=23.97 Sv**. Full-record gap vs. full array: 5.95 Sv (worse in absolute terms than `v2`'s old-definition 5.60 Sv, but not directly comparable -- different quantities). Coverage-threshold breakdown (`n_gaps_valid`, same diagnostic pattern as before) tells a different story than under the old definition:

| Filter | Days | Full array | Pilot | Diff |
|---|---|---|---|---|
| All days | 3382 | 18.0 Sv | 24.0 Sv | 5.95 |
| n_gaps≥5 | 2642 | 23.2 Sv | 22.7 Sv | **-0.41** |
| n_gaps=8 (perfect) | 1836 | 25.5 Sv | 23.0 Sv | -2.45 |

Correlation: 0.130 (all days) → **0.533** (n_gaps=8 only) -- a large improvement, though still well short of the near-1 agreement the paper's own pilot/full comparison implies. **Without the depth-mean removal, the 2017-2020 coverage outage (losing gap 3, D-P8, 71.3% of the total section width -- see the mechanism investigation above) now distorts the full array's RAW signal much more severely** (daily excursions down to -50 Sv, vs. only dipping toward 0 Sv under `v1`'s damped/residual definition) -- since there's no longer any depth-mean-removal step to partially compensate for the missing width. This means the full-record MEAN is dominated by these extreme excursions specific to the bad-coverage window; on good-coverage days the two methods are much more comparable (full array actually running slightly *higher* than Pilot, not lower). **Investigation continuing** -- 0.533 correlation on good days is a big improvement but not yet the near-identical agreement the paper reports, and not yet root-caused further.

**Reference-choice caveat for future work**: `moc_pilot_lnm_v1.m`'s finding that BPR and LNM references give identical results was specific to the abandoned depth-mean-removal definition (a per-day constant offset cancels exactly under that construction). Under this corrected (no depth-mean-removal) definition, a constant reference offset instead accumulates *linearly* with depth through the cumsum -- so BPR vs. LNM would generally give *different* `Psi_raw`/MOCup now. Not re-tested with `v3`; BPR remains the physically-preferred, user-chosen method regardless.

Outputs: `moc_streamfunction_v2_hovmoller.png`, `moc_streamfunction_v2_mean_profile.png`, `moc_streamfunction_v2_index_IES.png`, `moc_streamfunction_v2.mat`; `moc_pilot_v3_vs_full_IES.png`, `moc_pilot_v3_profile_comparison_IES.png`, `moc_pilot_v3.mat`. `v1`/`v2` (streamfunction and Pilot) are NOT deleted -- kept as a documented, legitimate-in-its-own-right RAPID-26°N-style "overturning strength" alternative, just not what Kersalé's Figure 2/5 "MOC" refers to.

## MOC shelf-transport fix: `moc_streamfunction_v3.m` + `moc_pilot_v4.m` (2026-08-15) -- found why the corrected full array (18.0 Sv) still ran high vs. the paper (17.3 Sv)

Follow-up: even with `v2`/`v3`'s corrected Kersalé definition, the full array's own MOCup over the paper's exact 2013-2017 window (25.75-25.88 Sv) ran ~8.5 Sv above the paper's reported 17.3 Sv -- too large a gap to be definitional. Root-caused by testing the corrected definition directly against **Marion's original, completely unmodified 2020-vintage Step E output** (`Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT.mat`, 1405 days, 2013-09-11 to 2017-07-17 -- predates every West/East source extension made in this session): `h_star=1310dbar` (paper: 1315dbar, near-exact) but `MOCup=22.86` Sv -- still 5.6 Sv above the paper, even with **zero reprocessing differences**. So roughly 3 Sv of `v2`'s gap came from this session's West/East source changes, but a larger ~5.6 Sv chunk was already present in Marion's own original pipeline output.

**Found the cause**: `Full_Depth_OverturningEstimate_Cst_for_MHT.m` (both the original and this session's extensions) loads `ecco/Wrk/ECCO_TPUD_ShelfBits.mat` (`West_TPUD`, `East_TPUD_new` -- static, time-invariant depth profiles, the paper's Section 2.1.3 "time-invariant estimate of the volume fluxes on the continental shelves/upper slopes inshore of the shallowest moorings," from ECCO2) but **never adds it to `Total_TPUD1..8`** -- loaded, pressure-trimmed alongside everything else, then silently unused for the rest of the MOC calculation. (It IS used downstream, in `mht_samba_marion_v2.m`'s `Q_shelfW`/`Q_shelfE` terms for heat transport -- this omission was specific to the volume/MOC side, which nobody had needed to compute in the paper's exact "basin-wide, coast-to-coast" convention until now.)

`West_TPUD` integrates to **-4.51 Sv** (a real, sizable southward shelf contribution); `East_TPUD_new` is negligible (+0.01 Sv). Adding both to Marion's original 2020 output before integrating: MOCup drops from 22.86 to **18.44 Sv** -- within ~1 Sv of the paper's 17.3 Sv, confirming this was the dominant missing piece. `h_star` is unaffected (1310dbar either way).

**Applied to both `moc_streamfunction_v3.m` (full array) and `moc_pilot_v4.m` (Pilot)** -- the Pilot needs the identical correction for a fair comparison (both methods exclude the same coastal shelf regions inshore of the A/P1 moorings in the same way; `West_TPUD`/`East_TPUD_new` are static so this is one depth profile added identically to every day in both scripts). Everything else unchanged from `v2`/`v3`.

**Result, over the paper's exact 2013-2017 window** (previously: 25.75/23.52 Sv, corr=0.502):

| | This work | Paper |
|---|---|---|
| Full array | 21.26 Sv | 17.3 Sv |
| Pilot | 19.07 Sv | 17.7 Sv |
| Correlation | 0.503 | 0.73 |

A substantial improvement (gap to the paper's full-array value shrunk from ~8.5 Sv to ~4 Sv) but not a perfect match. The residual ~2-4 Sv offset and the correlation gap (0.50 vs. 0.73) are plausibly attributable to real version/provenance differences already documented elsewhere in this repo -- the East array here is rebuilt from the published Daily_Tau SAMOC-SA dataset (not Marion's exact original calibration used in the 2021 paper), and GEM lookup tables may have been regenerated with updated CTD/Argo training data since 2021. **Investigation continuing** at the user's request -- this residual has not yet been further decomposed.

Coverage-threshold breakdown (full extended 2013-2022 record, same `n_gaps_valid` diagnostic as before) after the shelf fix: gap shrinks from 6.05 Sv (all days) to essentially flipping sign at `n_gaps≥5` (full=18.66, pilot=18.30, diff=-0.36), settling around -2.4 Sv at full coverage (full=21.00, pilot=18.60, corr=0.533) -- same qualitative pattern as `v2`/`v3` (coverage artifact dominates the full-record mean; good-coverage days show much better agreement), just with both curves now shifted down by the shelf correction.

Outputs: `moc_streamfunction_v3_hovmoller.png`, `moc_streamfunction_v3_index_IES.png`, `moc_streamfunction_v3_mean_profile.png`, `moc_streamfunction_v3.mat`; `moc_pilot_v4_vs_full_IES.png`, `moc_pilot_v4_profile_comparison_IES.png`, `moc_pilot_v4.mat`. (Superseded 2026-08-17 by `v4`(streamfunction)/`v5`(Pilot) below -- `v3`/`v4` still loaded Step E `v2`'s un-recalibrated output.)

### Decomposing the remaining ~4 Sv residual (2026-08-15, diagnostics only, no script changes)

After the shelf fix, the full array's MOCup over the paper's own 2013-2017 window (21.26 Sv) still ran ~4 Sv above the published 17.3 Sv. Split the gap into two independently-testable pieces:

1. **Marion's original pipeline vs. the published number**: her unmodified 2020 Step E output, with the shelf correction applied, gives 18.44 Sv vs. the paper's 17.3 Sv -- a **1.1 Sv** gap that involves none of this session's reprocessing at all. Likely irreducible version drift between her 2020 processing and whatever exact inputs went into the final 2021 publication (GEM table vintage, minor calibration revisions, etc.) -- not further investigated, diminishing returns expected.
2. **This session's reprocessing vs. Marion's original pipeline, identical 2013-2017 period**: 21.26 vs. 18.44 Sv -- a **2.82-2.89 Sv** gap (confirmed via two independent checks: whole-record MOCup difference, and a direct per-gap upper-layer (0-1300dbar) transport comparison, which totaled +2.89 Sv, consistent to within rounding). Broken down gap-by-gap (time-mean, 0-1300dbar, 1405 common days):

| Gap | Original (Sv) | Reprocessed (Sv) | Diff |
|---|---|---|---|
| A-C (West) | -4.50 | -2.82 | **+1.69** |
| C-D (West) | -2.09 | -2.44 | -0.35 |
| D-P8 (interior) | 14.31 | 15.04 | +0.73 |
| P8-P6 (East) | -6.50 | -7.76 | **-1.27** |
| P6-P5 (East) | -0.40 | 0.15 | +0.55 |
| P5-P4 (East) | 6.81 | 6.40 | -0.41 |
| P4-P2 (East) | 14.13 | 14.53 | +0.41 |
| P2-P1 (East) | 0.44 | 1.98 | **+1.54** |

**No single gap dominates** -- the three largest contributors (A-C/West +1.69, P2-P1/East +1.54, P8-P6/East -1.27) span both sides of the array, not concentrated in either the West source swap (`samba_w.mat`) or the East rebuild (Daily_Tau) alone. Consistent with genuine, expected calibration/provenance noise from switching underlying T/S/pressure sources -- already documented elsewhere in this repo that the East rebuild agrees with Marion's original East calibration only to ~0.0024s (tau1000) / 0.16dbar (pressure), a small per-measurement tolerance that compounds into Sv-scale differences once integrated across a whole gap and depth range. **Conclusion (superseded below, 2026-08-17): this "no dominant gap" reading turned out to be an artifact of an uncontrolled confound, not the real picture** -- see the next section. Session paused here (2026-08-15) at the user's request to think it over.

## West-vs-East isolation done properly (2026-08-17): the ECCO calibration-window confound, and a real GEM-vintage difference concentrated at site A/C -- the Brazil Current

Resumed with a direct test: rerun the West-vs-East decomposition using **hybrid** Step E runs (one side's source swapped, the other held at Marion's original) instead of inferring it indirectly from the full reprocessed output.

**Setup**: Marion's original West files (`ies/Wrk/IES_USA_SAM.mat`, `ies_profiles/Wrk/IES_Make_Profiles_USA_SAM.mat`) are still untouched on disk (dated 2020/2021, predate this whole session) -- but the original East files at `ies/Wrk/IES_FrSA_SAMBA_6PIES.mat`/`ies_profiles/Wrk/IES_Make_Profiles_FrSA_SAMBA_6PIES.mat` were overwritten in place by this session's Daily_Tau rebuild, with no backup (`.mat` files are gitignored, no git history either). **Before touching anything**, backed up the true-original Step E output to `Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_TRUEORIGINAL_backup.mat` (irreplaceable if lost). Then created `Full_Depth_OverturningEstimate_Cst_for_MHT_hybridB.m` -- an unmodified copy of the original script (only `addpath` + output-filename changed) -- which, run today, automatically uses Marion's original West (untouched) + this session's rebuilt East (since those files were overwritten at the same paths), giving a clean "Hybrid B" without needing her original East data at all.

**First pass showed something impossible**: comparing Hybrid B (orig West + new East) against `mine` (new West + new East) to isolate the West-alone effect, differences appeared even in gaps 4-8 (P8-P6 through P2-P1) which involve *only* East sites -- physically impossible if each gap's transport depends only on its own two endpoints' data. Root cause found: **the ECCO time-mean reference offset step calibrates against *this script's own date range's mean*** (`Offset<i>_mean=RefTPUD<i>_mean-nanmean(Absolute_TPUD<i>(pin,:))` at 1500dbar) -- `mine` (`_v2.m`) covers 2013-2022 and calibrates against that whole period's mean, while Hybrid B/the original script cover only 2013-2017. Verified directly on gap 8 (pure East): predicted shift from this alone (`mean(2013-2017 subset)-mean(full range)` at 1500dbar) = 980.65, closely matching the observed ~900-1000 unit shift at that depth. **Fixed by building `Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3.m`** (in `SAMBA_E_IES`, see that repo's `CLAUDE.md`): identical to `v2` except the offset calibration is restricted to the original 2013-09-11/2017-07-17 window (`cal_idx`) regardless of the script's own date range.

**Redone with matched calibration windows -- the confound vanishes and a clean signal emerges** (time-mean, upper 1300dbar, 1405 common days):

| Gap | Original | HybB (origW+newE) | v3 (matched cal.) | East-alone | West-alone |
|---|---|---|---|---|---|
| A-C (West) | -4.50 | -4.61 | -1.52 | -0.11 | **+3.10** |
| C-D (West) | -2.09 | -2.00 | -2.59 | +0.10 | -0.59 |
| D-P8 (interior) | 14.31 | 14.29 | 13.94 | -0.01 | -0.35 |
| P8-P6 (East) | -6.50 | -6.61 | -6.61 | -0.11 | **0.000** |
| P6-P5 (East) | -0.40 | -0.41 | -0.41 | -0.01 | **0.000** |
| P5-P4 (East) | 6.81 | 6.90 | 6.90 | +0.09 | **0.000** |
| P4-P2 (East) | 14.13 | 14.63 | 14.63 | +0.50 | **0.000** |
| P2-P1 (East) | 0.44 | 0.44 | 0.44 | 0.00 | **0.000** |

East-only gaps now show *exactly* zero West-alone effect (to rounding), confirming the fix. **East-alone total: +0.45 Sv** (small, as expected from the already-documented ~0.0024s/0.16dbar East-rebuild tolerance). **West-alone total: +2.16 Sv, concentrated almost entirely in A-C (+3.10 Sv), partially offset by C-D (-0.59) and D-P8 (-0.35)** -- this time a real, clean, well-localized signal, not noise spread evenly across the array as the earlier (confounded) test suggested.

**Traced to the source**: compared site A's and C's mean T/S profiles directly between `samba_w.mat` (new) and `IES_Make_Profiles_USA_SAM.mat` (old) over the same 1405-day window. Coverage is identical (100% both, not a data-availability issue). The values themselves differ systematically and *oppositely* at the two sites:
- **Site A: `samba_w.mat` is systematically warmer and saltier** than the original source in the upper ~800dbar (e.g. +0.63°C at 520dbar, +0.09 salinity at the surface).
- **Site C: `samba_w.mat` is systematically colder and slightly fresher** in the same layer (e.g. -0.26°C at 390dbar) -- the *opposite* sign from site A.

Because A and C shift in opposite directions, their dynamic-height difference (which drives the A-C gap's transport) doesn't partially cancel the way same-direction errors would -- it compounds, which is why this single gap absorbs most of the +2.16 Sv West-alone effect. This pattern (systematic, opposite-signed, largest at both sites simultaneously) is the signature of the two source pipelines using **different GEM (Gravest Empirical Mode) lookup table vintages** -- `samba_w.mat`'s GEM tables were very likely retrained against a different/more recent CTD+Argo dataset than whatever built `IES_Make_Profiles_USA_SAM.mat` (matches this repo's own note that the paper's GEM tables used CTD/Argo data "from January 1983 to December 2018," while `samba_w.mat` was built later and extends to 2022) -- not a bug in either pipeline, a genuine difference in training data/vintage between two independently-built empirical T/S estimators for the same site.

**Why this specific site matters more than the numbers alone suggest**: the A-C span is where the **Brazil Current** crosses the section -- the western boundary current that closes the subtropical gyre's northward Sverdrup transport (the return flow that balances the gyre's interior northward flow at this latitude). Because a western boundary current concentrates a large fraction of the total meridional transport into a narrow, swift jet, small errors in the T/S profile *at this one site* have an outsized effect on the total transport estimate compared to an equivalent error spread across a wide, slow interior gap -- exactly consistent with A-C alone absorbing +3.10 Sv out of the total +2.16 Sv West-alone effect (partially offset by C-D) while every wide interior/East gap shows a negligible, well-controlled discrepancy. This is a physical reason (not just a statistical one) to treat site A's/C's calibration and GEM-table provenance with more care than the rest of the array when refining this comparison further.

**Status: root cause confirmed, with an important refinement (2026-08-17).** Checked the original West GEM table's own training-data provenance (`gem/west/DATA/*.mat`, all dated Nov 2018 -- matches this repo's earlier note that the paper's GEM tables used CTD/Argo data "from January 1983 to December 2018") against `samba_w.mat` (no build script present in either repo -- provided externally by a collaborator, filename `IES_Make_Profiles_plusBrazil_plusDynHgt_Olga.mat`, dated Dec 2025). User confirmed `samba_w.mat`'s GEM table was retrained with more recent CTD/Argo data (needed to cover deployments through 2022) -- but flagged an important caveat on reflection: **standard practice is to use GEM training data *concomitant* with each PIES measurement's own era, not simply "the newest available."** Applying one GEM table (trained through ~2022) uniformly across the whole 2013-2022 record risks blending later hydrographic information into the interpretation of earlier (2013-2017) τ readings -- a look-ahead inconsistency, not just a beneficial recalibration.

**Tested this directly and confirmed it's a real effect, not just a theoretical concern**: the site A/C temperature discrepancy (`samba_w.mat` minus original, at 500dbar) is **not constant within the 2013-2017 window -- it grows over time**:

| Year | Days | ΔT at A | ΔT at C |
|---|---|---|---|
| 2013 | 112 | +0.41°C | -0.37°C |
| 2014 | 365 | +0.60°C | -0.14°C |
| 2015 | 365 | +0.56°C | -0.13°C |
| 2016 | 366 | +0.60°C | -0.16°C |
| 2017 | 197 | **+1.06°C** | **-0.69°C** |

Linear trend over the 1405-day window: A drifts +0.31°C, C drifts -0.42°C from start to end -- growing larger the closer the date gets to 2018-2022 (the extended training period). This is exactly the signature expected from a non-concomitant-training-data artifact (the retrained GEM table pulling earlier readings toward the extended period's character) rather than a uniform, stable recalibration -- confirms the user's practice-based concern with actual data, not just in principle.

**Practical implication**: the ~4 Sv residual between this repo's MOC estimate and Kersalé et al. (2021)'s published number is **not fully closable by matching methodology alone**, and the West-side (A-C, Brazil Current) contribution to it should NOT be read as "our number is simply more accurate/up to date" -- it's at least partly a genuine artifact of applying a GEM table trained on a wider period than the specific data being converted. This is a real methodological caveat to state explicitly in any write-up comparing this repo's 2013-2017 MOC values against the 2021 paper, and flags a possible concern for the FULL 2013-2022 record too (early-record years being interpreted with hydrographic information that postdates them). Not yet resolved -- would require either a GEM table built/validated separately per sub-period (concomitant with each), or explicitly quantifying/bounding this as a known source of secular drift in the West profiles.

### Does this also affect 2018-2022? Investigated 2026-08-17, inconclusive -- documented as an open limitation

Extended the same site A/C discrepancy check through the *full* old/new source overlap (2009-03-18 to 2018-04-30 -- the original West source, `IES_Make_Profiles_USA_SAM.mat`, actually runs to April 2018, not just to the paper's 2017-07-17 cutoff). Found **two overlaid effects**, not one:
1. **A stable seasonal cycle present from 2009 onward** (dA oscillates roughly 0.35-0.47°C mid-year vs. 0.72-1.08°C in Q1/Q2, every year, including years well before 2018) -- **confirmed below** to be a genuine methodology difference in seasonal correction between the two pipelines, not a training-window artifact.
2. **A secular drift that clearly accelerates from ~2013 onward and is still growing, not leveling off, at the last available comparison point (April 2018)** -- same-quarter comparisons: dA~0.37-0.46°C (2009-2012) → ~0.72-0.87°C (2013-2016) → ~1.03-1.14°C (2017-2018); dC shows an even sharper jump, from small positive values pre-2017 to consistently around -0.60 to -0.65°C by 2017-2018. Critically, this is **still accelerating right at the edge of available data**, giving no reassurance that it plateaus once "inside" the retrained table's own training window.

**Searched this machine thoroughly for independent (non-GEM-derived) CTD/Argo data to check 2018-2022 directly and found none usable**: confirmed the original West Argo training file (`gem/west/DATA/Reformat_Argo_HydroData_forGEM_fullres.mat`) covers exactly 2003-2018 (verified via the embedded year column, not just file dates) in the right region (lon -55 to -40, lat -38 to -32, covering A/C). Files that looked promising by name (`DATA/CTD/ctd_2019.mat`, `ctd_new.mat`) turned out to be **East array** servicing-cruise CTDs (`headerE`/`hydroE` variables), not West/Brazil-side casts. No West-side CTD or Argo collection extending past 2018 exists anywhere in either repo or the broader `~/research/sambar/` tree.

**Decision (user, 2026-08-17): document this as a known, open limitation rather than pursue further right now** (options considered: fetch live Argo profiles from an external database, or ask the collaborator who built `samba_w.mat` about the GEM table's exact training window/validation). **Conclusion to carry forward**: whether the A/C secular drift continues, plateaus, or reverses somewhere within 2018-2022 is **unverified** with data currently available in this project -- treat any MOC/Mov/MHT results for 2018-2022 with the same caution already established for 2013-2017's A-C-driven bias, and revisit if/when independent West-side hydrography for that period becomes available.

### The seasonal component, explained (2026-08-17): `samba_w.mat` very likely includes a seasonal correction the current MOC pipeline's old source lacks

Checked whether `IES_Make_Profiles_USA_SAM.m` (the script the *original*, non-extended Step E pipeline actually uses for the West source) applies any seasonal correction at all: it does not (`grep -i seascorr` on the script returns nothing). The seasonal correction is a **separate product**, only produced by the `_inclSeasGEM` variant script (`IES_Make_Profiles_USA_SAM_inclSeasGEM.m`), which is used elsewhere in this repo for MHT but **not** for the MOC/Step E pipeline. So the plain West source feeding Step E's `Total_TPUD` (and by extension `moc_streamfunction_v3.m`) is NOT seasonally corrected.

Tested directly whether `samba_w.mat` behaves as if it *does* include an equivalent correction: extracted the original pipeline's own seasonal-correction field (`gem/west/Wrk/seasCorr_temper_field_west.mat`, `SC_SMF_temperature` vs. `SCpresrange`/yearday, valid 0-300dbar only, built by `do_SeasCorr_temper_field_west.m` from the same 2003-2018 CTD/Argo training set) at 50dbar, and compared its shape (as a function of day-of-year) against the (`samba_w.mat` minus `IES_Make_Profiles_USA_SAM.mat`) discrepancy's own day-of-year climatology at site A, 50dbar, over the full 2009-2018 overlap:

| Yearday | Original seasonal-correction field | dA climatology |
|---|---|---|
| 1 | +0.30 | +0.40 |
| 61 (peak) | +1.86 | +0.50 |
| 241 (trough) | -1.66 | +0.08 |
| 361 | +0.12 | +0.40 |

**Correlation between the two shapes (both demeaned): r = 0.971** -- a near-perfect phase match, even though the original correction field's amplitude (peak-to-peak 3.64°C) is much larger than the observed discrepancy's (peak-to-peak 0.46°C, likely because the discrepancy also mixes in the unrelated secular-drift component documented above, and/or `samba_w.mat`'s own seasonal correction uses different vintage data or smoothing). This phase match is strong evidence that **`samba_w.mat`'s build applies a seasonal correction (or a very similar one) that the current MOC pipeline's plain West source does not** -- i.e. `samba_w.mat` is very likely the *more complete* product here (seasonal correction is a real, physically-motivated refinement -- the paper's own methodology treats it as standard practice for the upper 300dbar), not a bug or an inconsistency to fix.

**Net read on the two effects together**: the seasonal component of the A/C discrepancy (this section) is benign and expected -- an improvement in `samba_w.mat`, not a concern. The secular-drift component (previous section) remains a real, unresolved, non-concomitant-training-data concern, unaffected by this finding -- the two are independent and should not be conflated when interpreting the residual MOC gap.

## MOC rebuild with the calibration-window fix: `moc_streamfunction_v4.m` + `moc_pilot_v5.m` (2026-08-17)

User asked directly: did all the investigation above (ECCO calibration-window confound, West-vs-East isolation, GEM-vintage/Brazil Current, seasonal-component finding) actually change anything in the `moc_pilot_v4_vs_full_IES.png` figure already on disk? Checked the code: **yes, partially.** `moc_streamfunction_v3.m` (full array, blue curve) still loaded Step E `v2`'s output -- built and last run *before* Step E `v3` (the calibration-window fix, see "West-vs-East isolation" section above) existed -- so its absolute reference for every gap was still anchored to its own full 2013-2022 mean rather than the original 2013-2017 window. `moc_pilot_v4.m` (Pilot, red curve), by contrast, was **not** affected at all: its BPR-based reference is computed independently and never touches Step E's ECCO time-mean-reference step; the one Step E value it does borrow (`AREA_TOPO_1..8`, for the aggregate `AREA_TOPO_pilot(z)`) is a static, time-invariant bathymetric quantity, unaffected by the calibration-window issue either.

**Rebuilt only the full-array side**: `moc_streamfunction_v4.m` is identical to `v3.m` except it loads `Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3` (the calibration-fixed Step E output) instead of `_v2`. `moc_pilot_v5.m` is identical to `moc_pilot_v4.m`'s own Pilot computation (unchanged, per the above) -- only its comparison/plotting section now loads `moc_streamfunction_v4.mat` instead of `v3.mat`.

**Result**: a modest shift, not a reshuffle -- full-array mean 13.48→**12.97 Sv** (`h_star` 1210→1190dbar); Pilot unchanged at 19.53 Sv (as expected); correlation essentially unchanged (0.129→0.135). Confirms the calibration-window confound is mostly a near-constant per-gap offset, not something that changes the overall shape/story of the full-array-vs-Pilot comparison for the *whole 2013-2022 record* -- its main practical impact is specifically when comparing against the *paper's* 2013-2017-only numbers (see the earlier "West-vs-East isolation" section), less so for this internal, full-record comparison.

Outputs: `moc_streamfunction_v4_hovmoller.png`, `moc_streamfunction_v4_index_IES.png`, `moc_streamfunction_v4_mean_profile.png`, `moc_streamfunction_v4.mat`; `moc_pilot_v5_vs_full_IES.png`, `moc_pilot_v5_profile_comparison_IES.png`, `moc_pilot_v5.mat`. **`moc_streamfunction_v4.m`/`moc_pilot_v5.m` are now the current, most-correct MOC scripts** -- prefer these over `v3`/`v4` (or earlier) going forward. `v1`-`v3` (streamfunction) and `v1`-`v4` (Pilot) not deleted, kept as documented history.

## Ruled out BPR and AREA_TOPO as (additional) causes of the 2017-2019 full-array dip (2026-08-17, diagnostics only)

User asked directly whether the deep dip in the full-array MOCup curve (down toward/below 0, even to -50Sv on some days) during the 2017-08/2020-01 coverage-outage window could still be a BPR (bottom pressure) or `AREA_TOPO` problem, beyond the already-established "missing width" mechanism (gaps 3-6 mostly/fully absent, especially D-P8 at 71.3% of the total section width -- see the earlier mechanism-investigation section). Checked both directly:

- **`AREA_TOPO`, upper layer (0-1300dbar, where `h_star` lives): essentially 1.000 at every depth, for all 8 gaps, both the "surviving" gaps during the outage (1,2,7,8: A-C, C-D, P4-P2, P2-P1) and the "missing" ones (3,4,5,6)** -- no topographic blocking at all in the depth range that actually determines the MOCup index; it only starts to matter below ~1400dbar (e.g. P2-P1 drops to 0.833 there), well beneath where this diagnostic needs to look. Cannot explain or contribute to the dip.
- **BPR at the surviving sites (A, C, D, P4, P2, P1): no evidence of drift or malfunction specifically during the window** -- mean bottom pressure anomaly in vs. out of the 2017-2019 window is essentially identical at every site, and a day-to-day-jump check (>5x the outside-window std) found only a handful of isolated jumps (P2: 12/883 days, P1: 3/883, P4: 1/883) of unremarkable magnitude, not a systematic bias or an obvious sensor discontinuity.

**Conclusion: both ruled out.** The dip remains fully explained by the sampling/representativeness mechanism already documented -- losing ~78% of the section's nominal width leaves only a narrow ~21% sliver (gaps 1,2,7,8) which, regardless of how good ITS OWN data quality is, cannot represent the true basin-wide overturning. Not a data-quality artifact in BPR or a topographic-correction artifact in `AREA_TOPO`.

## Testing whether a fixed reduced-gap subset explains the dip, treating the Pilot as approximate ground truth (2026-08-17)

User's test design: since the Pilot (A-P1 only) never suffers the intermediate-site outages, treat it as approximate ground truth. If the full array is restricted to *only* the gaps/sites that have data during the 2017-2019 outage -- applied as a FIXED subset across the *entire* 2013-2022 record, not just during the outage -- and it still deviates from the Pilot outside that window too, that would mean the deviation isn't purely "too few gaps = unrepresentative", but something more general about using those particular gaps. If it tracks the Pilot fine outside the window, that supports the sampling/representativeness explanation specifically.

**First checked exactly which sites have data during 2017-08-01/2020-01-01** (per-site, not per-gap, to avoid conflating two sites' coverage): A 100%, C 100%, D 76.6%, P8 34.0%, P4 70.4%, P2 100%, P1 100% -- **only P5 and P6 are completely (0%) absent** throughout the window.

**Test 1 (`moc_streamfunction_test_no_P5P6.m`)**: restricted the full-array calculation to gaps 1,2,3,7,8 (A-C, C-D, D-P8, P4-P2, P2-P1 -- i.e. drop only 4,5,6, which all require P5 and/or P6) for every single day of the whole 2013-2022 record, and compared against the Pilot (`moc_pilot_v5.m`):

| | Full (5-gap restricted) | Pilot | corr |
|---|---|---|---|
| Outside 2017-2019 | 16.79 Sv | 18.53 Sv | 0.266 |
| Inside 2017-2019 | **-2.70 Sv** | 22.39 Sv | 0.059 |

**Result confirms the hypothesis, with a nuance.** Outside the outage window, the 5-gap-restricted full array's *mean* tracks the Pilot reasonably well (16.8 vs 18.5 Sv) -- consistent with "fewer gaps, when they all have data, still gives a roughly representative estimate." But its *variability* is much larger than the Pilot's (`std=21.3` Sv vs. `7.7` Sv) -- fewer gaps means more exposure to local/mesoscale noise not averaged out by having more independent measurements, a real but distinct effect from bias. Inside the outage window, the mean collapses to -2.70 Sv even though gaps 1,2,3,7,8 are nominally "included" -- because those same gaps are themselves only *partially* covered during 2017-2019 (D 76.6%, P8 34.0%, P4 70.4%, none at 100%), compounding on top of losing gaps 4-6 entirely. **So it's not simply "5 gaps vs. 8 gaps" -- it's the combination of losing 3 gaps outright AND partial coverage within the remaining 5 that produces the severe, biased dip, not just noise from having fewer gaps.**

Outputs: `moc_test_no_P5P6_vs_pilot.png`, `moc_test_no_P5P6.mat` (test script, `moc_streamfunction_test_no_P5P6.m`, kept for reference -- not a "vN" of the main line, a diagnostic).

### Test 2, stricter: only the 100%-covered sites (A, C, P2, P1 -- gaps 1 and 8 only)

Follow-up: restrict to *only* the sites with zero missing days anywhere in 2017-2019 -- A, C, P2, P1 -- i.e. gaps 1 (A-C) and 8 (P2-P1) only, dropping every other gap (2,3,4,5,6,7) entirely, again as a fixed subset for the whole 2013-2022 record (`moc_streamfunction_test_only100pct.m`).

**Result: degenerate, not just biased.** `h_star` collapses to 0dbar (the time-mean cumulative-transport profile never develops an interior maximum at all) and the MOCup index is flat at essentially 0 Sv throughout the *entire* record (mean -0.11 Sv, `moc_test_only100pct_vs_pilot.png` shows the blue curve pinned near zero against the Pilot's normal ~10-35 Sv variability) -- utterly unlike the Pilot, everywhere, not just during the outage.

**This is not a new problem -- it's a different construction than the Pilot, despite overlapping sites.** Gaps 1+8 alone only span ~1048km (A-C 366km + P2-P1 682km) -- 13.4% of the total ~7846km section width, entirely excluding the wide D-P8 interior gap (71.3% of the width) and everything else in between. The Pilot, by contrast, computes `Gpan_P1-Gpan_A` *directly* -- by the telescoping identity already proven exactly (`(Gpan_C-Gpan_A)+...+(Gpan_P1-Gpan_P2)=Gpan_P1-Gpan_A`), this implicitly integrates the density/dynamic-height structure across the *entire* basin width even though only 2 sites are physically measured. "Sum of gap1+gap8" and "the Pilot" are fundamentally different calculations that happen to share 2 of their 4 site names -- not two versions of the same thing.

### Consistency check: a "virtual single gap" A-to-P1, built via the Step-E-style V0/dx2 formula, exactly reproduces `moc_pilot_v5.m`

Direct request: build the A-to-P1 calculation using the *same* `V0`/`dx2`/`vel00` formula structure `moc_streamfunction_v4.m` and the gap-restriction tests use (not `moc_pilot_v5.m`'s own standalone code path), treating A-to-P1 as a single "virtual gap" with the full `dx_total` width and the aggregate `AREA_TOPO_pilot(z)`, then compare against `moc_pilot_v5.m`'s independently-coded result -- a from-scratch cross-check of the whole telescoping argument carried through the *entire* MOC calculation (BPR reference, `AREA_TOPO`, Ekman, shelf correction), not just the baroclinic `RelTPUD` piece already confirmed analytically days earlier.

**Result: identical to floating-point precision.** `h_star=1280dbar` both ways; mean `19.5323` Sv both ways; `max|diff|=1.4e-14` Sv over all 3382 common days. Confirms the telescoping identity holds through the full pipeline (not just the relative/baroclinic part), and that `moc_pilot_v5.m`'s standalone implementation and the Step-E-style `V0`/`dx2` formula pattern are algebraically the same calculation, just expressed two different ways -- a clean internal-consistency validation of the whole MOC codebase. (Not saved as a persistent script -- a one-off scratch verification; the derivation is fully documented here.)

**Together, these three tests give a complete, mechanistically-understood picture of the full-array-vs-Pilot comparison**: the baroclinic part telescopes exactly (proven twice now, algebraically and numerically); a reduced-but-still-wide gap subset (1,2,3,7,8) tracks the Pilot's mean reasonably outside outage periods but is noisier; a narrow-gap-only subset (1,8) is a fundamentally different, much-too-narrow construction that shouldn't be expected to match the Pilot at all; and the true full-array method's own remaining discrepancies (already documented in the sections above -- calibration-window, GEM-vintage, shelf correction) are the real, substantive open items, not something about "using few sites" per se.

## Quantifying the shelf correction's effect on the Pilot specifically (2026-08-17)

Both `moc_streamfunction_v4.m` and `moc_pilot_v4.m`/`v5.m` add the same static shelf-transport correction (`West_TPUD=-4.51`Sv, `East_TPUD_new≈0`, see the "MOC shelf-transport fix" section above) for a fair comparison -- but its effect had only been quantified for the full array directly. Since `moc_pilot_v3.mat` (no shelf) and `moc_pilot_v4.mat` (+shelf) differ *only* in this one correction (confirmed from `v4.m`'s own header: "adds the shelf correction" is the sole change from `v3`), compared them directly to isolate the effect on the Pilot:

| | `v3` (no shelf) | `v4` (+shelf) | Diff |
|---|---|---|---|
| `h_star` | 1270dbar | 1280dbar | +10dbar |
| Mean MOCup | 23.97 Sv | 19.53 Sv | **-4.44 Sv** |
| Std | 7.61 Sv | 7.67 Sv | ~unchanged |

**Nearly a pure constant offset**: day-to-day shift has mean -4.44 Sv but std only 0.069 Sv -- essentially the entire West shelf magnitude (-4.51 Sv), consistent with `West_TPUD`/`East_TPUD_new` being static/time-invariant. The time-mean profile shows *why* `h_star` shifts slightly despite the offset being nearly constant overall: the shelf correction has its own vertical structure, building up from ~0 at the surface to its full -4.494 Sv value by ~1200dbar, then staying exactly flat at every depth below that (physically sensible -- shelf/boundary transport is an upper-ocean phenomenon, no deep-reaching structure) -- so it's not *perfectly* uniform in the 0-1200dbar range specifically, which is enough to nudge where the cumulative profile peaks.

**Conclusion**: the shelf correction affects the Pilot the same systematic way it affects the full array -- a near-constant shift of magnitude matching the total shelf transport, no meaningful extra noise introduced. Expected and consistent, not a new concern.

## Restricting the time-mean profile to exactly Kersalé et al. (2021)'s own window: `moc_profile_comparison_2013_2017.m` (2026-08-17)

User's observation: every profile-comparison figure so far (`moc_pilot_v2/v3/v4/v5_profile_comparison_IES.png`) averages `Psi_raw` over the *full* 2013-2022 record -- not the paper's own 2013-09-11/2017-07-17 window specifically, which is what Figure 2a actually shows. Confirmed directly: none of the existing figures restrict to that period.

`moc_profile_comparison_2013_2017.m` reuses the already-saved `Psi_raw` `[depth x time]` fields from `moc_streamfunction_v4.mat` and `moc_pilot_v5.mat` directly (no pipeline rerun needed -- `Psi_raw` is a per-day quantity, identical regardless of which sub-period is later averaged over) and recomputes the time-mean profile and `h_star` using *only* the 1405 days in the paper's window.

**Result: a much closer match to the paper than the full-2013-2022 comparison.**

| | `h_star` | Peak |
|---|---|---|
| Full array (this repo, 2013-2017 only) | 1320dbar | 20.31 Sv |
| Full array (paper) | 1315dbar | 17.3 Sv |
| Pilot (this repo, 2013-2017 only) | 1270dbar | 18.77 Sv |
| Pilot (paper) | 1315dbar | 17.7 Sv |

Both `h_star` values now land within 45dbar of the paper's 1315dbar (vs. 1190/1280dbar for the full 2013-2022 comparison), and the Pilot's peak (18.77 Sv) is within ~1.1 Sv of the paper's 17.7 Sv -- much closer than any full-record comparison achieved. Visually (`moc_profile_comparison_2013_2017_IES.png`), the full-array and Pilot curves now closely track each other through most of the water column (0-3000dbar), only diverging more visibly below ~3500dbar -- structurally similar to what Figure 2a shows, unlike the full-2013-2022 comparison where the two curves separated much earlier and more severely.

**Confirms period-matching was a real, previously-uncontrolled factor** on top of the calibration-window/GEM-vintage/shelf issues already documented -- averaging over the extended 2013-2022 record (including the 2017-2019 coverage outage and whatever secular drift exists in the West GEM table) measurably degrades the match to the paper beyond what the methodology differences alone would predict.

Outputs: `moc_profile_comparison_2013_2017_IES.png`, `moc_profile_comparison_2013_2017.mat`.

## MAJOR FINDING: `~/matlab/stat/lowpass_filter.m` (shared utility, used throughout this project for every "45-day lowpass") does NOT preserve the mean -- confirmed root cause, 2026-08-17, project-wide implications not yet scoped

While chasing the remaining ~1.1 Sv gap between Marion's TRUE ORIGINAL 2020 pipeline output (same data, same general methodology as the paper -- the most apples-to-apples comparison available) and the published 17.3 Sv, applying this repo's correct Kersalé-definition + shelf correction directly to her untouched `Full_Depth_OverturningEstimate_Cst_for_MHT_TRUEORIGINAL_backup.mat`:

- **Mean of the RAW (unfiltered) daily MOCup series: 17.70 Sv** -- remarkably close to the paper's 17.3 Sv.
- **Mean of `lowpass_filter(MOCup_raw, 45)` (this project's standard "45-day lowpass" convention, used everywhere): 18.44 Sv** -- the number previously reported as the shelf-correction validation result.

**Root-caused exactly**, by replicating `lowpass_filter.m` step by step:

```matlab
function yg=lowpass_filter(var,npf);
um=ones(size(var));
inxs=(npf-1)/2;
myfilter=blackman(npf); myfilter=myfilter/sum(myfilter);
y=conv(var,myfilter); umf=conv(um,myfilter);
y=y(inxs+1:length(y)-inxs); umf=umf(inxs+1:length(umf)-inxs);
y=y./umf;                                    % <- proper, edge-corrected, MEAN-PRESERVING lowpass ends here
y=y*maxvar(var,y);                           % <- rescale #1: multiplies by a no-intercept regression gain
g=(y(:)'*y(:))\(y(:)'*var(:));               % <- rescale #2: another no-intercept regression gain
yg=y*g;                                      % <- final output
```

(`maxvar.m`: `g=(b(:)'*b(:))\(b(:)'*a(:));` -- "maximizes the amplitude of b in relation to a by least-squares fit.") The actual filtering (Blackman-window convolution + edge correction, ending at the `y=y./umf` line) is correctly mean-preserving on its own -- confirmed directly: mean after that step = 17.7001 Sv, matching the raw mean (17.7008) almost exactly. **The entire bias comes from the two subsequent amplitude-rescaling steps**, intended to "restore variance lost to smoothing" (a real, legitimate concern -- lowpass filtering does reduce a signal's variance) but implemented as a single **scalar multiplicative gain applied to the whole signal** -- which inflates the mean/DC component by the same factor as the AC/variance component, even though the mean needed no correction at all. Measured gain for this dataset: `g1=1.041878`, `g2≈1.000000` (the second pass is nearly a no-op once the first has already best-fit the scale) -- combined gain 1.0419, and `raw_mean * 1.0419 = 18.4421` Sv, matching the actual `lowpass_filter` output (18.4413 Sv) almost exactly. **Confirms the mechanism completely**, not just a suspicious coincidence.

**Implication -- likely explains a real chunk of the "unexplained" MOC residual**: since 17.70 Sv (unfiltered mean) is far closer to the paper's 17.3 Sv than 18.44 Sv (this project's usual lowpass-mean convention), it's plausible the paper's own reported statistics come from a mean of the raw/daily series (or an unbiased smoothing method), not something equivalent to this project's `lowpass_filter`-then-mean convention.

**Broader implication, NOT YET SCOPED**: `lowpass_filter(x,45)` is the standard smoothing convention used throughout this entire project -- `mov`, `Heat_total`/`Q_*`, and every MOC index reported in `CLAUDE.md` and in conversation this whole session used `nanmean(lowpass_filter(x,45))` (or an equivalent per-depth-level lowpass) as "the" reported mean. Every one of those means may carry a similar multiplicative inflation, proportional to how much variance-restoration gain (`g1*g2`) that specific signal happened to need -- the exact gain will differ per-signal (depends on how "peaky"/high-frequency the underlying daily series is relative to its 45-day-smoothed envelope), so this is NOT a single fixed correction factor across all reported numbers. **User's decision (2026-08-17): document this now, scope/re-audit which specific previously-reported means are affected in a later session** -- `lowpass_filter.m` itself has not been modified (shared utility outside both repos, in `~/matlab/stat/` -- same location as the earlier `sinfitb_tot.m` ill-conditioning fix, see `reference_sinfitb_tot_fix` memory). No fix applied yet; this section exists to make sure the finding isn't lost before the fuller audit happens.

## Scoping the lowpass bias within the MOC calculation, and reconstructing Kersalé et al. (2021) Figure 5a/5b (2026-08-17)

Two follow-ups to the finding above, both requested together: (1) scope the bias specifically for MOC before touching Mov/MHT, (2) reconstruct the paper's Figure 5 (temporal anomaly of MOCup, both configurations).

**`lowpass_filter_fixed.m`** (new, in `MOV_Project`, NOT touching the shared `~/matlab/stat/lowpass_filter.m`): identical to the original except it stops after the Blackman-window convolution + edge correction, omitting the two buggy amplitude-rescaling steps. Confirmed mean-preserving on its own.

**Scope for MOC, comparing raw / biased-lowpass / fixed-lowpass means**:

| | Raw | Biased (`lowpass_filter`) | Fixed (`lowpass_filter_fixed`) |
|---|---|---|---|
| Full array, full record (2013-2022) | 12.35 Sv | 12.97 Sv | 12.36 Sv |
| Pilot, full record | 19.19 Sv | 19.53 Sv | 19.19 Sv |
| Full array, 2013-2017 only | 20.11 Sv | 21.08 Sv | 20.08 Sv |
| Pilot, 2013-2017 only (paper: 17.7 Sv) | 18.77 Sv | 19.10 Sv | **18.76 Sv** |

**Correlation is unaffected** (0.135 full record, 0.514 for 2013-2017, identical whether biased or fixed filter is used) -- expected, since the bias is a pure multiplicative scalar applied uniformly to the whole signal, which cancels out of a correlation coefficient by construction. Confirms the bias only ever distorted *means*, never the *shape/timing* of variability already investigated in the gap-subset tests above.

**This meaningfully revises the earlier "Marion-original-vs-published, ~1.1 Sv unexplained" reading**: that number came from `nanmean(lowpass_filter(...))` too. Using the raw/unbiased mean instead (already computed in the previous section): 17.70 Sv vs. the paper's 17.3 Sv -- only ~0.4 Sv apart, not ~1.1 Sv. **Most of what looked like an unexplained baseline gap was actually this same lowpass bug**, not a separate mystery. The remaining, still-real, still-unexplained pieces are: (a) the West GEM-vintage/Brazil-Current effect (this repo's own reprocessing vs. Marion's original, ~2-3 Sv, concentrated in gap A-C -- unaffected by the lowpass finding, a genuine data-provenance difference), and (b) a small ~0.4 Sv baseline residual even with fully original data and an unbiased mean.

**Figure 5a/5b reconstruction** (`moc_anomaly_fig5_2013_2017.m`): temporal anomaly (each series minus its *own* 2013-2017 mean) for the full array and Pilot, using `lowpass_filter_fixed`, restricted to the paper's exact window. Sidesteps the whole absolute-level discrepancy (calibration-window/GEM-vintage/lowpass-bias) entirely -- isolates whether the two methods' *day-to-day variability* is consistent, which is what the paper's `r=0.73` actually measures.

**Result**: `corr(full anomaly, pilot anomaly) = 0.514` (paper: 0.73) -- visually (`moc_anomaly_fig5_2013_2017_IES.png`), the two anomaly curves track each other reasonably well through much of the record (peaks/troughs align around mid-2015, early/mid-2016), but the full array has much larger swings (anomaly `std=11.86` Sv vs. Pilot's `7.34` Sv) and one dramatic, Pilot-unmatched dip around January 2017 (full array down to -30Sv). **Not reconstructed**: panel 5c (`MOCab`, the abyssal cell -- a separate 3150-4300dbar depth-integrated quantity not yet computed in this project) and the paper's gray "estimated daily accuracy" shading (no equivalent error-budget methodology built here).

**Open question going forward**: why the correlation (0.51) still falls short of the paper's 0.73, now that the mean-level confounds are cleanly separated out -- likely candidates not yet tested: the remaining coverage-outage days even within 2013-2017 (some `n_gaps_valid<8` days exist even in this "good" window), the GEM-vintage secular drift affecting day-to-day shape (not just the mean level), or a genuine difference in how well `AREA_TOPO_pilot`'s aggregate approximation captures real day-to-day variability vs. per-gap-resolved topography.

Outputs: `lowpass_filter_fixed.m`, `moc_anomaly_fig5_2013_2017_IES.png`, `moc_anomaly_fig5_2013_2017.mat`.

### Daily (unfiltered) standard deviation matches the paper almost exactly for the full array

The paper's results text (Section 3.1) reports the DAILY (not lowpass-filtered) standard deviation of MOCup over 2013-2017: "the temporal standard deviation (σ) of the daily MOCup cell strength calculated over the continuous 4 years... with the full resolution array is equal to 15.4 Sv... This daily σ exceeds the previous estimate of ~9 Sv using the pilot array configuration." Computed the equivalent directly from `MOCup_raw` (no `lowpass_filter` involved at all, sidestepping the whole bias question above), restricted to the 2013-2017 window:

| | This repo | Paper |
|---|---|---|
| Full array, daily σ | **15.39 Sv** | **15.4 Sv** |
| Pilot, daily σ | 10.51 Sv | ~9 Sv |
| Ratio (full/pilot) | 1.46 | 1.71 |

**The full array's daily variability is an almost exact match** (15.39 vs. 15.4 Sv) -- strong evidence that this repo's day-to-day MOC calculation (formula, definition, physics) replicates the paper's methodology very precisely; the remaining discrepancy is specifically about *mean-level* calibration (GEM-vintage/Brazil-Current, already documented), not a structural problem with the approach. The Pilot's daily σ runs somewhat higher than the paper's (10.51 vs. ~9 Sv, +17%) -- not yet investigated further, but notably closer than the ratio of the two σ's, which undershoots the paper's more pronounced full-vs-pilot variability gap (1.46 vs. 1.71) -- consistent with the still-open "why is the correlation only 0.51 vs. 0.73" question from the section above, suggesting some of the full array's real high-frequency variability (which the paper's better-resolved 8-gap array captures more of, e.g. eddies/Agulhas Rings per the paper's own discussion) isn't as fully distinguished from the Pilot in this repo's implementation as in the original.

### Investigating the correlation shortfall (0.51 vs. paper's 0.73): coverage and lag ruled out, points toward excess noise specific to this repo's Pilot (2026-08-17)

Three candidate explanations tested directly:

1. **Residual coverage gaps within 2013-2017 -- ruled out.** Correlation is *identical* (0.514) regardless of the `n_gaps_valid` threshold applied (≥0, ≥5, ≥6, ≥7, ≥8) -- because coverage is already 8/8 for essentially the entire 2013-2017 window (the outage doesn't begin until August 2017, right at the window's edge). No coverage variation exists within this period to test against.
2. **Time lag between the two methods -- small effect, not the main driver.** Lagged cross-correlation (±30 days) peaks at +5/+6 days (`xcorr=0.542-0.543`) vs. `0.514` at zero lag -- a ~0.03 improvement, far short of closing the gap to 0.73.
3. **Excess noise concentrated in this repo's own Pilot reconstruction -- the likely remaining explanation, not yet further decomposed.** Consistent with the daily-σ finding above: the full array's σ matches the paper almost exactly (15.39 vs 15.4 Sv) while the Pilot's runs 17% high (10.51 vs ~9 Sv) -- suggesting some of the Pilot's variance in this repo is uncorrelated noise not present in the paper's own Pilot construction, diluting the correlation. Candidate sources not yet tested: the aggregate `AREA_TOPO_pilot(z)` approximation, GEM-vintage secular drift affecting day-to-day shape (not just the mean), or BPR calibration specifics.

**Before investigating further, stepping back to look at the raw `tau1000` data directly across all 9 sites** (A through P1) -- see next section.

### Looking at raw τ1000 directly: sites A and P1 (both used by the Pilot) are among the noisiest of the 9 -- likely explains the excess Pilot noise (2026-08-17)

Plotted raw `tau1000` (the fundamental PIES/CPIES acoustic travel-time measurement -- everything downstream, GEM-derived T/S, `Gpan`, BPR-referenced transport, is built on top of this) for all 9 sites, west-to-east, over the 2013-2017 window (`check_raw_tau1000_all_sites.m` → `raw_tau1000_all_sites_2013_2017.png`). Visually, the West sites (A, C, D) look noticeably more scattered/jagged point-to-point than most East sites (P8, P6, P5, P4, P2, P1) -- consistent with the already-documented fact that the East `Daily_Tau` source is pre-filtered (de-tided/de-drifted/de-spiked, 72h lowpass) while `samba_w.mat`'s West provenance doesn't carry the same documented treatment.

**Quantified directly** (`check_tau_noise.m`): day-to-day noise (`std(diff(tau1000))`) and a noise-to-signal ratio (`std(diff)/std(raw)`) per site, 2013-2017:

| Site | std(diff) | noise/signal (×1000) |
|---|---|---|
| A | 0.00082 | 271 |
| C | 0.00083 | 153 |
| D | 0.00084 | 248 |
| P8 | **0.00029** (smoothest overall) | 134 |
| P6 | 0.00051 | 162 |
| P5 | 0.00054 | 149 |
| P4 | 0.00082 (as noisy as West!) | 179 |
| P2 | 0.00063 | 325 |
| **P1** | 0.00052 | **315** (2nd highest) |

**Not a clean West-vs-East split** (P4 is as noisy as A/C/D in absolute terms) -- but directly relevant to the Pilot specifically: **A is tied for the noisiest site in absolute terms, and P1 has the 2nd-highest noise-to-signal ratio of all 9 sites.** Since the Pilot depends *entirely* on A and P1 (50% weight each, no other sites), both endpoints' elevated noise characteristics enter the calculation essentially undiluted. The full array, by contrast, averages across 8 gaps built from all 9 sites -- including P8 (the smoothest site by a wide margin) and several other lower-noise East sites -- naturally diluting any single site's noise contribution. **This is a well-supported, data-grounded explanation for both the Pilot's elevated daily σ (10.51 vs. the paper's ~9 Sv) and the correlation shortfall (0.51 vs. 0.73)** -- not a processing bug, but a real consequence of the Pilot method's inherent sensitivity to its two endpoint sites' individual data quality, which happen to include two relatively noisy measurements at 34.5°S specifically.

**Not yet resolved/quantified further**: whether this noise is itself instrumental (sensor/calibration precision differences between deployments) or partly real high-frequency oceanographic signal (e.g., local coastal-boundary variability at A and P1 specifically, both being edge/shelf-adjacent sites) that a 2-point method is simply more exposed to than an 8-gap spatial average. Distinguishing these would need either an independent noise-floor estimate per site (e.g., from published PIES instrument specs) or checking whether the "excess" noise at A/P1 has a spectral signature more consistent with sensor noise (white/high-frequency) or ocean variability (red/mesoscale).

Outputs: `check_raw_tau1000_all_sites.m`, `raw_tau1000_all_sites_2013_2017.png`, `check_tau_noise.m` (table reproduced above; script prints, doesn't save a `.mat`/figure).

**Extended to the full 2013-2022 record** (`check_raw_tau1000_all_sites_full.m` → `raw_tau1000_all_sites_2013_2022.png`), explicitly regridded onto a complete daily calendar axis (missing days set to `NaN` rather than letting the plot bridge across them, which would otherwise hide real gaps by connecting the nearest surrounding valid points). Visually confirms the already-documented coverage history cleanly: **P5 and P6** share one long simultaneous outage from roughly mid-2017 to late 2021; **P8** has a gap around 2018-2019; **D** has a shorter gap around 2018; **P4** has a shorter gap around 2019; **A, C, P2, P1** run continuously with no visible gaps the entire record. One thing noted: site **A** shows sharp, isolated downward spikes standing out from the surrounding signal -- investigated directly below.

### Site A's spikes: two are likely instrument glitches (not yet fixed), one is real ocean signal -- ⚠️ revisit tau processing with GEM later

Outlier-detected directly on `tau1000_A` (residual from a 31-day running median, 6σ threshold, `threshold=0.01116`):

| Date | `tau1000_A` | Residual |
|---|---|---|
| 2020-06-20 | 1.30906 | -0.0160 |
| 2022-10-17 | 1.30692 | -0.0189 |
| 2022-10-18 | 1.30838 | -0.0181 |

**These three points have the classic signature of an instrument glitch**: sharp, isolated, single/double-day V-shaped drops and immediate recovery, magnitude (~0.016-0.019s) roughly **20-25x** the typical day-to-day noise level already quantified for site A (`std(diff)≈0.0008s`, see the raw-`tau1000` noise section above) -- structurally the same *kind* of problem as the East `-9999` sentinel bug found and fixed earlier in this project (`SAMBA_E_IES`'s `read_daily_tau.m`), just not yet root-caused or fixed for this West source.

**By contrast, the dip visually noticed around late May/early June 2021 is NOT a glitch** -- checking the surrounding week shows a smooth ~10-day gradual descent and recovery (26-May to 5-Jun-2021), well below the 6σ outlier threshold -- consistent with real oceanographic variability (e.g., a mesoscale/submesoscale event), not an instrument artifact.

**Reinforces the earlier documented suspicion**: the East `Daily_Tau` dataset is explicitly documented as pre-processed ("de-tided/de-drifted/de-spiked," 72h lowpass) before this project ever sees it, while `samba_w.mat`'s West provenance has no equivalent documented QC step -- these three spurious, unfiltered points in `tau1000_A` are direct, concrete evidence supporting that asymmetry, not just an inference from noise-level statistics.

**⚠️ Not yet fixed -- explicit reminder for a later session**: the user wants to revisit the West τ processing pipeline using GEM (likely: apply a GEM-consistency-based de-spiking/QC step to `tau1000` before profile generation, analogous to what the published East `Daily_Tau` dataset already received, rather than a simple outlier-clipping fix applied after the fact). Not scoped further this session -- flagged here and in project memory so it isn't lost.

## Component breakdown of the 2017-2019 dip: relative/baroclinic is the component that "breaks" when gaps are missing -- not a physical Ekman/barotropic anomaly (2026-08-17)

Follow-up to "does the raw `tau1000` data itself show anything strange during 2017-2019, beyond the missing sites?" (user's own reading of the raw-data figures above: no). Tested directly using Step E v3's `_RefConst`/`_RelConst`/`_EKMANconst` variants -- Marion's own original design, built specifically "to see which is driving the MOC/MHT variations at different time scales" (Step E's header comment): each variant holds ONE component (reference/barotropic, relative/baroclinic, or Ekman) at its own long-term time-mean while letting the other two vary daily, using the *same* underlying per-gap data as the normal ("Full") calculation.

`moc_components_2017_2019.m` builds the Kersalé-definition MOCup index (same construction as `moc_streamfunction_v4.m`: `V0`/`dx2`, `h_star` from the variant's own time-mean, shelf correction) for all four versions (Full, RefConst, RelConst, EKMANconst) and plots them together, both zoomed on 2017-2020 (`moc_components_2017_2019_zoom.png`) and over the full record (`moc_components_full_record.png`).

**Result**: Full, `RefConst` (reference held constant), and `EKMANconst` (Ekman held constant) track each other almost exactly throughout the *entire* record, including the severe 2017-2019 dip (down to -30 to -55 Sv at various points in all three) -- holding either of those two components constant does **not** prevent or soften the dip. `RelConst` (relative/baroclinic held constant) is the one clear outlier -- visibly smoother and far less extreme throughout the whole record, not just 2017-2019, staying mostly in the 0-20 Sv range even when the other three plunge to -40/-55 Sv.

**Why**: `RelTPUD_i_Const` is built from `Gpan_XConst=nanmean(Gpan_X')'` -- a single, time-*invariant* profile (the whole-record mean), always available regardless of whether that specific gap has valid data on a given day. The normal ("Full") calculation instead needs each gap's *daily* `Gpan` (from that day's T/S), which goes to `NaN` outright whenever the site is missing. **This means the daily relative/baroclinic component is mechanically the one piece that cannot be computed at all without that day's real data** -- not because it's more "physically anomalous" during the outage, but because it's the only one of the three components with no fallback when data is missing. Reference (BPR) and Ekman, by contrast, don't get "rescued" the same way in their own `_Const` variants because -- as already found -- P5/P6's outage affects those components too on the same days (their own `_RefConst`/`_EKMANconst` variants still require the *other two* components' daily values, which are equally unavailable).

**Conclusion, directly answering the user's question**: the raw data itself shows nothing anomalous in Ekman or barotropic/BPR terms during 2017-2019 -- confirmed now from a third, independent angle (component decomposition, on top of the earlier coverage-percentage and raw-`tau1000` checks). The dip is entirely a **missing-relative-data** mechanical artifact of the coverage outage, localized specifically to the baroclinic/relative component's total dependence on that day's real T/S measurements, with no fallback -- not a physical signal, and not hiding in the other two components either.

Outputs: `moc_components_2017_2019.m`, `moc_components_2017_2019_zoom.png`, `moc_components_full_record.png`, `moc_components_2017_2019.mat`.

## Is the full-array-vs-Pilot correlation shortfall a phase/lag issue? `moc_lag_analysis_2013_2017.m` (2026-08-17)

User's direct observation: in Kersalé et al. (2021)'s Figure 5, the pilot and full-array anomaly curves look in-phase; in this repo's `moc_anomaly_fig5_2013_2017_IES.png` reconstruction, they visibly look out of phase in places (e.g. early 2014, early 2015).

**Time-varying lag on the 45-day-lowpass series** (6-month sliding window, ±30-day search): NOT a constant offset. Alternates between windows of near-zero lag with high correlation (Jul2015-Jan2016, Nov2016-Mar2017: `r=0.7-0.94`, matching or beating the paper's 0.73) and windows of unstable, large apparent lag with weak correlation. Extending the search to ±90 days destabilizes most of the "large lag" windows (several flip sign/magnitude entirely, e.g. Oct/Nov-14 goes from +30 to -72) -- a sign the lag estimate is unreliable there (tested lag becomes too large a fraction of the 180-day window for a robust estimate), not evidence of a genuine, larger phase shift.

**Same test on RAW (unfiltered) daily data -- the actual answer to the user's follow-up question ("if we skip the 45-day filter, are the curves more in phase?")**: **yes, clearly.** Whole-record best lag on raw data = **0 days** (already optimal -- no lag improves the raw correlation at all, `r=0.485`). The time-varying version is lag=0 in 76% of windows (vs. only 17% for the lowpass version) -- and the rare non-zero-lag raw windows land exactly on the same windows that also have near-zero raw correlation (Sep-Nov 2014), i.e. noise rather than a real, determinable lag there either.

**Conclusion**: the apparent phase offset visible in the lowpass-filtered anomaly plot is an **artifact of the 45-day lowpass filter itself** (smoothing broadens and can shift the apparent location of peaks/troughs), not a real property of the underlying daily data -- without the filter, the two methods are essentially in phase throughout. This means the remaining correlation shortfall (raw `r=0.485` vs. the paper's `0.73`) is *not* a phase/timing problem -- it's specifically about how well the *amplitude/shape* of fluctuations match during certain weaker windows (which, per the earlier `tau1000` noise findings, may trace back to sites A and P1's individually elevated raw measurement noise). Next natural step (not yet done): check whether the weak-correlation windows found here line up with the noisiest stretches of `tau1000_A`/`tau1000_P1`.

Outputs: `moc_lag_analysis_2013_2017.m`, `moc_lag_analysis_2013_2017_IES.png`.

### Following up: do the weak-correlation windows coincide with elevated A/P1 noise? Yes, especially for P1

Direct test (`moc_weak_windows_vs_noise.m`): for the same 6-month sliding windows as the lag analysis, computed the raw full-array-vs-Pilot correlation alongside each site's local noise level (`std(diff(tau1000))` within that window) for sites A and P1.

**Correlation between window-correlation and local noise**: site A `r=-0.114` (weak), site P1 `r=-0.343` (moderate, right-signed -- higher P1 noise → lower window correlation). Visually (`moc_weak_windows_vs_noise_IES.png`), P1's noise shows a clear long-term *declining* trend from 2014 through late 2016, reaching its minimum right around late 2016/early 2017 -- exactly when window correlation peaks (0.6-0.82, the best stretch of the whole 2013-2017 record, first identified in the lag-analysis section above). Site A's noise also declines somewhat over the same period but less monotonically and less tightly coupled to the correlation pattern.

**Conclusion**: a real, moderate, directionally-consistent relationship exists between P1's raw measurement noise and how well the Pilot tracks the full array, driven mainly by a *long-term trend* (P1 got quieter over 2014-2017, and Pilot-vs-full agreement improved over the same span) rather than noise explaining every individual weak window (e.g. the Sep-Nov 2014 dip doesn't coincide with P1's peak noise specifically). Site A's noise is not a strong driver of the window-to-window correlation pattern by this metric. This refines (doesn't just restate) the earlier `tau1000` noise finding: it's not only that A and P1 are noisier *sites* on average -- P1's noise *level itself varies over time*, and that time-variation measurably tracks the Pilot's reliability.

## Raw (unfiltered) Fig 5 reconstruction, and confirming r=0.73 is a 2013-2017-only, DAILY statistic: `moc_anomaly_fig5_raw_2013_2017.m` (2026-08-18)

Follow-up to the lowpass-filtered `moc_anomaly_fig5_2013_2017.m` (r=0.51) and the lag analysis above (raw lag=0, so the filtered version's apparent phase offset was a filtering artifact, not real): rebuilt the same Figure 5a/5b anomaly reconstruction on the RAW daily `MOCup_raw`/`MOCup_pilot_raw` series directly (no `lowpass_filter`/`lowpass_filter_fixed` involved at all), reusing the already-saved fields from `moc_streamfunction_v4.mat`/`moc_pilot_v5.mat` — no pipeline rerun.

**Verified directly in the paper's text first** (`pdftotext -layout kersale_jgr2021_samba.pdf`, Section 3.1), per the user's question of whether r=0.73 is a 2013-2017-only or whole-series statistic: *"the correlation coefficient... between the MOCup calculated over 2013–2017 using the pilot array configuration... and the estimate using the full resolution of moorings... is r = 0.73."* Confirms **2013-2017 only** — exactly the window this project has used throughout, ruling out a period mismatch as an explanation for the correlation gap. Also checked for any explicit lowpass/smoothing mention tied to Figure 5 or this correlation specifically (`grep -i "45-day\|low-pass\|lowpass\|smooth"`) — found none; the only low-pass filter mentioned anywhere in the paper is a 60-day Butterworth applied to a *different* figure (Figure 9, MHT reconstruction). Combined with the paper's own wording "daily MOCup cell strength" for the companion σ=15.4 Sv statistic in the same paragraph, this is good evidence **r=0.73 itself is a daily (unfiltered) correlation** — i.e. the raw reconstruction below, not the 45-day-lowpass version, is the more appropriate one to compare against it.

**Result** (1405 days, 2013-09-11 to 2017-07-17): full array anomaly std=15.39 Sv (paper: 15.4 Sv — near-exact, already established); Pilot anomaly std=10.51 Sv (paper: ~9 Sv); **`corr(full anomaly, pilot anomaly), RAW daily = 0.485`** (paper: r=0.73; the 45-day-lowpass version was r=0.51/0.514). Visually (`moc_anomaly_fig5_raw_2013_2017_IES.png`), the two curves track each other's broad swings but the full array's daily noise is visibly rougher/less smooth than the paper's Figure 5b — consistent with the still-open "excess Pilot/full-array noise" question from the sections above, now looked at without any filter-driven distortion.

Outputs: `moc_anomaly_fig5_raw_2013_2017.m`, `moc_anomaly_fig5_raw_2013_2017_IES.png`, `moc_anomaly_fig5_raw_2013_2017.mat`.

## Testing whether the Pilot's BPR reference is the noise source: `moc_pilot_component_test.m` (2026-08-18) — disproven, BPR is contributing real signal, not noise

Direct follow-up to "vamos focar no calculo do Pilot... vamos tentar achar a causa da diferença de fase e portanto da correlação, uma vez que a magnitude está correta": since sites A and P1's raw `tau1000` are already established as relatively noisy (see the noise sections above), and the Pilot's BPR (bottom-pressure) reference is one of the two components built directly from each site's own measurement (alongside the baroclinic/`Gpan` shear), tested whether specifically the BPR reference is the noise pathway degrading correlation with the full array.

Built two variants of the Pilot's own construction (mirroring `moc_pilot_v5.m`'s exact setup — same `samba_w.mat`/IES sources for A and P1, same aggregate `AREA_TOPO_pilot`, same Ekman and shelf correction), differing only in the reference/offset step:
- **Variant 1, "BPR reference" (normal, matches `moc_pilot_v5.m` exactly)**: `Absolute_TPUD = (RelTPUD_pilot + Offset_pilot).*AREA_TOPO_pilot`.
- **Variant 2, "baroclinic only"**: `Offset_pilot` dropped entirely (`Absolute_TPUD = RelTPUD_pilot.*AREA_TOPO_pilot`) — referenced to the sea surface instead of BPR, i.e. no absolute/barotropic correction of any kind.

Both go through the same Kersalé-definition streamfunction (own time-mean `h_star`, RAW daily, 2013-2017), compared against the full array's `MOCup_raw`:

| Variant | `h_star` | Mean | Std | Corr with full array (raw, 2013-2017) |
|---|---|---|---|---|
| BPR reference (normal) | 1280 dbar | 19.19 Sv | 11.45 Sv | **0.485** |
| Baroclinic only (no reference) | 60 dbar | 1.41 Sv | 4.83 Sv | **0.281** |

**Hypothesis disproven**: removing the BPR reference makes correlation with the full array WORSE (0.281 vs 0.485), not better. `AREA_TOPO_pilot(z)`'s own `h_star` search collapses to a shallow, near-meaningless 60dbar without the BPR-derived barotropic component to anchor the profile's real overturning shape — visually (`moc_pilot_component_test_IES.png`), the baroclinic-only curve (yellow) is much flatter/lower-amplitude than both the BPR-referenced Pilot (orange) and the full array (blue), which track each other's larger swings much more closely. This means the BPR reference is contributing real, necessary signal (the barotropic component, matching what the full array's own BPR-referenced gaps also carry) — not noise to be filtered out or downweighted. Rules out "exclude/downweight BPR" as a fix for the correlation shortfall.

**Where this leaves the open question**: with lag (ruled out, raw data already in phase) and BPR reference (ruled out here) both tested and disproven as causes, the remaining candidates not yet tested are more specifically about the baroclinic/`Gpan`-based signal itself — e.g. whether `Gpan_A`/`Gpan_P1`'s own day-to-day noise (as opposed to the reference/offset step) is where sites A/P1's elevated raw `tau1000` noise actually enters the calculation, or whether the aggregate `AREA_TOPO_pilot(z)` approximation itself distorts day-to-day shape (not just the mean level, already ruled out as a correlation driver since it's a static depth profile).

Outputs: `moc_pilot_component_test.m`, `moc_pilot_component_test_IES.png`.

## Investigating Gpan_A/Gpan_P1 directly: `moc_gpan_noise_test.m` (2026-08-18) — confirms the noise excess propagates from raw τ1000 through the GEM step into the actual baroclinic transport signal

Direct follow-up to the BPR-vs-baroclinic component test (BPR ruled out as the noise source) and the earlier raw-`tau1000` finding (sites A and P1 are among the 9 sites' noisiest, A tied-highest in absolute terms, P1 2nd-highest noise-to-signal ratio) — that finding was on the raw acoustic travel time, one step removed from the actual transport calculation. Tested directly on `Gpan` (the GEM-derived geopotential anomaly that `RelTPUD_pilot`/each gap's `RelTPUD_i` are built from) whether the same excess noise survives the GEM lookup-table conversion into the real baroclinic signal, or gets smoothed away by it.

Recomputed `Gpan` for all 9 sites directly from `concat_IESsamba.mat`'s per-site T/S (same `sw_gpan(...)` call as `moc_pilot_v5.m`/`mov_samba_marion*.m`), restricted to the paper's 2013-2017 window. For each of the full array's 8 gaps (A-C, C-D, D-P8, P8-P6, P6-P5, P5-P4, P4-P2, P2-P1) and the Pilot's own A-to-P1 span, built the depth-integral of `(Gpan_east-Gpan_west)` from the surface down to `h_star_pilot=1280dbar` (a proxy for `RelTPUD`'s own cumsum contribution to `MOCup`, without the `1/f` scaling, BPR reference, `AREA_TOPO`, Ekman, or shelf steps — isolating the pure baroclinic/Gpan signal, consistent with the component test above already showing baroclinic-only is the noisier piece), then computed each span's day-to-day noise/signal ratio (`std(diff(signal))/std(signal)`):

| Span | noise/signal (×1000) |
|---|---|
| A-C | 191.5 |
| C-D | 172.9 |
| D-P8 | 179.1 |
| P8-P6 | 153.7 |
| P6-P5 | 126.3 |
| P5-P4 | 158.4 |
| P4-P2 | 178.1 |
| P2-P1 | **350.4** (noisiest of all 9 spans) |
| **A-P1 (Pilot)** | **273.9** (2nd noisiest) |

**Confirms the mechanism directly on the transport-relevant quantity, not just the raw proxy**: mean across the 8 individual gaps = 188.8; the Pilot's A-P1 span (273.9) sits well above that average, and only P2-P1 — the short East gap that also uses site P1 — is noisier. This is the first test performed directly on `Gpan` rather than `tau1000`, and it shows the GEM lookup-table step does **not** damp the excess site-level noise already found in the raw measurement — it passes through essentially intact into the actual baroclinic signal that drives `MOCup`. Reinforces (with the real quantity, not a proxy) the standing explanation for the Pilot's elevated variance and correlation shortfall: the two-endpoint method has no spatial averaging to dilute a single noisy site's contribution, unlike the full array's 8-gap average, and P1 in particular (implicated in both the noisiest individual gap, P2-P1, and the 2nd-noisiest overall span, A-P1) is the main driver.

Outputs: `moc_gpan_noise_test.m`, `moc_gpan_noise_test_IES.png`, `moc_gpan_noise_test.mat`.

## Cross-spectral coherence, full array vs. Pilot: `moc_coherence_test.m` (2026-08-18) — weak, broadband decorrelation, not concentrated in a specific frequency band

Follow-up to the Gpan noise test above, to distinguish two competing explanations for the correlation shortfall: (a) **broadband instrumental/measurement noise** at sites A/P1, which would show up as weak coherence smeared across most of the spectrum (and could in principle be reduced by smoothing/de-spiking), vs. (b) a **specific-band physical difference** — the paper's own explanation for the full array's larger overall variance is that it resolves eddies/Agulhas Rings that a 2-point Pilot structurally cannot see, which would predict strong coherence everywhere *except* that particular mesoscale band (not fixable by better processing, an inherent limitation of only having 2 sites).

Computed magnitude-squared coherence (`mscohere`) and cross-spectral phase (`cpsd`) between the RAW daily full-array and Pilot MOCup anomalies (reusing `moc_anomaly_fig5_raw_2013_2017.mat`'s `full_anom`/`pilot_anom` directly — no NaNs, perfectly regular daily sampling, 1405 days), Welch-averaged with a 256-day Hamming window / 128-day overlap (9 segments, 95% significance threshold = 0.312 by the standard nd-segment formula). Also computed each series' own variance-preserving power spectrum (`pwelch`) to see where each one's variance is concentrated.

**Band-averaged coherence**:

| Band | Mean Cxy | 95% threshold |
|---|---|---|
| Long (>100d, 2 freq bins) | 0.368 | 0.312 |
| Medium (20-100d, 10 bins) | 0.375 | 0.312 |
| Short (2-20d, 116 bins) | 0.258 | 0.312 |

**Result: coherence is weak and mostly non-significant across essentially the whole spectrum**, oscillating above and below the 0.312 threshold almost continuously from ~150 days down to ~3 days, not settling into a clean "high here, low there" pattern. Only the very shortest resolvable periods (2-4 days, near the Nyquist limit) show sustained excursions up to 0.6-0.78 — treated with caution, since that end of the spectrum has the fewest effective degrees of freedom and the noisiest coherence estimates.

**Two supporting findings, both consistent with earlier work**:
1. **Phase lag is ~0 days across almost the entire resolvable spectrum** (only the lowest 1-2 frequency bins, essentially unresolved with 9 segments, show a large but unreliable apparent lag) — an independent, frequency-domain confirmation of `moc_lag_analysis_2013_2017.m`'s time-domain finding (raw whole-record best lag = 0 days).
2. **The power spectra show the full array carrying more variance than the Pilot specifically in the 20-150 day band** — matching the paper's own attribution of its larger variance to eddies/Agulhas Rings — but the two spectra converge in both magnitude and (already-shown) coherence at periods shorter than ~10 days.

**Interpretation**: if the shortfall were purely the "Pilot can't see mesoscale eddies" effect the paper describes, coherence should be strong everywhere *except* the 20-150 day band. Instead coherence is weak nearly everywhere, including well outside that band — more consistent with the Gpan-noise finding above (broadband, site-level measurement noise at A/P1 diluting agreement at essentially every timescale) than with a single specific-band physical limitation. This doesn't rule out the eddy-resolution effect entirely (the power-spectrum asymmetry in 20-150 days is real and expected) — but it's not the dominant, uniquely-identifiable driver of the coherence/correlation shortfall on its own.

Outputs: `moc_coherence_test.m`, `moc_coherence_test_IES.png`, `moc_coherence_test.mat`.

## Outputs

- `concat_IESsamba.mat`, `gpan_samba.mat`, `mov_samba_marion_v15.mat`, `mht_samba_marion_v1.mat` — gitignored (`.mat` files excluded generally; regenerate by running the corresponding script). `mov_samba_marion_v15.mat` (`dt`, `mov`, `mean_term`, `gyre`, `total_direct`, `n_gaps_valid`, `coverage_totalwidth`, `n_sites_valid`, `category`) is new as of `v13` — the save line existed but was commented out in every version before that; `v14` added `low_coverage` (superseded by `v15`'s `category`). `mht_samba_marion_v1.mat` (`dt`, `Heat_total` and its `_EkmanAnomaly`/`_RelativeAnomaly`/`_ReferenceAnomaly`/`_EKMANconst`/`_RelConst`/`_RefConst` variants, `n_sites_valid`, `category`) is new.
- `mov_IES.png` — plot output from `mov_samba_marion*.m`, committed (not gitignored). Since `v15`, shades reduced-coverage periods in 4 graded gray bands by how many of the 9 original sites were present each day.
- `mht_IES.png` — plot output from `mht_samba_marion*.m`, committed (not gitignored). Same graded coverage shading as `mov_IES.png`.
- `moc_streamfunction_v1.mat`, `moc_pilot_v1.mat`, `moc_pilot_lnm_v1.mat`, `moc_pilot_v2.mat`, `moc_streamfunction_v2.mat`, `moc_pilot_v3.mat`, `moc_streamfunction_v3.mat`, `moc_pilot_v4.mat`, `moc_streamfunction_v4.mat`, `moc_pilot_v5.mat` — gitignored. `moc_streamfunction_hovmoller.png`, `moc_streamfunction_mean_profile.png`, `moc_index_IES.png`, `moc_pilot_vs_full_IES.png`, `moc_profile_comparison_IES.png`, `moc_pilot_lnm_vs_bpr_IES.png`, `moc_pilot_v2_vs_full_IES.png`, `moc_pilot_v2_profile_comparison_IES.png`, `moc_streamfunction_v2_hovmoller.png`, `moc_streamfunction_v2_mean_profile.png`, `moc_streamfunction_v2_index_IES.png`, `moc_pilot_v3_vs_full_IES.png`, `moc_pilot_v3_profile_comparison_IES.png`, `moc_streamfunction_v3_hovmoller.png`, `moc_streamfunction_v3_index_IES.png`, `moc_streamfunction_v3_mean_profile.png`, `moc_pilot_v4_vs_full_IES.png`, `moc_pilot_v4_profile_comparison_IES.png`, `moc_streamfunction_v4_hovmoller.png`, `moc_streamfunction_v4_index_IES.png`, `moc_streamfunction_v4_mean_profile.png`, `moc_pilot_v5_vs_full_IES.png`, `moc_pilot_v5_profile_comparison_IES.png` — committed (not gitignored). **`moc_streamfunction_v4.m`/`moc_pilot_v5.m` are the current, most-correct MOC scripts (paper-correct definition + shelf transport + calibration-window fix) — prefer these over earlier versions going forward.**
- `moc_anomaly_fig5_raw_2013_2017.mat`, `moc_pilot_component_test` (no `.mat` saved), `moc_gpan_noise_test.mat`, `moc_coherence_test.mat` — gitignored/not saved. `moc_anomaly_fig5_raw_2013_2017_IES.png`, `moc_pilot_component_test_IES.png`, `moc_gpan_noise_test_IES.png`, `moc_coherence_test_IES.png` — committed (not gitignored).
