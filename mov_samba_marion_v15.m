% Step 2: MOV calculation using Marion's Total Velocities
% between each pair of sites. Her velocities are obtained at
% Marions_code/Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022.mat
%
% Load the Samba-W and Samba-E concatenated on Step 1
% Using only the 9 sites used by Marion: A, C, D, P1, P2 P4 P5 P6 P8
%
% Marion's paper link: https://doi.org/10.1029/2020JC016947

% Article with the definitions of MOV and MOC that I want to decompose: https://doi.org/10.1016/S0074-6142(01)80134-0, https://doi.org/10.1029/2023JC020558
%
% v2: adds the mean and gyre components to the freshwater transport
% decomposition (F = mean + Mov + gyre), following the standard
% decomposition v = V00 + v'(z) + v''(x,z), S = S0 + S'(z) + S''(x,z):
%   mean = net/barotropic transport (V00, S0 cancel algebraically)
%   Mov  = overturning component (already computed below, unchanged)
%   gyre = local (per-site-pair) deviation from the zonal-mean profiles
%
% v3: added a diagnostic (V0_masked) testing whether V0 and wsal using
% different weighting explains the mean+Mov+gyre vs. total_direct
% residual. Result: masking V0 the same way as wsal cut the mean
% residual from 0.3796 Sv to 0.0671 Sv, and the max from 54.4028 Sv to
% 5.7283 Sv -- confirms that was the dominant cause.
%
% v4: adopted the masked V0 as the primary calculation, and added a
% diagnostic checking whether dx1's validity mask varies in time.
% Result: dx1 is time-invariant (dx2 is an exact stand-in), so that
% wasn't the source of the remaining 0.0671/5.7283 Sv residual.
%
% v5: found and fixed the actual remaining cause. For mean+Mov+gyre to
% exactly reconstruct total_direct, the "primed" deviations must
% integrate to zero when weighted by dx2 across depth -- which requires
% V00 and S0 to themselves be the dx2-weighted depth-means of V0 and
% wsal. They weren't: V00 was a plain nanmean(V0) (unweighted across
% depth), and S0 was a double-unweighted mean (unweighted across sites
% via salm, then unweighted across depth), inconsistent with wsal's
% dx1/dx2 weighting. Since dx2 genuinely varies with depth, this left a
% real residual. Both are now defined as dx2-weighted depth-means of V0
% and wsal, which makes the decomposition close to numerical precision.
% Confirmed: residual mean=3.39e-16 Sv, max=8.527e-14 Sv (floating-
% point noise). NOTE: this changed mov's numeric value slightly
% relative to the original mov_samba_marion.m (used for OSM26), since
% mov divides by S0 and centers on V00.
%
% v6: robustness fixes ahead of adding more station data (more sites
% and/or a longer/different time window are both expected to break
% assumptions baked into v5):
%   - clear;clc; at the very start, so stale workspace variables from a
%     prior interactive run can't leak in (matches the convention
%     already used in concat_IES_v4.m/v5.m).
%   - assert(isequal(prew,pree)): pre=prew is used as the depth axis for
%     BOTH West and East sites throughout the script. If West/East
%     pressure grids ever differ, East site masking and the vertical
%     integration would silently misalign. Now fails loudly instead.
%   - The dx1 time-invariance check (v4) is kept permanently rather
%     than removed after one-time use, since dx2 silently goes stale if
%     new station data introduces real time-varying gaps.
%   - Commented the '*AA *B *CC' wildcard clear -- relies on '*B'
%     matching both 'B' and 'BB', which is easy to misread and would
%     silently over-match a future site name ending in 'B' or 'CC'.
%   - The monthly-means arrays (mov_mon/mean_mon/gyre_mon/dtm) are now
%     sized from the actual data span instead of a hardcoded 47, and
%     the first/last-year month boundaries (previously hardcoded
%     init=9/iend=7 tied to jj==1/jj==5) are derived from the real
%     first/last month in dt. Produces identical output on the current
%     2013-09-2017-07 dataset (47 months, same boundaries) but adapts
%     automatically if the time window changes.
%   NOT addressed here: the site list itself (9 sites, 8 gaps) is still
%   hardcoded throughout -- sal00's 3rd dimension, the "for ii=1:8"
%   gap loop, cpiesW([1 5 6])/cpiesE([1 2 4 5 6 8]) index selection,
%   and p_cpiesW/p_cpiesE cutoff arrays all assume exactly this site
%   set. Adding or removing a station requires editing all of those by
%   hand; none of the v6 fixes generalize that. Worth a separate pass
%   if the station list itself is changing, not just the time window.
%
% v7: the v6 assert(isequal(prew,pree)) fired for real -- prew and pree
% are NOT the same length (prew: 501 points, 0-5000 dbar; pree: 531
% points, 0-5300 dbar). Both grids have identical spacing and origin
% (10 dbar steps, starting at 0), so this isn't a real grid mismatch --
% pree is just uncut. `ind=pree<=5000` was already being used to trim
% sal_P*/tem_P* down to 501 points to match prew, but pree itself was
% never trimmed with that same ind, so it stayed at its original 531.
% Fixed by adding pree=pree(ind) right where ind is computed, before it
% gets reassigned to Pressure<=5000 for Marion's TPUD arrays. The
% assert should now pass.
%
% v8: concat_IESsamba.mat's date window was opened (concat_IES_v6.m,
% 2026-08-10) to 2013-09-06--2022-12-11, but Full_Depth_
% OverturningEstimate_Cst_for_MHT.mat (Total_TPUD1..8) still reflects
% the old 2013-2017 window -- README_MOC steps B-E haven't been
% rebuilt yet. Both files have their own 'dt' variable; previously the
% second load silently overwrote the first, which was harmless only
% because both files happened to cover identical dates. Now that they
% don't, this is fixed by keeping both dt's separately and explicitly
% intersecting them, subsetting every time-dependent array (T/S from
% concat_IESsamba, Total_TPUD1..8 from Marion's file) to the common
% dates via that intersection's index arrays -- not by positional
% assumption. Until steps B-E are rebuilt, the aligned dt/output will
% still be the narrower 2013-2017 window (Total_TPUD is the limiting
% side), but the script is now correct/safe for whenever that changes,
% instead of silently mismatching array sizes or crashing.
%
% v9: README_MOC steps B-E have now been rebuilt in Marions_code
% (2026-08-12) -- Total_TPUD1..8 live in a new file,
% Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022.mat, covering
% 2013-09-06 to 2022-12-11 (was 2013-2017). Updated the load path
% below accordingly -- West (samba_w.mat, via concat_IESsamba) is now
% the limiting side, so the aligned dt should come out close to
% concat_IESsamba's own 2013-09-06/2022-12-11 span. Also fixed a
% latent bug: pre_0/difpre (line ~241) hardcoded the day count to 1405
% (the old v8-era dataset size) instead of deriving it from the actual
% data -- harmless while difpre itself was unused downstream, but a
% landmine for whenever that changes. Now sized from n (the real
% number of aligned days).
%
% v10: with the wider 2013-2022 window, the dx1TimeVarying diagnostic
% (added in v6 as a safeguard, never expected to actually fire) fired
% for real -- 22.3% of depth/time/gap points now have a validity mask
% that differs from day 1's (CC has no instrument before 2019-06-27,
% P5/P6 are only ~58% valid throughout, etc.). This broke the
% mean+Mov+gyre vs. total_direct residual check that v5 had closed to
% ~1e-16 Sv: residual mean=-2.24 Sv, max=219.9 Sv. Root cause: dx2 (the
% width profile mov/mean_term are weighted by) was built from day 1's
% coverage mask only and silently reused as if constant across all
% days, while gyre/total_direct use the real per-day mask via dx1
% directly -- fine when coverage genuinely was constant (2013-2017),
% wrong once it isn't. Fixed by making dx2 itself time-resolved
% (dx2=squeeze(nansum(dx1,3)), i.e. actually sum dx1 over gaps at
% every depth AND time, instead of just at t=1) so mov/mean_term are
% weighted by the same real per-day width as gyre/total_direct. The
% dx1TimeVarying warning block is replaced by a plain informational
% print, since coverage varying in time is no longer a source of bias
% now that dx2 tracks it exactly.
%
% v11: v10's fix alone wasn't enough -- residual only dropped to
% mean=0.34 Sv, max=171.3 Sv, still far from the ~1e-16 Sv baseline.
% Two further bugs found and fixed, both exposed by v10 making dx2
% time-resolved rather than the depth-only column vector it always was
% before:
%  1. dx1 was still built from salinity's NaN pattern only. Diagnosed
%     via a debug mask comparison: 1.27M of 13.55M (depth,time,gap)
%     points have valid salinity but NaN velocity (Total_TPUD has its
%     own real coverage gaps, independent of salinity's -- e.g. P5/P6's
%     ~58% coverage). Those points kept their full width in dx2/wsal's
%     denominator while contributing zero to V0's numerator (geo's own
%     NaNs), biasing V0 relative to wsal. Fixed by moving the geo/
%     vel00 construction before dx1, then building dx1 from the UNION
%     of both NaN patterns (isnan(sal00_mid) | isnan(geo)).
%  2. The real culprit: S0=nansum(wsal.*dx2)/nansum(dx2) and
%     V00=nansum(V0.*dx2)/nansum(dx2) used '/' (matrix right division),
%     not './' (elementwise). While dx2 was a plain [depth x 1] column
%     (pre-v10), nansum(dx2) was always a true scalar, where '/' and
%     './' agree -- harmless. Once dx2 became [depth x time] (v10),
%     nansum(dx2) became a [1 x time] row too, and MATLAB's A/B between
%     two same-size non-scalar arrays is a least-squares matrix
%     division, not elementwise -- this silently collapsed S0 (and
%     V00) to a SINGLE best-fit scalar for the entire record (confirmed
%     via size(S0) == [1 1], not [1 3381]) instead of one value per
%     day. Every downstream quantity (sal_prime, vel_prime, mov,
%     mean_term) was computed against this one wrong constant. Fixed by
%     changing both to './'. Also fixed mean_term's sum(dx2) ->
%     nansum(dx2) while in here: this repo's custom nansum.m (functions/
%     nan/) returns NaN, not 0, when an entire reduced dimension is
%     NaN, so a plain sum() would have propagated NaN across a whole
%     day the moment any single depth level had zero valid gaps --
%     didn't fire on this dataset (0 NaN days either way) but was
%     inconsistent with every other reduction in the script, all of
%     which already use nansum.
% Confirmed: residual mean=-1.062e-14 Sv, max=3.126e-13 Sv -- back to
% the v5-era floating-point-noise baseline, now on the full 2013-2022,
% 3381-day dataset.
%
% v12: fixes the AREA_TOPO/salinity-cutoff inconsistency flagged in
% CLAUDE.md. Total_TPUD (used through v11 for vel00/geo) already has
% Step E's per-gap AREA_TOPO_i topographic open-area fraction baked in
% (Total_TPUD = true_velocity*dx(gap)*AREA_TOPO(z)), but dx1 (the width
% V0/wsal/dx2 are weighted by) was always the flat, nominal dx(gap) --
% so at any depth where a gap is topographically ~blocked, its full
% nominal width still counted toward the zonal average's denominator
% while contributing ~nothing to the numerator, diluting V0 toward
% zero. Quantified per gap ("phantom width", fraction of the nominal
% width within each gap's salinity-valid depth range that's actually
% topographically blocked): P2-P1 64.8%, P4-P2 21.9%, A-C 24.7%, D-P8
% 15.6%; C-D/P8-P6/P6-P5/P5-P4 negligible (0-4%).
%
% Fixed by switching to Absolute_TPUD1..8_preTopo (Step E v2's new
% output: the velocity BEFORE the AREA_TOPO multiplication) for vel00/
% geo, and building dxA = dx1 .* AREA_TOPO_gap(z) -- the REAL per-gap
% open width -- used everywhere dx1 previously weighted an average
% (dx2, wsal, V0, aux_gyre, aux_total). dx1 itself is kept only for
% mask construction and the coverage-variability print. Re-derived the
% v5 cross-term-cancellation algebra with these new definitions from
% scratch (not just assumed it still holds) -- it does, unchanged in
% structure, since it only depends on dx2=sum_gap(dxA) and V0/wsal
% being *weighted-mean-of-something-times-dxA*, not on what dxA
% physically represents. A sanity check recomputing vel00.*dxA and
% comparing against the loaded Total_TPUD (which should match exactly,
% since Total_TPUD = true_velocity*dx(gap)*AREA_TOPO(z) = vel00*dxA by
% construction) is printed below as a safeguard against the coding
% getting this backwards.
%
% Two bugs found and fixed while validating v12 against that sanity
% check and the residual:
%  1. V0 was `nansum(geo,3)./nansum(dxA,3)` -- geo alone
%     (Absolute_TPUD_preTopo) has no AREA_TOPO factor at all, so this
%     mismatched a non-topo-weighted numerator against a
%     topo-weighted denominator. Residual blew back up to mean=-14.49
%     Sv, max=81.68 Sv. Fixed to `nansum(vel00.*dxA,3)./nansum(dxA,3)`
%     -- consistent with how the sanity check and aux_total/aux_gyre
%     already used vel00.*dxA.
%  2. Absolute_TPUD_preTopo is captured in Step E *before* the
%     AREA_TOPO multiplication, which is itself before Ekman is added
%     -- so it excludes Ekman entirely, at every row, not just where
%     topography matters. Confirmed AREA_TOPO_i==1 at all 7 Ekman rows
%     (0-60dbar) for all 8 gaps, so fixed by having Step E v2 also
%     save Ekman_TPUD1..8 and adding it straight into geo/vel00 at
%     those 7 rows here (no extra topo-scaling needed there).
% Confirmed: residual mean=-8.826e-16 Sv, max=1.705e-13 Sv -- matches
% v11's floating-point-noise baseline, sanity check ~0 at every row
% including the Ekman ones.
%
% v12 vs v11 impact on the actual output: mean_term essentially
% unchanged (diff ~1e-15 Sv); mov changes modestly (mean|diff|=0.007
% Sv, max|diff|=0.12 Sv on 2018-09-27, corr=0.997 against v11); gyre
% barely moves (mean|diff|=0.001 Sv, corr=0.9998). The fix is a real,
% necessary correction (removes a genuine bias source, see CLAUDE.md's
% AREA_TOPO note), but its numerical footprint on the published mov
% series is modest rather than dramatic.
%
% v13: addresses the OTHER open item in CLAUDE.md -- the 2017-08/
% 2019-10 window where mov's variance drops ~40% because P5/P6/P4/P8
% outages left 3 of 8 gaps at 0% coverage (see CLAUDE.md's "Known
% limitation" section for the full analysis). Not fixable (the data
% simply isn't there), but now flagged automatically: two new per-day
% diagnostics, n_gaps_valid (how many of the 8 gaps have any valid
% depth that day) and coverage_totalwidth (nansum(dx2), the same total
% open width mean_term is weighted by), saved alongside mov/mean_term/
% gyre/total_direct -- previously this script computed everything but
% never actually saved it (the save line was commented out through
% v12). Downstream users can now screen low-coverage days
% programmatically instead of re-deriving the 2017-2019 window by hand
% each time.
%
% v14 (2026-08-13): visualizes the low-coverage window on the mov_IES
% plot directly, instead of either (a) leaving it unmarked (v13 -- the
% period just looks quiet, easy to misread as a real signal) or (b)
% masking it to NaN (tried as a one-off diagnostic script,
% mov_samba_marion_v13_9pies.m, not kept -- requiring all 9 original
% sites simultaneously present leaves a blank ~4-year hole, 2017-08/
% 2021-09, and throws away the ~54% of days that DO have partial,
% still-informative coverage). Adds a semi-transparent gray shaded
% background behind the mov line wherever a 30-day rolling mean of
% coverage_totalwidth drops below 50% of the record's median -- a
% single, stable stretch (~Oct 2017-Oct 2019, insensitive to the exact
% threshold in the 30-60% range) that's narrower than the full P5/P6
% outage (2017-08/2021-09) because coverage_totalwidth recovers enough
% once P4/P8 return in Nov 2019, even though P5/P6 themselves don't
% return until 2021. mov itself is NOT modified or masked -- every
% value is the same real, fully-computed estimate as v13, just visually
% flagged as lower-confidence where the array was reduced.
%
% v15 (2026-08-13): replaces v14's single binary low-coverage flag with
% graded shading bands by n_sites_valid (how many of the 9 original
% sites are present each day) -- user's follow-up question: is there a
% minimum site count that gives a more informative picture than just
% "low" vs. "not low"? A one-off test (test_minsites.m, not kept) swept
% thresholds of >=6/7/8/9 sites and found >=7 almost exactly reproduces
% v14's coverage_totalwidth-based stretch (a good cross-check), while
% >=8 and >=9 are nearly identical to each other (P5 and P6 are usually
% both missing or both present) and >=6 keeps 90.7% of days with only
% two short gaps. A second test (test_colorband.m, not kept) turned
% this into 4 shading bands (9/7-8/6/<6 sites) instead of a single
% threshold, revealing a previously unflagged patch of reduced (7-8/9)
% coverage in 2020-2022 that v14's single stretch missed. mov itself is
% still NOT modified or masked, same as v14 -- only the shading changed.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')

load concat_IESsamba
dt_TS=dt; clear dt

% Keep only A, C, D (drop AA, B, BB, CC) -- relies on '*B' matching
% both 'B' and 'BB' (both end in the letter B), and '*CC' matching only
% CC (leaving C alone, since C doesn't end in CC). Fragile: a future
% site name ending in 'B' or 'CC' would be swept up too.
clear *AA *B *CC

% --------------------------------------------------------
% These calculations are for the OSM26 presentation and
% use the velocity fields estimated by Marion's paper in 2021.

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2']);
dt_TPUD=dt; clear dt

% --------------------------------------------------------
% Align the T/S dates (concat_IESsamba) with Total_TPUD's dates
% (Marion's file) by explicit intersection rather than assuming the
% two arrays are already the same length/order. See v8 note above.
[dt,i_TS,i_TPUD]=intersect(dt_TS,dt_TPUD);
fprintf('Aligned dt: %d common days (%s to %s) out of %d (T/S) and %d (Total_TPUD)\n', ...
    numel(dt),datestr(dt(1)),datestr(dt(end)),numel(dt_TS),numel(dt_TPUD));

sal_A=sal_A(:,i_TS);   tem_A=tem_A(:,i_TS);
sal_C=sal_C(:,i_TS);   tem_C=tem_C(:,i_TS);
sal_D=sal_D(:,i_TS);   tem_D=tem_D(:,i_TS);
sal_P1=sal_P1(:,i_TS); tem_P1=tem_P1(:,i_TS);
sal_P2=sal_P2(:,i_TS); tem_P2=tem_P2(:,i_TS);
sal_P4=sal_P4(:,i_TS); tem_P4=tem_P4(:,i_TS);
sal_P5=sal_P5(:,i_TS); tem_P5=tem_P5(:,i_TS);
sal_P6=sal_P6(:,i_TS); tem_P6=tem_P6(:,i_TS);
sal_P8=sal_P8(:,i_TS); tem_P8=tem_P8(:,i_TS);

% v12: align both Total_TPUD (kept only for the vel00.*dxA sanity check
% below) and Absolute_TPUD*_preTopo (the actual input to vel00/geo now).
Total_TPUD1=Total_TPUD1(:,i_TPUD);
Total_TPUD2=Total_TPUD2(:,i_TPUD);
Total_TPUD3=Total_TPUD3(:,i_TPUD);
Total_TPUD4=Total_TPUD4(:,i_TPUD);
Total_TPUD5=Total_TPUD5(:,i_TPUD);
Total_TPUD6=Total_TPUD6(:,i_TPUD);
Total_TPUD7=Total_TPUD7(:,i_TPUD);
Total_TPUD8=Total_TPUD8(:,i_TPUD);

Absolute_TPUD1_preTopo=Absolute_TPUD1_preTopo(:,i_TPUD);
Absolute_TPUD2_preTopo=Absolute_TPUD2_preTopo(:,i_TPUD);
Absolute_TPUD3_preTopo=Absolute_TPUD3_preTopo(:,i_TPUD);
Absolute_TPUD4_preTopo=Absolute_TPUD4_preTopo(:,i_TPUD);
Absolute_TPUD5_preTopo=Absolute_TPUD5_preTopo(:,i_TPUD);
Absolute_TPUD6_preTopo=Absolute_TPUD6_preTopo(:,i_TPUD);
Absolute_TPUD7_preTopo=Absolute_TPUD7_preTopo(:,i_TPUD);
Absolute_TPUD8_preTopo=Absolute_TPUD8_preTopo(:,i_TPUD);

% v12: Ekman_TPUD1..8 (7 rows, 0-60dbar -- no pressure trim needed,
% already that size) -- needed to add Ekman back into geo/vel00 below,
% since Absolute_TPUD_preTopo excludes it (see v12 header note).
Ekman_TPUD1=Ekman_TPUD1(:,i_TPUD);
Ekman_TPUD2=Ekman_TPUD2(:,i_TPUD);
Ekman_TPUD3=Ekman_TPUD3(:,i_TPUD);
Ekman_TPUD4=Ekman_TPUD4(:,i_TPUD);
Ekman_TPUD5=Ekman_TPUD5(:,i_TPUD);
Ekman_TPUD6=Ekman_TPUD6(:,i_TPUD);
Ekman_TPUD7=Ekman_TPUD7(:,i_TPUD);
Ekman_TPUD8=Ekman_TPUD8(:,i_TPUD);

clear dt_TS dt_TPUD i_TS i_TPUD

% --------------------------------------------------------
% --------------------------------------------------------
% Building the salinity matrices
% Cut the pressure at 5000 m (using Samba_W limit)

ind=pree<=5000;
pree=pree(ind);
sal_P1=sal_P1(ind,:);
sal_P2=sal_P2(ind,:);
sal_P4=sal_P4(ind,:);
sal_P5=sal_P5(ind,:);
sal_P6=sal_P6(ind,:);
sal_P8=sal_P8(ind,:);

tem_P1=tem_P1(ind,:);
tem_P2=tem_P2(ind,:);
tem_P4=tem_P4(ind,:);
tem_P5=tem_P5(ind,:);
tem_P6=tem_P6(ind,:);
tem_P8=tem_P8(ind,:);
% --------------------------------------------------------
% These are Marion's tranport per unit depth (TPUD) which is
% the absolute velocity times dx.
ind=Pressure<=5000;
Total_TPUD1=Total_TPUD1(ind,:);
Total_TPUD2=Total_TPUD2(ind,:);
Total_TPUD3=Total_TPUD3(ind,:);
Total_TPUD4=Total_TPUD4(ind,:);
Total_TPUD5=Total_TPUD5(ind,:);
Total_TPUD6=Total_TPUD6(ind,:);
Total_TPUD7=Total_TPUD7(ind,:);
Total_TPUD8=Total_TPUD8(ind,:);

% v12: same <=5000dbar trim for Absolute_TPUD*_preTopo and for
% AREA_TOPO_i (static per-gap profile, no time dimension).
Absolute_TPUD1_preTopo=Absolute_TPUD1_preTopo(ind,:);
Absolute_TPUD2_preTopo=Absolute_TPUD2_preTopo(ind,:);
Absolute_TPUD3_preTopo=Absolute_TPUD3_preTopo(ind,:);
Absolute_TPUD4_preTopo=Absolute_TPUD4_preTopo(ind,:);
Absolute_TPUD5_preTopo=Absolute_TPUD5_preTopo(ind,:);
Absolute_TPUD6_preTopo=Absolute_TPUD6_preTopo(ind,:);
Absolute_TPUD7_preTopo=Absolute_TPUD7_preTopo(ind,:);
Absolute_TPUD8_preTopo=Absolute_TPUD8_preTopo(ind,:);

AREA_TOPO_1=AREA_TOPO_1(ind); AREA_TOPO_2=AREA_TOPO_2(ind);
AREA_TOPO_3=AREA_TOPO_3(ind); AREA_TOPO_4=AREA_TOPO_4(ind);
AREA_TOPO_5=AREA_TOPO_5(ind); AREA_TOPO_6=AREA_TOPO_6(ind);
AREA_TOPO_7=AREA_TOPO_7(ind); AREA_TOPO_8=AREA_TOPO_8(ind);
% --------------------------------------------------------
pre=prew;
assert(isequal(prew,pree), ...
    'prew and pree (West/East pressure grids) differ -- pre=prew is used as the depth axis for both, so a mismatch would silently misalign East site masking and the vertical integration.');
clear prew Pressure ind
% --------------------------------------------------------
% Only consider the values up to each site's depth: values used by Marion.
p_cpiesW=[1370,4620,4850];
p_cpiesE=[1280,2150,4560,5060,5290,4690];

sal_A(find(pre>p_cpiesW(1)),:)=nan;
sal_C(find(pre>p_cpiesW(2)),:)=nan;
sal_D(find(pre>p_cpiesW(3)),:)=nan;
sal_P8(find(pre>p_cpiesE(1)),:)=nan;
sal_P6(find(pre>p_cpiesE(2)),:)=nan;
sal_P5(find(pre>p_cpiesE(3)),:)=nan;
sal_P4(find(pre>p_cpiesE(4)),:)=nan;
sal_P2(find(pre>p_cpiesE(5)),:)=nan;
sal_P1(find(pre>p_cpiesE(6)),:)=nan;

% --------------------------------------------------------
% Estimate the section mean salinity (So). Making a continuous array along
% the basin, from 1 to 9, relative to A, C, D, P8, P6, P5, P4, P2, P1.
[m,n]=size(sal_A);
sal00=nan*ones([m n 9]);
sal00(:,:,1)=sal_A;
sal00(:,:,2)=sal_C;
sal00(:,:,3)=sal_D;
sal00(:,:,4)=sal_P8;
sal00(:,:,5)=sal_P6;
sal00(:,:,6)=sal_P5;
sal00(:,:,7)=sal_P4;
sal00(:,:,8)=sal_P2;
sal00(:,:,9)=sal_P1;

% v15: per-day count of how many of the 9 original sites
% (A,C,D,P8,P6,P5,P4,P2,P1) have at least one valid depth -- used below
% to shade the plot in graded bands (9/7-8/6/<6 sites) instead of v14's
% single binary low-coverage flag.
n_sites_valid=reshape(sum(squeeze(any(~isnan(sal00),1)),2),1,n);

clear sal_* tem_*
% --------------------------------------------------------
% --------------------------------------------------------
% Estimate the distance between sites to divide the TPUD and get the velocity
position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);
% --------------------------------------------------------
% Salinity at mid-point

col01=squeeze(sal00(:,:,1:8));
col02=squeeze(sal00(:,:,2:9));
mcol(:,:,:,1)=col01;
mcol(:,:,:,2)=col02;
sal00_mid=nanmean(mcol,4);

% mean sal in dz
salm=(sal00_mid(1:end-1,:,:)+sal00_mid(2:end,:,:))/2;
pre_0=repmat(pre,[1 n 9]);
difpre=diff(pre_0,1);
[l,m,n]=size(sal00_mid);

% Do the same for velocity. Find the barotropic component first. Built
% here, BEFORE dx1, so its own NaN pattern (real coverage gaps) can be
% folded into dx1's mask below (v11 fix, see note).
% v12: geo/vel00 now come from Absolute_TPUD*_preTopo (Step E v2's
% pre-AREA_TOPO output), not Total_TPUD -- vel00 is now the TRUE
% (topo-unattenuated) velocity; the AREA_TOPO correction is applied
% explicitly below via dxA instead of being silently baked into geo.
vel00=nan*ones([l m n]);
geo=nan*ones([l m n]);
for ii=1:8
str=['vel00(:,:,',num2str(ii),')=Absolute_TPUD',num2str(ii),'_preTopo/dx(',num2str(ii),');'];
eval(str)
str=['geo(:,:,',num2str(ii),')=Absolute_TPUD',num2str(ii),'_preTopo;'];
eval(str)
end

% v12: add Ekman back in at the 7 surface rows (0-60dbar) -- excluded
% from Absolute_TPUD_preTopo by construction (added even later than
% the AREA_TOPO step in Step E). AREA_TOPO_i==1 at these rows for all
% 8 gaps (confirmed directly against topo_corr_msm60.mat), so adding
% Ekman_TPUD straight into geo (transport) and geo/dx(ii) into vel00
% (velocity) needs no extra topo-scaling to exactly reproduce what
% Total_TPUD has there -- verified by the sanity check further below.
for ii=1:8
str=['geo(1:7,:,',num2str(ii),')=geo(1:7,:,',num2str(ii),')+Ekman_TPUD',num2str(ii),';'];
eval(str)
str=['vel00(1:7,:,',num2str(ii),')=geo(1:7,:,',num2str(ii),')/dx(',num2str(ii),');'];
eval(str)
end
clear Ekman_TPUD1 Ekman_TPUD2 Ekman_TPUD3 Ekman_TPUD4 Ekman_TPUD5 Ekman_TPUD6 Ekman_TPUD7 Ekman_TPUD8

clear Absolute_TPUD*_preTopo

% Make a 3D dx to multiply by the salinity'
dx1=repmat(reshape(dx,1,1,[]),l,m,1);

% v12: AREA_TOPO_i (static per-gap open-area fraction, [pre x 1] each)
% stacked into a 3D array matching dx1's shape, broadcast across time.
AREA_TOPO_3d=nan*ones([l 1 n]);
for ii=1:8
str=['AREA_TOPO_3d(:,1,',num2str(ii),')=AREA_TOPO_',num2str(ii),'(:);'];
eval(str)
end
AREA_TOPO_3d=repmat(AREA_TOPO_3d,[1 m 1]);
clear AREA_TOPO_1 AREA_TOPO_2 AREA_TOPO_3 AREA_TOPO_4 AREA_TOPO_5 AREA_TOPO_6 AREA_TOPO_7 AREA_TOPO_8

% v11: dx1's validity mask must be the UNION of salinity's NaN pattern
% AND velocity's (geo's) -- v10 only used sal00_mid's, then
% one-directionally re-masked geo by it (geo(isnan(dx1))=nan below).
% That's fine only if velocity is NaN exclusively where salinity
% already is. Diagnosed on the 2013-2022 dataset: 1.27M of 13.55M
% (depth,time,gap) points have valid salinity but NaN velocity (real,
% independent coverage gaps, e.g. P5/P6's ~58% coverage) -- those
% points still counted their full width in dx2/wsal's denominator while
% contributing nothing to V0's numerator, biasing V0 relative to wsal
% and breaking the mean+Mov+gyre vs. total_direct identity (residual up
% to 219.9 Sv even after v10's dx2 time-resolution fix alone). Fixed by
% combining both NaN patterns before anything downstream (wsal, dx2,
% V0) is computed from dx1. v12 additionally folds in AREA_TOPO_3d's
% own NaN pattern (should be empty in practice, but not assumed).
ind=find(isnan(sal00_mid) | isnan(geo) | isnan(AREA_TOPO_3d));
dx1(ind)=nan;
geo(isnan(dx1))=nan;

% v12: dxA is the REAL per-gap open width (dx(gap)*AREA_TOPO(z)),
% used everywhere dx1 previously weighted an average (dx2, wsal, V0,
% aux_gyre, aux_total below) -- see v12 header note for the full
% motivation and a from-scratch re-derivation confirming this preserves
% the mean+Mov+gyre==total_direct identity exactly. dx1 itself is kept
% only for the coverage-variability print below.
dxA=dx1.*AREA_TOPO_3d;

% dx2: total valid OPEN width profile, summed across the 9 gaps -- now
% time-resolved ([depth x time], v10) AND topography-weighted (v12),
% built directly from dxA's real per-day mask rather than a day-1
% stand-in or a flat nominal width. S0/V00/mov/mean_term below are all
% weighted by dx2, and gyre/total_direct are weighted by dxA directly
% -- both sides now use the exact same real per-day, per-gap open width.
dx2=squeeze(nansum(dxA,3));

% Informational only (not a bias risk since v10/v11 -- dx2 now tracks
% dx1's real, combined per-day mask exactly): how much the coverage
% mask actually varies across days, for context.
mask=isnan(dx1);
mask_t1=repmat(mask(:,1,:),[1 size(mask,2) 1]);
frac_diff=mean(mask(:)~=mask_t1(:));
fprintf('Site/depth coverage mask varies across days for %.4f%% of depth,time,gap points (dx2 is time-resolved, so this is not a source of bias).\n', frac_diff*100);

% Distance weighted salinity mean for each depth
% v12: dx1 -> dxA (real open width, not nominal).
wsal=nansum(sal00_mid.*dxA,3)./nansum(dxA,3);

% Mean salinity time series: dx2-weighted depth-mean of wsal, so that
% sal_prime=wsal-S0 integrates to zero when weighted by dx2 across
% depth -- required for mean+Mov+gyre to reconstruct total_direct
% exactly (v5 fix).
% v11: / -> ./. With dx2 now [depth x time] (v10) instead of the old
% [depth x 1], nansum(dx2) is itself [1 x time], not a scalar -- and
% MATLAB's A/B between two same-size non-scalar arrays is matrix right
% division (a least-squares fit), not elementwise. That silently
% collapsed S0 to a single best-fit scalar for the whole record
% instead of one value per day (confirmed: size(S0) was [1 1], not
% [1 3381]) -- the actual cause of the still-large residual after v10
% alone. Harmless before v10 only because dx2 being a column vector
% always made nansum(dx2) a true scalar, where / and ./ agree.
S0=nansum(wsal.*dx2)./nansum(dx2);

% --------------------------------------------------------
% Estimate Sal'=sal00-S0, subtract S0 and zonally mean to get the
% the <S(z)>.
sal_prime=nanmean(wsal-S0,3);
% --------------------------------------------------------
% V0: mean zonal velocity as function of depth, masked/weighted the
% same way as wsal (nansum(dxA,3) instead of the constant sum(dx)).
% V00: dx2-weighted depth-mean of V0, for the same reason S0 was
% redefined above (v5 fix). vel_prime is the baroclinic (overturning)
% component.
% v12: dx1 -> dxA, AND geo -> vel00.*dxA. geo alone (Absolute_TPUD_
% preTopo = true_velocity*dx(gap)) has no AREA_TOPO factor at all, so
% nansum(geo,3)/nansum(dxA,3) would mismatch a topo-weighted
% denominator against a NON-topo-weighted numerator. vel00.*dxA =
% true_velocity*dx(gap)*AREA_TOPO(z), the actual topo-corrected
% transport (equal to Total_TPUD, verified by the sanity check below)
% -- summing THAT across gaps, divided by nansum(dxA,3), gives a V0
% consistently weighted by the real open width on both sides.
V0=nansum(vel00.*dxA,3)./nansum(dxA,3);
V00=nansum(V0.*dx2)./nansum(dx2);  % v11: / -> ./, same reason as S0 above
vel_prime=V0-V00;

% Integrating <v(z)><S(z)> vertically
aux=vel_prime.*sal_prime*mean(diff(pre));

% as the ocean width varies with depth (and, since v10, with time too):
mov=-nansum(aux.*dx2)/1e6./S0;

% --------------------------------------------------------
% Mean (net/barotropic) component.
% F_mean = -(1/S0)*V00*S0*A = -V00*A -- S0 cancels algebraically, so it
% does not appear below. Uses the same dx2 width as mov (now
% depth-and-time-resolved, v10), for consistency with that term.
% v11: sum(dx2) -> nansum(dx2). This repo's custom nansum.m (functions/
% nan/) returns NaN, not 0, when an ENTIRE reduced dimension is NaN --
% so at any depth level where all 8 gaps are invalid that day, dx2(z,t)
% itself is NaN (not 0). Every other reduction in this script already
% uses nansum for exactly this reason (S0, V00, mov, gyre, total_direct
% all do); this one used plain sum, so it would silently turn the
% WHOLE day's mean_term into NaN as soon as even one depth level had
% zero valid gaps -- increasingly likely now that P5/P6 (~58% coverage)
% and West's depth-cutoff-narrowed gaps overlap at the deepest levels.
mean_term=-V00.*nansum(dx2)*mean(diff(pre))/1e6;

% --------------------------------------------------------
% Gyre component: local (per-site-pair) deviation from the zonal-mean
% profiles V0(z,t) and wsal(z,t). Needs per-gap resolution, so it uses
% dxA (per-gap, per-time, topography-weighted width, v12) rather than
% dx2.
vel_pp=vel00-V0;          % v''_i(z,t) = per-gap velocity - zonal-mean profile
sal_pp=sal00_mid-wsal;    % S''_i(z,t) = per-gap salinity - zonal-mean profile
aux_gyre=vel_pp.*sal_pp.*dxA*mean(diff(pre));  % v12: dx1 -> dxA
gyre=-nansum(nansum(aux_gyre,3),1)/1e6./S0;

% --------------------------------------------------------
% Consistency check: a directly-computed total (no decomposition, using
% the real per-gap v and S, weighted by dxA, v12) vs. the sum of the
% three components above. Should be close to zero (v5 fix).
aux_total=vel00.*sal00_mid.*dxA*mean(diff(pre));  % v12: dx1 -> dxA
total_direct=-nansum(nansum(aux_total,3),1)/1e6./S0;
residual=total_direct-(mean_term+mov+gyre);
fprintf('Residual (total_direct - [mean+Mov+gyre]): mean=%.4g Sv, max|.|=%.4g Sv\n', ...
        nanmean(residual), nanmax(abs(residual)));

% v12 sanity check: vel00.*dxA should exactly reconstruct Total_TPUD at
% EVERY row now, including the 7 Ekman rows (0-60dbar) -- geo/vel00
% were augmented with Ekman_TPUD there above, and AREA_TOPO_i==1 at
% those rows for all 8 gaps, so no extra topo-scaling is needed for the
% reconstruction to hold. A mismatch anywhere would mean the vel00/dxA/
% Ekman plumbing above has a bug.
recon_err=nan(1,8);
for ii=1:8
    recon=squeeze(vel00(:,:,ii)).*squeeze(dxA(:,:,ii));
    actual=eval(['Total_TPUD',num2str(ii)]);
    recon_err(ii)=nanmax(abs(recon(:)-actual(:)));
end
fprintf('Sanity check (vel00.*dxA vs Total_TPUD, should be ~0 at ALL rows): max abs diff per gap = %s\n', mat2str(recon_err,4));
clear Total_TPUD1 Total_TPUD2 Total_TPUD3 Total_TPUD4 Total_TPUD5 Total_TPUD6 Total_TPUD7 Total_TPUD8 recon recon_err ii

% --------------------------------------------------------
% Coverage diagnostics (v13): per-day summaries of how much of the
% array actually had data that day, so downstream users can screen for
% low-coverage periods (like 2017-08/2019-10, when P5/P6/P8/P4 outages
% dropped 3 of 8 gaps to 0% coverage and suppressed mov's variance by
% ~40% -- see CLAUDE.md's "Known limitation" section) automatically
% instead of re-deriving it by hand each time.
%   n_gaps_valid(t)        -- how many of the 8 gaps have at least one
%                              valid depth on day t (0-8).
%   coverage_totalwidth(t) -- total open cross-sectional width (dxA,
%                              already topo-weighted, v12) summed over
%                              depth on day t, in meters. The same
%                              quantity mean_term is weighted by
%                              (nansum(dx2)); low values flag days
%                              where mov/mean_term/gyre are being
%                              computed from a much-reduced section.
% Explicit reshape to [1 x time], matching mov/mean_term/gyre's own
% row-vector convention (rather than relying on squeeze's orientation,
% which is ambiguous for a size-1-in-two-dims input).
n_gaps_valid=reshape(sum(any(~isnan(dxA),1),3),1,m);
coverage_totalwidth=reshape(nansum(dx2,1),1,m);

% --------------------------------------------------------
% v15: replaces v14's single binary low-coverage flag with graded
% shading bands by n_sites_valid (how many of the 9 original sites are
% present that day) -- a follow-up test (test_minsites.m, not kept)
% found that >=7/9 sites almost exactly reproduces v14's
% coverage_totalwidth-based low-coverage stretch (a good cross-check),
% and that useful additional structure exists below that (6/9 vs <6/9)
% that a single binary flag couldn't show, including a previously
% unflagged patch of reduced (7-8/9) coverage in 2020-2022. Bands
% (from that test's actual day counts): 9/9 (54.4%, unshaded), 7-8/9
% (23.7%, light gray), 6/9 (12.5%, medium gray), <6/9 (9.3%, dark gray).
% mov itself is still NOT modified or masked, same as v14 -- only the
% shading changed.
band_colors={[1 1 1],[0.88 0.88 0.88],[0.72 0.72 0.72],[0.5 0.5 0.5]};
band_labels={'9/9 sites','7-8/9 sites','6/9 sites','<6/9 sites'};

category=nan(size(n_sites_valid));
category(n_sites_valid==9)=1;
category(n_sites_valid>=7 & n_sites_valid<=8)=2;
category(n_sites_valid==6)=3;
category(n_sites_valid<6)=4;

fprintf('Coverage categories:\n');
for cc=1:4
    fprintf('  %s: %d of %d days (%.1f%%)\n', band_labels{cc}, sum(category==cc), numel(category), 100*mean(category==cc));
end

% --------------------------------------------------------
% Some plottings
mov_low=lowpass_filter(mov,45);
[nmov,coef]=sinfitb_tot(dt,mov,365.25);
slope=coef(1)+coef(2)*dt;

figure
hold on
ylim0=[min(mov) max(mov)];
ylim0=ylim0+[-1 1]*0.05*diff(ylim0);
legend_h=[]; legend_str={};
for cc=2:4  % skip category 1 (==9 sites, unshaded)
    iscat=(category==cc);
    d=diff([0 iscat 0]);
    starts=find(d==1); ends=find(d==-1)-1;
    for kk=1:numel(starts)
        xpatch=[dt(starts(kk)) dt(ends(kk)) dt(ends(kk)) dt(starts(kk))];
        ypatch=[ylim0(1) ylim0(1) ylim0(2) ylim0(2)];
        hh=fill(xpatch,ypatch,band_colors{cc},'EdgeColor','none');
        if kk==1
            legend_h(end+1)=hh; %#ok<SAGROW>
            legend_str{end+1}=band_labels{cc}; %#ok<SAGROW>
        end
    end
end
p=plot(dt,mov,dt,slope,'k',dt,mov_low,'r');
set(p(2),'linewidth',1.5)
set(p(1),'linewidth',1.5)
set(p(3),'linewidth',1.5)
title('Freshwater transport (Mov) at 34.5^oS from IES','fontsize',14)
ylabel('Mov (Sv)')
datetick('x',12);grid;axis('tight')
ylim(ylim0)
set(gca,'linewidth',1,'fontsize',14)
box on
legend(legend_h,legend_str,'Location','southoutside','Orientation','horizontal','fontsize',8)

print -dpng -r150  mov_IES

% v13: enabled (was commented out through v12) -- now includes the
% coverage diagnostics above so downstream users can screen for
% low-coverage periods without re-deriving n_gaps_valid/
% coverage_totalwidth by hand. v15 adds n_sites_valid and category (the
% shading bands used in the plot above) to the saved output too.
save('mov_samba_marion_v15.mat','dt','mov','mean_term','gyre','total_direct', ...
    'n_gaps_valid','coverage_totalwidth','n_sites_valid','category')

% --------------------------------------------------------
% Estimate the monthly means. Array lengths and the first/last-year
% month boundaries are derived from the actual span of dt instead of
% hardcoded (v6) -- adapts automatically to whatever the aligned dt
% turns out to be (v8/v9).

[yy,mm,dd]=datevec(dt);
year=unico(yy);
nyears=length(year);

first_month=mm(find(yy==year(1),1,'first'));
last_month=mm(find(yy==year(end),1,'last'));

nmonths_per_year=12*ones(nyears,1);
if nyears==1
    nmonths_per_year(1)=last_month-first_month+1;
else
    nmonths_per_year(1)=12-first_month+1;
    nmonths_per_year(end)=last_month;
end
total_months=sum(nmonths_per_year);

mov_mon=nan*ones(total_months,1);
mean_mon=nan*ones(total_months,1);
gyre_mon=nan*ones(total_months,1);
dtm=nan*ones(total_months,1);
count=1;
for jj=1:nyears
if jj==1
init=first_month;
else
init=1;
end
if jj==nyears
iend=last_month;
else
iend=12;
end
for kk=init:iend
ind=yy==year(jj)&mm==kk;
mov_mon(count)=mean(mov(ind));
mean_mon(count)=mean(mean_term(ind));
gyre_mon(count)=mean(gyre(ind));
dtm(count)=datenum(year(jj),kk,15);
count=count+1;
end
end
