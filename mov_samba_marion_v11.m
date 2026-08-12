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
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022']);
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

Total_TPUD1=Total_TPUD1(:,i_TPUD);
Total_TPUD2=Total_TPUD2(:,i_TPUD);
Total_TPUD3=Total_TPUD3(:,i_TPUD);
Total_TPUD4=Total_TPUD4(:,i_TPUD);
Total_TPUD5=Total_TPUD5(:,i_TPUD);
Total_TPUD6=Total_TPUD6(:,i_TPUD);
Total_TPUD7=Total_TPUD7(:,i_TPUD);
Total_TPUD8=Total_TPUD8(:,i_TPUD);

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
% here, BEFORE dx1, so its own NaN pattern (real Total_TPUD coverage
% gaps) can be folded into dx1's mask below (v11 fix, see note).
vel00=nan*ones([l m n]);
geo=nan*ones([l m n]);
for ii=1:8
str=['vel00(:,:,',num2str(ii),')=Total_TPUD',num2str(ii),'/dx(',num2str(ii),');'];
eval(str)
str=['geo(:,:,',num2str(ii),')=Total_TPUD',num2str(ii),';'];
eval(str)
end

clear Total_*

% Make a 3D dx to multiply by the salinity'
dx1=repmat(reshape(dx,1,1,[]),l,m,1);

% v11: dx1's validity mask must be the UNION of salinity's NaN pattern
% AND velocity's (geo's) -- v10 only used sal00_mid's, then
% one-directionally re-masked geo by it (geo(isnan(dx1))=nan below).
% That's fine only if velocity is NaN exclusively where salinity
% already is. Diagnosed on the 2013-2022 dataset: 1.27M of 13.55M
% (depth,time,gap) points have valid salinity but NaN velocity (real,
% independent Total_TPUD coverage gaps, e.g. P5/P6's ~58% coverage) --
% those points still counted their full width in dx2/wsal's
% denominator while contributing nothing to V0's numerator, biasing V0
% relative to wsal and breaking the mean+Mov+gyre vs. total_direct
% identity (residual up to 219.9 Sv even after v10's dx2 time-
% resolution fix alone). Fixed by combining both NaN patterns before
% anything downstream (wsal, dx2, V0) is computed from dx1.
ind=find(isnan(sal00_mid) | isnan(geo));
dx1(ind)=nan;
geo(isnan(dx1))=nan;

% dx2: total valid width profile, summed across the 9 gaps -- now
% time-resolved ([depth x time], v10), built directly from dx1's real
% per-day mask rather than a day-1 stand-in. Necessary once coverage
% actually varies over time (v10 note above): S0/V00/mov/mean_term
% below are all weighted by dx2, and gyre/total_direct are weighted by
% dx1 directly -- both sides now use the exact same real per-day width.
dx2=squeeze(nansum(dx1,3));

% Informational only (not a bias risk since v10/v11 -- dx2 now tracks
% dx1's real, combined per-day mask exactly): how much the coverage
% mask actually varies across days, for context.
mask=isnan(dx1);
mask_t1=repmat(mask(:,1,:),[1 size(mask,2) 1]);
frac_diff=mean(mask(:)~=mask_t1(:));
fprintf('Site/depth coverage mask varies across days for %.4f%% of depth,time,gap points (dx2 is time-resolved, so this is not a source of bias).\n', frac_diff*100);

% Distance weighted salinity mean for each depth
wsal=nansum(sal00_mid.*dx1,3)./nansum(dx1,3);

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
% same way as wsal (nansum(dx1,3) instead of the constant sum(dx)).
% V00: dx2-weighted depth-mean of V0, for the same reason S0 was
% redefined above (v5 fix). vel_prime is the baroclinic (overturning)
% component.
V0=nansum(geo,3)./nansum(dx1,3);
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
% dx1 (per-gap, per-time, topography-aware width) rather than dx2.
vel_pp=vel00-V0;          % v''_i(z,t) = per-gap velocity - zonal-mean profile
sal_pp=sal00_mid-wsal;    % S''_i(z,t) = per-gap salinity - zonal-mean profile
aux_gyre=vel_pp.*sal_pp.*dx1*mean(diff(pre));
gyre=-nansum(nansum(aux_gyre,3),1)/1e6./S0;

% --------------------------------------------------------
% Consistency check: a directly-computed total (no decomposition, using
% the real per-gap v and S, weighted by dx1) vs. the sum of the three
% components above. Should be close to zero (v5 fix).
aux_total=vel00.*sal00_mid.*dx1*mean(diff(pre));
total_direct=-nansum(nansum(aux_total,3),1)/1e6./S0;
residual=total_direct-(mean_term+mov+gyre);
fprintf('Residual (total_direct - [mean+Mov+gyre]): mean=%.4g Sv, max|.|=%.4g Sv\n', ...
        nanmean(residual), nanmax(abs(residual)));

% --------------------------------------------------------
% Some plottings
mov_low=lowpass_filter(mov,45);
[nmov,coef]=sinfitb_tot(dt,mov,365.25);
slope=coef(1)+coef(2)*dt;

p=plot(dt,mov,dt,slope,'k',dt,mov_low,'r');
set(p(2),'linewidth',1.5)
set(p(1),'linewidth',1.5)
set(p(3),'linewidth',1.5)
title('Freshwater transport (Mov) at 34.5^oS from IES','fontsize',14)
ylabel('Mov (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)

print -dpng -r150  mov_IES

%save mov_samba_marion dt mov mean_term gyre total_direct

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
