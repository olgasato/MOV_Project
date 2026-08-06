% Step 2: MOV calculation using Marion's Total Velocities
% between each pair of sites. Her velocities are obtained at
% Marions_code/Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT.mat
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

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')

load concat_IESsamba

% Keep only A, C, D (drop AA, B, BB, CC) -- relies on '*B' matching
% both 'B' and 'BB' (both end in the letter B), and '*CC' matching only
% CC (leaving C alone, since C doesn't end in CC). Fragile: a future
% site name ending in 'B' or 'CC' would be swept up too.
clear *AA *B *CC

% --------------------------------------------------------
% These calculations are for the OSM26 presentation and
% use the velocity fields estimated by Marion's paper in 2021.

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT']);

% --------------------------------------------------------
% --------------------------------------------------------
% Building the salinity matrices
% Cut the pressure at 5000 m (using Samba_W limit)

ind=pree<=5000;
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
pre_0=repmat(pre,[1 1405 9]);
difpre=diff(pre_0,1);
[l,m,n]=size(sal00_mid);

% Make a 3D dx to multiply by the salinity'
dx1=repmat(reshape(dx,1,1,[]),l,m,1);
% Make the points where the salinity is NaN (limit of the topography).
ind=find(isnan(sal00_mid));
dx1(ind)=nan;
dx2=nansum(squeeze(dx1(:,1,:)),2);

% Diagnostic (kept as a standing safeguard, not a one-time check): dx2
% assumes the site-validity mask in dx1 is time-invariant, i.e. fully
% explained by the fixed p_cpies* depth cutoffs rather than real data
% gaps. Confirmed true for the 2013-2017/9-site dataset (v4). As more
% stations/years are added this could silently change.
mask=isnan(dx1);
mask_t1=repmat(mask(:,1,:),[1 size(mask,2) 1]);
frac_diff=mean(mask(:)~=mask_t1(:));
if frac_diff>0
    warning('mov_samba_marion_v6:dx1TimeVarying', ...
        'dx1 validity mask VARIES in time (%.4f%% of depth,time,gap points differ from the t=1 pattern used by dx2) -- dx2/mov/mean_term may be biased relative to gyre/total_direct, which use the full time-resolved dx1.', frac_diff*100);
else
    fprintf('dx1 validity mask is time-invariant -- dx2 is an exact stand-in for the full dx1.\n');
end

% Distance weighted salinity mean for each depth
wsal=nansum(sal00_mid.*dx1,3)./nansum(dx1,3);

% Mean salinity time series: dx2-weighted depth-mean of wsal, so that
% sal_prime=wsal-S0 integrates to zero when weighted by dx2 across
% depth -- required for mean+Mov+gyre to reconstruct total_direct
% exactly (v5 fix).
S0=nansum(wsal.*dx2)/nansum(dx2);

% --------------------------------------------------------
% Estimate Sal'=sal00-S0, subtract S0 and zonally mean to get the
% the <S(z)>.
sal_prime=nanmean(wsal-S0,3);
% --------------------------------------------------------
% Do the same for velocity. Find the barotropic component first
vel00=nan*ones([l m n]);
geo=nan*ones([l m n]);
for ii=1:8
str=['vel00(:,:,',num2str(ii),')=Total_TPUD',num2str(ii),'/dx(',num2str(ii),');'];
eval(str)
str=['geo(:,:,',num2str(ii),')=Total_TPUD',num2str(ii),';'];
eval(str)
end

clear Total_*
% --------------------------------------------------------
% Mask geo the same way wsal/dx1 are masked (each site's cutoff depth),
% so V0 and wsal are computed over the same set of valid gaps (v3 fix).
geo(isnan(dx1))=nan;

% --------------------------------------------------------
% V0: mean zonal velocity as function of depth, masked/weighted the
% same way as wsal (nansum(dx1,3) instead of the constant sum(dx)).
% V00: dx2-weighted depth-mean of V0, for the same reason S0 was
% redefined above (v5 fix). vel_prime is the baroclinic (overturning)
% component.
V0=nansum(geo,3)./nansum(dx1,3);
V00=nansum(V0.*dx2)/nansum(dx2);
vel_prime=V0-V00;

% Integrating <v(z)><S(z)> vertically
aux=vel_prime.*sal_prime*mean(diff(pre));

% as the ocean width varies with depth:
mov=-nansum(aux.*dx2)/1e6./S0;

% --------------------------------------------------------
% Mean (net/barotropic) component.
% F_mean = -(1/S0)*V00*S0*A = -V00*A -- S0 cancels algebraically, so it
% does not appear below. Uses the same dx2 (depth-only) width as mov,
% for consistency with that term.
mean_term=-V00.*sum(dx2)*mean(diff(pre))/1e6;

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
% hardcoded (v6) -- produces the same 47 months/same boundaries on the
% current 2013-09-2017-07 dataset, but adapts if the time window
% changes (e.g. once more station data extends the record).

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


