% Step 2 (heat), phase 3: mean+overturning+gyre decomposition of MHT,
% analogous to mov_samba_marion*.m's freshwater decomposition, but
% derived from scratch for heat rather than ported directly -- heat
% has no natural reference-quantity analog to freshwater's S0 (the
% "reference temperature problem" in the MHT literature). See the
% derivation note below for exactly where that shows up.
%
% Citations: see mov_samba_marion*.m's header -- same underlying
% method/array, Marion's 2021 paper (https://doi.org/10.1029/2020JC016947).
%
% v2 (2026-08-14): builds on mht_samba_marion_v1.m (phase 1-2: faithful
% reproduction of Marion's MHT_Estimate_constituents.m, then extended
% to 2013-2022) but restructures the computation into a mean/
% overturning/gyre split instead of the total+Ekman/relative/reference-
% sensitivity style v1 used. Also switches the velocity side from
% Total_TPUD (already AREA_TOPO-multiplied) to Absolute_TPUD*_preTopo +
% explicit dxA=dx(gap)*AREA_TOPO(z) weighting, same as
% mov_samba_marion_v12's fix -- the same phantom-width bias that
% motivated that fix for freshwater applies here too, since it's the
% same array/velocity field.
%
% --- Derivation (why T0 is needed explicitly here, unlike S0) ---
% Freshwater: F = -(1/S0)*integral(v*S dA), and S0 cancels
% algebraically out of the "mean" term (F_mean=-(1/S0)*V00*S0*Area=
% -V00*Area) because freshwater transport is DEFINED as a salt-flux-to-
% freshwater-flux conversion using a reference salinity -- the 1/S0
% factor is not optional, it's what makes the quantity "freshwater
% transport" rather than "salt transport".
%
% Heat has no equivalent conversion: Q=rho*Cp*integral(v*T dA) is
% already in physical units (Watts) with no natural reference-
% quantity to divide by. Decomposing T=T0+T'(z)+T''(x,z) the same way
% S is decomposed, and re-deriving the cross-term cancellation from
% scratch (T0 defined as the dx2-weighted section-mean temperature,
% analogous to S0):
%   Q_overturning = -rho*Cp*dz*sum_z[(V0-V00)*(wtemp-T0)*dx2]
%   Q_mean        = -rho*Cp*dz*V00*T0*sum_z(dx2)
% (sum_z is sum-over-depth; wtemp is the dxA-weighted zonal-mean
% temperature profile, analogous to wsal). Confirmed algebraically that
% Q_mean+Q_overturning+Q_gyre==Q_total_direct exactly regardless of the
% T0 choice (same cross-term-cancellation argument as mov's v5/v12
% fixes), BUT: Q_overturning and Q_gyre are individually INVARIANT to
% the T0 reference choice (verified: shifting T0 by a constant adds
% V0-weighted-by-dx2 integrated over depth to Q_overturning, which is
% zero by construction since V00 IS that weighted mean) while Q_mean is
% NOT invariant -- it scales directly with whatever T0 is chosen. This
% is exactly "the reference temperature problem": the overturning
% (baroclinic) heat transport is a robust, reference-independent
% physical quantity, but the net/mean (barotropic) heat transport is
% only meaningful relative to an explicit, stated reference temperature.
% T0 here = dx2-weighted section-mean temperature (same convention as
% S0), NOT T0=0C -- stated explicitly so Q_mean's absolute value can be
% interpreted/recomputed under a different convention if needed.
%
% Static shelf (Q_shelfW/E, from ECCO) and deep (Q_deep, from WOA,
% assuming a fixed 6.6 Sv deep transport) contributions are NOT
% included in the mean/overturning/gyre split below -- they're external
% additions to the CPIES-array-measured baroclinic structure, not part
% of what this decomposition can meaningfully attribute. Reported
% separately (Q_total_full) for comparison against v1's Heat_total.
%
% Same deliberate simplification as v1 for gap 3 (D-P8): direct
% Temperature_D/Temperature_P8 average instead of the (not yet
% extended) altimetry-based "interior" estimate -- see v1's header and
% CLAUDE.md for the quantified impact.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Shelf and deep heat transports (static, climatological -- kept
%%% separate from the mean/overturning/gyre split, see header note)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';

load([prefix_marion,'ecco/Wrk/ECCO_Heat_ShelfBits.mat'],'West_Heat_Total','East_Heat_Total');
Q_shelfW=West_Heat_Total;clear West_Heat_Total
Q_shelfE=East_Heat_Total;clear East_Heat_Total

load([prefix_marion,'woa/Wrk/WOA_Heat_Deep.mat'],'Deep_Heat_Total');
Q_deep=Deep_Heat_Total;clear Deep_Heat_Total

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% West temperature (samba_w.mat + seasonal correction, 2009-2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ies_profiles/Wrk/IES_Make_Profiles_USA_SAM_SeasCorr_2013_2022.mat'], ...
    'Temperature_A','Temperature_C','Temperature_D','dt','presrange')
dt_SAM=dt; clear dt
Pressure_SAM=presrange; clear presrange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% East temperature (Daily_Tau rebuild + seasonal correction, 2013-2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ies_profiles/Wrk/IES_Make_Profiles_FrSA_SAMBA_6PIES_SeasCorr_2013_2022.mat'], ...
    'Temperature_P1','Temperature_P2','Temperature_P4','Temperature_P5', ...
    'Temperature_P6','Temperature_P8','dt','presrange')
dt_GH=dt; clear dt
Pressure_GH=presrange; clear presrange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Velocity: Absolute_TPUD*_preTopo + AREA_TOPO + Ekman_TPUD from
%%% Step E v2 (AREA_TOPO-consistent, same as mov_samba_marion_v12+)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat'])
dt_TPUD=dt; clear dt
clear Total_TPUD*_RefConst Total_TPUD*_RelConst Total_TPUD*_EKMANconst  % only need the plain Total_TPUD1..8 for the sanity check below

% AREA_TOPO_i are saved as [1 x 516] row vectors -- force to columns so
% the >4700dbar row-deletion below (bad,:)=[] works consistently with
% the other (already-column) variables it's applied alongside.
AREA_TOPO_1=AREA_TOPO_1(:); AREA_TOPO_2=AREA_TOPO_2(:); AREA_TOPO_3=AREA_TOPO_3(:); AREA_TOPO_4=AREA_TOPO_4(:);
AREA_TOPO_5=AREA_TOPO_5(:); AREA_TOPO_6=AREA_TOPO_6(:); AREA_TOPO_7=AREA_TOPO_7(:); AREA_TOPO_8=AREA_TOPO_8(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Align on the common dt (3-way intersect, same pattern as v1 and
%%%  everywhere else in this project)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dt_tmp,i_SAM,i_GH]=intersect(dt_SAM,dt_GH);
[dt,i_tmp,i_TPUD]=intersect(dt_tmp,dt_TPUD);
i_SAM=i_SAM(i_tmp); i_GH=i_GH(i_tmp);
clear dt_tmp i_tmp

fprintf('Aligned dt: %d common days (%s to %s) out of %d (West), %d (East), %d (Total_TPUD)\n', ...
    numel(dt),datestr(dt(1)),datestr(dt(end)),numel(dt_SAM),numel(dt_GH),numel(dt_TPUD));

Temperature_A=Temperature_A(:,i_SAM); Temperature_C=Temperature_C(:,i_SAM); Temperature_D=Temperature_D(:,i_SAM);
Temperature_P1=Temperature_P1(:,i_GH); Temperature_P2=Temperature_P2(:,i_GH); Temperature_P4=Temperature_P4(:,i_GH);
Temperature_P5=Temperature_P5(:,i_GH); Temperature_P6=Temperature_P6(:,i_GH); Temperature_P8=Temperature_P8(:,i_GH);

for vv={'Total_TPUD1','Total_TPUD2','Total_TPUD3','Total_TPUD4','Total_TPUD5','Total_TPUD6','Total_TPUD7','Total_TPUD8', ...
        'Absolute_TPUD1_preTopo','Absolute_TPUD2_preTopo','Absolute_TPUD3_preTopo','Absolute_TPUD4_preTopo', ...
        'Absolute_TPUD5_preTopo','Absolute_TPUD6_preTopo','Absolute_TPUD7_preTopo','Absolute_TPUD8_preTopo', ...
        'Ekman_TPUD1','Ekman_TPUD2','Ekman_TPUD3','Ekman_TPUD4','Ekman_TPUD5','Ekman_TPUD6','Ekman_TPUD7','Ekman_TPUD8'}
    eval([vv{1},'=',vv{1},'(:,i_TPUD);']);
end
clear vv i_SAM i_GH i_TPUD dt_SAM dt_GH dt_TPUD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Throw out all data below 4700 dbar at all sites (same cutoff as
%%%  v1/the original MHT script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bad=find(Pressure_SAM>4700);
Temperature_A(bad,:)=[]; Temperature_C(bad,:)=[]; Temperature_D(bad,:)=[];
Pressure_SAM(bad)=[];
clear bad

bad=find(Pressure_GH>4700);
Temperature_P1(bad,:)=[]; Temperature_P2(bad,:)=[]; Temperature_P4(bad,:)=[];
Temperature_P5(bad,:)=[]; Temperature_P6(bad,:)=[]; Temperature_P8(bad,:)=[];
Pressure_GH(bad)=[];
clear bad

bad=find(Pressure>4700);
for vv={'Total_TPUD1','Total_TPUD2','Total_TPUD3','Total_TPUD4','Total_TPUD5','Total_TPUD6','Total_TPUD7','Total_TPUD8', ...
        'Absolute_TPUD1_preTopo','Absolute_TPUD2_preTopo','Absolute_TPUD3_preTopo','Absolute_TPUD4_preTopo', ...
        'Absolute_TPUD5_preTopo','Absolute_TPUD6_preTopo','Absolute_TPUD7_preTopo','Absolute_TPUD8_preTopo', ...
        'AREA_TOPO_1','AREA_TOPO_2','AREA_TOPO_3','AREA_TOPO_4','AREA_TOPO_5','AREA_TOPO_6','AREA_TOPO_7','AREA_TOPO_8'}
    eval([vv{1},'(bad,:)=[];']);
end
Pressure(bad)=[];
clear bad

assert(isequal(Pressure_SAM,Pressure_GH) && isequal(Pressure_SAM,Pressure), ...
    'Pressure_SAM/Pressure_GH/Pressure differ after the >4700dbar trim.');
Pressure=Pressure_SAM;
clear Pressure_SAM Pressure_GH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Build the 8 gap-midpoint temperatures (Temp_3 uses the D/P8
%%%  average simplification, see header note) and the per-day site
%%%  coverage diagnostic (before the Temperature_* variables are
%%%  cleared).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sal00_like=cat(3,Temperature_A,Temperature_C,Temperature_D,Temperature_P8, ...
    Temperature_P6,Temperature_P5,Temperature_P4,Temperature_P2,Temperature_P1);
n_sites_valid=reshape(sum(squeeze(any(~isnan(sal00_like),1)),2),1,numel(dt));
clear sal00_like

Temp(:,:,1)=(Temperature_A+Temperature_C)./2;
Temp(:,:,2)=(Temperature_C+Temperature_D)./2;
Temp(:,:,3)=(Temperature_D+Temperature_P8)./2;
Temp(:,:,4)=(Temperature_P8+Temperature_P6)./2;
Temp(:,:,5)=(Temperature_P6+Temperature_P5)./2;
Temp(:,:,6)=(Temperature_P5+Temperature_P4)./2;
Temp(:,:,7)=(Temperature_P4+Temperature_P2)./2;
Temp(:,:,8)=(Temperature_P2+Temperature_P1)./2;
clear Temperature_A Temperature_C Temperature_D Temperature_P1 Temperature_P2 Temperature_P4 Temperature_P5 Temperature_P6 Temperature_P8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distances between sites (dx, per gap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build vel00/geo from Absolute_TPUD*_preTopo, add Ekman back at the
%%% surface rows, build dxA -- identical logic to mov_samba_marion_v12+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,m]=size(Total_TPUD1);
n=8;
vel00=nan*ones([l m n]);
geo=nan*ones([l m n]);
for ii=1:8
str=['vel00(:,:,',num2str(ii),')=Absolute_TPUD',num2str(ii),'_preTopo/dx(',num2str(ii),');'];
eval(str)
str=['geo(:,:,',num2str(ii),')=Absolute_TPUD',num2str(ii),'_preTopo;'];
eval(str)
end
clear Absolute_TPUD*_preTopo

for ii=1:8
str=['geo(1:7,:,',num2str(ii),')=geo(1:7,:,',num2str(ii),')+Ekman_TPUD',num2str(ii),';'];
eval(str)
str=['vel00(1:7,:,',num2str(ii),')=geo(1:7,:,',num2str(ii),')/dx(',num2str(ii),');'];
eval(str)
end
clear Ekman_TPUD1 Ekman_TPUD2 Ekman_TPUD3 Ekman_TPUD4 Ekman_TPUD5 Ekman_TPUD6 Ekman_TPUD7 Ekman_TPUD8 ii

AREA_TOPO_3d=nan*ones([l 1 n]);
for ii=1:8
str=['AREA_TOPO_3d(:,1,',num2str(ii),')=AREA_TOPO_',num2str(ii),'(:);'];
eval(str)
end
AREA_TOPO_3d=repmat(AREA_TOPO_3d,[1 m 1]);
clear AREA_TOPO_1 AREA_TOPO_2 AREA_TOPO_3 AREA_TOPO_4 AREA_TOPO_5 AREA_TOPO_6 AREA_TOPO_7 AREA_TOPO_8 ii

dx1=repmat(reshape(dx,1,1,[]),l,m,1);
ind=find(isnan(Temp) | isnan(geo) | isnan(AREA_TOPO_3d));
dx1(ind)=nan;
geo(isnan(dx1))=nan;
dxA=dx1.*AREA_TOPO_3d;
dx2=squeeze(nansum(dxA,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% wtemp/T0/V0/V00 -- same structure as mov_samba_marion's
%%% wsal/S0/V0/V00, T instead of S.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wtemp=nansum(Temp.*dxA,3)./nansum(dxA,3);
T0=nansum(wtemp.*dx2)./nansum(dx2);
temp_prime=wtemp-T0;

V0=nansum(vel00.*dxA,3)./nansum(dxA,3);
V00=nansum(V0.*dx2)./nansum(dx2);
vel_prime=V0-V00;

dz=mean(diff(Pressure));
rho=1025; %kg/m3
cp=4000; %J/(kg C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Q_mean/Q_overturning/Q_gyre -- see header derivation. Q_mean is the
%%% ONLY term that depends on the T0 reference choice.
%%%
%%% NOTE: no leading minus sign here, unlike mov_samba_marion's
%%% mean_term/mov/gyre. That "-" is specific to freshwater's
%%% salt-flux-to-freshwater-flux conversion (F=-(1/S0)*salt_flux) --
%%% not a general sign convention. Heat transport is
%%% Q=rho*Cp*integral(v*T dA) directly, matching v1/the original
%%% MHT_Estimate_constituents.m's Heat_total (also no leading "-").
%%% First draft of this script copied the "-" over by reflex and got
%%% Q_total_direct=-0.57 PW (wrong sign AND wrong magnitude vs. v1's
%%% +0.33 PW for the same period) -- caught by comparing against v1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_mean=rho.*cp.*V00.*T0.*nansum(dx2)*dz;

aux=vel_prime.*temp_prime*dz;
Q_overturning=rho.*cp.*nansum(aux.*dx2,1);

temp_pp=Temp-wtemp;
vel_pp=vel00-V0;
aux_gyre=rho.*cp.*vel_pp.*temp_pp.*dxA*dz;
Q_gyre=nansum(nansum(aux_gyre,3),1);

aux_total=rho.*cp.*vel00.*Temp.*dxA*dz;
Q_total_direct=nansum(nansum(aux_total,3),1);

residual=Q_total_direct-(Q_mean+Q_overturning+Q_gyre);
fprintf('Residual (Q_total_direct - [Q_mean+Q_overturning+Q_gyre]): mean=%.4g W, max|.|=%.4g W\n', ...
    nanmean(residual), nanmax(abs(residual)));

Q_total_full=Q_total_direct+Q_deep+Q_shelfE+Q_shelfW;

fprintf('mean(Q_mean)=%.4f PW, mean(Q_overturning)=%.4f PW, mean(Q_gyre)=%.4f PW, mean(Q_total_direct)=%.4f PW, mean(Q_total_full)=%.4f PW\n', ...
    nanmean(Q_mean)*1e-15, nanmean(Q_overturning)*1e-15, nanmean(Q_gyre)*1e-15, nanmean(Q_total_direct)*1e-15, nanmean(Q_total_full)*1e-15);
fprintf('std(Q_mean)=%.4f PW, std(Q_overturning)=%.4f PW, std(Q_gyre)=%.4f PW\n', ...
    nanstd(Q_mean)*1e-15, nanstd(Q_overturning)*1e-15, nanstd(Q_gyre)*1e-15);

% v12-style sanity check: vel00.*dxA should reconstruct Total_TPUD
recon_err=nan(1,8);
for ii=1:8
    recon=squeeze(vel00(:,:,ii)).*squeeze(dxA(:,:,ii));
    actual=eval(['Total_TPUD',num2str(ii)]);
    recon_err(ii)=nanmax(abs(recon(:)-actual(:)));
end
fprintf('Sanity check (vel00.*dxA vs Total_TPUD, should be ~0 at ALL rows): max abs diff per gap = %s\n', mat2str(recon_err,4));
clear Total_TPUD1 Total_TPUD2 Total_TPUD3 Total_TPUD4 Total_TPUD5 Total_TPUD6 Total_TPUD7 Total_TPUD8 recon recon_err ii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coverage shading (same n_sites_valid/graded-band approach as
%%% mov_samba_marion_v15/mht_samba_marion_v1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band_colors={[1 1 1],[0.88 0.88 0.88],[0.72 0.72 0.72],[0.5 0.5 0.5]};
band_labels={'9/9 sites','7-8/9 sites','6/9 sites','<6/9 sites'};
category=nan(size(n_sites_valid));
category(n_sites_valid==9)=1;
category(n_sites_valid>=7 & n_sites_valid<=8)=2;
category(n_sites_valid==6)=3;
category(n_sites_valid<6)=4;

figure
hold on
ylim0=[min(Q_overturning*1e-15) max(Q_overturning*1e-15)];
ylim0=ylim0+[-1 1]*0.05*diff(ylim0);
legend_h=[]; legend_str={};
for cc=2:4
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
p=plot(dt,Q_overturning*1e-15);
set(p,'linewidth',1)
title('MHT overturning component at 34.5^oS from IES (2013-2022)','fontsize',14)
ylabel('Q_{overturning} (PW)')
datetick('x',12);grid;axis('tight')
ylim(ylim0)
set(gca,'linewidth',1,'fontsize',14)
box on
legend([p(:);legend_h(:)],[{'Q overturning'},legend_str], ...
    'Location','southoutside','Orientation','horizontal','fontsize',8)
print -dpng -r150 mht_overturning_IES

save('mht_samba_marion_v2.mat','dt','Q_mean','Q_overturning','Q_gyre','Q_total_direct','Q_total_full', ...
    'T0','n_sites_valid','category')
