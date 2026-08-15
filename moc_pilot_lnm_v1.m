% MOC "Pilot" method, LEVEL-OF-NO-MOTION (LNM) variant, for comparison
% against moc_pilot_v1.m's BPR-absolute-referenced version. This is
% option (1) from the original 3-way clarifying question the user was
% asked before building moc_pilot_v1.m ("baroclinic-only shear with a
% level-of-no-motion assumption") -- revisited afterward: "quais eram
% as outras opcoes... sera que vale a pena investigar?", user then
% asked to implement it for comparison.
%
% Purpose: quantify how much the BPR absolute reference actually
% matters for the two-endpoint Pilot method, by computing the SAME
% A-to-P1 baroclinic shear (RelTPUD_pilot) but referencing it to ZERO
% velocity at the deepest common level both profiles resolve, instead
% of to the two sites' real bottom-pressure records. No pres_A/pres_P1/
% rhob_A/rhob_P1 needed at all -- this variant only needs Gpan_A/
% Gpan_P1, unlike moc_pilot_v1.m.
%
% Reference level: deepest pressure row where NEITHER Gpan_A nor
% Gpan_P1 is entirely NaN across the whole record (A's real profile is
% capped at 5000dbar by samba_w.mat -- see Step E v2's header note --
% and padded with NaN beyond that, so this lands at 5000dbar in
% practice). Found programmatically rather than hardcoded, so it stays
% correct if either source's depth coverage ever changes.
%
% v1 (2026-08-15): first cut. Shares the RelTPUD_pilot derivation
% (Gpan/dt-alignment/Ekman/distance) with moc_pilot_v1.m verbatim --
% see that script's header for the reasoning behind the East-minus-West
% sign convention, the A/P1 shallow-site identification, etc.

clear;clc;

prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';
addpath([prefix_marion,'functions/seawater/']);
addpath([prefix_marion,'functions/positions_pies']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ekman (whole-basin A-to-P1 total)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ccmp/Wrk/Ekman_transports_9pies.mat'],'ekmanN_AtoZ','dt');
dt_ekman=dt; clear dt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% West: site A T/S, from samba_w.mat (no pres_A needed for LNM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'MOV/samba_w.mat'],'dt','Salinity_A','Temperature_A','presrange')
dt_SAM=dt; clear dt

n_pad=531-numel(presrange);
presrange=[presrange; (presrange(end)+10:10:presrange(end)+10*n_pad)'];
Salinity_A=[Salinity_A; NaN(n_pad,size(Salinity_A,2))];
Temperature_A=[Temperature_A; NaN(n_pad,size(Temperature_A,2))];
clear n_pad
Pressure_SAM=presrange; clear presrange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% East: site P1 T/S (no pres_P1 needed for LNM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'dt')
dt_GH=dt; clear dt

load([prefix_marion,'ies_profiles/Wrk/IES_Make_Profiles_FrSA_SAMBA_6PIES.mat'], ...
    'Temperature_P1','Salinity_P1','presrange')
Pressure_GH=presrange; clear presrange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Positions (for A-P1 distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position_West
position_East
lon=[cpiesW(1).lon cpiesE(1).lon];
lat=[cpiesW(1).lat cpiesE(1).lat];
dx_pilot=gsw_distance(lon,lat,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Geopotential anomaly (relative to sea surface)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gpan_A=sw_gpan(Salinity_A,Temperature_A,Pressure_SAM(:));
Gpan_P1=sw_gpan(Salinity_P1,Temperature_P1,Pressure_GH(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Align West/East/Ekman on their common dt (same pattern as
%%% moc_pilot_v1.m/Step E v2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dt_tmp,i_SAM,i_GH]=intersect(dt_SAM,dt_GH);
[dt,i_tmp,i_ekman]=intersect(dt_tmp,dt_ekman);
i_SAM=i_SAM(i_tmp); i_GH=i_GH(i_tmp);
clear dt_tmp i_tmp

fprintf('Aligned dt: %d common days (%s to %s) out of %d (West), %d (East), %d (Ekman)\n', ...
    numel(dt),datestr(dt(1)),datestr(dt(end)),numel(dt_SAM),numel(dt_GH),numel(dt_ekman));

Gpan_A=Gpan_A(:,i_SAM);
Gpan_P1=Gpan_P1(:,i_GH);
ekmanN_AtoZ=ekmanN_AtoZ(i_ekman);
clear dt_SAM dt_GH dt_ekman i_SAM i_GH i_ekman

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Throw out data below 5150dbar (same cutoff as Step E v2/moc_pilot_v1.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bad=find(Pressure_SAM>5150);
Gpan_A(bad,:)=[]; Pressure_SAM(bad)=[];
clear bad
bad=find(Pressure_GH>5150);
Gpan_P1(bad,:)=[]; Pressure_GH(bad)=[];
clear bad
assert(isequal(Pressure_SAM,Pressure_GH),'Pressure_SAM/Pressure_GH differ after the >5150dbar trim.');
Pressure=Pressure_SAM;
clear Pressure_SAM Pressure_GH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relative (baroclinic, surface-referenced) transport per unit depth
%%% -- identical to moc_pilot_v1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RelTPUD_pilot=(1./sw_f(-34.5)).*(Gpan_P1 - Gpan_A);
RelTPUD_pilot=-RelTPUD_pilot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Level-of-no-motion reference: deepest pressure row where NEITHER
%%% Gpan_A nor Gpan_P1 is entirely NaN across the whole record.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref_idx=find(~all(isnan(Gpan_A),2) & ~all(isnan(Gpan_P1),2), 1, 'last');
fprintf('LNM reference level: %d dbar (index %d of %d)\n',Pressure(ref_idx),ref_idx,numel(Pressure));

Offset_LNM=-RelTPUD_pilot(ref_idx,:);
Absolute_TPUD_pilot_LNM=NaN*ones(size(RelTPUD_pilot));
for i=1:length(Pressure)
    Absolute_TPUD_pilot_LNM(i,:)=RelTPUD_pilot(i,:)+Offset_LNM;
end
clear i ref_idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add Ekman at the 7 surface rows (0-60dbar) -- same construction as
%%% moc_pilot_v1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pressure_Ekman=[0:10:60]';
Total_TPUD_pilot_LNM=Absolute_TPUD_pilot_LNM;
for i=1:length(Pressure_Ekman)
    pin=find(Pressure==Pressure_Ekman(i));
    Total_TPUD_pilot_LNM(pin,:)=Absolute_TPUD_pilot_LNM(pin,:)+ekmanN_AtoZ(:)'./60;
end
clear i pin Absolute_TPUD_pilot_LNM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Streamfunction: identical construction to moc_pilot_v1.m (remove
%%% the depth-mean/barotropic component, cumsum, 45-day lowpass,
%%% MOC(t)=max_z Psi(z,t))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vel_pilot=Total_TPUD_pilot_LNM/dx_pilot;
V00_pilot=nanmean(vel_pilot,1);
vel_prime_pilot=vel_pilot-V00_pilot;
transport_profile_pilot=vel_prime_pilot*dx_pilot;

pre=Pressure; clear Pressure
dz=mean(diff(pre));

Psi_pilot_LNM=cumsum(transport_profile_pilot,1,'omitnan')*dz/1e6;
fprintf('Sanity check: mean(Psi_pilot_LNM at max depth)=%.4g Sv (should be ~0 by construction)\n', nanmean(Psi_pilot_LNM(end,:)));

Psi_pilot_LNM_low=nan(size(Psi_pilot_LNM));
for zz=1:size(Psi_pilot_LNM,1)
    Psi_pilot_LNM_low(zz,:)=lowpass_filter(Psi_pilot_LNM(zz,:),45);
end
[MOC_pilot_LNM_low,imax]=max(Psi_pilot_LNM_low,[],1);
MOC_pilot_LNM_depth_low=pre(imax)';
fprintf('LOWPASS-FILTERED (45-day) Pilot-LNM MOC: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOC_pilot_LNM_low), nanstd(MOC_pilot_LNM_low));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare against the BPR-referenced Pilot and the full array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_full=load('moc_streamfunction_v1.mat','dt','MOC_low');
s_bpr=load('moc_pilot_v1.mat','dt','MOC_pilot_low');

[dt_common,ia,ib]=intersect(dt,s_full.dt);
[~,~,ic]=intersect(dt_common,s_bpr.dt);
fprintf('Overlay common dt: %d days (%s to %s)\n',numel(dt_common),datestr(dt_common(1)),datestr(dt_common(end)));

figure('Position',[100 100 900 500])
plot(dt_common,s_full.MOC_low(ib),'linewidth',1.5); hold on
plot(dt_common,s_bpr.MOC_pilot_low(ic),'r','linewidth',1.5)
plot(dt_common,MOC_pilot_LNM_low(ia),'color',[0 0.6 0],'linewidth',1.5)
legend('Full array (8 gaps)','Pilot, BPR reference','Pilot, LNM reference')
title('MOC index at 34.5^oS: full array vs. Pilot (BPR vs. LNM) (45-day lowpass)','fontsize',14)
ylabel('MOC (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_pilot_lnm_vs_bpr_IES

fprintf('\nMeans over common period: Full=%.2f Sv, Pilot-BPR=%.2f Sv, Pilot-LNM=%.2f Sv\n', ...
    nanmean(s_full.MOC_low(ib)), nanmean(s_bpr.MOC_pilot_low(ic)), nanmean(MOC_pilot_LNM_low(ia)));
fprintf('corr(BPR,LNM)=%.3f, corr(full,LNM)=%.3f\n', ...
    corr2_nan(s_bpr.MOC_pilot_low(ic),MOC_pilot_LNM_low(ia)), corr2_nan(s_full.MOC_low(ib),MOC_pilot_LNM_low(ia)));

save('moc_pilot_lnm_v1.mat','dt','pre','Psi_pilot_LNM','Psi_pilot_LNM_low','MOC_pilot_LNM_low','MOC_pilot_LNM_depth_low','dx_pilot')

function r=corr2_nan(a,b)
ok=~isnan(a)&~isnan(b);
cc=corrcoef(a(ok),b(ok));
r=cc(1,2);
end
