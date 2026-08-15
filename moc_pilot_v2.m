% MOC "Pilot" method, v2: adds an aggregate AREA_TOPO correction,
% following the root-cause investigation (see CLAUDE.md, "Root cause of
% the full-array-vs-Pilot MOC intensity gap") which found AREA_TOPO --
% present per-gap in the full 8-gap method, absent from v1's Pilot -- is
% the dominant real cause of the intensity/shape discrepancy between the
% two methods (up to 10,900 m^2/s by ~4400dbar, versus the BPR-vs-LNM
% reference choice and the telescoping/baroclinic part, both confirmed
% NOT to matter). User's question: "nao tem um jeito de considerarmos
% essa AREA_TOPO para o Pilot?"
%
% No independent bathymetry exists for a single A-to-P1 span (Step E's
% topo_corr_msm60.mat only has per-adjacent-gap AREA_TOPO_1..8). Rather
% than acquire new bathymetry data, this builds an AGGREGATE, width-
% weighted average of the 8 REAL per-gap AREA_TOPO_i(z) profiles
% (themselves static/depth-only, no time dimension -- confirmed in Step
% E's own header notes) as a data-driven proxy for "what fraction of
% the total A-P1 width is topographically blocked at depth z":
%
%   AREA_TOPO_pilot(z) = sum_i[ dx_i * AREA_TOPO_i(z) ] / sum_i[ dx_i ]
%
% Deliberately built from the STATIC AREA_TOPO_i(z) profiles and STATIC
% dx_i weights only -- not from the full array's own day-varying dx2(z,t)
% (which already reflects which of the 8 gaps had valid data that day).
% Using dx2(z,t) directly would have silently imported the full array's
% own coverage outages (e.g. the 2017-2019 P5/P6/P4/P8 dropouts) INTO
% the Pilot, defeating the entire point of having an independent
% cross-check unaffected by those specific instrument gaps.
%
% dx_i (8-gap path) vs dx_pilot (direct A-P1 great-circle, v1's width):
% sum(dx_i)=7846km vs dx_pilot=6197km, a 26.6% difference -- NOT a bug,
% but a real geometric fact: two points on the same latitude (except the
% equator) are connected by a great circle that bows toward the pole,
% shorter than the constant-latitude arc the array physically follows.
% v1's dx_pilot (direct great-circle) never affected v1's Psi/MOC output
% (it's used to convert TPUD->velocity and then multiplied back by the
% identical value, canceling exactly -- same "any consistent depth-
% constant factor cancels via the depth-mean removal" mechanism already
% confirmed for the BPR-vs-LNM reference choice). But AREA_TOPO_pilot's
% weighted average DOES need the physically-correct along-section width,
% so this version uses dx_total=sum(dx_i) throughout instead.
%
% v2 (2026-08-15): adds the AREA_TOPO_pilot correction and switches to
% dx_total=sum(dx_i) for the width. Everything else (Gpan/RefTPUD/
% RelTPUD/Offset derivation, Ekman, streamfunction construction) is
% identical to v1 -- see that script's header for the full reasoning.

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
%%% West: site A pressure + T/S, from samba_w.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'MOV/samba_w.mat'],'pres_A','dt','Salinity_A','Temperature_A','presrange')
dt_SAM=dt; clear dt

n_pad=531-numel(presrange);
presrange=[presrange; (presrange(end)+10:10:presrange(end)+10*n_pad)'];
Salinity_A=[Salinity_A; NaN(n_pad,size(Salinity_A,2))];
Temperature_A=[Temperature_A; NaN(n_pad,size(Temperature_A,2))];
clear n_pad
Pressure_SAM=presrange; clear presrange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% East: site P1 pressure (My Fields) + T/S (Daily_Tau rebuild)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'pres_P1','dt')
dt_GH=dt; clear dt

load([prefix_marion,'ies_profiles/Wrk/IES_Make_Profiles_FrSA_SAMBA_6PIES.mat'], ...
    'Temperature_P1','Salinity_P1','presrange')
Pressure_GH=presrange; clear presrange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Positions and along-section width. dx_i: same 8 adjacent-site
%%% distances moc_streamfunction_v1.m uses; dx_total=sum(dx_i) is the
%%% physically-correct along-34.5S-parallel width (see header note on
%%% why the direct great-circle A-P1 distance is NOT the same thing).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position_West
position_East
lon8=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat8=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx_i=gsw_distance(lon8,lat8,0);
dx_total=sum(dx_i);
fprintf('sum(dx_i), 8-gap path (used as the Pilot width): %.1f km\n',dx_total/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aggregate AREA_TOPO_pilot(z): width-weighted average of the 8 real,
%%% static per-gap topographic-width profiles from Step E v2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat'], ...
    'AREA_TOPO_1','AREA_TOPO_2','AREA_TOPO_3','AREA_TOPO_4','AREA_TOPO_5','AREA_TOPO_6','AREA_TOPO_7','AREA_TOPO_8','Pressure');
AT_grid=Pressure; clear Pressure
AT=[AREA_TOPO_1(:) AREA_TOPO_2(:) AREA_TOPO_3(:) AREA_TOPO_4(:) AREA_TOPO_5(:) AREA_TOPO_6(:) AREA_TOPO_7(:) AREA_TOPO_8(:)];
AREA_TOPO_pilot_grid=(AT*dx_i(:))/dx_total;
clear AREA_TOPO_1 AREA_TOPO_2 AREA_TOPO_3 AREA_TOPO_4 AREA_TOPO_5 AREA_TOPO_6 AREA_TOPO_7 AREA_TOPO_8 AT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Geopotential anomaly (relative to sea surface)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gpan_A=sw_gpan(Salinity_A,Temperature_A,Pressure_SAM(:));
Gpan_P1=sw_gpan(Salinity_P1,Temperature_P1,Pressure_GH(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bottom pressure indices/densities (same as v1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ind_A=find(Pressure_SAM==1370);
p_ind_P1=find(Pressure_GH==1280);

rhob_A=nanmean(sw_dens(Salinity_A(p_ind_A,:),Temperature_A(p_ind_A,:),Pressure_SAM(p_ind_A)));
rhob_P1=nanmean(sw_dens(Salinity_P1(p_ind_P1,:),Temperature_P1(p_ind_P1,:),Pressure_GH(p_ind_P1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BPR reference correction (same as v1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gpan_AtoP1_corr=Gpan_A(p_ind_A,:)-Gpan_A(p_ind_P1,:);
Gpan_AtoP1_corr=Gpan_AtoP1_corr-nanmean(Gpan_AtoP1_corr);

pres_A=pres_A'*1e+04;
pres_P1=pres_P1*1e+04;
pres_A=pres_A-nanmean(pres_A);
pres_P1=pres_P1-nanmean(pres_P1);

pres_AtoP1=pres_A-(Gpan_AtoP1_corr.*rhob_A);
clear Gpan_AtoP1_corr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Align West/East/Ekman on their common dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dt_tmp,i_SAM,i_GH]=intersect(dt_SAM,dt_GH);
[dt,i_tmp,i_ekman]=intersect(dt_tmp,dt_ekman);
i_SAM=i_SAM(i_tmp); i_GH=i_GH(i_tmp);
clear dt_tmp i_tmp

fprintf('Aligned dt: %d common days (%s to %s)\n',numel(dt),datestr(dt(1)),datestr(dt(end)));

Gpan_A=Gpan_A(:,i_SAM); pres_A=pres_A(i_SAM); pres_AtoP1=pres_AtoP1(i_SAM);
Gpan_P1=Gpan_P1(:,i_GH); pres_P1=pres_P1(i_GH);
ekmanN_AtoZ=ekmanN_AtoZ(i_ekman);
clear dt_SAM dt_GH dt_ekman i_SAM i_GH i_ekman

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Throw out data below 5150dbar (same cutoff as Step E v2)
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

assert(isequal(Pressure,AT_grid(AT_grid<=5150)),'AREA_TOPO_pilot grid does not match Pressure after trim.');
AREA_TOPO_pilot=AREA_TOPO_pilot_grid(AT_grid<=5150);
clear AT_grid AREA_TOPO_pilot_grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relative + reference, same construction as v1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RelTPUD_pilot=(1./sw_f(-34.5)).*(Gpan_P1 - Gpan_A);
RelTPUD_pilot=-RelTPUD_pilot;

RefTPUD_pilot=((pres_AtoP1 - pres_P1))./(sw_f(-34.5).*rhob_P1);
RefTPUD_pilot=RefTPUD_pilot-nanmean(RefTPUD_pilot);

Offset_pilot=RefTPUD_pilot-RelTPUD_pilot(p_ind_P1,:);

Absolute_TPUD_pilot=NaN*ones(size(RelTPUD_pilot));
for i=1:length(Pressure)
    Absolute_TPUD_pilot(i,:)=RelTPUD_pilot(i,:)+Offset_pilot;
end
clear i p_ind_A p_ind_P1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply the aggregate AREA_TOPO_pilot(z) correction -- same order as
%%% Step E v2 (topo multiplication BEFORE Ekman is added)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Absolute_TPUD_pilot=Absolute_TPUD_pilot.*AREA_TOPO_pilot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add Ekman at the 7 surface rows (0-60dbar), unscaled (AREA_TOPO_i==1
%%% there for all 8 real gaps, so AREA_TOPO_pilot==1 there too)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pressure_Ekman=[0:10:60]';
Total_TPUD_pilot=Absolute_TPUD_pilot;
for i=1:length(Pressure_Ekman)
    pin=find(Pressure==Pressure_Ekman(i));
    Total_TPUD_pilot(pin,:)=Absolute_TPUD_pilot(pin,:)+ekmanN_AtoZ(:)'./60;
end
clear i pin Absolute_TPUD_pilot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Streamfunction, same construction as v1/moc_streamfunction_v1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vel_pilot=Total_TPUD_pilot/dx_total;
V00_pilot=nanmean(vel_pilot,1);
vel_prime_pilot=vel_pilot-V00_pilot;
transport_profile_pilot=vel_prime_pilot*dx_total;

pre=Pressure; clear Pressure
dz=mean(diff(pre));

Psi_pilot=cumsum(transport_profile_pilot,1,'omitnan')*dz/1e6;
fprintf('Sanity check: mean(Psi_pilot at max depth)=%.4g Sv (should be ~0 by construction)\n', nanmean(Psi_pilot(end,:)));

Psi_pilot_low=nan(size(Psi_pilot));
for zz=1:size(Psi_pilot,1)
    Psi_pilot_low(zz,:)=lowpass_filter(Psi_pilot(zz,:),45);
end
[MOC_pilot_low,imax_pilot_low]=max(Psi_pilot_low,[],1);
MOC_pilot_depth_low=pre(imax_pilot_low)';
fprintf('LOWPASS-FILTERED (45-day) Pilot v2 (with AREA_TOPO) MOC: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOC_pilot_low), nanstd(MOC_pilot_low));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 3-way comparison: full array, Pilot v1 (no topo), Pilot v2 (topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_full=load('moc_streamfunction_v1.mat','dt','MOC_low');
s_v1=load('moc_pilot_v1.mat','dt','MOC_pilot_low');

[dt_common,ia,ib]=intersect(dt,s_full.dt);
[~,~,ic]=intersect(dt_common,s_v1.dt);
fprintf('Overlay common dt: %d days\n',numel(dt_common));

figure('Position',[100 100 900 500])
plot(dt_common,s_full.MOC_low(ib),'linewidth',1.5); hold on
plot(dt_common,s_v1.MOC_pilot_low(ic),'color',[0.6 0.6 0.6],'linewidth',1)
plot(dt_common,MOC_pilot_low(ia),'r','linewidth',1.5)
legend('Full array (8 gaps)','Pilot v1, no AREA\_TOPO','Pilot v2, with AREA\_TOPO')
title('MOC index at 34.5^oS: full array vs. Pilot, with/without AREA\_TOPO (45-day lowpass)','fontsize',13)
ylabel('MOC (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_pilot_v2_vs_full_IES

fprintf('\nMeans over common period: Full=%.2f Sv, Pilot-v1(no topo)=%.2f Sv, Pilot-v2(topo)=%.2f Sv\n', ...
    nanmean(s_full.MOC_low(ib)), nanmean(s_v1.MOC_pilot_low(ic)), nanmean(MOC_pilot_low(ia)));
fprintf('corr(full,pilot-v2)=%.3f (was %.3f for v1)\n', ...
    corr2_nan(s_full.MOC_low(ib),MOC_pilot_low(ia)), corr2_nan(s_full.MOC_low(ib),s_v1.MOC_pilot_low(ic)));

% Time-mean profile comparison (Kersale Figure 2 style)
s_full2=load('moc_streamfunction_v1.mat','pre','Psi_mean','Psi_std');
assert(isequal(s_full2.pre,pre),'pressure grids differ for profile comparison');
Psi_pilot_mean=nanmean(Psi_pilot_low,2);
Psi_pilot_std=nanstd(Psi_pilot_low,0,2);

figure('Position',[500 100 650 650])
hold on
fill([s_full2.Psi_mean-s_full2.Psi_std;flipud(s_full2.Psi_mean+s_full2.Psi_std)],[pre;flipud(pre)],[0.8 0.85 1],'EdgeColor','none','FaceAlpha',0.6)
fill([Psi_pilot_mean-Psi_pilot_std;flipud(Psi_pilot_mean+Psi_pilot_std)],[pre;flipud(pre)],[1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.6)
p1=plot(s_full2.Psi_mean,pre,'b','linewidth',2);
p2=plot(Psi_pilot_mean,pre,'r','linewidth',2);
plot([0 0],[0 max(pre)],'k--')
set(gca,'YDir','reverse')
xlabel('\Psi (Sv)'); ylabel('Pressure (dbar)')
title({'Time-mean MOC streamfunction vs. depth at 34.5^oS','Pilot v2 with AREA\_TOPO (shading = +/-1 std)'},'fontsize',12)
legend([p1 p2],{'Full array (8 gaps)','Pilot v2 (A-P1, with AREA\_TOPO)'},'Location','southeast')
grid on
print -dpng -r150 moc_pilot_v2_profile_comparison_IES

fprintf('Full array:  peak %.1f Sv\nPilot v2:    peak %.1f Sv\n', max(s_full2.Psi_mean), max(Psi_pilot_mean));

save('moc_pilot_v2.mat','dt','pre','Psi_pilot','Psi_pilot_low','MOC_pilot_low','MOC_pilot_depth_low','dx_total','AREA_TOPO_pilot')

function r=corr2_nan(a,b)
ok=~isnan(a)&~isnan(b);
cc=corrcoef(a(ok),b(ok));
r=cc(1,2);
end
