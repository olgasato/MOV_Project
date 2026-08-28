% MOC "Pilot" method, v6: fixes two real bugs found by comparing this
% script line-by-line against Kersale/Marion's own ORIGINAL 2-CPIES
% script (Marions_code/Full_Depth_MOC/MOC_Estimate_2cpies_constituents.m,
% authored by marion.kersale, predates this whole project) --
% discovered while testing whether feeding OUR pipeline with HER
% original published raw data (~/research/sambar/kersale/data_A.txt,
% data_Z.txt -- the calibrated tau1000/BPR supplementary data package
% for Kersale et al. 2020, Science Advances) could reproduce her own
% saved output. It could, almost exactly (r=0.999 against her saved
% MOC_2cpies_FixeP.mat, r=0.712 against her saved full-array output --
% matching the paper's r=0.73) -- ONLY once both bugs below were fixed.
% Using OUR sign/missing-step (i.e. v5's own formula) on HER data
% reproduced our own r=0.485-ish shortfall almost exactly instead,
% proving the correlation gap chased across many prior sessions
% (lag analysis, BPR-vs-baroclinic component test, Gpan noise,
% cross-spectral coherence, P1 glitch check -- all in CLAUDE.md) was
% predominantly these two coding bugs, not measurement noise, phase,
% or GEM-vintage drift as previously suspected.
%
% BUG 1 -- SIGN: v5's RefTPUD_pilot = (pres_AtoP1 - pres_P1)/(f*rhob_P1).
% Her original: RefTPUD1 = (pres_P1 - pres_AtoZ)/(f*rhob_P1) -- the
% exact negative. Fixed below.
%
% BUG 2 -- MISSING STEP: her original script adds a time-mean ABSOLUTE
% reference velocity from ECCO at 1500dbar (DATA/ECCO/VVEL_Mean_34p5N.mat,
% averaged over the ECCO longitude band between sites A and P1) on top
% of the demeaned-BPR reference -- exactly analogous to what Step E
% already does for the FULL array (Step E's own "Add the time-mean
% velocity from OFES" step, already documented in
% moc_streamfunction_v2.m's header). v5.m (and v1-v4 before it) never
% had an equivalent step -- its BPR reference was purely relative/
% demeaned with no absolute anchor at all. Added below, inserted
% BEFORE the AREA_TOPO_pilot multiplication (matching Step E's own
% ordering: absolute reference -> preTopo save point -> AREA_TOPO ->
% Ekman -> shelf).
%
% NOTE on distance units: gsw_distance returns METERS (not km, unlike
% sw_dist used in her original script with an explicit *1000 km->m
% conversion) -- confirmed directly from gsw_distance.m's own header
% ("Note. The output is in m not km."). Used directly below with no
% extra unit conversion.
%
% Everything else (data sources, AREA_TOPO_pilot aggregate correction,
% Ekman, shelf correction, Kersale streamfunction definition) is
% UNCHANGED from v5.m -- this version isolates exactly the two fixes
% above, on our own current data sources (not Kersale's raw data,
% which was only used for the isolated diagnostic in
% moc_pilot_kersale_original_data_test.m).
%
% v6 (2026-08-28): first cut with both fixes applied.

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
%%% Positions and along-section width (v2's dx_total=sum(dx_i)),
%%% plus the direct A-to-P1 distance needed for the ECCO-mean step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position_West
position_East
lon8=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat8=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx_i=gsw_distance(lon8,lat8,0);
dx_total=sum(dx_i);
dx_pilot=gsw_distance([cpiesW(1).lon cpiesE(1).lon],[cpiesW(1).lat cpiesE(1).lat],0); % meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aggregate AREA_TOPO_pilot(z) -- identical construction to v2-v5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_AT=load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat'], ...
    'AREA_TOPO_1','AREA_TOPO_2','AREA_TOPO_3','AREA_TOPO_4','AREA_TOPO_5','AREA_TOPO_6','AREA_TOPO_7','AREA_TOPO_8','Pressure');
AT_grid=S_AT.Pressure;
AT=[S_AT.AREA_TOPO_1(:) S_AT.AREA_TOPO_2(:) S_AT.AREA_TOPO_3(:) S_AT.AREA_TOPO_4(:) S_AT.AREA_TOPO_5(:) S_AT.AREA_TOPO_6(:) S_AT.AREA_TOPO_7(:) S_AT.AREA_TOPO_8(:)];
AREA_TOPO_pilot_grid=(AT*dx_i(:))/dx_total;
clear S_AT AT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Geopotential anomaly (relative to sea surface)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gpan_A=sw_gpan(Salinity_A,Temperature_A,Pressure_SAM(:));
Gpan_P1=sw_gpan(Salinity_P1,Temperature_P1,Pressure_GH(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bottom pressure indices/densities (A=1370dbar, P1=1280dbar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ind_A=find(Pressure_SAM==1370);
p_ind_P1=find(Pressure_GH==1280);

rhob_A=nanmean(sw_dens(Salinity_A(p_ind_A,:),Temperature_A(p_ind_A,:),Pressure_SAM(p_ind_A)));
rhob_P1=nanmean(sw_dens(Salinity_P1(p_ind_P1,:),Temperature_P1(p_ind_P1,:),Pressure_GH(p_ind_P1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BPR reference correction (same as v1-v5)
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
%%% Relative (unchanged) + reference (BUG 1 FIX: sign flipped to
%%% match Kersale's original RefTPUD1=(pres_P1-pres_AtoZ)/(f*rhob_P1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RelTPUD_pilot=(1./sw_f(-34.5)).*(Gpan_P1 - Gpan_A);
RelTPUD_pilot=-RelTPUD_pilot;

RefTPUD_pilot=((pres_P1 - pres_AtoP1))./(sw_f(-34.5).*rhob_P1); % BUG 1 FIX (was pres_AtoP1-pres_P1)
RefTPUD_pilot=RefTPUD_pilot-nanmean(RefTPUD_pilot);

Offset_pilot=RefTPUD_pilot-RelTPUD_pilot(p_ind_P1,:);

Absolute_TPUD_pilot=NaN*ones(size(RelTPUD_pilot));
for i=1:length(Pressure)
    Absolute_TPUD_pilot(i,:)=RelTPUD_pilot(i,:)+Offset_pilot;
end
clear i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUG 2 FIX: add the time-mean ABSOLUTE reference velocity from
%%% ECCO at 1500dbar -- entirely missing from v1-v5. Mirrors Step E's
%%% own "Add the time-mean velocity from OFES/ECCO" step for the full
%%% array, and Kersale's original 2-CPIES script's identical step
%%% (MOC_Estimate_2cpies_constituents.m, "Add the time-mean velocity
%%% from OFES" section). Inserted BEFORE AREA_TOPO_pilot, matching
%%% Step E's own ordering (absolute reference -> preTopo -> AREA_TOPO).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'DATA/ECCO/VVEL_Mean_34p5N.mat'],'lon','V','dpt');
Pressure_ECCO=sw_prs(dpt',-34.5);
vmean2=interp2(Pressure_ECCO,lon,V',Pressure,lon); warning('off','all');
vmean2=vmean2';

pref=1500;
p_ind_1500=find(Pressure==pref);
lon_ind1=find(lon>=cpiesW(1).lon & lon<=cpiesE(1).lon);
vel_ecco1=nanmean(vmean2(p_ind_1500,lon_ind1));
RefTPUD_mean=vel_ecco1.*dx_pilot; % dx_pilot already in meters (gsw_distance)
fprintf('ECCO time-mean reference: vel=%.6f m/s, dist=%.1f m, RefTPUD_mean=%.2f m^2/s\n',vel_ecco1,dx_pilot,RefTPUD_mean);
clear vmean2 Pressure_ECCO lon V dpt lon_ind1

Offset_mean=RefTPUD_mean-nanmean(Absolute_TPUD_pilot(p_ind_1500,:));
Absolute_TPUD_pilot=Absolute_TPUD_pilot+Offset_mean;
clear p_ind_A p_ind_P1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply the aggregate AREA_TOPO_pilot(z) correction, then Ekman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Absolute_TPUD_pilot=Absolute_TPUD_pilot.*AREA_TOPO_pilot;

Pressure_Ekman=[0:10:60]';
Total_TPUD_pilot=Absolute_TPUD_pilot;
for i=1:length(Pressure_Ekman)
    pin=find(Pressure==Pressure_Ekman(i));
    Total_TPUD_pilot(pin,:)=Absolute_TPUD_pilot(pin,:)+ekmanN_AtoZ(:)'./60;
end
clear i pin Absolute_TPUD_pilot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Shelf transport correction, same static profile the full array uses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ecco/Wrk/ECCO_TPUD_ShelfBits.mat'],'West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(Pressure,pres_shelf_grid);
assert(isequal(ia(:),(1:numel(Pressure))'),'Pressure is not a subset of the shelf grid as expected');
shelf_profile=West_TPUD(ib)+East_TPUD_new(ib);
fprintf('Shelf correction: West=%.2f Sv, East=%.2f Sv (static, added to every day)\n', ...
    nansum(West_TPUD)*10/1e6, nansum(East_TPUD_new)*10/1e6);
clear West_TPUD East_TPUD_new pres_shelf_grid ia ib

Total_TPUD_pilot=Total_TPUD_pilot+shelf_profile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kersale definition: RAW cumulative transport (now +shelf), no
%%% depth-mean removal, evaluated at h_star found from the Pilot's OWN
%%% time-mean profile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pre=Pressure; clear Pressure
dz=mean(diff(pre));

Psi_raw=cumsum(Total_TPUD_pilot,1,'omitnan')*dz/1e6;
Psi_raw_mean=nanmean(Psi_raw,2);
[~,imax_mean]=max(Psi_raw_mean);
h_star_pilot=pre(imax_mean);
fprintf('h_star (Pilot, time-mean transition depth): %d dbar\n',h_star_pilot);

MOCup_pilot_raw=Psi_raw(imax_mean,:);
MOCup_pilot_low=lowpass_filter(MOCup_pilot_raw,45);
fprintf('LOWPASS (45-day) Pilot v6 MOCup at h_star: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOCup_pilot_low), nanstd(MOCup_pilot_low));

Psi_low=nan(size(Psi_raw));
for zz=1:size(Psi_raw,1)
    Psi_low(zz,:)=lowpass_filter(Psi_raw(zz,:),45);
end
Psi_std=nanstd(Psi_low,0,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare against the calibration-window-fixed full-array method (v4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_full=load('moc_streamfunction_v4.mat');

[dt_common,ia,ib]=intersect(dt,s_full.dt);
fprintf('Overlay common dt: %d days\n',numel(dt_common));

w1=datenum(2013,9,11); w2=datenum(2017,7,17);
win=dt_common>=w1 & dt_common<=w2;
full_raw_win=s_full.MOCup_raw(ib(win));
pilot_raw_win=MOCup_pilot_raw(ia(win));
cc_raw=corrcoef(full_raw_win,pilot_raw_win,'rows','complete');
fprintf('\ncorr(full,pilot), RAW daily, 2013-2017 (paper: r=0.73): %.3f\n',cc_raw(1,2));
fprintf('Means, 2013-2017: Full=%.2f Sv, Pilot=%.2f Sv\n',nanmean(full_raw_win),nanmean(pilot_raw_win));
fprintf('Stds,  2013-2017: Full=%.2f Sv, Pilot=%.2f Sv\n',nanstd(full_raw_win),nanstd(pilot_raw_win));

figure('Position',[100 100 900 500])
plot(dt_common,s_full.MOCup_low(ib),'linewidth',1.5); hold on
plot(dt_common,MOCup_pilot_low(ia),'r','linewidth',1.5)
legend('Full array (8 gaps)','Pilot v6 (A-P1, BPR+ECCO-mean+AREA\_TOPO)')
title('MOC_{up} index at 34.5^oS, calib-fixed + shelf: full array vs. Pilot v6 (45-day lowpass)','fontsize',13)
ylabel('MOC_{up} (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_pilot_v6_vs_full_IES

fprintf('\nMeans over common period (whole record, lowpass): Full=%.2f Sv (h*=%d), Pilot=%.2f Sv (h*=%d)\n', ...
    nanmean(s_full.MOCup_low(ib)), s_full.h_star, nanmean(MOCup_pilot_low(ia)), h_star_pilot);
fprintf('corr(full,pilot), lowpass, whole record=%.3f\n', corr2_nan(s_full.MOCup_low(ib),MOCup_pilot_low(ia)));

figure('Position',[500 100 650 650])
hold on
fill([s_full.Psi_raw_mean-s_full.Psi_std;flipud(s_full.Psi_raw_mean+s_full.Psi_std)],[s_full.pre;flipud(s_full.pre)],[0.8 0.85 1],'EdgeColor','none','FaceAlpha',0.6)
fill([Psi_raw_mean-Psi_std;flipud(Psi_raw_mean+Psi_std)],[pre;flipud(pre)],[1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.6)
p1=plot(s_full.Psi_raw_mean,s_full.pre,'b','linewidth',2);
p2=plot(Psi_raw_mean,pre,'r','linewidth',2);
plot([0 0],[0 max(pre)],'k--')
set(gca,'YDir','reverse')
xlabel('\Psi (Sv)'); ylabel('Pressure (dbar)')
title({'Time-mean MOC cumulative transport vs. depth at 34.5^oS','Calibration-window-fixed + shelf, Pilot v6 (shading = +/-1 std)'},'fontsize',12)
legend([p1 p2],{'Full array (8 gaps)','Pilot v6 (A-P1)'},'Location','southeast')
grid on
print -dpng -r150 moc_pilot_v6_profile_comparison_IES

save('moc_pilot_v6.mat','dt','pre','Psi_raw','Psi_low','h_star_pilot','MOCup_pilot_raw','MOCup_pilot_low','Psi_raw_mean','Psi_std')

function r=corr2_nan(a,b)
ok=~isnan(a)&~isnan(b);
cc=corrcoef(a(ok),b(ok));
r=cc(1,2);
end
