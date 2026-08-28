% Runs OUR Pilot pipeline fed with KERSALE'S OWN ORIGINAL raw data
% (calibrated tau1000 + bottom pressure, sites A and Z/P1, from the
% published supplementary data package at
% ~/research/sambar/kersale/data_A.txt and data_Z.txt, 2013-09-11 to
% 2017-07-16) instead of our own reprocessed samba_w.mat/Daily_Tau
% sources -- to separate "our own calculation differs from hers" from
% "our data/calibration provenance differs from hers".
%
% Motivation: found that Kersale/Marion's own ORIGINAL saved outputs
% (MOC_Estimate_2cpies_constituents.m's MOC_2cpies_FixeP.mat, and the
% original full-array Full_Depth_OvertuningEstimates.mat) correlate at
% r=0.712 with each other -- essentially the paper's r=0.73 -- while
% OUR reimplementation (moc_pilot_v5.m vs moc_streamfunction_v4.m)
% only reaches r=0.485. Comparing MOC_Estimate_2cpies_constituents.m
% line-by-line against moc_pilot_v5.m found TWO real, structural
% differences (not just noise):
%   1. SIGN: her RefTPUD1=(pres_P1-pres_AtoZ)/(f*rhob_P1); ours
%      (moc_pilot_v5.m) has RefTPUD_pilot=(pres_AtoP1-pres_P1)/(f*
%      rhob_P1) -- these are exact negatives of each other, both
%      before AND after demeaning.
%   2. MISSING STEP: her script adds a time-mean ABSOLUTE reference
%      velocity from ECCO at 1500dbar (VVEL_Mean_34p5N.mat, averaged
%      over the ECCO longitude band between sites A and P1) on top of
%      the demeaned-BPR reference -- exactly analogous to what Step E
%      already does for the FULL array (already documented in
%      CLAUDE.md's "Corrected construction" section: "Step E's Add the
%      time-mean velocity from OFES step"). moc_pilot_v5.m never had
%      an equivalent step at all -- its BPR reference is purely
%      relative/demeaned, with no absolute anchor.
%
% Tests four variants, all built from her raw tau1000/BPR run through
% the ORIGINAL-era GEM tables (gem/west|east/Wrk/do_temper_field.mat +
% do_sal_field.mat, dated March 2021 -- NOT samba_w.mat's later
% retrained vintage), so the GEM-vintage confound documented elsewhere
% in CLAUDE.md doesn't contaminate this specific test:
%   A: faithful replica of her script (her sign, +ECCO-mean, no AREA_TOPO)
%      -- sanity check: should closely reproduce MOC_2cpies_FixeP.mat.
%   B: her sign, +ECCO-mean, +our AREA_TOPO_pilot correction
%   C: OUR sign, no ECCO-mean, no AREA_TOPO -- isolates the two bugs together
%   D: OUR sign, no ECCO-mean, +our AREA_TOPO_pilot -- literally our
%      current moc_pilot_v5.m formula, just fed her data instead of ours
%
% v1 (2026-08-27): first cut.

clear;clc;
prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';
prefix_kersale='/Users/olga/research/sambar/kersale/';
addpath([prefix_marion,'functions/seawater/']);
addpath([prefix_marion,'functions/positions_pies']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load Kersale's own raw calibrated data (site A, site Z=P1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dt_A,tau_A,bpr_A]=read_kersale_site([prefix_kersale,'data_A.txt']);
[dt_P1,tau_P1,bpr_P1]=read_kersale_site([prefix_kersale,'data_Z.txt']);
assert(isequal(dt_A,dt_P1),'Site A and Z date vectors are not identical -- unexpected for the published package.');
dt=dt_A; clear dt_P1
fprintf('Kersale original data: %d days, %s to %s\n',numel(dt),datestr(dt(1)),datestr(dt(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GEM lookup, ORIGINAL-era tables (March 2021 -- predates
%%% samba_w.mat's later retrained vintage), same interp2 call her own
%%% IES_Make_Profiles_*.m scripts use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gW=load([prefix_marion,'gem/west/Wrk/do_temper_field.mat'],'taurange','presrange','SMF_temperature');
gWs=load([prefix_marion,'gem/west/Wrk/do_sal_field.mat'],'SMF_salinity');
gE=load([prefix_marion,'gem/east/Wrk/do_temper_field.mat'],'taurange','presrange','SMF_temperature');
gEs=load([prefix_marion,'gem/east/Wrk/do_sal_field.mat'],'SMF_salinity');

Pressure_SAM=gW.presrange(:);
Temperature_A=nan(numel(Pressure_SAM),numel(dt));
Salinity_A=nan(numel(Pressure_SAM),numel(dt));
for i=1:numel(dt)
    if ~isnan(tau_A(i))
        Temperature_A(:,i)=interp2(gW.taurange,Pressure_SAM,gW.SMF_temperature,tau_A(i),Pressure_SAM);
        Salinity_A(:,i)=interp2(gW.taurange,Pressure_SAM,gWs.SMF_salinity,tau_A(i),Pressure_SAM);
    end
end

Pressure_GH=gE.presrange(:);
Temperature_P1=nan(numel(Pressure_GH),numel(dt));
Salinity_P1=nan(numel(Pressure_GH),numel(dt));
for i=1:numel(dt)
    if ~isnan(tau_P1(i))
        Temperature_P1(:,i)=interp2(gE.taurange,Pressure_GH,gE.SMF_temperature,tau_P1(i),Pressure_GH);
        Salinity_P1(:,i)=interp2(gE.taurange,Pressure_GH,gEs.SMF_salinity,tau_P1(i),Pressure_GH);
    end
end
clear gW gWs gE gEs i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Common pressure grid (both West/East GEM tables use the same
%%% 10dbar spacing from 0 -- trim to the shorter one, West to 5000dbar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmin=min(numel(Pressure_SAM),numel(Pressure_GH));
assert(isequal(Pressure_SAM(1:nmin),Pressure_GH(1:nmin)),'West/East GEM pressure grids do not share a common origin/spacing.');
Pressure=Pressure_SAM(1:nmin);
Temperature_A=Temperature_A(1:nmin,:); Salinity_A=Salinity_A(1:nmin,:);
Temperature_P1=Temperature_P1(1:nmin,:); Salinity_P1=Salinity_P1(1:nmin,:);
clear Pressure_SAM Pressure_GH nmin

% Trim to <=5150dbar to match the AREA_TOPO grid's cap (same cutoff
% moc_pilot_v5.m uses for the same reason).
bad=find(Pressure>5150);
Pressure(bad)=[];
Temperature_A(bad,:)=[]; Salinity_A(bad,:)=[];
Temperature_P1(bad,:)=[]; Salinity_P1(bad,:)=[];
clear bad

Gpan_A=sw_gpan(Salinity_A,Temperature_A,Pressure(:));
Gpan_P1=sw_gpan(Salinity_P1,Temperature_P1,Pressure(:));

p_ind_A=find(Pressure==1370);
p_ind_P1=find(Pressure==1280);
rhob_A=nanmean(sw_dens(Salinity_A(p_ind_A,:),Temperature_A(p_ind_A,:),Pressure(p_ind_A)));
rhob_P1=nanmean(sw_dens(Salinity_P1(p_ind_P1,:),Temperature_P1(p_ind_P1,:),Pressure(p_ind_P1)));

Gpan_AtoP1_corr=Gpan_A(p_ind_A,:)-Gpan_A(p_ind_P1,:);
Gpan_AtoP1_corr=Gpan_AtoP1_corr-nanmean(Gpan_AtoP1_corr);

pres_A=bpr_A(:)'*1e+04; pres_A=pres_A-nanmean(pres_A);
pres_P1=bpr_P1(:)'*1e+04; pres_P1=pres_P1-nanmean(pres_P1);
pres_AtoP1=pres_A-(Gpan_AtoP1_corr.*rhob_A);

RelTPUD=(1./sw_f(-34.5)).*(Gpan_P1 - Gpan_A);
RelTPUD=-RelTPUD;

RefTPUD_ours=((pres_AtoP1 - pres_P1))./(sw_f(-34.5).*rhob_P1);
RefTPUD_ours=RefTPUD_ours-nanmean(RefTPUD_ours);
RefTPUD_hers=((pres_P1 - pres_AtoP1))./(sw_f(-34.5).*rhob_P1);
RefTPUD_hers=RefTPUD_hers-nanmean(RefTPUD_hers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AREA_TOPO_pilot (our v2+ addition -- not in her original script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position_West
position_East
lon8=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat8=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx_i=gsw_distance(lon8,lat8,0);
dx_total=sum(dx_i);
dx_pilot=gsw_distance([cpiesW(1).lon cpiesE(1).lon],[cpiesW(1).lat cpiesE(1).lat],0);

S_AT=load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat'], ...
    'AREA_TOPO_1','AREA_TOPO_2','AREA_TOPO_3','AREA_TOPO_4','AREA_TOPO_5','AREA_TOPO_6','AREA_TOPO_7','AREA_TOPO_8','Pressure');
AT_grid=S_AT.Pressure; % NOT our own Pressure -- loaded into a struct to avoid shadowing it
AT=[S_AT.AREA_TOPO_1(:) S_AT.AREA_TOPO_2(:) S_AT.AREA_TOPO_3(:) S_AT.AREA_TOPO_4(:) S_AT.AREA_TOPO_5(:) S_AT.AREA_TOPO_6(:) S_AT.AREA_TOPO_7(:) S_AT.AREA_TOPO_8(:)];
AREA_TOPO_pilot_grid=(AT*dx_i(:))/dx_total;
clear S_AT AT

[~,ia,ib]=intersect(Pressure,AT_grid);
assert(isequal(ia(:),(1:numel(Pressure))'),'Pressure is not a subset of the AREA_TOPO grid as expected');
AREA_TOPO_pilot=AREA_TOPO_pilot_grid(ib);
clear AT_grid AREA_TOPO_pilot_grid ia ib

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ekman, shelf -- same static/basin-wide inputs both versions use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ccmp/Wrk/Ekman_transports_9pies.mat'],'ekmanN_AtoZ','dt');
dt_ekman=dt; clear dt
load([prefix_marion,'ecco/Wrk/ECCO_TPUD_ShelfBits.mat'],'West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(Pressure,pres_shelf_grid);
assert(isequal(ia(:),(1:numel(Pressure))'),'Pressure is not a subset of the shelf grid as expected');
shelf_profile=West_TPUD(ib)+East_TPUD_new(ib);
clear West_TPUD East_TPUD_new pres_shelf_grid ia ib

Pressure_Ekman=[0:10:60]';
ekmanN_full=nan(1,numel(dt_A));
[~,ia2,ib2]=intersect(dt_A,dt_ekman);
ekmanN_full(ia2)=ekmanN_AtoZ(ib2);
ekmanN_AtoZ=ekmanN_full;
assert(~any(isnan(ekmanN_AtoZ)),'Ekman record does not fully cover the Kersale data window.');
clear ekmanN_full ia2 ib2 dt_ekman

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ECCO time-mean absolute reference at 1500dbar (her original
%%% script's step -- entirely missing from moc_pilot_v5.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'DATA/ECCO/VVEL_Mean_34p5N.mat'],'lon','V','dpt');
Pressure_ECCO=sw_prs(dpt',-34.5);
vmean2=interp2(Pressure_ECCO,lon,V',Pressure,lon); warning('off','all');
vmean2=vmean2';
pref=1500;
p_ind_1500=find(Pressure==pref);
lon_ind1=find(lon>=cpiesW(1).lon & lon<=cpiesE(1).lon);
vel_ecco1=nanmean(vmean2(p_ind_1500,lon_ind1));
dumdist1=dx_pilot; % gsw_distance already returns meters, unlike sw_dist(...,'km')
RefTPUD_mean=vel_ecco1.*dumdist1;
clear vmean2 Pressure_ECCO lon V dpt lon_ind1 dumdist1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build all 4 variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w1=datenum(2013,9,11); w2=datenum(2017,7,17);
dz=mean(diff(Pressure));

variants=struct('name',{},'sign',{},'ecco',{},'topo',{});
variants(1)=struct('name','A: faithful replica (her sign, +ECCO, no topo)','sign','hers','ecco',true,'topo',false);
variants(2)=struct('name','B: her sign, +ECCO, +our AREA_TOPO','sign','hers','ecco',true,'topo',true);
variants(3)=struct('name','C: OUR sign, no ECCO, no topo','sign','ours','ecco',false,'topo',false);
variants(4)=struct('name','D: OUR sign, no ECCO, +our AREA_TOPO (=our v5.m formula on her data)','sign','ours','ecco',false,'topo',true);

fprintf('\n%-58s %8s %8s %8s\n','Variant','h_star','mean','std');
results=struct();
for v=1:numel(variants)
    if strcmp(variants(v).sign,'hers')
        RefTPUD=RefTPUD_hers;
    else
        RefTPUD=RefTPUD_ours;
    end
    Offset=RefTPUD-RelTPUD(p_ind_P1,:);
    Absolute_TPUD=NaN*ones(size(RelTPUD));
    for i=1:numel(Pressure)
        Absolute_TPUD(i,:)=RelTPUD(i,:)+Offset;
    end

    if variants(v).ecco
        Offset_mean=RefTPUD_mean-nanmean(Absolute_TPUD(p_ind_1500,:));
        Absolute_TPUD=Absolute_TPUD+Offset_mean;
    end

    if variants(v).topo
        Absolute_TPUD=Absolute_TPUD.*AREA_TOPO_pilot;
    end

    Total_TPUD=Absolute_TPUD;
    for i=1:numel(Pressure_Ekman)
        pin=find(Pressure==Pressure_Ekman(i));
        Total_TPUD(pin,:)=Absolute_TPUD(pin,:)+ekmanN_AtoZ./60;
    end
    Total_TPUD=Total_TPUD+shelf_profile;

    Psi_raw=cumsum(Total_TPUD,1,'omitnan')*dz/1e6;
    Psi_raw_mean=nanmean(Psi_raw,2);
    [~,imax]=max(Psi_raw_mean);
    h_star=Pressure(imax);
    MOCup=Psi_raw(imax,:);

    win=dt_A>=w1 & dt_A<=w2;
    fprintf('%-58s %6d dbar %6.2f Sv %6.2f Sv\n',variants(v).name,h_star,nanmean(MOCup(win)),nanstd(MOCup(win)));

    results(v).name=variants(v).name;
    results(v).dt=dt_A; results(v).MOCup=MOCup; results(v).h_star=h_star;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare against Kersale/Marion's OWN saved original outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_orig=load([prefix_marion,'Full_Depth_MOC/MOC_2cpies_FixeP.mat']);
f_orig=load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OvertuningEstimates.mat']);

fprintf('\nReference: her original pilot (MOC_2cpies_FixeP.mat), 2013-2017: mean=%.2f Sv, std=%.2f Sv\n', ...
    nanmean(p_orig.MOC_Total(p_orig.dt>=w1&p_orig.dt<=w2)),nanstd(p_orig.MOC_Total(p_orig.dt>=w1&p_orig.dt<=w2)));
fprintf('Reference: her original full array (Full_Depth_OvertuningEstimates.mat): mean=%.2f Sv, std=%.2f Sv\n', ...
    nanmean(f_orig.Upper_Transport)/1e6,nanstd(f_orig.Upper_Transport)/1e6);

fprintf('\n%-58s %14s %14s\n','Variant','corr vs HER pilot','corr vs HER full array');
for v=1:numel(variants)
    [dtc,ia,ib]=intersect(results(v).dt,p_orig.dt);
    win=dtc>=w1&dtc<=w2;
    cc1=corrcoef(results(v).MOCup(ia(win)),p_orig.MOC_Total(ib(win)),'rows','complete');

    [dtc2,ia2,ib2]=intersect(results(v).dt,f_orig.dt);
    win2=dtc2>=w1&dtc2<=w2;
    cc2=corrcoef(results(v).MOCup(ia2(win2)),f_orig.Upper_Transport(ib2(win2))/1e6,'rows','complete');

    fprintf('%-58s %14.3f %14.3f\n',variants(v).name,cc1(1,2),cc2(1,2));
end

save('moc_pilot_kersale_original_data_test.mat','results','variants','dt_A')

function [dt,tau,bpr]=read_kersale_site(fname)
fid=fopen(fname,'r');
raw=textscan(fid,'%f%f%f%f%f%f%f','CommentStyle','%');
fclose(fid);
dt=datenum(raw{1},raw{2},raw{3},raw{4},raw{5},0);
tau=raw{6};
bpr=raw{7};
end
