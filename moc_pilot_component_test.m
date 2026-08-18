% Focuses specifically on the Pilot's OWN construction to find the
% source of the phase/correlation shortfall vs. the full array
% (already confirmed: magnitude is essentially correct -- full array
% daily sigma matches the paper almost exactly, 15.39 vs 15.4 Sv; raw
% data has ~zero lag; the remaining gap is in point-to-point
% amplitude/shape agreement, r=0.485 vs the paper's 0.73).
%
% Tests whether the BPR (bottom pressure, absolute reference) part of
% the Pilot's own construction is a source of extra, uncorrelated
% noise not present in the baroclinic (relative/Gpan-based) part alone
% -- by building a "baroclinic-only" Pilot variant (RelTPUD_pilot with
% NO reference correction at all, i.e. referenced to the sea surface)
% alongside the normal BPR-referenced Pilot (moc_pilot_v5.m), both
% under the identical Kersale-definition/shelf-correction construction,
% and comparing each one's correlation with the full array
% (moc_streamfunction_v4.m).
%
% Note: this is NOT the same test as the old moc_pilot_lnm_v1.m
% (BPR-vs-level-of-no-motion), which was built under the ABANDONED
% depth-mean-removal MOC definition, where BOTH BPR and LNM references
% were proven mathematically identical (any per-day constant offset
% cancels under that definition -- see CLAUDE.md's MOC Pilot LNM
% section). Under the CURRENT Kersale definition (no depth-mean
% removal), a constant reference offset does NOT cancel -- it
% accumulates linearly with depth through the cumsum -- so the
% reference choice can genuinely matter here, unlike before.
%
% v1 (2026-08-18): first cut.

clear;clc;

prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';
addpath([prefix_marion,'functions/seawater/']);
addpath([prefix_marion,'functions/positions_pies']);
addpath('/Users/olga/matlab/stat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Same loading/setup as moc_pilot_v5.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'ccmp/Wrk/Ekman_transports_9pies.mat'],'ekmanN_AtoZ','dt');
dt_ekman=dt; clear dt

load([prefix_marion,'MOV/samba_w.mat'],'pres_A','dt','Salinity_A','Temperature_A','presrange')
dt_SAM=dt; clear dt
n_pad=531-numel(presrange);
presrange=[presrange; (presrange(end)+10:10:presrange(end)+10*n_pad)'];
Salinity_A=[Salinity_A; NaN(n_pad,size(Salinity_A,2))];
Temperature_A=[Temperature_A; NaN(n_pad,size(Temperature_A,2))];
Pressure_SAM=presrange; clear presrange n_pad

load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'pres_P1','dt')
dt_GH=dt; clear dt
load([prefix_marion,'ies_profiles/Wrk/IES_Make_Profiles_FrSA_SAMBA_6PIES.mat'],'Temperature_P1','Salinity_P1','presrange')
Pressure_GH=presrange; clear presrange

position_West
position_East
lon8=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat8=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx_i=gsw_distance(lon8,lat8,0);
dx_total=sum(dx_i);

load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat'], ...
    'AREA_TOPO_1','AREA_TOPO_2','AREA_TOPO_3','AREA_TOPO_4','AREA_TOPO_5','AREA_TOPO_6','AREA_TOPO_7','AREA_TOPO_8','Pressure');
AT_grid=Pressure; clear Pressure
AT=[AREA_TOPO_1(:) AREA_TOPO_2(:) AREA_TOPO_3(:) AREA_TOPO_4(:) AREA_TOPO_5(:) AREA_TOPO_6(:) AREA_TOPO_7(:) AREA_TOPO_8(:)];
AREA_TOPO_pilot_grid=(AT*dx_i(:))/dx_total;
clear AREA_TOPO_1 AREA_TOPO_2 AREA_TOPO_3 AREA_TOPO_4 AREA_TOPO_5 AREA_TOPO_6 AREA_TOPO_7 AREA_TOPO_8 AT

Gpan_A=sw_gpan(Salinity_A,Temperature_A,Pressure_SAM(:));
Gpan_P1=sw_gpan(Salinity_P1,Temperature_P1,Pressure_GH(:));

p_ind_A=find(Pressure_SAM==1370);
p_ind_P1=find(Pressure_GH==1280);
rhob_A=nanmean(sw_dens(Salinity_A(p_ind_A,:),Temperature_A(p_ind_A,:),Pressure_SAM(p_ind_A)));
rhob_P1=nanmean(sw_dens(Salinity_P1(p_ind_P1,:),Temperature_P1(p_ind_P1,:),Pressure_GH(p_ind_P1)));

Gpan_AtoP1_corr=Gpan_A(p_ind_A,:)-Gpan_A(p_ind_P1,:);
Gpan_AtoP1_corr=Gpan_AtoP1_corr-nanmean(Gpan_AtoP1_corr);
pres_A=pres_A'*1e+04; pres_P1=pres_P1*1e+04;
pres_A=pres_A-nanmean(pres_A); pres_P1=pres_P1-nanmean(pres_P1);
pres_AtoP1=pres_A-(Gpan_AtoP1_corr.*rhob_A);

[dt_tmp,i_SAM,i_GH]=intersect(dt_SAM,dt_GH);
[dt,i_tmp,i_ekman]=intersect(dt_tmp,dt_ekman);
i_SAM=i_SAM(i_tmp); i_GH=i_GH(i_tmp);

Gpan_A=Gpan_A(:,i_SAM); pres_A=pres_A(i_SAM); pres_AtoP1=pres_AtoP1(i_SAM);
Gpan_P1=Gpan_P1(:,i_GH); pres_P1=pres_P1(i_GH);
ekmanN_AtoZ=ekmanN_AtoZ(i_ekman);

bad=find(Pressure_SAM>5150);
Gpan_A(bad,:)=[]; Pressure_SAM(bad)=[];
bad=find(Pressure_GH>5150);
Gpan_P1(bad,:)=[]; Pressure_GH(bad)=[];
Pressure=Pressure_SAM;
AREA_TOPO_pilot=AREA_TOPO_pilot_grid(AT_grid<=5150);

RelTPUD_pilot=(1./sw_f(-34.5)).*(Gpan_P1 - Gpan_A);
RelTPUD_pilot=-RelTPUD_pilot;

RefTPUD_pilot=((pres_AtoP1 - pres_P1))./(sw_f(-34.5).*rhob_P1);
RefTPUD_pilot=RefTPUD_pilot-nanmean(RefTPUD_pilot);
Offset_pilot=RefTPUD_pilot-RelTPUD_pilot(p_ind_P1,:);

pre=Pressure; clear Pressure
dz=mean(diff(pre));

load([prefix_marion,'ecco/Wrk/ECCO_TPUD_ShelfBits.mat'],'West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(pre,pres_shelf_grid);
shelf_profile=zeros(size(pre));
shelf_profile(ia)=West_TPUD(ib)+East_TPUD_new(ib);

Pressure_Ekman=(0:10:60)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variant 1: FULL (BPR reference), matches moc_pilot_v5.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Absolute_TPUD_bpr=(RelTPUD_pilot+Offset_pilot).*AREA_TOPO_pilot;
Total_TPUD_bpr=Absolute_TPUD_bpr;
for i=1:numel(Pressure_Ekman)
    pin=find(pre==Pressure_Ekman(i));
    Total_TPUD_bpr(pin,:)=Absolute_TPUD_bpr(pin,:)+ekmanN_AtoZ(:)'./60;
end
Total_TPUD_bpr=Total_TPUD_bpr+shelf_profile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variant 2: BAROCLINIC ONLY (no reference correction at all,
%%% referenced to the sea surface -- Offset=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Absolute_TPUD_baroclinic=RelTPUD_pilot.*AREA_TOPO_pilot;
Total_TPUD_baroclinic=Absolute_TPUD_baroclinic;
for i=1:numel(Pressure_Ekman)
    pin=find(pre==Pressure_Ekman(i));
    Total_TPUD_baroclinic(pin,:)=Absolute_TPUD_baroclinic(pin,:)+ekmanN_AtoZ(:)'./60;
end
Total_TPUD_baroclinic=Total_TPUD_baroclinic+shelf_profile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kersale definition applied to each variant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MOCup_raw=struct();
labels={'BPR reference (v5, normal)','Baroclinic only (no reference)'};
variants={Total_TPUD_bpr,Total_TPUD_baroclinic};
for vv=1:2
    Psi=cumsum(variants{vv},1,'omitnan')*dz/1e6;
    Psi_mean=nanmean(Psi,2);
    [~,imax]=max(Psi_mean);
    MOCup_raw.(sprintf('v%d',vv))=Psi(imax,:);
    fprintf('%-32s h_star=%d dbar, mean=%.2f Sv, std=%.2f Sv\n',labels{vv},pre(imax),nanmean(Psi(imax,:)),nanstd(Psi(imax,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare each against the full array, 2013-2017, RAW daily
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_full=load('moc_streamfunction_v4.mat','dt','MOCup_raw');
[dt_common,ia,ib]=intersect(dt,s_full.dt);
w1=datenum(2013,9,11); w2=datenum(2017,7,17);
win=dt_common>=w1&dt_common<=w2;

full_win=s_full.MOCup_raw(ib(win));
fprintf('\n--- RAW daily correlation with full array, 2013-2017 ---\n');
for vv=1:2
    pv=MOCup_raw.(sprintf('v%d',vv));
    pv_win=pv(ia(win));
    cc=corrcoef(full_win,pv_win,'rows','complete');
    fprintf('%-32s r=%.3f\n',labels{vv},cc(1,2));
end

figure('Position',[100 100 950 500])
plot(dt_common(win),full_win-nanmean(full_win),'linewidth',0.8); hold on
for vv=1:2
    pv=MOCup_raw.(sprintf('v%d',vv));
    pv_win=pv(ia(win));
    plot(dt_common(win),pv_win-nanmean(pv_win),'linewidth',0.8)
end
legend(['Full array' labels])
title('Pilot component test: BPR reference vs. baroclinic-only, vs. full array','fontsize',12)
ylabel('MOC_{up} anomaly (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',12)
print -dpng -r150 moc_pilot_component_test_IES
