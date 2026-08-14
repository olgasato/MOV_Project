% Step 2 (heat): MHT calculation using Marion's Total Velocities,
% analogous to mov_samba_marion*.m but for heat transport instead of
% freshwater transport. Adapted from Marion's
% Marions_code/MHT/MHT_Estimate_constituents.m (the "final step" of
% README_MHT), which computes total heat transport (not yet decomposed
% into mean+overturning+gyre) plus a sensitivity/attribution analysis
% (how much of the variability comes from Ekman vs. the relative/
% baroclinic vs. the reference/barotropic velocity component).
%
% Citations: see mov_samba_marion*.m's header -- same underlying
% method/array, Marion's 2021 paper (https://doi.org/10.1029/2020JC016947).
%
% v1 (2026-08-13): first extends README_MHT's original pipeline from
% 2013-2017 to 2013-2022, rather than jumping straight to a mean/
% overturning/gyre decomposition (that's a planned v2). Verified first
% (see CLAUDE.md) that running MHT_Estimate_constituents.m completely
% unmodified (just an addpath fix) exactly reproduces the existing
% reference output (Marions_code/MHT/Wrk/MHT_estimates_Constituents.mat,
% mean(Heat_total)=0.5466 PW, 2013-09-11 to 2017-07-16) before touching
% anything -- confidence that the method itself is understood correctly
% before extending it.
%
% Three source updates for the 2013-2022 extension:
%  1. West temperature: samba_w.mat (2009-2022) run through the same
%     yearday-indexed seasonal correction as the original West source
%     (Marions_code/ies_profiles/IES_Make_Profiles_USA_SAM_inclSeasGEM_2013_2022.m),
%     instead of the old IES_Make_Profiles_USA_SAM_SeasCorr.mat (capped
%     2018-04-30).
%  2. East temperature: IES_Make_Profiles_FrSA_SAMBA_6PIES.mat
%     (2013-2023, the Daily_Tau-based rebuild) run through the same
%     East seasonal correction
%     (IES_Make_Profiles_FrSA_SAMBA_6PIES_inclSeasGEM_2013_2022.m),
%     instead of the old version (capped 2017-07-21).
%  3. Total_TPUD1..8 (+ RefConst/RelConst/EKMANconst variants): from
%     Step E v2 (Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat,
%     2013-2022), already Ekman/topo-corrected the same way the original
%     MHT script expects.
%
% Deliberate simplification, NOT yet resolved (flagged explicitly in
% CLAUDE.md, discussed with the user before proceeding): gap 3 (D-P8,
% the widest gap, spanning the open-ocean interior between the West and
% East boundary arrays) originally used a SEPARATE, more sophisticated
% "interior" temperature estimate derived from satellite altimetry
% (ALTI_Profiles_SAMBA_Int_SeasCorr.mat) rather than a simple D/P8
% average -- that altimetry pipeline is its own multi-step sub-pipeline
% (README_MHT steps A-E: AVISO SLA, WOA/Argo hydrography, GEM dynamic-
% height-temperature fields) and is itself capped at 2017-08-01, not
% yet extended. Until that's tackled, this version uses
% Temp_3=(Temperature_D+Temperature_P8)/2 instead -- the same simple
% adjacent-site-average convention mov_samba_marion*.m already uses for
% every other gap. This is a real methodological simplification
% relative to Marion's original approach for this specific gap, not
% just a data-currency gap; treat Heat_total's magnitude/trend with
% that in mind until the altimetry pipeline is revisited.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First up load the shelf and deep heat transports (static,
%%% climatological -- no dt, no rebuild needed, same as MOV's Steps C/D)
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
%%% Absolute TPUD (Step E v2, 2013-2022, AREA_TOPO/Ekman consistent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prefix_marion,'Full_Depth_MOC/Wrk/Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2.mat'])
dt_TPUD=dt; clear dt
clear Absolute_TPUD*_preTopo AREA_TOPO_* Ekman_TPUD*  % not needed here, only Total_TPUD*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Align the three sources on their common dt (replaces the original
%%%  script's hardcoded 2013-09-11/2017-07-17 trim -- same 3-way
%%%  intersect() pattern used in concat_IES_v6/mov_samba_marion_v8+/
%%%  Step E v2).
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

varnames={'Total_TPUD1','Total_TPUD2','Total_TPUD3','Total_TPUD4','Total_TPUD5','Total_TPUD6','Total_TPUD7','Total_TPUD8', ...
    'Total_TPUD1_RefConst','Total_TPUD2_RefConst','Total_TPUD3_RefConst','Total_TPUD4_RefConst','Total_TPUD5_RefConst','Total_TPUD6_RefConst','Total_TPUD7_RefConst','Total_TPUD8_RefConst', ...
    'Total_TPUD1_RelConst','Total_TPUD2_RelConst','Total_TPUD3_RelConst','Total_TPUD4_RelConst','Total_TPUD5_RelConst','Total_TPUD6_RelConst','Total_TPUD7_RelConst','Total_TPUD8_RelConst', ...
    'Total_TPUD1_EKMANconst','Total_TPUD2_EKMANconst','Total_TPUD3_EKMANconst','Total_TPUD4_EKMANconst','Total_TPUD5_EKMANconst','Total_TPUD6_EKMANconst','Total_TPUD7_EKMANconst','Total_TPUD8_EKMANconst'};
for vv=1:numel(varnames)
    eval([varnames{vv},'=',varnames{vv},'(:,i_TPUD);']);
end
clear vv varnames dt_SAM dt_GH dt_TPUD i_SAM i_GH i_TPUD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Throw out all data below 4700 dbar at all sites (matches the
%%%  original script's cutoff exactly)
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
        'Total_TPUD1_RefConst','Total_TPUD2_RefConst','Total_TPUD3_RefConst','Total_TPUD4_RefConst','Total_TPUD5_RefConst','Total_TPUD6_RefConst','Total_TPUD7_RefConst','Total_TPUD8_RefConst', ...
        'Total_TPUD1_RelConst','Total_TPUD2_RelConst','Total_TPUD3_RelConst','Total_TPUD4_RelConst','Total_TPUD5_RelConst','Total_TPUD6_RelConst','Total_TPUD7_RelConst','Total_TPUD8_RelConst', ...
        'Total_TPUD1_EKMANconst','Total_TPUD2_EKMANconst','Total_TPUD3_EKMANconst','Total_TPUD4_EKMANconst','Total_TPUD5_EKMANconst','Total_TPUD6_EKMANconst','Total_TPUD7_EKMANconst','Total_TPUD8_EKMANconst'}
    eval([vv{1},'(bad,:)=[];']);
end
Pressure(bad)=[];
clear bad vv

assert(isequal(Pressure_SAM,Pressure_GH) && isequal(Pressure_SAM,Pressure), ...
    'Pressure_SAM/Pressure_GH/Pressure differ after the >4700dbar trim -- they must share the same grid for Temp_3''s D/P8 average and the Q1..Q8 integration below to line up depth-for-depth.');
Pressure=Pressure_SAM;
clear Pressure_SAM Pressure_GH

% Per-day count of how many of the 9 original sites have at least one
% valid depth -- same coverage diagnostic as mov_samba_marion_v13+
% (n_sites_valid), used below to shade the plot the same way
% mov_samba_marion_v15 does (this is the same array/outages, so the
% same P5/P6/P4/P8 gaps apply here too).
sal00_like=cat(3,Temperature_A,Temperature_C,Temperature_D,Temperature_P8, ...
    Temperature_P6,Temperature_P5,Temperature_P4,Temperature_P2,Temperature_P1);
n_sites_valid=reshape(sum(squeeze(any(~isnan(sal00_like),1)),2),1,numel(dt));
clear sal00_like

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Build the 8 gap-midpoint temperatures. Temp_3 (D-P8, the widest
%%%  gap) uses a direct D/P8 average here instead of the original
%%%  altimetry-based "interior" estimate -- see the header note above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp_1=(Temperature_A+Temperature_C)./2;
Temp_2=(Temperature_C+Temperature_D)./2;
Temp_3=(Temperature_D+Temperature_P8)./2;
Temp_4=(Temperature_P8+Temperature_P6)./2;
Temp_5=(Temperature_P6+Temperature_P5)./2;
Temp_6=(Temperature_P5+Temperature_P4)./2;
Temp_7=(Temperature_P4+Temperature_P2)./2;
Temp_8=(Temperature_P2+Temperature_P1)./2;
clear Temperature_A Temperature_C Temperature_D Temperature_P1 Temperature_P2 Temperature_P4 Temperature_P5 Temperature_P6 Temperature_P8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Compute MHT: same math as the original MHT_Estimate_constituents.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=1025; %kg/m3
cp=4000; %J/(kg C)

z=sw_dpth(Pressure,-34.5);
diffz=diff(z);

Temp_all={Temp_1,Temp_2,Temp_3,Temp_4,Temp_5,Temp_6,Temp_7,Temp_8};

for gg=1:8
    Temp_g=Temp_all{gg};
    TPUD=eval(['Total_TPUD',num2str(gg)]);
    TPUD_Ek=eval(['Total_TPUD',num2str(gg),'_EKMANconst']);
    TPUD_Rel=eval(['Total_TPUD',num2str(gg),'_RelConst']);
    TPUD_Ref=eval(['Total_TPUD',num2str(gg),'_RefConst']);

    Temp_mid=(Temp_g(1:end-1,:)+Temp_g(2:end,:))/2;
    TPUD_mid=(TPUD(1:end-1,:)+TPUD(2:end,:))/2;
    TPUD_Ek_mid=(TPUD_Ek(1:end-1,:)+TPUD_Ek(2:end,:))/2;
    TPUD_Rel_mid=(TPUD_Rel(1:end-1,:)+TPUD_Rel(2:end,:))/2;
    TPUD_Ref_mid=(TPUD_Ref(1:end-1,:)+TPUD_Ref(2:end,:))/2;

    Q(gg,:)=nansum(rho.*cp.*TPUD_mid.*Temp_mid.*diffz(:),1); %#ok<SAGROW>
    Q_Ek(gg,:)=nansum(rho.*cp.*TPUD_Ek_mid.*Temp_mid.*diffz(:),1); %#ok<SAGROW>
    Q_Rel(gg,:)=nansum(rho.*cp.*TPUD_Rel_mid.*Temp_mid.*diffz(:),1); %#ok<SAGROW>
    Q_Ref(gg,:)=nansum(rho.*cp.*TPUD_Ref_mid.*Temp_mid.*diffz(:),1); %#ok<SAGROW>
end
clear gg Temp_g TPUD TPUD_Ek TPUD_Rel TPUD_Ref Temp_mid TPUD_mid TPUD_Ek_mid TPUD_Rel_mid TPUD_Ref_mid
clear Temp_1 Temp_2 Temp_3 Temp_4 Temp_5 Temp_6 Temp_7 Temp_8 Temp_all
clear Total_TPUD1 Total_TPUD2 Total_TPUD3 Total_TPUD4 Total_TPUD5 Total_TPUD6 Total_TPUD7 Total_TPUD8
clear Total_TPUD1_RefConst Total_TPUD2_RefConst Total_TPUD3_RefConst Total_TPUD4_RefConst Total_TPUD5_RefConst Total_TPUD6_RefConst Total_TPUD7_RefConst Total_TPUD8_RefConst
clear Total_TPUD1_RelConst Total_TPUD2_RelConst Total_TPUD3_RelConst Total_TPUD4_RelConst Total_TPUD5_RelConst Total_TPUD6_RelConst Total_TPUD7_RelConst Total_TPUD8_RelConst
clear Total_TPUD1_EKMANconst Total_TPUD2_EKMANconst Total_TPUD3_EKMANconst Total_TPUD4_EKMANconst Total_TPUD5_EKMANconst Total_TPUD6_EKMANconst Total_TPUD7_EKMANconst Total_TPUD8_EKMANconst

Heat_total=sum(Q,1)+Q_deep+Q_shelfE+Q_shelfW;
Heat_total_RelConst=sum(Q_Rel,1)+Q_deep+Q_shelfE+Q_shelfW;
Heat_total_RefConst=sum(Q_Ref,1)+Q_deep+Q_shelfE+Q_shelfW;
Heat_total_EKMANconst=sum(Q_Ek,1)+Q_deep+Q_shelfE+Q_shelfW;

Heat_total_EkmanAnomaly=Heat_total-Heat_total_EKMANconst;
Heat_total_RelativeAnomaly=Heat_total-Heat_total_RelConst;
Heat_total_ReferenceAnomaly=Heat_total-Heat_total_RefConst;

fprintf('mean(Heat_total) = %.4f PW, std = %.4f PW, %d days (%s to %s)\n', ...
    nanmean(Heat_total)*1e-15, nanstd(Heat_total)*1e-15, numel(dt), datestr(dt(1)), datestr(dt(end)));

disp('Corr. Ekman')
[c2]=corrcoef(Heat_total,Heat_total_EkmanAnomaly)
disp('Corr. Barocline')
[c2]=corrcoef(Heat_total,Heat_total_RelativeAnomaly)
disp('Corr. Barotrope')
[c2]=corrcoef(Heat_total,Heat_total_ReferenceAnomaly)

% Graded coverage shading, same bands/colors as mov_samba_marion_v15
% (same array, same P5/P6/P4/P8 outages) -- Heat_total etc. are NOT
% modified or masked, only the plot background is shaded.
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

figure
hold on
ylim0=[-17 6];
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
p=plot(dt,Heat_total*1e-15,dt,Heat_total_EkmanAnomaly*1e-15-5, ...
    dt,Heat_total_RelativeAnomaly*1e-15-10,dt,Heat_total_ReferenceAnomaly*1e-15-15);
legend([p(:);legend_h(:)],[{'Heat total','Ekman anomaly - 5','Relative(baroclinic) anomaly - 10','Reference(barotropic) anomaly - 15'},legend_str], ...
    'Location','southoutside','Orientation','horizontal','fontsize',7)
title('MHT at 34.5^oS from IES (2013-2022)','fontsize',14)
ylabel('PW')
datetick('x',12);grid;axis('tight')
ylim(ylim0)
set(gca,'linewidth',1,'fontsize',14)
box on
print -dpng -r150 mht_IES

save('mht_samba_marion_v1.mat','dt','Heat_total','Heat_total_EkmanAnomaly', ...
    'Heat_total_ReferenceAnomaly','Heat_total_RelativeAnomaly','Heat_total_EKMANconst', ...
    'Heat_total_RelConst','Heat_total_RefConst','n_sites_valid','category')
