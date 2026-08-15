% MOC "Pilot" method: a two-endpoint (site A West edge, site P1 East
% edge) dynamic-height-difference estimate of the overturning
% streamfunction, replicating the red "Pilot" curve in Kersale et al.
% (2021) Figure 2a for comparison against moc_streamfunction_v1.m's
% full 8-gap streamfunction. Requested after seeing that figure: "tem
% uma curva com a estimativa da MOC que eles chamam de 'Pilot'. Esse
% metodo considera a diferenca da altura dinamica entre as duas bordas
% da bacia... ao inves de uma integral considerando todos os pontos
% disponiveis." User's methodological choice, given a 3-way clarifying
% question: use BPR (bottom pressure) absolute referencing at BOTH A
% and P1 (not just baroclinic/relative shear with a level-of-no-motion
% assumption).
%
% This is NOT built from mov_samba_marion/mht_samba_marion/
% moc_streamfunction_v1's Step E v2 output directly -- that .mat only
% saves per-gap Total_TPUD/Absolute_TPUD_preTopo (already the
% BPR-absolute-referenced 8-gap decomposition), not the intermediate
% Gpan/pres/rhob quantities needed to build a NEW, single "gap"
% spanning the whole A-to-P1 width in one hop. So this script reloads
% the same raw sources Step E v2 does (samba_w.mat for A;
% ies/Wrk/IES_FrSA_SAMBA_6PIES.mat + ies_profiles/Wrk/
% IES_Make_Profiles_FrSA_SAMBA_6PIES.mat for P1;
% Ekman_transports_9pies.mat for the whole-basin ekmanN_AtoZ) and
% reruns just the A/P1 slice of Step E v2's own Gpan/RefTPUD/RelTPUD/
% Offset machinery.
%
% Formula choice: mirrors gap 1 (A-C) of Step E v2 EXACTLY, substituting
% P1 for C, rather than deriving a fresh formula or guessing which of
% the 8 gaps' (non-uniform, direction-dependent) sign/ordering patterns
% to copy. This is a safe structural analogy because A and P1 have the
% SAME shallow/deep relationship as A and C do: p_cpiesW(1)=1370dbar
% (A) and p_cpiesE(1)=1280dbar (P1) -- P1 is in fact the SHALLOWEST
% East site (not the deepest, despite being the far-East edge; the
% array's bottom depths are non-monotonic along the section), so P1
% plays the "shallow, reference-index" role A plays in gap 1, and A
% (1370dbar, deeper than P1) plays the "deep, own-profile-for-
% correction" role C plays in gap 1:
%   Gpan_AtoP1_corr = Gpan_A(p_ind_A,:) - Gpan_A(p_ind_P1,:)   [A's own
%       profile evaluated at both bottom depths, demeaned]
%   pres_AtoP1 = pres_A - Gpan_AtoP1_corr.*rhob_A               [project
%       A's bottom pressure up to P1's shallower depth]
%   RefTPUD_pilot = (pres_AtoP1 - pres_P1) / (f*rhob_P1)
%   RelTPUD_pilot = -(1/f)*(Gpan_P1 - Gpan_A)                  [same
%       east-minus-west convention + sign flip used for all 8 gaps]
%   Offset_pilot = RefTPUD_pilot - RelTPUD_pilot(p_ind_P1,:)    [applied
%       at the shallower site's index, P1]
%
% No AREA_TOPO correction: there is no pre-computed topographic-width
% profile for a single A-to-P1 span (topo_corr_msm60.mat only has
% per-adjacent-gap AREA_TOPO_1..8). Using the plain geographic distance
% (no fine-bathymetry open-width reduction) is a genuine simplification
% versus the full 8-gap method, but is consistent with "Pilot" being
% the simpler/coarser of the two methods being compared -- not an
% oversight.
%
% Streamfunction: same construction as moc_streamfunction_v1.m --
% remove the depth-mean (barotropic) component so the profile
% integrates to ~0 at the deepest common level, cumsum from the
% surface down, 45-day lowpass, MOC(t)=max_z Psi(z,t). Since there is
% only one "gap" here (no zonal averaging across multiple gaps), the
% depth-mean removal is a simple nanmean over depth of the A-P1
% velocity profile itself (dx_pilot has no depth dependence, unlike the
% full method's AREA_TOPO-weighted dx2).
%
% v1 (2026-08-15): first cut.

clear;clc;

prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';
addpath([prefix_marion,'functions/seawater/']);
addpath([prefix_marion,'functions/positions_pies']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ekman (whole-basin A-to-P1 total, ekmanN_AtoZ, already computed by
%%% Ekman_transports_9pies_2009_2022.m alongside the per-gap ekmanN_i)
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
%%% Positions (for A-P1 distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position_West
position_East
lon=[cpiesW(1).lon cpiesE(1).lon];
lat=[cpiesW(1).lat cpiesE(1).lat];
dx_pilot=gsw_distance(lon,lat,0);
fprintf('A-to-P1 direct distance: %.1f km\n',dx_pilot/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Geopotential anomaly (relative to sea surface)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gpan_A=sw_gpan(Salinity_A,Temperature_A,Pressure_SAM(:));
Gpan_P1=sw_gpan(Salinity_P1,Temperature_P1,Pressure_GH(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bottom pressure indices/densities -- A=1370dbar (West), P1=1280dbar
%%% (East, the SHALLOWEST East site despite being the far-East edge;
%%% confirmed from Step E v2's p_cpiesE=[1280,2150,2890,4560,5060,
%%% 5290,5130,4690] with p_ind_P1=find(Pressure_GH==p_cpiesE(1))).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ind_A=find(Pressure_SAM==1370);
p_ind_P1=find(Pressure_GH==1280);

rhob_A=nanmean(sw_dens(Salinity_A(p_ind_A,:),Temperature_A(p_ind_A,:),Pressure_SAM(p_ind_A)));
rhob_P1=nanmean(sw_dens(Salinity_P1(p_ind_P1,:),Temperature_P1(p_ind_P1,:),Pressure_GH(p_ind_P1)));
fprintf('rhob_A=%.2f, rhob_P1=%.2f\n',rhob_A,rhob_P1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correct A's (the deeper site's) own pressure record up to P1's
%%% shallower depth, using A's own Gpan profile -- exact mirror of
%%% gap 1's Gpan_CtoA_corr/pres_CtoA.
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
%%% Align West/East/Ekman on their common dt (same 3-way intersect
%%% pattern as Step E v2 -- gives the same dt Step E v2/
%%% moc_streamfunction_v1.m end up with, since it's the identical
%%% West/East/Ekman sources and identical intersect logic).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dt_tmp,i_SAM,i_GH]=intersect(dt_SAM,dt_GH);
[dt,i_tmp,i_ekman]=intersect(dt_tmp,dt_ekman);
i_SAM=i_SAM(i_tmp); i_GH=i_GH(i_tmp);
clear dt_tmp i_tmp

fprintf('Aligned dt: %d common days (%s to %s) out of %d (West), %d (East), %d (Ekman)\n', ...
    numel(dt),datestr(dt(1)),datestr(dt(end)),numel(dt_SAM),numel(dt_GH),numel(dt_ekman));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relative (baroclinic, surface-referenced) transport per unit depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RelTPUD_pilot=(1./sw_f(-34.5)).*(Gpan_P1 - Gpan_A);
RelTPUD_pilot=-RelTPUD_pilot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reference (BPR-absolute) transport per unit depth, applied at P1's
%%% (the shallower site's) index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RefTPUD_pilot=((pres_AtoP1 - pres_P1))./(sw_f(-34.5).*rhob_P1);
RefTPUD_pilot=RefTPUD_pilot-nanmean(RefTPUD_pilot);

Offset_pilot=RefTPUD_pilot-RelTPUD_pilot(p_ind_P1,:);

Absolute_TPUD_pilot=NaN*ones(size(RelTPUD_pilot));
for i=1:length(Pressure)
    Absolute_TPUD_pilot(i,:)=RelTPUD_pilot(i,:)+Offset_pilot;
end
clear i p_ind_A p_ind_P1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add Ekman at the 7 surface rows (0-60dbar), same construction as
%%% Step E v2's Ekman_TPUD_i (ekmanN/60 spread across Pressure_Ekman)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pressure_Ekman=[0:10:60]';
Total_TPUD_pilot=Absolute_TPUD_pilot;
for i=1:length(Pressure_Ekman)
    pin=find(Pressure==Pressure_Ekman(i));
    Total_TPUD_pilot(pin,:)=Absolute_TPUD_pilot(pin,:)+ekmanN_AtoZ(:)'./60;
end
clear i pin Absolute_TPUD_pilot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Streamfunction: remove the depth-mean (barotropic) component so
%%% the profile integrates to ~0 at depth, matching
%%% moc_streamfunction_v1.m's vel_prime construction -- here there's
%%% only one "gap" so the depth-mean is a plain nanmean over depth of
%%% the A-P1 velocity (dx_pilot has no depth dependence, no AREA_TOPO).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vel_pilot=Total_TPUD_pilot/dx_pilot;
V00_pilot=nanmean(vel_pilot,1);
vel_prime_pilot=vel_pilot-V00_pilot;
transport_profile_pilot=vel_prime_pilot*dx_pilot;

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
fprintf('LOWPASS-FILTERED (45-day) Pilot MOC: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOC_pilot_low), nanstd(MOC_pilot_low));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Overlay against the full 8-gap streamfunction MOC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=load('moc_streamfunction_v1.mat','dt','MOC_low');
[dt_common,ia,ib]=intersect(dt,s.dt);
fprintf('Overlay common dt: %d days (%s to %s)\n',numel(dt_common),datestr(dt_common(1)),datestr(dt_common(end)));

figure('Position',[100 100 900 500])
plot(dt_common,s.MOC_low(ib),'linewidth',1.5); hold on
plot(dt_common,MOC_pilot_low(ia),'r','linewidth',1.5)
legend('Full array (8 gaps)','Pilot (A-P1 only)')
title('MOC index at 34.5^oS: full array vs. Pilot method (45-day lowpass)','fontsize',14)
ylabel('MOC (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_pilot_vs_full_IES

fprintf('corr(full,pilot)=%.3f over the common period\n', ...
    corr2_nan(s.MOC_low(ib),MOC_pilot_low(ia)));

save('moc_pilot_v1.mat','dt','pre','Psi_pilot','Psi_pilot_low','MOC_pilot_low','MOC_pilot_depth_low','dx_pilot')

function r=corr2_nan(a,b)
ok=~isnan(a)&~isnan(b);
cc=corrcoef(a(ok),b(ok));
r=cc(1,2);
end
