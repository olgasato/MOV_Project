% Step 2 (volume), take 5: points at Step E `v3` (SAMBA_E_IES) instead
% of `v2`, to fix the ECCO calibration-window confound found while
% isolating West-vs-East reprocessing effects (see CLAUDE.md's "West-vs-
% East isolation done properly" section).
%
% `v3`'s own ECCO time-mean reference offset step calibrates each gap
% against THIS SCRIPT'S OWN date range's mean at 1500dbar -- since `v2`
% spans 2013-2022, every gap's absolute reference was silently anchored
% to that whole-record mean, not to the 2013-09-11/2017-07-17 window
% Marion's original script (and the published Kersale et al. 2021
% numbers) actually used. Step E `v3`
% (Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3.m, in
% SAMBA_E_IES) fixes this by restricting that calibration to the
% original window regardless of the script's own date range -- so this
% script's full-array output stays comparable to the paper/Marion's
% original, no matter how far the record has been extended.
%
% Note: moc_pilot_v4.m's OWN Pilot calculation is NOT affected by this
% confound at all -- its BPR-based reference is computed independently,
% with no dependency on Step E's ECCO offset step. Only the full-array
% side needed rebuilding; moc_pilot_v4.m itself doesn't need a new
% version, just a rerun of its comparison section against this script's
% new output (see moc_pilot_v5.m).
%
% v4 (2026-08-17): identical to v3 except loading Step E `v3`'s output
% instead of `v2`'s -- see v2's header for the full Kersale-definition
% derivation (h_star, raw cumulative transport, no depth-mean removal)
% and v3's header for the shelf-transport fix, both unchanged here.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/')

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3']);

fprintf('dt range: %s to %s (%d days)\n', datestr(dt(1)), datestr(dt(end)), numel(dt));

AREA_TOPO_1=AREA_TOPO_1(:); AREA_TOPO_2=AREA_TOPO_2(:); AREA_TOPO_3=AREA_TOPO_3(:); AREA_TOPO_4=AREA_TOPO_4(:);
AREA_TOPO_5=AREA_TOPO_5(:); AREA_TOPO_6=AREA_TOPO_6(:); AREA_TOPO_7=AREA_TOPO_7(:); AREA_TOPO_8=AREA_TOPO_8(:);

% --------------------------------------------------------
% Distance between sites (dx, per gap) -- identical to v1/v2/v3
position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);

% --------------------------------------------------------
% vel00/geo/dxA/dx2/V0 -- identical construction to v1/v2/v3
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
ind=find(isnan(geo) | isnan(AREA_TOPO_3d));
dx1(ind)=nan;
geo(isnan(dx1))=nan;
dxA=dx1.*AREA_TOPO_3d;
dx2=squeeze(nansum(dxA,3));

n_gaps_valid=reshape(sum(any(~isnan(dxA),1),3),1,m);

V0=nansum(vel00.*dxA,3)./nansum(dxA,3);

pre=Pressure; clear Pressure
dz=mean(diff(pre));

% --------------------------------------------------------
% Shelf transport correction (West_TPUD+East_TPUD_new, static,
% time-invariant -- broadcasts identically across every day). Matched
% to `pre` by pressure VALUE (shelf profile is on 0:10:5300dbar, 531
% levels; `pre` is already trimmed to <=5150dbar, 516 levels).
load('/Users/olga/research/sambar/renellys_sent/Marions_code/ecco/Wrk/ECCO_TPUD_ShelfBits.mat','West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(pre,pres_shelf_grid);
assert(isequal(ia(:),(1:numel(pre))'),'pre is not a subset of the shelf grid as expected');
shelf_profile=West_TPUD(ib)+East_TPUD_new(ib);
fprintf('Shelf correction: West=%.2f Sv, East=%.2f Sv (static, added to every day)\n', ...
    nansum(West_TPUD)*10/1e6, nansum(East_TPUD_new)*10/1e6);
clear West_TPUD East_TPUD_new pres_shelf_grid ia ib

% --------------------------------------------------------
% Kersale's actual definition: RAW cumulative transport (incl. shelf),
% no depth-mean removal, evaluated at a FIXED depth found from the
% time-mean profile.
transport_profile=V0.*dx2+shelf_profile; % raw absolute transport-per-unit-depth, all 8 gaps + shelves
Psi_raw=cumsum(transport_profile,1,'omitnan')*dz/1e6; % [depth x time], Sv

Psi_raw_mean=nanmean(Psi_raw,2);
[~,imax_mean]=max(Psi_raw_mean);
h_star=pre(imax_mean);
fprintf('h_star (time-mean transition depth, MOCup base): %d dbar (paper reports 1315 dbar for their 2013-2017 window)\n',h_star);

MOCup_raw=Psi_raw(imax_mean,:);
MOCup_low=lowpass_filter(MOCup_raw,45);
fprintf('RAW (unfiltered) MOCup at h_star: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOCup_raw), nanstd(MOCup_raw));
fprintf('LOWPASS (45-day) MOCup at h_star: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOCup_low), nanstd(MOCup_low));

% --------------------------------------------------------
% Hovmoller of the raw (no depth-mean-removal) time-mean-lowpassed Psi.
Psi_low=nan(size(Psi_raw));
for zz=1:size(Psi_raw,1)
    Psi_low(zz,:)=lowpass_filter(Psi_raw(zz,:),45);
end

figure('Position',[100 100 900 500])
imagesc(dt,pre,Psi_low)
set(gca,'YDir','reverse')
datetick('x',12,'keeplimits')
colormap(redblue_or_default())
cb=colorbar; ylabel(cb,'Sv')
xlabel('Date'); ylabel('Pressure (dbar)')
hold on
plot(get(gca,'xlim'),[h_star h_star],'k--','linewidth',1)
title('MOC cumulative transport \Psi(z,t) at 34.5^oS (raw+shelf, calib-fixed, 45-day lowpass)','fontsize',13)
clim_val=prctile(abs(Psi_low(:)),99);
caxis([-clim_val clim_val])
print -dpng -r150 moc_streamfunction_v4_hovmoller

figure
plot(dt,MOCup_raw,dt,MOCup_low,'r','linewidth',1.5)
legend('raw daily','45-day lowpass')
title('MOCup index (transport at fixed h\_star, +shelf, calib-fixed) at 34.5^oS','fontsize',14)
ylabel('MOC_{up} (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_streamfunction_v4_index_IES

% --------------------------------------------------------
% Time-mean Psi(z) vs. depth -- Kersale Figure 2a style.
Psi_std=nanstd(Psi_low,0,2);

figure('Position',[500 100 650 650])
hold on
fill([Psi_raw_mean-Psi_std;flipud(Psi_raw_mean+Psi_std)],[pre;flipud(pre)],[1 0.8 0.8],'EdgeColor','none')
plot(Psi_raw_mean,pre,'r','linewidth',2)
plot([0 0],[0 max(pre)],'k--')
plot(xlim,[h_star h_star],'k:','linewidth',1)
text(0,h_star-100,sprintf('h_{star}=%d dbar',h_star),'fontsize',10)
set(gca,'YDir','reverse')
xlabel('\Psi (Sv)'); ylabel('Pressure (dbar)')
title({'Time-mean MOC cumulative transport vs. depth at 34.5^oS','Calibration-window-fixed + shelf (shading = +/-1 std)'},'fontsize',12)
grid on
print -dpng -r150 moc_streamfunction_v4_mean_profile

save('moc_streamfunction_v4.mat','dt','pre','Psi_raw','Psi_low','h_star','MOCup_raw','MOCup_low','n_gaps_valid','Psi_raw_mean','Psi_std')

function cmap=redblue_or_default()
if exist('redblue','file')
    cmap=redblue(64);
else
    cmap=jet(64);
end
end
