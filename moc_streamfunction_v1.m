% Step 2 (volume), take 2: MOC via the classical overturning
% STREAMFUNCTION definition (RAPID-26degN style), after
% moc_samba_marion_v1.m's "reuse mov_samba_marion, drop the tracer"
% approach turned out to be mathematically guaranteed to give zero for
% both the overturning and gyre terms -- see CLAUDE.md's MOC section
% for the full derivation of why. This script instead computes:
%
%   Psi(z,t) = -cumsum_{z'=0}^{z} [ V0(z',t) * dx2(z',t) ] dz'  / 1e6
%
% i.e. the cumulative (surface-downward) volume transport streamfunction
% across the WHOLE section, and MOC(t) = max_z Psi(z,t) -- the depth of
% maximum cumulative transport, the standard "MOC strength" diagnostic
% (e.g. "17 Sv at 26N" from the RAPID array). V0(z,t) is the same
% AREA_TOPO/dxA-weighted zonal-mean velocity profile
% mov_samba_marion_v12+/mht_samba_marion_v2 already compute; V0.*dx2 is
% the total (all 8 gaps) volume-transport-per-unit-depth profile, and
% its cumulative sum from the surface down is exactly what a
% streamfunction is.
%
% Sign convention: matches this project's established "-V00*width" sign
% (mov_samba_marion's mean_term, mht_samba_marion's Q_mean) for
% consistency, applied here to the whole depth-resolved profile rather
% than a single depth-integrated scalar.
%
% v1 (2026-08-14): first cut. Loads Step E v2 exactly like
% mov_samba_marion_v12+/mht_samba_marion_v2 (Absolute_TPUD*_preTopo +
% AREA_TOPO + Ekman_TPUD, no salinity/temperature needed at all -- see
% CLAUDE.md's abandoned-MOC section for why this is structurally
% simpler than mov_samba_marion). Plots a depth-time Hovmoller of
% Psi(z,t) first, to visually confirm the vertical structure (single
% extremum vs. more complex shape) before committing to max(Psi) as
% "the" MOC index -- standard in the literature but worth checking
% against this specific section/array rather than assuming.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/')

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2']);

fprintf('dt range: %s to %s (%d days)\n', datestr(dt(1)), datestr(dt(end)), numel(dt));

% AREA_TOPO_i are [1x516] row vectors -- force to columns (same fix as
% mht_samba_marion_v2.m needed).
AREA_TOPO_1=AREA_TOPO_1(:); AREA_TOPO_2=AREA_TOPO_2(:); AREA_TOPO_3=AREA_TOPO_3(:); AREA_TOPO_4=AREA_TOPO_4(:);
AREA_TOPO_5=AREA_TOPO_5(:); AREA_TOPO_6=AREA_TOPO_6(:); AREA_TOPO_7=AREA_TOPO_7(:); AREA_TOPO_8=AREA_TOPO_8(:);

% --------------------------------------------------------
% Distance between sites (dx, per gap)
position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);

% --------------------------------------------------------
% vel00/geo from Absolute_TPUD*_preTopo, Ekman restored at the surface
% rows, dxA (AREA_TOPO-weighted real open width) -- identical
% construction to mov_samba_marion_v12+/mht_samba_marion_v2.
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
V00=nansum(V0.*dx2)./nansum(dx2);
vel_prime=V0-V00;

pre=Pressure; clear Pressure
dz=mean(diff(pre));

% --------------------------------------------------------
% Streamfunction. v1's FIRST DRAFT used V0 (the full zonal-mean
% velocity, including the barotropic/net component) here -- Psi(bottom)
% then exactly equals mean_term/moc_mean by construction (verified:
% mean(Psi(end,:))=6.285 Sv, matching the abandoned moc_samba_marion_v1
% exactly), which sounds like a good sanity check, but it means every
% individual day's Psi(z) profile is riding on top of that day's own
% (often large -- mean_term's std is ~40 Sv vs. a ~6 Sv mean, already
% documented from the 2017-2019 coverage-artifact investigation)
% barotropic offset, which swamps the baroclinic "hump" structure and
% never reliably closes back toward zero at depth on any GIVEN day
% (e.g. 2018-09-23 monotonically climbed to 177 Sv). A 45-day lowpass
% on Psi(z,:) barely helped (still ~24 Sv std) -- the barotropic
% component is simply too energetic to filter away after the fact.
%
% Fixed by building the streamfunction from vel_prime=V0-V00 (the
% baroclinic/overturning deviation mov_samba_marion*.m already computes
% this same way) instead of V0. Since V00 is *by definition* the
% dx2-weighted depth-mean of V0, sum_z(vel_prime.*dx2)=0 EXACTLY (same
% identity that made the abandoned moc_samba_marion_v1's overturning
% term trivially zero -- see CLAUDE.md) -- but applied here to just the
% BOTTOM boundary condition of a cumulative sum, not the whole profile,
% this is exactly right: it guarantees Psi_baroclinic(deepest depth)=0
% by construction, cleanly isolating the overturning signal from the
% noisy barotropic offset. This is, not coincidentally, much closer to
% how MOC streamfunctions are actually defined in the literature
% (relative to a reference/no-net-flow level).
%
% Sign: NO leading minus here, unlike mean_term/mov/Q_mean elsewhere in
% this project. Plotted both ways (moc_streamfunction_mean_profile.png)
% at the user's prompt ("a camada superior deveria dar positiva... sera
% que e problema de convencao de sinal?") -- with a leading "-" (this
% project's usual convention for V00-derived quantities), the mean
% profile is negative through most of the water column, peaking around
% -20 Sv near 1000-1500dbar; WITHOUT it, the profile matches the
% expected/standard MOC shape (e.g. Kersale et al. 2021's Figure 2):
% ~0 at the surface, growing to +20 Sv around 1000-1500dbar (upper
% cell), tapering back down through the NADW layer, small negative dip
% near the very bottom (hint of an abyssal cell). Adopted the
% no-leading-minus sign as correct on that basis.
transport_profile=vel_prime.*dx2; % baroclinic transport-per-unit-depth, all 8 gaps
Psi=cumsum(transport_profile,1,'omitnan')*dz/1e6; % [depth x time], Sv

fprintf('Sanity check: mean(Psi at max depth)=%.4g Sv (should be ~0 by construction)\n', nanmean(Psi(end,:)));

[MOC_raw,imax]=max(Psi,[],1);
MOC_depth_raw=pre(imax)';
fprintf('RAW (unfiltered) daily MOC: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOC_raw), nanstd(MOC_raw));

% Still apply the same 45-day lowpass this project uses elsewhere
% (mov_low=lowpass_filter(mov,45)) -- daily snapshots are still noisy
% (mesoscale/synoptic variability in the baroclinic field itself), just
% no longer dominated by the barotropic component specifically.
Psi_low=nan(size(Psi));
for zz=1:size(Psi,1)
    Psi_low(zz,:)=lowpass_filter(Psi(zz,:),45);
end
[MOC_low,imax_low]=max(Psi_low,[],1);
MOC_depth_low=pre(imax_low)';
fprintf('LOWPASS-FILTERED (45-day) daily MOC: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOC_low), nanstd(MOC_low));
fprintf('LOWPASS MOC_depth: mean=%.1f dbar, median=%.1f dbar, range=[%.0f,%.0f]\n', ...
    nanmean(MOC_depth_low), nanmedian(MOC_depth_low), min(MOC_depth_low), max(MOC_depth_low));

% --------------------------------------------------------
% Hovmoller of the LOWPASS-FILTERED Psi(z,t). YDir set to 'reverse' --
% v1's first draft used 'normal', which (since pre runs 0->5150
% ascending) put the SURFACE at the BOTTOM of the figure and the seafloor
% at the TOP, backwards from the standard oceanographic convention
% (depth increasing downward) and from how Kersale et al. (2021)-style
% figures are drawn. Caught while building the mean-profile plot below,
% at the user's prompt ("sera que e problema de convencao de sinal da
% profundidade?").
figure('Position',[100 100 900 500])
imagesc(dt,pre,Psi_low)
set(gca,'YDir','reverse')
datetick('x',12,'keeplimits')
colormap(redblue_or_default())
cb=colorbar; ylabel(cb,'Sv')
xlabel('Date'); ylabel('Pressure (dbar)')
title('MOC streamfunction \Psi(z,t) at 34.5^oS from IES (45-day lowpass)','fontsize',14)
clim_val=prctile(abs(Psi_low(:)),99);
caxis([-clim_val clim_val])
print -dpng -r150 moc_streamfunction_hovmoller

figure
plot(dt,MOC_raw,dt,MOC_low,'r','linewidth',1.5)
legend('raw daily','45-day lowpass')
title('MOC index (depth of max \Psi) at 34.5^oS from IES','fontsize',14)
ylabel('MOC (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_index_IES

% --------------------------------------------------------
% Time-mean Psi(z) vs. depth -- the "vertically integrated transport in
% the basin as a function of depth" plot requested, analogous to
% Kersale et al. (2021)'s Figure 2. Depth on the y-axis, increasing
% downward (oceanographic convention); shaded band = +/-1 std across
% the (lowpass-filtered) daily record.
Psi_mean=nanmean(Psi_low,2);
Psi_std=nanstd(Psi_low,0,2);

figure('Position',[500 100 650 650])
hold on
fill([Psi_mean-Psi_std;flipud(Psi_mean+Psi_std)],[pre;flipud(pre)],[1 0.8 0.8],'EdgeColor','none')
plot(Psi_mean,pre,'r','linewidth',2)
plot([0 0],[0 max(pre)],'k--')
set(gca,'YDir','reverse')
xlabel('\Psi (Sv)'); ylabel('Pressure (dbar)')
title({'Time-mean MOC streamfunction vs. depth at 34.5^oS','(shading = +/-1 std)'},'fontsize',12)
grid on
print -dpng -r150 moc_streamfunction_mean_profile

save('moc_streamfunction_v1.mat','dt','pre','Psi','Psi_low','MOC_raw','MOC_depth_raw','MOC_low','MOC_depth_low','n_gaps_valid','Psi_mean','Psi_std')

function cmap=redblue_or_default()
% Small helper so this runs even if a redblue colormap isn't on the
% path -- falls back to MATLAB's builtin 'jet' otherwise. Not
% critical, purely cosmetic for the Hovmoller.
if exist('redblue','file')
    cmap=redblue(64);
else
    cmap=jet(64);
end
end
