% Step 2 (volume), take 3: MOC using Kersale et al. (2021)'s ACTUAL
% definition, after confirming (by reading kersale_jgr2021_samba.pdf
% directly, section 2.1.1/2.1.4) that v1's construction was NOT what
% the paper does:
%
%   "Because the present study does not use a residual method, and
%   obtains the reference (barotropic) velocity variability directly
%   from data (bottom pressure gradients), no 'residual' zero net
%   volume transport assumption is made here. In other words, the
%   velocities are not uniformly adjusted to ensure zero net volume
%   transport..."
%
%   "the strength of the MOCup flow is defined as the basin-wide
%   transport integrated from the surface down to the time-mean
%   pressure interface where the zonally integrated meridional flow
%   changes from northward to southward (1,315 dbar here at 34.5S for
%   both [pilot and full array] configurations)."
%
% v1 (moc_streamfunction_v1.m) instead subtracted V00 (the depth-mean
% of V0) before integrating -- a RAPID-26N-style "residual/zero-net"
% construction -- then took a PER-DAY max over depth. Both of those are
% exactly what Kersale's paper explicitly does NOT do. Rebuilt here:
%
%   1) Psi_raw(z,t) = cumsum_z[ V0(z,t)*dx2(z,t) ] dz -- the RAW
%      absolute-velocity cumulative transport, NO depth-mean removal.
%      V0 already contains the real BPR-derived barotropic component
%      AND the real ECCO-based time-mean reference (Step E's "Add the
%      time-mean velocity from OFES for a time mean" step, already
%      baked into Absolute_TPUD*_preTopo) -- exactly the ingredients
%      the paper describes combining, no extra data needed.
%   2) h_star = argmax_z[ mean_t(Psi_raw(z,t)) ] -- a SINGLE, FIXED
%      depth found once from the time-mean profile (where the
%      zonally-integrated time-mean flow crosses from northward to
%      southward), NOT a per-day search.
%   3) MOCup(t) = Psi_raw(h_star, t) for every day (same fixed depth),
%      45-day lowpass applied to this single time series (not to the
%      whole depth-time field, since only one depth level is needed
%      now).
%
% Why this matters: in the paper, the Pilot (2-endpoint) and full-array
% MOCup time-means come out at 17.7 vs 17.3 Sv -- essentially
% identical, well within error bars -- vs. the 5-9 Sv gap v1's
% (incorrect) definition produced between our own two implementations.
% Confirms the user's telescoping-identity intuition was right all
% along; v1 was just measuring something else (a RAPID-style
% "overturning strength" that legitimately excludes the barotropic
% component by design -- a fine quantity in its own right, just not
% what Kersale's Figure 2/5 "MOC" refers to).
%
% Bonus cross-check found while re-reading the paper: Kersale reports
% "a time-mean southward transport of 6.6 +/- 2.7 Sv integrated between
% the surface and [4700dbar]" -- essentially identical to the abandoned
% moc_samba_marion_v1.m's moc_mean=6.285 Sv (deleted at the time because
% its overturning/gyre terms were trivially zero -- see CLAUDE.md's MOC
% section). That "mean" term itself was actually a correct, useful,
% independently-validated quantity all along.
%
% v2 (2026-08-15): loads identically to v1 (same Step E v2 source, same
% vel00/geo/dxA/dx2/V0 construction -- that part was never wrong, only
% the streamfunction/MOC-index definition built on top of it changes.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/')

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v2']);

fprintf('dt range: %s to %s (%d days)\n', datestr(dt(1)), datestr(dt(end)), numel(dt));

AREA_TOPO_1=AREA_TOPO_1(:); AREA_TOPO_2=AREA_TOPO_2(:); AREA_TOPO_3=AREA_TOPO_3(:); AREA_TOPO_4=AREA_TOPO_4(:);
AREA_TOPO_5=AREA_TOPO_5(:); AREA_TOPO_6=AREA_TOPO_6(:); AREA_TOPO_7=AREA_TOPO_7(:); AREA_TOPO_8=AREA_TOPO_8(:);

% --------------------------------------------------------
% Distance between sites (dx, per gap) -- identical to v1
position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);

% --------------------------------------------------------
% vel00/geo/dxA/dx2/V0 -- identical construction to v1 (this part was
% never in question, only what's built from V0 afterward)
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
% Kersale's actual definition: RAW cumulative transport, no depth-mean
% removal, evaluated at a FIXED depth found from the time-mean profile.
transport_profile=V0.*dx2; % raw absolute transport-per-unit-depth, all 8 gaps
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
% Hovmoller of the raw (no depth-mean-removal) time-mean-lowpassed Psi,
% for visual reference -- same lowpass-per-depth-level construction as
% v1, just without the vel_prime step.
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
title('MOC cumulative transport \Psi(z,t) at 34.5^oS (raw, Kersale definition, 45-day lowpass)','fontsize',13)
clim_val=prctile(abs(Psi_low(:)),99);
caxis([-clim_val clim_val])
print -dpng -r150 moc_streamfunction_v2_hovmoller

figure
plot(dt,MOCup_raw,dt,MOCup_low,'r','linewidth',1.5)
legend('raw daily','45-day lowpass')
title('MOCup index (transport at fixed h\_star) at 34.5^oS, Kersale definition','fontsize',14)
ylabel('MOC_{up} (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',14)
print -dpng -r150 moc_streamfunction_v2_index_IES

% --------------------------------------------------------
% Time-mean Psi(z) vs. depth -- Kersale Figure 2a style, raw (no
% depth-mean removal), with h_star marked.
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
title({'Time-mean MOC cumulative transport vs. depth at 34.5^oS','Kersale definition (shading = +/-1 std)'},'fontsize',12)
grid on
print -dpng -r150 moc_streamfunction_v2_mean_profile

save('moc_streamfunction_v2.mat','dt','pre','Psi_raw','Psi_low','h_star','MOCup_raw','MOCup_low','n_gaps_valid','Psi_raw_mean','Psi_std')

function cmap=redblue_or_default()
if exist('redblue','file')
    cmap=redblue(64);
else
    cmap=jet(64);
end
end
