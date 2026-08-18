% Investigates whether the Pilot's own baroclinic signal -- the
% depth-integrated Gpan_P1-Gpan_A difference feeding RelTPUD_pilot,
% cumsum'd 0-to-h_star into MOCup_pilot -- is intrinsically noisier
% (day-to-day) than each of the full array's 8 individual gap
% differences, which is what the earlier raw-tau1000 noise finding
% (site A tied-highest, P1 2nd-highest of 9 sites) predicts should be
% true, but hadn't been tested on the actual GEM-derived Gpan quantity
% that enters the transport calculation -- only on the underlying raw
% tau1000. Also follows the moc_pilot_component_test.m finding that the
% BPR reference is NOT the noise source (removing it makes correlation
% worse) -- redirecting attention to the baroclinic/Gpan signal itself.
%
% Recomputes Gpan for all 9 sites directly from concat_IESsamba.mat's
% per-site T/S (same sw_gpan(...,...,Pressure) call as moc_pilot_v5.m
% and mov_samba_marion*.m), restricted to the paper's 2013-09-11 to
% 2017-07-17 window. For each of the 8 full-array gaps (A-C, C-D,
% D-P8, P8-P6, P6-P5, P5-P4, P4-P2, P2-P1) and the Pilot's own A-to-P1
% span, builds the depth-integral of (Gpan_east-Gpan_west) from the
% surface down to h_star_pilot=1280dbar (from moc_pilot_v5.mat) -- a
% proxy for RelTPUD's own cumsum contribution to MOCup, without the
% 1/f scaling (a constant, irrelevant to noise/signal ratios) or the
% BPR/AREA_TOPO/Ekman/shelf steps (already tested separately).
%
% v1 (2026-08-18): first cut.

clear;clc;
prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';
addpath([prefix_marion,'functions/seawater/']);

load('concat_IESsamba.mat');

w1=datenum(2013,9,11); w2=datenum(2017,7,17);
win=dt(:)>=w1 & dt(:)<=w2;
dt_win=dt(win);
fprintf('%d days in the Kersale window\n',numel(dt_win));

% Common depth grid: prew is 0:10:5000 (501 levels), pree is
% 0:10:5300 (531 levels) -- same spacing/origin, so pree(1:501) lines
% up exactly with prew.
assert(isequal(prew,pree(1:numel(prew))),'prew is not a prefix of pree as expected.');
pre=prew;

sites_W=struct('name',{'A','C','D'},'tem',{tem_A,tem_C,tem_D},'sal',{sal_A,sal_C,sal_D});
sites_E=struct('name',{'P1','P2','P4','P5','P6','P8'},'tem',{tem_P1,tem_P2,tem_P4,tem_P5,tem_P6,tem_P8},'sal',{sal_P1,sal_P2,sal_P4,sal_P5,sal_P6,sal_P8});

Gpan=struct();
for s=1:numel(sites_W)
    nm=sites_W(s).name;
    Gpan.(nm)=sw_gpan(sites_W(s).sal(:,win),sites_W(s).tem(:,win),pre(:));
end
for s=1:numel(sites_E)
    nm=sites_E(s).name;
    g=sw_gpan(sites_E(s).sal(1:numel(pre),win),sites_E(s).tem(1:numel(pre),win),pre(:));
    Gpan.(nm)=g;
end
clear tem_A tem_AA tem_B tem_BB tem_C tem_CC tem_D tem_P1 tem_P2 tem_P4 tem_P5 tem_P6 tem_P8
clear sal_A sal_AA sal_B sal_BB sal_C sal_CC sal_D sal_P1 sal_P2 sal_P4 sal_P5 sal_P6 sal_P8

h_star_pilot=1280; % from moc_pilot_v5.mat
iz=find(pre==h_star_pilot);
dz=mean(diff(pre));

gaps={'A','C';'C','D';'D','P8';'P8','P6';'P6','P5';'P5','P4';'P4','P2';'P2','P1'};
gap_labels={'A-C','C-D','D-P8','P8-P6','P6-P5','P5-P4','P4-P2','P2-P1'};

fprintf('\nDepth-integrated (surface to h_star=%d dbar) Gpan(east)-Gpan(west), day-to-day noise:\n',h_star_pilot);
fprintf('%-10s %10s %10s %10s\n','Span','std(diff)','std(sig)','noise/signal (x1000)');

results=nan(numel(gap_labels)+1,3);
for g=1:size(gaps,1)
    w_series=sum(Gpan.(gaps{g,1})(1:iz,:),1)*dz/1e6; % same units as Psi_raw (1e6 = Sv scaling used elsewhere)
    e_series=sum(Gpan.(gaps{g,2})(1:iz,:),1)*dz/1e6;
    sig=e_series-w_series;
    n=nanstd(diff(sig)); s=nanstd(sig);
    results(g,:)=[n s n/s*1000];
    fprintf('%-10s %10.4f %10.4f %10.1f\n',gap_labels{g},n,s,n/s*1000);
end

% Pilot: direct A-to-P1 span
w_series=sum(Gpan.A(1:iz,:),1)*dz/1e6;
e_series=sum(Gpan.P1(1:iz,:),1)*dz/1e6;
sig_pilot=e_series-w_series;
n_pilot=nanstd(diff(sig_pilot)); s_pilot=nanstd(sig_pilot);
results(end,:)=[n_pilot s_pilot n_pilot/s_pilot*1000];
fprintf('%-10s %10.4f %10.4f %10.1f  <-- Pilot (A-to-P1 direct)\n','A-P1',n_pilot,s_pilot,n_pilot/s_pilot*1000);

fprintf('\nMean noise/signal across the 8 gaps: %.1f (x1000)\n',mean(results(1:8,3)));
fprintf('Pilot (A-P1) noise/signal: %.1f (x1000)\n',results(end,3));
fprintf('Pilot rank among the 9 spans (1=noisiest): %d\n',sum(results(:,3)>=results(end,3)));

figure('Position',[100 100 800 450])
bar(results(:,3))
set(gca,'xticklabel',[gap_labels {'A-P1 (Pilot)'}])
ylabel('noise/signal ratio (x1000)')
title({'Day-to-day noise/signal of depth-integrated (0-1280dbar) Gpan difference,', ...
    'each full-array gap vs. the Pilot''s direct A-to-P1 span, 2013-2017'},'fontsize',11)
grid on
print -dpng -r150 moc_gpan_noise_test_IES

save('moc_gpan_noise_test.mat','gap_labels','results','dt_win','sig_pilot')
