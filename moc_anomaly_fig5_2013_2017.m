% Reconstructs Kersale et al. (2021) Figure 5a/5b: temporal ANOMALIES
% (each series minus its OWN time-mean) of MOCup for the Pilot and the
% full array, restricted to exactly the paper's 2013-09-11/2017-07-17
% window, using the mean-preserving lowpass_filter_fixed.m (not the
% buggy shared ~/matlab/stat/lowpass_filter.m -- see CLAUDE.md's "MAJOR
% FINDING" section).
%
% Why anomalies specifically: subtracting each series' own mean
% sidesteps the whole absolute-level discrepancy this session has been
% chasing (calibration-window, GEM-vintage/Brazil Current, lowpass
% bias) and isolates whether the two methods' DAY-TO-DAY VARIABILITY is
% consistent -- which is what the paper's own r=0.73 correlation
% actually measures. Reuses the already-saved MOCup_raw series from
% moc_streamfunction_v4.mat/moc_pilot_v5.mat -- no pipeline rerun.
%
% NOT reconstructed: panel 5c (MOCab, the abyssal cell -- a separate
% depth-integrated quantity, 3150-4300dbar, not yet computed in this
% project) and the paper's gray "estimated daily accuracy" shading (no
% equivalent error-budget methodology built here yet). This script
% covers 5a/5b only (MOCup, both configurations), overlaid on one axis
% rather than as separate subplots, for direct comparison.
%
% v1 (2026-08-17): first cut.

clear;clc;
addpath('/Users/olga/matlab/stat'); % for blackman() via lowpass_filter_fixed

s_full=load('moc_streamfunction_v4.mat','dt','MOCup_raw');
s_pilot=load('moc_pilot_v5.mat','dt','MOCup_pilot_raw');

w1=datenum(2013,9,11); w2=datenum(2017,7,17);

full_fixed=lowpass_filter_fixed(s_full.MOCup_raw,45);
pilot_fixed=lowpass_filter_fixed(s_pilot.MOCup_pilot_raw,45);

[dt_common,ia,ib]=intersect(s_full.dt,s_pilot.dt);
win=dt_common>=w1 & dt_common<=w2;
dt_win=dt_common(win);

full_win=full_fixed(ia(win));
pilot_win=pilot_fixed(ib(win));

full_anom=full_win-nanmean(full_win);
pilot_anom=pilot_win-nanmean(pilot_win);

fprintf('%d days in the Kersale window\n',numel(dt_win));
fprintf('Full array: mean=%.2f Sv (removed), anomaly std=%.2f Sv\n',nanmean(full_win),nanstd(full_anom));
fprintf('Pilot:      mean=%.2f Sv (removed), anomaly std=%.2f Sv\n',nanmean(pilot_win),nanstd(pilot_anom));

cc=corrcoef(full_anom,pilot_anom,'rows','complete');
fprintf('corr(full anomaly, pilot anomaly) = %.3f (paper: r=0.73)\n',cc(1,2));

figure('Position',[100 100 950 500])
plot(dt_win,full_anom,'linewidth',1.3); hold on
plot(dt_win,pilot_anom,'r','linewidth',1.3)
plot(xlim,[0 0],'k--')
legend('Full array (8 gaps)','Pilot (A-P1)')
title({'MOC_{up} temporal anomaly at 34.5^oS, 2013-2017 (Kersale et al. 2021''s window)', ...
    sprintf('mean-preserving 45-day lowpass, r=%.2f (paper: r=0.73)',cc(1,2))},'fontsize',12)
ylabel('MOC_{up} anomaly (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',13)
print -dpng -r150 moc_anomaly_fig5_2013_2017_IES

save('moc_anomaly_fig5_2013_2017.mat','dt_win','full_win','pilot_win','full_anom','pilot_anom')
