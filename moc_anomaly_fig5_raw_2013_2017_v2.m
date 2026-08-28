% Reconstructs Kersale et al. (2021) Figure 5a/5b using the RAW
% (unfiltered) daily MOCup series, same construction as
% moc_anomaly_fig5_raw_2013_2017.m -- but now loads moc_pilot_v6.mat
% instead of moc_pilot_v5.mat, i.e. the Pilot AFTER fixing the two
% bugs found by comparing against Kersale's own original 2-CPIES
% script (BPR reference sign + missing ECCO time-mean absolute
% reference -- see CLAUDE.md's "MAJOR FINDING" section). v1's r=0.485
% should now come out close to the paper's r=0.73 (moc_pilot_v6.m
% itself already reported r=0.726 for this same comparison; this
% script exists to regenerate the actual Figure-5-style plot with the
% corrected Pilot, not just reprint the number).
%
% v2 (2026-08-28): moc_pilot_v5.mat -> moc_pilot_v6.mat, otherwise
% identical to v1.

clear;clc;

s_full=load('moc_streamfunction_v4.mat','dt','MOCup_raw');
s_pilot=load('moc_pilot_v6.mat','dt','MOCup_pilot_raw');

w1=datenum(2013,9,11); w2=datenum(2017,7,17);

[dt_common,ia,ib]=intersect(s_full.dt,s_pilot.dt);
win=dt_common>=w1 & dt_common<=w2;
dt_win=dt_common(win);

full_win=s_full.MOCup_raw(ia(win));
pilot_win=s_pilot.MOCup_pilot_raw(ib(win));

full_anom=full_win-nanmean(full_win);
pilot_anom=pilot_win-nanmean(pilot_win);

fprintf('%d days in the Kersale window\n',numel(dt_win));
fprintf('Full array: mean=%.2f Sv (removed), anomaly std=%.2f Sv (paper daily sigma: 15.4 Sv)\n',nanmean(full_win),nanstd(full_anom));
fprintf('Pilot (v6, fixed): mean=%.2f Sv (removed), anomaly std=%.2f Sv (paper daily sigma: ~9 Sv)\n',nanmean(pilot_win),nanstd(pilot_anom));

cc=corrcoef(full_anom,pilot_anom,'rows','complete');
fprintf('corr(full anomaly, pilot anomaly), RAW daily = %.3f (paper: r=0.73; pre-fix v1 was r=0.485)\n',cc(1,2));

figure('Position',[100 100 950 500])
plot(dt_win,full_anom,'linewidth',0.8); hold on
plot(dt_win,pilot_anom,'r','linewidth',0.8)
plot(xlim,[0 0],'k--')
legend('Full array (8 gaps)','Pilot v6 (A-P1, fixed)')
title({'MOC_{up} temporal anomaly at 34.5^oS, 2013-2017 (Kersale et al. 2021''s window)', ...
    sprintf('RAW DAILY, no lowpass filter, Pilot v6 (bug-fixed), r=%.2f (paper: r=0.73)',cc(1,2))},'fontsize',12)
ylabel('MOC_{up} anomaly (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',13)
print -dpng -r150 moc_anomaly_fig5_raw_2013_2017_v2_IES

save('moc_anomaly_fig5_raw_2013_2017_v2.mat','dt_win','full_win','pilot_win','full_anom','pilot_anom')
