% Restricts the time-mean MOC cumulative-transport profile (and h_star)
% to EXACTLY Kersale et al. (2021)'s own 2013-09-11/2017-07-17 window,
% for a direct, apples-to-apples comparison against their Figure 2a --
% every moc_streamfunction_v*/moc_pilot_v* profile-comparison figure so
% far (including v5, the current best) averaged over the FULL extended
% 2013-2022 record instead, which is a real difference from what the
% paper actually plots (asked directly by the user, confirmed: no
% existing version restricts to the paper's own period).
%
% Reuses the ALREADY-SAVED Psi_raw [depth x time] fields from
% moc_streamfunction_v4.mat (full array) and moc_pilot_v5.mat (Pilot)
% directly -- no need to rerun the full pipeline, since Psi_raw itself
% is a per-day quantity computed identically regardless of which
% sub-period you later average over. Only the time-mean/h_star step
% needs restricting.
%
% v1 (2026-08-17): first cut.

clear;clc;

s_full=load('moc_streamfunction_v4.mat','dt','pre','Psi_raw');
s_pilot=load('moc_pilot_v5.mat','dt','pre','Psi_raw');

assert(isequal(s_full.pre,s_pilot.pre),'Pressure grids differ between full array and Pilot.');
pre=s_full.pre;

w1=datenum(2013,9,11); w2=datenum(2017,7,17);

win_full=s_full.dt>=w1 & s_full.dt<=w2;
win_pilot=s_pilot.dt>=w1 & s_pilot.dt<=w2;
fprintf('Full array: %d days in the Kersale window (%s to %s)\n',sum(win_full),datestr(s_full.dt(find(win_full,1))),datestr(s_full.dt(find(win_full,1,'last'))));
fprintf('Pilot:      %d days in the Kersale window (%s to %s)\n',sum(win_pilot),datestr(s_pilot.dt(find(win_pilot,1))),datestr(s_pilot.dt(find(win_pilot,1,'last'))));

Psi_full_mean=nanmean(s_full.Psi_raw(:,win_full),2);
Psi_full_std=nanstd(s_full.Psi_raw(:,win_full),0,2);
[~,imax_full]=max(Psi_full_mean);
h_star_full=pre(imax_full);

Psi_pilot_mean=nanmean(s_pilot.Psi_raw(:,win_pilot),2);
Psi_pilot_std=nanstd(s_pilot.Psi_raw(:,win_pilot),0,2);
[~,imax_pilot]=max(Psi_pilot_mean);
h_star_pilot=pre(imax_pilot);

fprintf('\n2013-2017 window only:\n');
fprintf('Full array: h_star=%d dbar, peak=%.2f Sv (paper: 1315dbar, 17.3 Sv)\n',h_star_full,Psi_full_mean(imax_full));
fprintf('Pilot:      h_star=%d dbar, peak=%.2f Sv (paper: 1315dbar, 17.7 Sv)\n',h_star_pilot,Psi_pilot_mean(imax_pilot));

figure('Position',[500 100 650 650])
hold on
fill([Psi_full_mean-Psi_full_std;flipud(Psi_full_mean+Psi_full_std)],[pre;flipud(pre)],[0.8 0.85 1],'EdgeColor','none','FaceAlpha',0.6)
fill([Psi_pilot_mean-Psi_pilot_std;flipud(Psi_pilot_mean+Psi_pilot_std)],[pre;flipud(pre)],[1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.6)
p1=plot(Psi_full_mean,pre,'b','linewidth',2);
p2=plot(Psi_pilot_mean,pre,'r','linewidth',2);
plot([0 0],[0 max(pre)],'k--')
set(gca,'YDir','reverse')
xlabel('\Psi (Sv)'); ylabel('Pressure (dbar)')
title({'Time-mean MOC cumulative transport vs. depth at 34.5^oS','2013-09-11 to 2017-07-17 ONLY (Kersale et al. 2021''s own window)'},'fontsize',12)
legend([p1 p2],{'Full array (8 gaps)','Pilot (A-P1)'},'Location','southeast')
grid on
print -dpng -r150 moc_profile_comparison_2013_2017_IES

save('moc_profile_comparison_2013_2017.mat','pre','Psi_full_mean','Psi_full_std','h_star_full','Psi_pilot_mean','Psi_pilot_std','h_star_pilot')
