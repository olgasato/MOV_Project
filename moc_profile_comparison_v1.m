% Kersale et al. (2021) Figure 2-style plot: time-mean MOC streamfunction
% vs. depth, full array (8 gaps, moc_streamfunction_v1.m) overlaid with
% the Pilot method (A-P1 only, moc_pilot_v1.m), each shaded +/-1 std.
% Reuses both scripts' saved .mat output directly -- no recomputation.
%
% v1 (2026-08-15): first cut.

clear;clc;

s_full=load('moc_streamfunction_v1.mat','pre','Psi_mean','Psi_std');
s_pilot=load('moc_pilot_v1.mat','pre','Psi_pilot_low');

assert(isequal(s_full.pre,s_pilot.pre),'Full-array and Pilot pressure grids differ -- cannot overlay directly.');
pre=s_full.pre;

Psi_pilot_mean=nanmean(s_pilot.Psi_pilot_low,2);
Psi_pilot_std=nanstd(s_pilot.Psi_pilot_low,0,2);

figure('Position',[500 100 650 650])
hold on
fill([s_full.Psi_mean-s_full.Psi_std;flipud(s_full.Psi_mean+s_full.Psi_std)],[pre;flipud(pre)],[0.8 0.85 1],'EdgeColor','none','FaceAlpha',0.6)
fill([Psi_pilot_mean-Psi_pilot_std;flipud(Psi_pilot_mean+Psi_pilot_std)],[pre;flipud(pre)],[1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.6)
p1=plot(s_full.Psi_mean,pre,'b','linewidth',2);
p2=plot(Psi_pilot_mean,pre,'r','linewidth',2);
plot([0 0],[0 max(pre)],'k--')
set(gca,'YDir','reverse')
xlabel('\Psi (Sv)'); ylabel('Pressure (dbar)')
title({'Time-mean MOC streamfunction vs. depth at 34.5^oS','(shading = +/-1 std)'},'fontsize',12)
legend([p1 p2],{'Full array (8 gaps)','Pilot (A-P1)'},'Location','southeast')
grid on
print -dpng -r150 moc_profile_comparison_IES

fprintf('Full array:  peak %.1f Sv at %d dbar\n', max(s_full.Psi_mean), pre(find(s_full.Psi_mean==max(s_full.Psi_mean),1)));
fprintf('Pilot:       peak %.1f Sv at %d dbar\n', max(Psi_pilot_mean), pre(find(Psi_pilot_mean==max(Psi_pilot_mean),1)));
