% Tests whether the weak-correlation windows found in
% moc_lag_analysis_2013_2017.m coincide with elevated raw tau1000 noise
% at sites A and/or P1 (the Pilot's two endpoints, already found to be
% relatively noisy sites -- see the raw-tau1000 investigation above).
% Same 6-month sliding window as the lag analysis; for each window,
% computes the raw (unfiltered) full-array-vs-Pilot correlation
% alongside each site's local noise level (std of day-to-day
% tau1000 difference within that window).
%
% v1 (2026-08-17): first cut.

clear;clc;
addpath('/Users/olga/matlab/stat');
prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';

% MOC series
s_full=load('/Users/olga/research/claude/MOV_Project/moc_streamfunction_v4.mat','dt','MOCup_raw');
s_pilot=load('/Users/olga/research/claude/MOV_Project/moc_pilot_v5.mat','dt','MOCup_pilot_raw');
[dt_common,ia,ib]=intersect(s_full.dt,s_pilot.dt);
w1=datenum(2013,9,11); w2=datenum(2017,7,17);
win=dt_common>=w1&dt_common<=w2;
dt_win=dt_common(win);
full_raw=s_full.MOCup_raw(ia(win));
pilot_raw=s_pilot.MOCup_pilot_raw(ib(win));

% raw tau1000 for A and P1
load([prefix_marion,'MOV/samba_w.mat'],'dt','tau1000_A');
dt_A=dt(:); tauA=tau1000_A(:); clear dt tau1000_A
load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'dt','tau1000_P1');
dt_P1=dt(:); tauP1=tau1000_P1(:); clear dt tau1000_P1

win_days=180; step=30;
n=numel(dt_win);

fprintf('%-10s %8s %10s %10s %10s\n','Mid','corr','noiseA','noiseP1','noiseA+P1');
results=[];
for start=1:step:(n-win_days)
    idx=start:(start+win_days-1);
    mid=dt_win(idx(round(win_days/2)));
    t1=dt_win(idx(1)); t2=dt_win(idx(end));

    a=full_raw(idx); b=pilot_raw(idx);
    cc=corrcoef(a,b,'rows','complete');
    r=cc(1,2);

    iA=dt_A>=t1&dt_A<=t2;
    iP1=dt_P1>=t1&dt_P1<=t2;
    nA=nanstd(diff(tauA(iA)));
    nP1=nanstd(diff(tauP1(iP1)));

    fprintf('%-10s %8.3f %10.5f %10.5f %10.5f\n',datestr(mid,'mmm-yy'),r,nA,nP1,nA+nP1);
    results=[results; mid r nA nP1];
end

fprintf('\nCorrelation between window-correlation and noiseA: %.3f\n',corr(results(:,2),results(:,3)));
fprintf('Correlation between window-correlation and noiseP1: %.3f\n',corr(results(:,2),results(:,4)));
fprintf('Correlation between window-correlation and (noiseA+noiseP1): %.3f\n',corr(results(:,2),results(:,3)+results(:,4)));

figure('Position',[100 100 900 500])
subplot(2,1,1)
plot(results(:,1),results(:,2),'o-','linewidth',1.5); hold on
plot(xlim,[0 0],'k--')
datetick('x',12,'keeplimits')
ylabel('Raw corr (full,pilot)')
title('Window correlation (top) vs. tau1000_A/P1 local noise (bottom)','fontsize',12)
grid on

subplot(2,1,2)
plot(results(:,1),results(:,3),'o-','linewidth',1.5); hold on
plot(results(:,1),results(:,4),'s-','linewidth',1.5)
datetick('x',12,'keeplimits')
ylabel('std(diff(tau1000))')
legend('Site A','Site P1')
grid on
print -dpng -r150 moc_weak_windows_vs_noise_IES
exit
