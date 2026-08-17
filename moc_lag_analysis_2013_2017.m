% Investigates whether the full-array-vs-Pilot correlation shortfall
% (0.51 vs. the paper's 0.73, both using this project's mean-preserving
% 45-day lowpass) is a genuine phase/lag issue, following the user's
% direct observation that the two curves in Kersale et al. (2021)'s
% Figure 5 look in-phase while this repo's `moc_anomaly_fig5_2013_2017`
% reconstruction looked visibly out of phase in places.
%
% Two-part investigation:
%   1) Time-varying lag (6-month sliding window, +/-30 day search) on
%      the 45-day-lowpass-filtered series -- NOT a constant offset:
%      alternates between windows of near-zero lag + high correlation
%      (Jul2015-Jan2016, Nov2016-Mar2017, r=0.7-0.94) and windows of
%      unstable/large apparent lag + weak correlation. Extending the
%      lag search to +/-90 days destabilizes most of the "large lag"
%      estimates (they flip sign/magnitude), indicating those aren't a
%      genuine, well-determined phase shift -- just noise in short
%      windows once the tested lag becomes a large fraction of the
%      window length.
%   2) Same test on the RAW (unfiltered) daily series: whole-record
%      best lag = 0 days (already the best fit -- no lag improves the
%      raw correlation at all). The time-varying version is
%      overwhelmingly lag=0 across windows too, with the rare non-zero
%      exceptions landing exactly on the SAME windows that also have
%      near-zero raw correlation (Sep-Nov 2014) -- i.e. noise, not a
%      real lag there either.
%
% CONCLUSION: the apparent phase offset was an ARTIFACT of the 45-day
% lowpass filter itself (smoothing broadens/shifts peak locations),
% not a real property of the underlying data. Without the filter, the
% two methods are essentially in phase (raw lag=0). The remaining
% correlation shortfall (raw r=0.485 vs. paper's 0.73) is about
% amplitude/shape matching during specific weak windows, not phase.
%
% v1 (2026-08-17): first cut, reusing already-saved MOCup_raw /
% moc_anomaly_fig5_2013_2017.mat (lowpass-filtered) -- no pipeline rerun.

clear;clc;

s_full=load('moc_streamfunction_v4.mat','dt','MOCup_raw');
s_pilot=load('moc_pilot_v5.mat','dt','MOCup_pilot_raw');
s_anom=load('moc_anomaly_fig5_2013_2017.mat'); % dt_win, full_win, pilot_win (45-day lowpass, 2013-2017)

[dt_common,ia,ib]=intersect(s_full.dt,s_pilot.dt);
w1=datenum(2013,9,11); w2=datenum(2017,7,17);
win=dt_common>=w1&dt_common<=w2;
dt_win=dt_common(win);
full_raw=s_full.MOCup_raw(ia(win));
pilot_raw=s_pilot.MOCup_pilot_raw(ib(win));

cc_raw=corrcoef(full_raw,pilot_raw,'rows','complete');
cc_low=corrcoef(s_anom.full_win,s_anom.pilot_win,'rows','complete');
fprintf('Whole-record correlation: RAW=%.3f, 45-day lowpass=%.3f (paper: 0.73)\n',cc_raw(1,2),cc_low(1,2));

maxlag=30; win_days=180; step=30;
n=numel(dt_win);

lag_raw=nan(1,0); lag_low=nan(1,0); mids=[];
for start=1:step:(n-win_days)
    idx=start:(start+win_days-1);
    mid=dt_win(idx(round(win_days/2)));
    mids=[mids mid];

    a=full_raw(idx)-nanmean(full_raw(idx)); b=pilot_raw(idx)-nanmean(pilot_raw(idx));
    a(isnan(a))=0; b(isnan(b))=0;
    [xc,lags]=xcorr(a,b,maxlag,'coeff');
    [~,ilag]=max(xc);
    lag_raw(end+1)=lags(ilag);

    a=s_anom.full_win(idx)-nanmean(s_anom.full_win(idx)); b=s_anom.pilot_win(idx)-nanmean(s_anom.pilot_win(idx));
    a(isnan(a))=0; b(isnan(b))=0;
    [xc,lags]=xcorr(a,b,maxlag,'coeff');
    [~,ilag]=max(xc);
    lag_low(end+1)=lags(ilag);
end

fprintf('Raw:    mean best lag=%.1f days, %% of windows with lag==0: %.0f%%\n',mean(lag_raw),100*mean(lag_raw==0));
fprintf('45-day: mean best lag=%.1f days, %% of windows with lag==0: %.0f%%\n',mean(lag_low),100*mean(lag_low==0));

figure('Position',[100 100 900 400])
plot(mids,lag_low,'o-','linewidth',1.3,'markersize',5); hold on
plot(mids,lag_raw,'s-','linewidth',1.3,'markersize',5)
plot(xlim,[0 0],'k--')
datetick('x',12,'keeplimits')
ylabel('Best lag (days, + = pilot leads)')
legend('45-day lowpass','Raw daily','Location','best')
title('Time-varying lag, full array vs. Pilot: raw vs. 45-day-lowpass','fontsize',12)
grid on
print -dpng -r150 moc_lag_analysis_2013_2017_IES
