% Cross-spectral coherence between the full-array and Pilot MOCup
% anomalies, RAW daily, 2013-2017 -- follow-up to the Gpan noise test
% above. Goal: distinguish whether the correlation shortfall (raw
% r=0.485 vs. the paper's 0.73) is broadband (consistent with the
% already-confirmed excess instrumental/measurement noise at sites
% A/P1 propagating through Gpan) or concentrated in a specific
% frequency band (consistent with the paper's own explanation that the
% full array's larger variance comes from resolving eddies/Agulhas
% Rings the 2-point Pilot structurally cannot see, i.e. real ocean
% signal, not noise). These have different practical implications: a
% broadband noise floor could in principle be reduced (smoothing/
% de-spiking A and P1's raw tau1000); a specific-band physical
% difference cannot be fixed by better processing, only by adding more
% moorings.
%
% Reuses moc_anomaly_fig5_raw_2013_2017.mat's full_anom/pilot_anom
% directly (no NaNs, perfectly regular 1-day sampling, 1405 days) --
% no pipeline rerun.
%
% v1 (2026-08-18): first cut.

clear;clc;
load('moc_anomaly_fig5_raw_2013_2017.mat','dt_win','full_anom','pilot_anom');

fs=1; % cycles/day
N=numel(dt_win);
winlen=256; noverlap=128;
win=hamming(winlen);

[Cxy,f]=mscohere(full_anom,pilot_anom,win,noverlap,[],fs);
[Pxy,~]=cpsd(full_anom,pilot_anom,win,noverlap,[],fs);
[Pxx,~]=pwelch(full_anom,win,noverlap,[],fs);
[Pyy,~]=pwelch(pilot_anom,win,noverlap,[],fs);

period=1./f; % days

% 95% significance threshold for coherence (Welch-averaged), following
% the standard nd-segment formula.
nseg=floor((N-noverlap)/(winlen-noverlap));
alpha=0.05;
coh_thresh=1-alpha^(1/(nseg-1));
fprintf('Welch segments used: %d, 95%% coherence significance threshold: %.3f\n',nseg,coh_thresh);

% Phase -> time lag (days) at each frequency (skip f=0)
lag_days=nan(size(f));
lag_days(2:end)=-angle(Pxy(2:end))./(2*pi*f(2:end));

% Band-averaged coherence, for a quick summary
bands={[100 inf],[20 100],[2 20]};
band_labels={'long (>100d)','medium (20-100d)','short (2-20d)'};
fprintf('\nBand-averaged coherence:\n');
for b=1:numel(bands)
    idx=period>=bands{b}(1) & period<bands{b}(2);
    fprintf('  %-18s: mean Cxy=%.3f (n=%d freq bins)\n',band_labels{b},mean(Cxy(idx)),sum(idx));
end

figure('Position',[100 100 900 750])

subplot(3,1,1)
semilogx(period(2:end),Cxy(2:end),'linewidth',1.5); hold on
plot(xlim,[coh_thresh coh_thresh],'k--')
set(gca,'XDir','reverse')
ylabel('Coherence (Cxy)')
title('Full array vs. Pilot MOCup anomaly: magnitude-squared coherence, RAW daily, 2013-2017','fontsize',12)
legend('Coherence','95% significance','Location','best')
grid on; ylim([0 1])

subplot(3,1,2)
semilogx(period(2:end),lag_days(2:end),'linewidth',1.5); hold on
plot(xlim,[0 0],'k--')
set(gca,'XDir','reverse')
ylabel('Lag (days, + = pilot leads)')
title('Frequency-dependent phase lag (only meaningful where coherence is significant)','fontsize',11)
grid on

subplot(3,1,3)
loglog(period(2:end),f(2:end).*Pxx(2:end),'linewidth',1.5); hold on
loglog(period(2:end),f(2:end).*Pyy(2:end),'linewidth',1.5)
set(gca,'XDir','reverse')
xlabel('Period (days)')
ylabel('Variance-preserving power (f*S(f))')
legend('Full array','Pilot','Location','best')
title('Power spectra (variance-preserving form) -- where each series'' variance is concentrated','fontsize',11)
grid on

print -dpng -r150 moc_coherence_test_IES
save('moc_coherence_test.mat','f','period','Cxy','lag_days','Pxx','Pyy','coh_thresh','nseg')
