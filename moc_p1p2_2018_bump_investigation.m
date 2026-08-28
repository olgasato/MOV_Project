% Investigates the time-varying component (on top of the now-confirmed
% PERSISTENT PIES-vs-CTD tau1000 offset, see CLAUDE.md's "fourth
% cruise" section) that appeared slightly elevated at P1/P2 specifically
% around October 2018, based on only 3-4 sparse CTD snapshots. Uses the
% FULL, continuous tau1000_P1/tau1000_P2 daily record instead (no CTD
% needed) to pin down the timing/shape of any 2018-specific anomaly
% much more precisely than 3-4 cruise dates can.
%
% Method: fit a smooth seasonal climatology (annual + semi-annual
% harmonics, least squares) to each site's tau1000 using the WHOLE
% record, then look at the residual (raw minus seasonal fit) over time
% -- if there's a real, temporary 2018-specific excursion on top of the
% persistent offset, it should show up as a distinct bump in this
% residual, not explained by the regular seasonal cycle.
%
% v1 (2026-08-28): first cut.

clear;clc;
prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';

load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'tau1000_P1','tau1000_P2','dt');
dt=dt(:);
tauP1=tau1000_P1(:); tauP2=tau1000_P2(:);

t_frac=(dt-dt(1))/365.25;
doy_frac=mod(t_frac,1);
Xdes=[ones(size(t_frac)) cos(2*pi*doy_frac) sin(2*pi*doy_frac) cos(4*pi*doy_frac) sin(4*pi*doy_frac)];

function [resid,fitted]=deseasonalize(y,X)
    ok=~isnan(y);
    coef=X(ok,:)\y(ok);
    fitted=X*coef;
    resid=y-fitted;
end

[residP1,fitP1]=deseasonalize(tauP1,Xdes);
[residP2,fitP2]=deseasonalize(tauP2,Xdes);

% 30-day running mean of the residual, to see the slow-varying part
smoothP1=movmean(residP1,30,'omitnan');
smoothP2=movmean(residP2,30,'omitnan');

fprintf('Residual (raw - seasonal fit) stats, whole record:\n');
fprintf('  P1: mean=%.5f, std=%.5f\n',nanmean(residP1),nanstd(residP1));
fprintf('  P2: mean=%.5f, std=%.5f\n',nanmean(residP2),nanstd(residP2));

fprintf('\n30-day smoothed residual, monthly snapshots, 2016-2020:\n');
fprintf('%-10s %10s %10s\n','Month','P1','P2');
for y=2016:2020
    for m=1:12
        d=datenum(y,m,15);
        idx=find(abs(dt-d)<15,1);
        if isempty(idx), continue; end
        fprintf('%-10s %10.5f %10.5f\n',datestr(d,'mmm-yyyy'),smoothP1(idx),smoothP2(idx));
    end
end

figure('Position',[100 100 950 600])
subplot(2,1,1)
plot(dt,residP1,'color',[0.7 0.7 1]); hold on
plot(dt,smoothP1,'b','linewidth',1.8)
plot(xlim,[0 0],'k--')
datetick('x',12,'keeplimits')
ylabel('P1 residual (s)')
title('tau1000_{P1}: deseasonalized residual (thin) and 30-day smoothed (thick)','fontsize',12)
grid on

subplot(2,1,2)
plot(dt,residP2,'color',[1 0.7 0.7]); hold on
plot(dt,smoothP2,'r','linewidth',1.8)
plot(xlim,[0 0],'k--')
datetick('x',12,'keeplimits')
ylabel('P2 residual (s)')
title('tau1000_{P2}: deseasonalized residual (thin) and 30-day smoothed (thick)','fontsize',12)
grid on
print -dpng -r150 moc_p1p2_2018_bump_investigation_IES

% Zoom on 2016-2020
figure('Position',[100 100 950 500])
plot(dt,smoothP1,'b','linewidth',1.8); hold on
plot(dt,smoothP2,'r','linewidth',1.8)
plot(xlim,[0 0],'k--')
xlim([datenum(2016,1,1) datenum(2021,1,1)])
datetick('x',12,'keeplimits')
ylabel('30-day smoothed residual (s)')
legend('P1','P2','Location','best')
title('Zoom 2016-2020: does a distinct 2018 bump stand out from the persistent baseline?','fontsize',12)
grid on
print -dpng -r150 moc_p1p2_2018_bump_zoom_IES

save('moc_p1p2_2018_bump_investigation.mat','dt','residP1','residP2','smoothP1','smoothP2')
