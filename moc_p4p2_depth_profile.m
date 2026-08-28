% Depth profile of gap P4-P2's own transport-per-unit-depth bias
% during the 2017-2020 outage, follow-up to the component-decomposition
% test above (which isolated the bias to the relative/baroclinic
% component) and the CTD ground-truth checks (which found a ~7ms
% PIES-vs-CTD tau1000 offset at P4/P2/P1). Compares the time-mean
% Total_TPUD7 profile (normal days vs. outage days, using only days
% gap 7 actually has data) across the FULL depth range, not just the
% 0-h_star integral already computed -- to see whether the bias is
% concentrated at a specific depth (favoring a localized real signal)
% or has the shape of the gap's own dominant GEM mode, just scaled
% (favoring a tau1000-driven, GEM-mediated effect).
%
% v1 (2026-08-28): first cut.

clear;clc;
prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
S=load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3'], ...
    'Total_TPUD7','Pressure','dt');
pre=S.Pressure(:); dt=S.dt(:);
outage1=datenum(2017,8,1); outage2=datenum(2020,1,1);
normal=~(dt>=outage1 & dt<=outage2);
outage=~normal;

prof_n=nanmean(S.Total_TPUD7(:,normal),2);
prof_o=nanmean(S.Total_TPUD7(:,outage),2);
diffp=prof_o-prof_n;

fprintf('%-8s %10s %10s %10s\n','Pressure','normal','outage','diff');
for pp=[0 50 100 200 300 500 700 900 1190 1500 2000 2500 3000 3500 4000 4500 5000]
    iz=find(pre==pp);
    if isempty(iz), continue; end
    fprintf('%-8d %10.1f %10.1f %10.1f\n',pp,prof_n(iz),prof_o(iz),diffp(iz));
end

figure('Position',[100 100 650 700])
plot(prof_n,pre,'b','linewidth',1.8); hold on
plot(prof_o,pre,'r','linewidth',1.8);
plot([0 0],[0 max(pre)],'k--')
set(gca,'YDir','reverse')
xlabel('Total\_TPUD7 (P4-P2), m^2/s'); ylabel('Pressure (dbar)')
legend('Normal-period mean','Outage-period mean','Location','southeast')
title({'Gap P4-P2: time-mean transport-per-unit-depth profile,','normal vs. 2017-2020 outage period'},'fontsize',12)
grid on
print -dpng -r150 moc_p4p2_depth_profile_IES
