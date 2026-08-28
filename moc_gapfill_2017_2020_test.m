% Quantifies, via direct gap-filling, how much of the full array's mean
% MOCup shortfall during the 2017-08/2020-01 coverage outage (gaps
% 3-6 -- D-P8, P8-P6, P6-P5, P5-P4 -- missing or badly degraded, see
% "Root cause of the full-array-vs-Pilot MOC intensity gap" and
% "Component breakdown of the 2017-2019 dip" above) is mechanically
% explained by those gaps' missing contribution to the width-integral,
% versus something else.
%
% Follow-up to the user's specific question after the Pilot v6 fix:
% now that the Pilot (ground truth-ish, unaffected by this outage) and
% the full array are both known-correct, why does the MEAN specifically
% diverge so much during 2017-2020 while it matches well before/after?
%
% Method: moc_streamfunction_v4.m's own derivation shows
% transport_profile = nansum_i(Total_TPUD_i) + shelf_profile exactly
% (confirmed algebraically from its vel00/dxA/dx2/V0 construction) --
% i.e. each day's basin-wide transport-per-unit-depth is a literal SUM
% across the 8 gaps' own Total_TPUD_i, where a gap with no data that
% day contributes exactly ZERO (via nansum), not "unknown". This
% script builds a GAP-FILLED alternative: whenever a gap's Total_TPUD_i
% is NaN on a given day, substitute that gap's own long-term
% climatological (time-mean, computed from "normal" -- outside the
% outage -- days only) profile instead of letting it silently
% contribute zero. If the resulting reconstructed MOCup mean during
% 2017-2020 comes close to the Pilot's mean over the same days, that
% directly confirms the missing MEAN (not just the missing variance)
% is explained by this literal missing-contribution mechanism.
%
% v1 (2026-08-28): first cut.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/')

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3'], ...
    'Total_TPUD1','Total_TPUD2','Total_TPUD3','Total_TPUD4','Total_TPUD5','Total_TPUD6','Total_TPUD7','Total_TPUD8', ...
    'Pressure','dt');
pre=Pressure(:); clear Pressure
dz=mean(diff(pre));

outage1=datenum(2017,8,1); outage2=datenum(2020,1,1);
normal=~(dt(:)>=outage1 & dt(:)<=outage2);
outage=~normal;
fprintf('%d "normal" days, %d "outage" days (2017-08-01 to 2020-01-01)\n',sum(normal),sum(outage));

names={'A-C','C-D','D-P8','P8-P6','P6-P5','P5-P4','P4-P2','P2-P1'};
Total_all=cat(3,Total_TPUD1,Total_TPUD2,Total_TPUD3,Total_TPUD4,Total_TPUD5,Total_TPUD6,Total_TPUD7,Total_TPUD8);
clear Total_TPUD1 Total_TPUD2 Total_TPUD3 Total_TPUD4 Total_TPUD5 Total_TPUD6 Total_TPUD7 Total_TPUD8

% Each gap's own climatological (time-mean, normal days only) profile
clim=nanmean(Total_all(:,normal,:),2); % [depth x 1 x 8]
fprintf('\nFraction of outage days each gap is missing, and its normal-period mean transport at h_star-ish depth (500dbar):\n');
iz500=find(pre==500);
for i=1:8
    miss=sum(isnan(Total_all(iz500,outage,i)))/sum(outage);
    fprintf('  Gap %-6s: %5.1f%% missing during outage, normal-period mean(500dbar)=%.1f m^2/s\n',names{i},100*miss,clim(iz500,1,i));
end

% Gap-filled version: substitute the gap's own climatology wherever
% NaN, ONLY during the outage window (leave normal days untouched --
% irrelevant there anyway since coverage is ~complete).
Total_filled=Total_all;
for i=1:8
    sl=Total_filled(:,:,i);
    badmask=isnan(sl) & repmat(outage(:)',size(sl,1),1);
    climrep=repmat(clim(:,1,i),1,size(sl,2));
    sl(badmask)=climrep(badmask);
    Total_filled(:,:,i)=sl;
end

transport_actual=nansum(Total_all,3);
% nansum (not plain sum): Total_filled has no NaNs left WITHIN the
% outage window by construction, but stray missing days can still
% exist elsewhere (e.g. the smaller, separately-documented 2020-2022
% reduced-coverage patch) -- plain sum would propagate those to a
% whole-day NaN even where 7 of 8 gaps are fine, inconsistent with how
% transport_actual itself is built. nansum treats both the same way.
transport_filled=nansum(Total_filled,3);

load('/Users/olga/research/sambar/renellys_sent/Marions_code/ecco/Wrk/ECCO_TPUD_ShelfBits.mat','West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(pre,pres_shelf_grid);
shelf_profile=West_TPUD(ib)+East_TPUD_new(ib);
clear West_TPUD East_TPUD_new pres_shelf_grid ia ib

transport_actual=transport_actual+shelf_profile;
transport_filled=transport_filled+shelf_profile;

Psi_actual=cumsum(transport_actual,1,'omitnan')*dz/1e6;
Psi_filled=cumsum(transport_filled,1,'omitnan')*dz/1e6;

% Use the SAME h_star as the real, unfilled full array (moc_streamfunction_v4.mat)
s_full=load('moc_streamfunction_v4.mat','h_star','dt','MOCup_raw');
h_star=s_full.h_star;
iz=find(pre==h_star);
fprintf('\nUsing h_star=%d dbar (from moc_streamfunction_v4.mat)\n',h_star);

MOCup_actual=Psi_actual(iz,:);
MOCup_filled=Psi_filled(iz,:);

% Sanity check: actual reconstruction here should match v4's own saved
% MOCup_raw closely (same formula, just recombined from Total_TPUD_i
% directly instead of via the V0/dx2 route -- confirms the algebraic
% equivalence claimed in the header).
[~,ia,ib]=intersect(dt,s_full.dt);
d=MOCup_actual(ia)-s_full.MOCup_raw(ib);
fprintf('Sanity check (actual reconstruction vs. v4''s saved MOCup_raw): max|diff|=%.4f Sv\n',max(abs(d)));

s_pilot=load('moc_pilot_v6.mat','dt','MOCup_pilot_raw');
[dtc,ip,ipp]=intersect(dt,s_pilot.dt);

fprintf('\n%-30s %10s %10s\n','Period','mean (Sv)','std (Sv)');
for w=1:2
    if w==1
        win=outage(ip); lbl='Outage (2017-08/2020-01)';
    else
        win=normal(ip); lbl='Normal (rest of record)';
    end
    act=MOCup_actual(ip(win)); fil=MOCup_filled(ip(win)); pil=s_pilot.MOCup_pilot_raw(ipp(win));
    fprintf('%s:\n',lbl);
    fprintf('  %-28s %10.2f %10.2f\n','Full array, ACTUAL (as-is)',nanmean(act),nanstd(act));
    fprintf('  %-28s %10.2f %10.2f\n','Full array, GAP-FILLED',nanmean(fil),nanstd(fil));
    fprintf('  %-28s %10.2f %10.2f\n','Pilot v6',nanmean(pil),nanstd(pil));
end

win=outage(ip);
gap_actual_vs_pilot=nanmean(s_pilot.MOCup_pilot_raw(ipp(win)))-nanmean(MOCup_actual(ip(win)));
gap_filled_vs_pilot=nanmean(s_pilot.MOCup_pilot_raw(ipp(win)))-nanmean(MOCup_filled(ip(win)));
pct_closed=100*(1-gap_filled_vs_pilot/gap_actual_vs_pilot);
fprintf('\nDuring the outage: Pilot-minus-ACTUAL gap = %.2f Sv; Pilot-minus-FILLED gap = %.2f Sv\n',gap_actual_vs_pilot,gap_filled_vs_pilot);
fprintf('Gap-filling closes %.0f%% of the mean discrepancy.\n',pct_closed);

figure('Position',[100 100 950 500])
plot(dt,MOCup_actual,'b','linewidth',0.8); hold on
plot(dt,MOCup_filled,'color',[0 0.6 0],'linewidth',1.1);
plot(dtc,s_pilot.MOCup_pilot_raw(ipp),'r','linewidth',0.8);
plot([outage1 outage1],ylim,'k--'); plot([outage2 outage2],ylim,'k--')
legend('Full array, ACTUAL','Full array, GAP-FILLED (missing gaps -> own climatology)','Pilot v6','Location','best')
title('Does gap-filling the missing 2017-2020 sites recover the Pilot''s mean level?','fontsize',12)
ylabel('MOC_{up} (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',12)
print -dpng -r150 moc_gapfill_2017_2020_test_IES

save('moc_gapfill_2017_2020_test.mat','dt','MOCup_actual','MOCup_filled','h_star','outage','normal')
