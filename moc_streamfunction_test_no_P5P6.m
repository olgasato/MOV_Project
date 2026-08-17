% TEST script (not a new "vN" of the main line -- delete if inconclusive,
% keep+rename if it reveals something worth preserving): does the
% full-array MOCup still deviate from the Pilot even when restricted,
% for the WHOLE 2013-2022 record, to only the gaps that have at least
% some data during the 2017-2019 coverage outage?
%
% User's hypothesis: treating the Pilot (moc_pilot_v5.m, unaffected by
% intermediate-site outages since it only uses A/P1) as approximate
% "ground truth", if the full array is restricted to ONLY the gaps that
% survive the 2017-2019 outage -- applied to the ENTIRE record, not
% just that window -- and it STILL deviates from the Pilot outside that
% window too, that would suggest the deviation isn't purely a "narrow
% sliver of gaps = unrepresentative sample" effect specific to 2017-
% 2019, but something more general about using few/these particular
% gaps.
%
% Checked exactly which SITES have ANY data during 2017-08-01 to
% 2020-01-01 (not just gap-level, which conflates two sites): A 100%,
% C 100%, D 76.6%, P8 34.0%, P4 70.4%, P2 100%, P1 100% -- only P5 and
% P6 are completely (100%) absent throughout the window. Per user's
% choice, restricting to gaps 1,2,3,7,8 (A-C, C-D, D-P8, P4-P2, P2-P1)
% -- i.e. dropping only gaps 4,5,6 (P8-P6, P6-P5, P5-P4, all of which
% require P5 and/or P6) -- applied as a FIXED subset for every single
% day in the 2013-2022 record, not day-varying.
%
% Otherwise identical to moc_streamfunction_v4.m (Step E v3 source,
% Kersale definition, shelf correction) -- only the gap subset changes.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/')

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3']);

fprintf('dt range: %s to %s (%d days)\n', datestr(dt(1)), datestr(dt(end)), numel(dt));

AREA_TOPO_1=AREA_TOPO_1(:); AREA_TOPO_2=AREA_TOPO_2(:); AREA_TOPO_3=AREA_TOPO_3(:); AREA_TOPO_4=AREA_TOPO_4(:);
AREA_TOPO_5=AREA_TOPO_5(:); AREA_TOPO_6=AREA_TOPO_6(:); AREA_TOPO_7=AREA_TOPO_7(:); AREA_TOPO_8=AREA_TOPO_8(:);

% --------------------------------------------------------
% Distance between sites (dx, per gap) -- identical to v4
position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);

% --------------------------------------------------------
% vel00/geo/dxA/dx2/V0 -- identical construction to v4
[l,m]=size(Total_TPUD1);
n=8;
vel00=nan*ones([l m n]);
geo=nan*ones([l m n]);
for ii=1:8
str=['vel00(:,:,',num2str(ii),')=Absolute_TPUD',num2str(ii),'_preTopo/dx(',num2str(ii),');'];
eval(str)
str=['geo(:,:,',num2str(ii),')=Absolute_TPUD',num2str(ii),'_preTopo;'];
eval(str)
end
clear Absolute_TPUD*_preTopo

for ii=1:8
str=['geo(1:7,:,',num2str(ii),')=geo(1:7,:,',num2str(ii),')+Ekman_TPUD',num2str(ii),';'];
eval(str)
str=['vel00(1:7,:,',num2str(ii),')=geo(1:7,:,',num2str(ii),')/dx(',num2str(ii),');'];
eval(str)
end
clear Ekman_TPUD1 Ekman_TPUD2 Ekman_TPUD3 Ekman_TPUD4 Ekman_TPUD5 Ekman_TPUD6 Ekman_TPUD7 Ekman_TPUD8 ii

AREA_TOPO_3d=nan*ones([l 1 n]);
for ii=1:8
str=['AREA_TOPO_3d(:,1,',num2str(ii),')=AREA_TOPO_',num2str(ii),'(:);'];
eval(str)
end
AREA_TOPO_3d=repmat(AREA_TOPO_3d,[1 m 1]);
clear AREA_TOPO_1 AREA_TOPO_2 AREA_TOPO_3 AREA_TOPO_4 AREA_TOPO_5 AREA_TOPO_6 AREA_TOPO_7 AREA_TOPO_8 ii

dx1=repmat(reshape(dx,1,1,[]),l,m,1);
ind=find(isnan(geo) | isnan(AREA_TOPO_3d));
dx1(ind)=nan;
geo(isnan(dx1))=nan;
dxA=dx1.*AREA_TOPO_3d;

% --------------------------------------------------------
% NEW: restrict to gaps 1,2,3,7,8 for EVERY day (fixed subset, not
% day-varying) -- drop 4,5,6 (P8-P6,P6-P5,P5-P4) entirely by forcing
% their dxA to NaN everywhere.
keep_gaps=[1 2 3 7 8];
drop_gaps=setdiff(1:8,keep_gaps);
dxA(:,:,drop_gaps)=NaN;
vel00(:,:,drop_gaps)=NaN;
fprintf('Restricted to gaps %s (dropped %s) for the WHOLE record.\n',mat2str(keep_gaps),mat2str(drop_gaps));

dx2=squeeze(nansum(dxA,3));
n_gaps_valid=reshape(sum(any(~isnan(dxA),1),3),1,m);

V0=nansum(vel00.*dxA,3)./nansum(dxA,3);

pre=Pressure; clear Pressure
dz=mean(diff(pre));

% --------------------------------------------------------
% Shelf transport correction (same as v3/v4)
load('/Users/olga/research/sambar/renellys_sent/Marions_code/ecco/Wrk/ECCO_TPUD_ShelfBits.mat','West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(pre,pres_shelf_grid);
assert(isequal(ia(:),(1:numel(pre))'),'pre is not a subset of the shelf grid as expected');
shelf_profile=West_TPUD(ib)+East_TPUD_new(ib);
clear West_TPUD East_TPUD_new pres_shelf_grid ia ib

% --------------------------------------------------------
% Kersale definition
transport_profile=V0.*dx2+shelf_profile;
Psi_raw=cumsum(transport_profile,1,'omitnan')*dz/1e6;

Psi_raw_mean=nanmean(Psi_raw,2);
[~,imax_mean]=max(Psi_raw_mean);
h_star=pre(imax_mean);
fprintf('h_star: %d dbar\n',h_star);

MOCup_raw=Psi_raw(imax_mean,:);
MOCup_low=lowpass_filter(MOCup_raw,45);
fprintf('LOWPASS (45-day) MOCup at h_star: mean=%.4f Sv, std=%.4f Sv\n', nanmean(MOCup_low), nanstd(MOCup_low));

% --------------------------------------------------------
% Compare against the Pilot (v5, unaffected by this restriction)
s_pilot=load('moc_pilot_v5.mat','dt','MOCup_pilot_low');
[dt_common,ia,ib]=intersect(dt,s_pilot.dt);

figure('Position',[100 100 900 500])
hold on
patch([datenum(2017,8,1) datenum(2020,1,1) datenum(2020,1,1) datenum(2017,8,1)],[-100 -100 100 100],[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off')
p1=plot(dt_common,MOCup_low(ia),'b','linewidth',1.5);
p2=plot(dt_common,s_pilot.MOCup_pilot_low(ib),'r','linewidth',1.5);
ylim([min([MOCup_low(ia) s_pilot.MOCup_pilot_low(ib)])-2 max([MOCup_low(ia) s_pilot.MOCup_pilot_low(ib)])+2])
legend([p1 p2],{'Full array, gaps 1,2,3,7,8 ONLY (fixed, whole record)','Pilot (A-P1)'})
title('Test: full array restricted to 2017-2019-surviving gaps, vs. Pilot','fontsize',13)
ylabel('MOC_{up} (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',13)
print -dpng -r150 moc_test_no_P5P6_vs_pilot

fprintf('\nMeans: Full(restricted)=%.2f Sv, Pilot=%.2f Sv\n',nanmean(MOCup_low(ia)),nanmean(s_pilot.MOCup_pilot_low(ib)));
fprintf('corr(full_restricted,pilot), ALL days = %.3f\n', corr2_nan(MOCup_low(ia),s_pilot.MOCup_pilot_low(ib)));

win=dt_common>=datenum(2017,8,1)&dt_common<=datenum(2020,1,1);
fprintf('\n--- Inside 2017-2019 (n=%d) ---\n',sum(win));
fprintf('Full(restricted)=%.2f Sv, Pilot=%.2f Sv, corr=%.3f\n', ...
    nanmean(MOCup_low(ia(win))),nanmean(s_pilot.MOCup_pilot_low(ib(win))),corr2_nan(MOCup_low(ia(win)),s_pilot.MOCup_pilot_low(ib(win))));
fprintf('--- Outside 2017-2019 (n=%d) ---\n',sum(~win));
fprintf('Full(restricted)=%.2f Sv, Pilot=%.2f Sv, corr=%.3f\n', ...
    nanmean(MOCup_low(ia(~win))),nanmean(s_pilot.MOCup_pilot_low(ib(~win))),corr2_nan(MOCup_low(ia(~win)),s_pilot.MOCup_pilot_low(ib(~win))));

% For reference, also compare against the FULL (all-8-gap) v4 result
s_full=load('moc_streamfunction_v4.mat','dt','MOCup_low');
[~,ia2,ib2]=intersect(dt,s_full.dt);
fprintf('\nFor reference, all-8-gap full array (v4) mean over same common dt: %.2f Sv\n',nanmean(s_full.MOCup_low(ib2)));

save('moc_test_no_P5P6.mat','dt','pre','Psi_raw','h_star','MOCup_raw','MOCup_low','n_gaps_valid')

function r=corr2_nan(a,b)
ok=~isnan(a)&~isnan(b);
cc=corrcoef(a(ok),b(ok));
r=cc(1,2);
end
