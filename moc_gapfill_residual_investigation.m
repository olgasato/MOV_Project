% Investigates the ~39% (10.14 Sv) of the 2017-2020 mean shortfall
% left unexplained by moc_gapfill_2017_2020_test.m's simple
% whole-period-climatology gap-fill. Two candidate explanations tested:
%
%   (a) The SURVIVING gaps (A-C, partial C-D/P4-P2, P2-P1 -- the ones
%       that DO have data during the outage) might carry their own
%       systematic BIAS during 2017-2020 specifically, not just their
%       usual noise level -- i.e. the days that survive might not be a
%       representative sample of "normal" conditions at those gaps.
%   (b) A single whole-period time-mean climatology may be a poor
%       proxy for the missing gaps' actual transport during 2017-2020
%       specifically (seasonal mismatch) -- tested by rebuilding the
%       gap-fill with a DAY-OF-YEAR climatology instead of a flat mean.
%
% Reuses the same Total_TPUD1-8/shelf/h_star setup as
% moc_gapfill_2017_2020_test.m -- no need to touch moc_pilot_v6.m or
% moc_streamfunction_v4.m.
%
% v1 (2026-08-28): first cut.

clear;clc;

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3'], ...
    'Total_TPUD1','Total_TPUD2','Total_TPUD3','Total_TPUD4','Total_TPUD5','Total_TPUD6','Total_TPUD7','Total_TPUD8', ...
    'Pressure','dt');
pre=Pressure(:); clear Pressure
dz=mean(diff(pre));

outage1=datenum(2017,8,1); outage2=datenum(2020,1,1);
normal=~(dt(:)>=outage1 & dt(:)<=outage2);
outage=~normal;

names={'A-C','C-D','D-P8','P8-P6','P6-P5','P5-P4','P4-P2','P2-P1'};
Total_all=cat(3,Total_TPUD1,Total_TPUD2,Total_TPUD3,Total_TPUD4,Total_TPUD5,Total_TPUD6,Total_TPUD7,Total_TPUD8);
clear Total_TPUD1 Total_TPUD2 Total_TPUD3 Total_TPUD4 Total_TPUD5 Total_TPUD6 Total_TPUD7 Total_TPUD8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (a) Do the SURVIVING gaps carry their own bias during the outage?
%%% Compare each gap's mean transport (where it HAS data) during the
%%% outage vs. its own normal-period mean, at 500dbar and depth-
%%% integrated (0-1190dbar, matching h_star).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iz500=find(pre==500);
h_star=1190; izh=find(pre<=h_star);

fprintf('%-8s %12s %12s %10s %14s %14s %10s\n','Gap','norm@500', 'outage@500','diff@500','norm 0-hstar','outage 0-hstar','diff int');
for i=1:8
    n500=nanmean(Total_all(iz500,normal,i));
    o500=nanmean(Total_all(iz500,outage,i));

    % depth-integrated (0-h_star) mean transport, Sv, using only days
    % with data at that gap (nanmean over time handles missing days
    % automatically; integrate the resulting time-mean PROFILE, not
    % integrate-then-average, to avoid double-counting partial days
    % oddly -- same convention as the streamfunction's own Psi_raw_mean)
    prof_n=nanmean(Total_all(izh,normal,i),2);
    prof_o=nanmean(Total_all(izh,outage,i),2);
    int_n=sum(prof_n)*dz/1e6;
    int_o=sum(prof_o)*dz/1e6;

    fprintf('%-8s %12.1f %12.1f %10.1f %14.2f %14.2f %10.2f\n',names{i},n500,o500,o500-n500,int_n,int_o,int_o-int_n);
end

fprintf('\n(Gaps 3-6 above are computed only from their scarce/zero remaining outage-days -- not very meaningful; the ones that matter for "surviving gap bias" are 1,2,7,8: A-C, C-D, P4-P2, P2-P1.)\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (b) Day-of-year (seasonal) climatology instead of a flat mean,
%%% for the gap-fill substitution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doy=day(datetime(dt,'convertfrom','datenum'),'dayofyear');
% Build a smooth (harmonic) seasonal climatology per gap/depth rather
% than a raw day-bin mean (which would be noisy with only ~7 normal
% years of data per calendar day). Fit annual + semi-annual harmonics
% by least squares at each depth, using normal days only.
t_frac=doy/365.25;
Xdes=[ones(size(t_frac)) cos(2*pi*t_frac) sin(2*pi*t_frac) cos(4*pi*t_frac) sin(4*pi*t_frac)];

clim_flat=nanmean(Total_all(:,normal,:),2); % [depth x 1 x 8], from part 1 above

Total_filled_seasonal=Total_all;
for i=1:8
    sl=Total_all(:,:,i);
    for zz=1:size(sl,1)
        y=sl(zz,normal)';
        ok=~isnan(y);
        if sum(ok)<10, continue; end
        b=Xdes(normal,:); b=b(ok,:);
        coef=b\y(ok);
        pred_all=Xdes*coef; % seasonal prediction for every day in the record
        row=Total_filled_seasonal(zz,:,i);
        badmask=isnan(row) & outage(:)';
        row(badmask)=pred_all(badmask)';
        Total_filled_seasonal(zz,:,i)=row;
    end
end

transport_actual=nansum(Total_all,3);
Total_filled_flat=Total_all;
for i=1:8
    sl=Total_filled_flat(:,:,i);
    badmask=isnan(sl) & repmat(outage(:)',size(sl,1),1);
    climrep=repmat(clim_flat(:,1,i),1,size(sl,2));
    sl(badmask)=climrep(badmask);
    Total_filled_flat(:,:,i)=sl;
end
transport_filled_flat=nansum(Total_filled_flat,3);
transport_filled_seasonal=nansum(Total_filled_seasonal,3);

load('/Users/olga/research/sambar/renellys_sent/Marions_code/ecco/Wrk/ECCO_TPUD_ShelfBits.mat','West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(pre,pres_shelf_grid);
shelf_profile=West_TPUD(ib)+East_TPUD_new(ib);
clear West_TPUD East_TPUD_new pres_shelf_grid ia ib

transport_actual=transport_actual+shelf_profile;
transport_filled_flat=transport_filled_flat+shelf_profile;
transport_filled_seasonal=transport_filled_seasonal+shelf_profile;

Psi_actual=cumsum(transport_actual,1,'omitnan')*dz/1e6;
Psi_flat=cumsum(transport_filled_flat,1,'omitnan')*dz/1e6;
Psi_seasonal=cumsum(transport_filled_seasonal,1,'omitnan')*dz/1e6;

iz=find(pre==h_star);
MOCup_actual=Psi_actual(iz,:);
MOCup_flat=Psi_flat(iz,:);
MOCup_seasonal=Psi_seasonal(iz,:);

s_pilot=load('moc_pilot_v6.mat','dt','MOCup_pilot_raw');
[dtc,ip,ipp]=intersect(dt,s_pilot.dt);
win=outage(ip);

pilot_mean=nanmean(s_pilot.MOCup_pilot_raw(ipp(win)));
actual_mean=nanmean(MOCup_actual(ip(win)));
flat_mean=nanmean(MOCup_flat(ip(win)));
seasonal_mean=nanmean(MOCup_seasonal(ip(win)));

fprintf('\nDuring the outage (2017-08/2020-01):\n');
fprintf('  Pilot v6:                    %.2f Sv\n',pilot_mean);
fprintf('  Full array, ACTUAL:          %.2f Sv\n',actual_mean);
fprintf('  Full array, FLAT climatology:%.2f Sv (closes %.0f%%)\n',flat_mean,100*(1-(pilot_mean-flat_mean)/(pilot_mean-actual_mean)));
fprintf('  Full array, SEASONAL clim.:  %.2f Sv (closes %.0f%%)\n',seasonal_mean,100*(1-(pilot_mean-seasonal_mean)/(pilot_mean-actual_mean)));

save('moc_gapfill_residual_investigation.mat','dt','MOCup_actual','MOCup_flat','MOCup_seasonal','outage')
