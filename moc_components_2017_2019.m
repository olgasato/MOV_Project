% Component-breakdown test: does the 2017-2019 dip show up in ALL of
% the baroclinic (relative), barotropic (reference/BPR), and Ekman
% components individually, or is it specific to just one?
%
% Uses Step E v3's already-computed "held constant" variants (Marion's
% own original design, per her Step E header comment: "calculates the
% MOC holding alternately the reference, relative or Ekman contribution
% constant in order to see which is driving the MOC/MHT variations at
% different time scales"):
%   Total_TPUD_i            -- everything varies daily (the normal case)
%   Total_TPUD_i_RefConst   -- reference/barotropic held at its own time-
%                              mean; relative+Ekman still vary daily
%   Total_TPUD_i_RelConst   -- relative/baroclinic held at its own time-
%                              mean; reference+Ekman still vary daily
%   Total_TPUD_i_EKMANconst -- Ekman held at its own time-mean;
%                              relative+reference still vary daily
%
% Same Kersale-definition/h_star/shelf-correction construction as
% moc_streamfunction_v4.m, applied to all four in parallel, using Step E
% v3 (calibration-window-fixed) as the source.
%
% Expectation if the dip is a "missing width" effect (as already
% established -- gaps 3-6 fully/mostly absent 2017-2019): since ALL
% FOUR variants are built from the SAME underlying per-gap data (only
% the CONSTANT-vs-VARYING treatment of one component differs -- missing
% gaps are missing in ALL FOUR), all four should show a similarly
% severe dip. If instead only ONE variant's dip disappears when a
% specific component is held constant, that component would be the
% real driver -- a genuinely different, more specific finding.
%
% v1 (2026-08-17): first cut.

clear;clc;

addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/positions_pies')
addpath('/Users/olga/research/sambar/renellys_sent/Marions_code/functions/seawater/')
addpath('/Users/olga/matlab/stat');
addpath('/Users/olga/research/claude/MOV_Project');

prefix='/Users/olga/research/sambar/renellys_sent/Marions_code/Full_Depth_MOC/Wrk/';
load([prefix,'Full_Depth_OverturningEstimate_Cst_for_MHT_2013_2022_v3']);

AREA_TOPO_1=AREA_TOPO_1(:); AREA_TOPO_2=AREA_TOPO_2(:); AREA_TOPO_3=AREA_TOPO_3(:); AREA_TOPO_4=AREA_TOPO_4(:);
AREA_TOPO_5=AREA_TOPO_5(:); AREA_TOPO_6=AREA_TOPO_6(:); AREA_TOPO_7=AREA_TOPO_7(:); AREA_TOPO_8=AREA_TOPO_8(:);

position_West
position_East
lon=[cpiesW([1 5 6]).lon cpiesE([1 2 4 5 6 8]).lon];
lat=[cpiesW([1 5 6]).lat cpiesE([1 2 4 5 6 8]).lat];
dx=gsw_distance(lon,lat,0);

[l,m]=size(Total_TPUD1);
n=8;

AREA_TOPO_3d=nan*ones([l 1 n]);
for ii=1:8
str=['AREA_TOPO_3d(:,1,',num2str(ii),')=AREA_TOPO_',num2str(ii),'(:);'];
eval(str)
end
AREA_TOPO_3d=repmat(AREA_TOPO_3d,[1 m 1]);

dx1=repmat(reshape(dx,1,1,[]),l,m,1);

pre=Pressure;
dz=mean(diff(pre));

load('/Users/olga/research/sambar/renellys_sent/Marions_code/ecco/Wrk/ECCO_TPUD_ShelfBits.mat','West_TPUD','East_TPUD_new');
pres_shelf_grid=(0:10:5300)';
[~,ia,ib]=intersect(pre,pres_shelf_grid);
shelf_profile=zeros(size(pre));
shelf_profile(ia)=West_TPUD(ib)+East_TPUD_new(ib);

variants={'','_RefConst','_RelConst','_EKMANconst'};
labels={'Full (all varying)','Reference held constant','Relative held constant','Ekman held constant'};
colors={[0 0.447 0.741],[0.85 0.325 0.098],[0.929 0.694 0.125],[0.494 0.184 0.556]};

MOCup_all=nan(numel(variants),m);
for vv=1:numel(variants)
    geo=nan*ones([l m n]);
    for ii=1:8
        str=['geo(:,:,',num2str(ii),')=Total_TPUD',num2str(ii),variants{vv},';'];
        eval(str)
    end
    % geo here already includes AREA_TOPO+Ekman (it's Total_TPUD, the
    % final per-gap product) -- just need dx1 for the NaN mask (same
    % mask for all four, from the plain geo/AREA_TOPO NaN pattern)
    dx1v=dx1; dx1v(isnan(geo)|isnan(AREA_TOPO_3d))=NaN;
    dxAv=dx1v.*AREA_TOPO_3d;
    dx2v=squeeze(nansum(dxAv,3));
    transport_profile=nansum(geo,3)+shelf_profile; % geo already = Total_TPUD (topo+ekman applied per gap), sum over gaps directly
    Psi=cumsum(transport_profile,1,'omitnan')*dz/1e6;
    Psi_mean=nanmean(Psi,2);
    [~,imax]=max(Psi_mean);
    MOCup_raw=Psi(imax,:);
    MOCup_all(vv,:)=lowpass_filter_fixed(MOCup_raw,45);
    fprintf('%-28s h_star=%d dbar, mean=%.2f Sv\n',labels{vv},pre(imax),nanmean(MOCup_all(vv,:)));
end

figure('Position',[100 100 950 550])
hold on
w1=datenum(2017,1,1); w2=datenum(2020,6,1);
for vv=1:numel(variants)
    plot(dt,MOCup_all(vv,:),'color',colors{vv},'linewidth',1.3)
end
patch([datenum(2017,8,1) datenum(2020,1,1) datenum(2020,1,1) datenum(2017,8,1)],[-60 -60 60 60],[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off')
uistack(findobj(gca,'type','patch'),'bottom')
legend(labels,'Location','southoutside','Orientation','horizontal')
xlim([w1 w2])
title('MOC_{up}, holding each component constant in turn (gray = 2017-2019 outage)','fontsize',13)
ylabel('MOC_{up} (Sv)')
datetick('x',12,'keeplimits');grid
set(gca,'linewidth',1,'fontsize',12)
print -dpng -r150 moc_components_2017_2019_zoom

figure('Position',[100 100 950 550])
hold on
for vv=1:numel(variants)
    plot(dt,MOCup_all(vv,:),'color',colors{vv},'linewidth',1)
end
legend(labels,'Location','southoutside','Orientation','horizontal')
title('MOC_{up}, holding each component constant in turn -- full record','fontsize',13)
ylabel('MOC_{up} (Sv)')
datetick('x',12);grid;axis('tight')
set(gca,'linewidth',1,'fontsize',12)
print -dpng -r150 moc_components_full_record

save('moc_components_2017_2019.mat','dt','pre','MOCup_all','labels')
