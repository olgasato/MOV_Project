% Plots raw tau1000 (acoustic travel time, the fundamental PIES/CPIES
% measurement -- everything else, GEM-derived T/S, Gpan, BPR-referenced
% transport, is built on top of this) for all 9 sites A through P1, in
% west-to-east geographic order, restricted to the 2013-2017 window
% currently under investigation (correlation shortfall between the
% full array and Pilot methods). Purpose: look at the raw data directly
% before investigating further -- does anything look visually odd
% (jumps, dropouts, drift) at any specific site that might explain the
% excess Pilot noise / correlation shortfall already documented in
% CLAUDE.md?

clear;clc;

prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';

load([prefix_marion,'MOV/samba_w.mat'],'dt','tau1000_A','tau1000_C','tau1000_D');
dt_W=dt; clear dt

load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'dt','tau1000_P1','tau1000_P2','tau1000_P4','tau1000_P5','tau1000_P6','tau1000_P8');
dt_E=dt(:); clear dt

w1=datenum(2013,9,11); w2=datenum(2017,7,17);

names={'A','C','D','P8','P6','P5','P4','P2','P1'};
data={tau1000_A,tau1000_C,tau1000_D,tau1000_P8,tau1000_P6,tau1000_P5,tau1000_P4,tau1000_P2,tau1000_P1};
dts={dt_W,dt_W,dt_W,dt_E,dt_E,dt_E,dt_E,dt_E,dt_E};

figure('Position',[50 50 1000 1400])
for i=1:9
    subplot(9,1,i)
    d=dts{i}; v=data{i};
    idx=d>=w1&d<=w2;
    plot(d(idx),v(idx),'.','markersize',3)
    ylabel(names{i},'fontsize',11,'fontweight','bold')
    xlim([w1 w2])
    grid on
    set(gca,'fontsize',9)
    if i<9
        set(gca,'xticklabel',[])
    else
        datetick('x',12,'keeplimits')
        xlabel('Date')
    end
end
sgtitle('Raw \tau_{1000} (s), all 9 sites A\rightarrowP1, 2013-09-11 to 2017-07-17','fontsize',13,'fontweight','bold')
print -dpng -r150 raw_tau1000_all_sites_2013_2017
