% Same as check_raw_tau1000_all_sites.m but for the FULL 2013-2022
% record, explicitly regridded onto a complete daily calendar axis so
% missing days show as real gaps (NaN) in the plotted line, rather than
% being invisibly bridged by connecting the nearest surrounding valid
% points (which is what would happen if dt simply skips missing days
% without a NaN placeholder).

clear;clc;

prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';

load([prefix_marion,'MOV/samba_w.mat'],'dt','tau1000_A','tau1000_C','tau1000_D');
dt_W=dt(:); clear dt

load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'dt','tau1000_P1','tau1000_P2','tau1000_P4','tau1000_P5','tau1000_P6','tau1000_P8');
dt_E=dt(:); clear dt

names={'A','C','D','P8','P6','P5','P4','P2','P1'};
data={tau1000_A,tau1000_C,tau1000_D,tau1000_P8,tau1000_P6,tau1000_P5,tau1000_P4,tau1000_P2,tau1000_P1};
dts={dt_W,dt_W,dt_W,dt_E,dt_E,dt_E,dt_E,dt_E,dt_E};

w1=datenum(2013,9,6); w2=datenum(2022,12,11);
dt_full=(w1:w2)'; % complete daily calendar grid, no gaps

figure('Position',[50 50 1000 1400])
for i=1:9
    d=dts{i}(:); v=data{i}(:);
    v_full=NaN(size(dt_full));
    [tf,loc]=ismember(round(d),dt_full);
    v_full(loc(tf))=v(tf);

    subplot(9,1,i)
    plot(dt_full,v_full,'-','linewidth',0.5)
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
sgtitle('Raw \tau_{1000} (s), all 9 sites A\rightarrowP1, 2013-09-06 to 2022-12-11 (gaps = NaN, not bridged)','fontsize',13,'fontweight','bold')
print -dpng -r150 raw_tau1000_all_sites_2013_2022
