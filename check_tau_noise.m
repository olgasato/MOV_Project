load('/Users/olga/research/sambar/renellys_sent/Marions_code/MOV/samba_w.mat','dt','tau1000_A','tau1000_C','tau1000_D');
dt_W=dt; clear dt
load('/Users/olga/research/sambar/renellys_sent/Marions_code/ies/Wrk/IES_FrSA_SAMBA_6PIES.mat','dt','tau1000_P1','tau1000_P2','tau1000_P4','tau1000_P5','tau1000_P6','tau1000_P8');
dt_E=dt(:); clear dt

w1=datenum(2013,9,11); w2=datenum(2017,7,17);
names={'A','C','D','P8','P6','P5','P4','P2','P1'};
data={tau1000_A,tau1000_C,tau1000_D,tau1000_P8,tau1000_P6,tau1000_P5,tau1000_P4,tau1000_P2,tau1000_P1};
dts={dt_W,dt_W,dt_W,dt_E,dt_E,dt_E,dt_E,dt_E,dt_E};

fprintf('%-6s %12s %12s %12s\n','Site','std(diff)','std(raw)','ratio diff/raw (x1000)');
for i=1:9
    d=dts{i}; v=data{i};
    idx=d>=w1&d<=w2;
    vv=v(idx);
    dd=diff(vv);
    fprintf('%-6s %12.5f %12.5f %12.3f\n',names{i},nanstd(dd),nanstd(vv),1000*nanstd(dd)/nanstd(vv));
end
exit
