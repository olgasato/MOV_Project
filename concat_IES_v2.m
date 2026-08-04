% Step 1 for MOV calculation
% Concatenate and select the TS files from the Samba_West and Samba_East.
%

t0=datenum(2013,9,11);
t1=datenum(2017,7,17);

% ------------------------------------------------------------------
load samba_w
ind=find(dt>=t0 & dt<t1);

% Selecting the period: 11-Sep-2013 to 17-Jul-2017 (to match samba_e)
% ------------------------------------------------------------------
sitesW={'A','AA','B','BB','C','CC','D'};
for k=1:numel(sitesW)
    eval(['tem_' sitesW{k} ' = Temperature_' sitesW{k} '(:,ind);']);
    eval(['sal_' sitesW{k} ' = Salinity_' sitesW{k} '(:,ind);']);
end

% This is the bottom pressure
%pr_A=pres_A(ind);
%pr_AA=pres_AA(ind);
%pr_B=pres_B(ind);
%pr_BB=pres_BB(ind);
%pr_C=pres_C(ind);
%pr_CC=pres_CC(ind);
%pr_D=pres_D(ind);

% Pressure vector and depth of samba_w: A, AA,  to D
prew=presrange;
depw=[1350 2902 3510 4173 4558 4730 4756];
lonw=sort([lon_SAM lon_Brazil]);

dt1=dt(ind);
% ------------------------------------------------------------------
load samba_e

ind=find(dt>=t0 & dt<t1);

sitesE={'P1','P2','P4','P5','P6','P8'};
for k=1:numel(sitesE)
    eval(['tem_' sitesE{k} ' = Temperature_' sitesE{k} '(:,ind);']);
    eval(['sal_' sitesE{k} ' = Salinity_' sitesE{k} '(:,ind);']);
end

% Pressure vector and depth of samba_e: P1:8 (3 and 7 missing)
pree=presrange;
depe=[1266 2129 4482 4969 5185 4608];
lone=[17.5576 17.3006 15.0027 11.2031 7.4505 0.0];

dt=dt(ind);

clear Tempera* Salin* DynHgt* pres_* tau100* Err* ind* presrange lon_* dt1 sitesW sitesE k t0 t1

%save concat_IESsamba
