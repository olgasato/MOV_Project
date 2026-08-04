% Step 1 for MOV calculation
% Concatenate and select the TS files from the Samba_West and Samba_East.
% 

load samba_w
ind=find(dt>=datenum(2013,9,11) & dt<datenum(2017,7,17));

% Selecting the period: 09/06/2013 to 21/07/2017 (to match samba_e)
% ------------------------------------------------------------------
tem_A=Temperature_A(:,ind);
tem_AA=Temperature_AA(:,ind);
tem_B=Temperature_B(:,ind);
tem_BB=Temperature_BB(:,ind);
tem_C=Temperature_C(:,ind);
tem_CC=Temperature_CC(:,ind);
tem_D=Temperature_D(:,ind);

sal_A=Salinity_A(:,ind);
sal_AA=Salinity_AA(:,ind);
sal_B=Salinity_B(:,ind);
sal_BB=Salinity_BB(:,ind);
sal_C=Salinity_C(:,ind);
sal_CC=Salinity_CC(:,ind);
sal_D=Salinity_D(:,ind);

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

ind1=6;
ind2=1411;

dt1=dt(ind);
% ------------------------------------------------------------------
load samba_e

ind=find(dt>=datenum(2013,9,11) & dt<datenum(2017,7,17));

tem_P1=Temperature_P1(:,ind);
tem_P2=Temperature_P2(:,ind);
tem_P4=Temperature_P4(:,ind);
tem_P5=Temperature_P5(:,ind);
tem_P6=Temperature_P6(:,ind);
tem_P8=Temperature_P8(:,ind);

sal_P1=Salinity_P1(:,ind);
sal_P2=Salinity_P2(:,ind);
sal_P4=Salinity_P4(:,ind);
sal_P5=Salinity_P5(:,ind);
sal_P6=Salinity_P6(:,ind);
sal_P8=Salinity_P8(:,ind);

% Pressure vector and depth of samba_e: P1:8 (3 and 7 missing)
pree=presrange;
depe=[1266 2129 4482 4969 5185 4608];
lone=[17.5576 17.3006 15.0027 11.2031 7.4505 0.0];

dt=dt(ind);
	 

clear Tempera* Salin* DynHgt* pres_* tau100* Err* ind* presrange lon_* dt1

%save concat_IESsamba 

