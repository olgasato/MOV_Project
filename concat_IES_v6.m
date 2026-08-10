%% Step 1 for MOV calculation: opens the date window to the full overlap
% available after the 2026-08 rebuild of samba_e.mat (see CLAUDE.md).
%
% v5 hardcoded 11-Sep-2013 to 17-Jul-2017, matching the East array's
% availability at the time. Since then:
%   - samba_w.mat already covered 2009-03-18 to 2022-12-11 the whole
%     time (it just wasn't being used past 2017 due to this window).
%   - samba_e.mat was rebuilt from Marions_code's Daily_Tau-based
%     pipeline (a published, already-calibrated dataset) and now
%     covers 2013-09-06 to 2023-09-24, instead of stopping at 2017.
% The real overlap between the two is now bounded by samba_w.mat's
% end date (2022-12-11), NOT samba_e.mat's -- West is the limiting
% side now, the opposite of the v5-era situation.
%
% This script extracts temperature and salinity data from Samba moorings
% for the common time period: 06-Sep-2013 to 11-Dec-2022
%
% Temperature/salinity are kept as separate per-mooring variables
% (tem_A, tem_AA, ..., sal_P1, sal_P2, ...) rather than concatenated
% into one matrix, since mov_samba.m computes dynamic height anomaly
% per mooring against that mooring's own pressure grid.

clear; clc;

%% Define time period (common to both datasets)
t0 = datenum(2013, 9, 6);
t1 = datenum(2022, 12, 11);

%% ========================================================================
% Process Samba_West
% Reminder: The same as IES_Make_Profiles_plusBrazil_plusDynHgt_Olga.mat
% on folder: /Users/olga/research/sambar/renellys_sent/IES profile data
% ========================================================================
fprintf('Processing Samba_West...\n');
load('samba_w.mat');

% Select time period
ind = find(dt >= t0 & dt < t1);

% Per-mooring temperature/salinity
sitesW = {'A','AA','B','BB','C','CC','D'};
for k = 1:numel(sitesW)
    eval(['tem_' sitesW{k} ' = Temperature_' sitesW{k} '(:,ind);']);
    eval(['sal_' sitesW{k} ' = Salinity_' sitesW{k} '(:,ind);']);
end

% Optional: Bottom pressure (uncomment if needed)
% pr_A=pres_A(ind); pr_AA=pres_AA(ind); pr_B=pres_B(ind);
% pr_BB=pres_BB(ind); pr_C=pres_C(ind); pr_CC=pres_CC(ind); pr_D=pres_D(ind);

% Metadata: pressure vector and depth of samba_w, sites A, AA, ... to D
prew = presrange;
depw = [1350, 2902, 3510, 4173, 4558, 4730, 4756];
lonw = sort([lon_SAM, lon_Brazil]);
dt1 = dt(ind);

%% ========================================================================
% Process Samba_East
% ========================================================================
fprintf('Processing Samba_East...\n');
load('samba_e.mat');

% Select time period
ind = find(dt >= t0 & dt < t1);

% Per-mooring temperature/salinity (P3 and P7 are missing)
sitesE = {'P1','P2','P4','P5','P6','P8'};
for k = 1:numel(sitesE)
    eval(['tem_' sitesE{k} ' = Temperature_' sitesE{k} '(:,ind);']);
    eval(['sal_' sitesE{k} ' = Salinity_' sitesE{k} '(:,ind);']);
end

% Optional: Bottom pressure (uncomment if needed)
% pr_P1=pres_P1(ind); pr_P2=pres_P2(ind); pr_P4=pres_P4(ind);
% pr_P5=pres_P5(ind); pr_P6=pres_P6(ind); pr_P8=pres_P8(ind);

% Metadata: pressure vector and depth of samba_e, sites P1:8 (3 and 7 missing)
pree = presrange;
depe = [1266, 2129, 4482, 4969, 5185, 4608];
lone = [17.5576, 17.3006, 15.0027, 11.2031, 7.4505, 0.0];
dt = dt(ind);

%% ========================================================================
% Quick summary (optional)
% ========================================================================
fprintf('\nSummary:\n');
fprintf('Samba_West: %d moorings (%s)\n', numel(sitesW), strjoin(sitesW, ', '));
fprintf('Samba_East: %d moorings (%s)\n', numel(sitesE), strjoin(sitesE, ', '));
fprintf('Time period: %s to %s\n', datestr(t0), datestr(t1));

%% ========================================================================
% Clean up workspace (remove raw loaded variables no longer needed)
% ========================================================================
clear Tempera* Salin* DynHgt* pres_* tau100* Err* ind* presrange lon_* dt1 sitesW sitesE k t0 t1

%% ========================================================================
% Save processed data
% ========================================================================
fprintf('Saving data to concat_IESsamba.mat...\n');
save('concat_IESsamba.mat', '-v7.3');
fprintf('Done! Data saved successfully.\n');
