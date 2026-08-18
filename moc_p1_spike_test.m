% Checks whether site P1's raw tau1000 has instrument-glitch spikes
% like the three already found at site A (2020-06-20, 2022-10-17,
% 2022-10-18 -- see "Site A's spikes" in CLAUDE.md), using the exact
% same methodology: residual from a 31-day running median, flagged as
% an outlier past a 6-sigma threshold. Re-runs the same test on A too
% (for a clean side-by-side, same session/thresholding) plus P1,
% over the FULL available record for each site, then flags which (if
% any) land inside the paper's 2013-09-11/2017-07-17 window -- the
% only period that matters for the r=0.73/raw-correlation comparison
% this whole investigation has been chasing.
%
% v1 (2026-08-18): first cut.

clear;clc;
prefix_marion='/Users/olga/research/sambar/renellys_sent/Marions_code/';

load([prefix_marion,'MOV/samba_w.mat'],'dt','tau1000_A');
dt_A=dt(:); tauA=tau1000_A(:); clear dt tau1000_A

load([prefix_marion,'ies/Wrk/IES_FrSA_SAMBA_6PIES.mat'],'dt','tau1000_P1');
dt_P1=dt(:); tauP1=tau1000_P1(:); clear dt tau1000_P1

w1=datenum(2013,9,11); w2=datenum(2017,7,17);

sites={'A','P1'};
dts={dt_A,dt_P1};
taus={tauA,tauP1};

for s=1:2
    dt_s=dts{s}; tau_s=taus{s};
    med=movmedian(tau_s,31,'omitnan');
    resid=tau_s-med;
    sig=nanstd(resid);
    thresh=6*sig;
    out=find(abs(resid)>thresh);

    fprintf('\n=== Site %s: %d days (%s to %s), residual sigma=%.5f, 6-sigma threshold=%.5f ===\n', ...
        sites{s},numel(dt_s),datestr(dt_s(1)),datestr(dt_s(end)),sig,thresh);
    fprintf('%d outlier(s) found over the full record:\n',numel(out));
    for k=1:numel(out)
        i=out(k);
        inwin=dt_s(i)>=w1 & dt_s(i)<=w2;
        fprintf('  %s: tau1000=%.5f, residual=%.5f (%.1f sigma)%s\n', ...
            datestr(dt_s(i)),tau_s(i),resid(i),abs(resid(i))/sig, ...
            ternary(inwin,'  <-- INSIDE 2013-2017 window','  (outside 2013-2017 window)'));
    end

    n_inwin=sum(dt_s(out)>=w1 & dt_s(out)<=w2);
    fprintf('Of these, %d fall inside the 2013-2017 window.\n',n_inwin);
end

function r=ternary(cond,a,b)
if cond, r=a; else, r=b; end
end
