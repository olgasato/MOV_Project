function yg=lowpass_filter_fixed(var,npf)
% Mean-preserving 45-day lowpass filter -- identical to
% ~/matlab/stat/lowpass_filter.m EXCEPT it omits that function's two
% post-hoc amplitude-rescaling steps (via maxvar.m / an inline
% no-intercept least-squares gain), which were found (2026-08-17) to
% inflate the MEAN along with the variance whenever the underlying
% signal has a nonzero mean -- see MOV_Project/CLAUDE.md's "MAJOR
% FINDING" section for the full root-cause derivation. Confirmed there
% that the Blackman-window convolution + edge correction below is, on
% its own, exactly mean-preserving (17.7001 Sv vs. raw 17.7008 Sv for
% the test case that exposed the bug).
%
% Kept as a SEPARATE function (not editing the shared
% ~/matlab/stat/lowpass_filter.m in place) -- that's a shared utility
% used throughout other parts of this project too, and the user chose
% to scope/fix it project-wide in a later session rather than change
% it now. This function is for the MOC-specific investigation only.

um=ones(size(var));
inxs=(npf-1)/2;
myfilter=blackman(npf);
myfilter=myfilter/sum(myfilter);
y=conv(var,myfilter);
umf=conv(um,myfilter);
y=y(inxs+1:length(y)-inxs);
umf=umf(inxs+1:length(umf)-inxs);
yg=y./umf;
end
