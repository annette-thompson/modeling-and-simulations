%All units are uM and sec.

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(336,1);

% init_cond(3) = 500;%s3 (Acetyl-CoA)
% init_cond(4) = 10;%s6 (holo ACP)
% init_cond(5) = 1000;%s7 (NADPH)
% init_cond(6) = 1000;%s8 (NADH)
% init_cond(8) = 500;%p2 (malonyl-CoA)

init_cond(3) = 300;%s3 (Acetyl-CoA)
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1300;%s7 (NADPH)
init_cond(6) = 1300;%s8 (NADH)
init_cond(8) = 1500;%p2 (malonyl-CoA)

compart=1E-3*6.7e-16;
init_cond(309) = 9.61452e-19/compart;%s3 (FtsH)
init_cond(310) = 6.39307e-19/compart;%s6 (LpxC)
init_cond(311) = 2.54062e-19/compart;%s8 (KdtA)

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)

%Specify time range to solve system (seconds)
time_range = [0 3600];

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
% enz_conc = [0  1    1    1    1    1    10    1    1    1]; %Reference
enz_conc = [0  1    10   1    1    1    10   0.1    1    1];

[total_FA,F_weighted,frac_FA] = CPS_Base_solv(init_cond,enz_conc,time_range);

%% Plotting can be done below.
bar(frac_FA)


