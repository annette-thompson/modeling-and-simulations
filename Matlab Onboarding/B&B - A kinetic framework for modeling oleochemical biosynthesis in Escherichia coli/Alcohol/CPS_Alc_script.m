%All units are uM and sec.

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(360,1);

init_cond(1) = 1000; %s1 ATP
init_cond(3) = 500;%s3 (Acetyl-CoA)
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1000;%s7 (NADP)
init_cond(6) = 1000;%s8 (NADPH)
init_cond(8) = 500;%p2 (malonyl-CoA)
compart=1E-3*6.7e-16;
init_cond(309) = 9.61452e-19/compart;%s3 (FtsH)
init_cond(310) = 6.39307e-19/compart;%s6 (LpxC)
init_cond(311) = 2.54062e-19/compart;%s8 (KdtA)

%Specify time range to solve system (seconds)
time_range = [0 720];

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    2.59018648976186	1   10 0.497550254250183	0.0394246097521289	1.52877593179746];

%Run the solver file
FA_alc = CPS_Alc_solv(init_cond,enz_conc,time_range);

%Plotting and Data Handling
exp = [0.42 0.28  0.05 0.01 0.1  0.12  0.15  0.04];
FA_ratio_fig = FA_alc(5:12)./sum(FA_alc(5:12));

figure()
barnames = {'C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
bar([exp; FA_ratio_fig]')
set(gca,'xticklabel',barnames)
legend('exp','model')
ylabel('Fraction Fatty Alcohol')
title('Fatty Alcohol, Figure 2A')
