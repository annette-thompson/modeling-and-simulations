%All units are uM and sec.

%optimization vector
p_vec_k = [0.0224653079869399	0.321909503820145	1.29949048145055];

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(345,1);
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

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)

%Specify time range to solve system (seconds)
time_range = [0 720];

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    2.59018648976186	1   0    0.497550254250183	   0.0394246097521289	   1.52877593179746];

%run the solver file
FA_alc = CPS_Alc_solv(init_cond,enz_conc,time_range,p_vec_k);

%Plotting and Data Handling
exp = [0 1.04309037 1.543565773 0.655563998 21.08202425 131.1736761];
exp_FA_titer = 0.7/256.4*10^6; %uM - reported g/L / MW of C16 * scaing factor
f = 33.85/exp_FA_titer;
exp = f.*exp;

FA_ratio_all = FA_alc./sum(FA_alc);
FA_fig = [FA_alc(1:3) FA_alc(4)+FA_alc(7) FA_alc(5)+FA_alc(8) FA_alc(6)+FA_alc(9)];

figure()
barnames = {'C6','C8','C10','C12','C14','C16'};
bar([exp; FA_fig]')
set(gca,'xticklabel',barnames)
legend('exp','model')
ylabel('Fatty Alcohol Titer')
title('Fatty Alcohol through acyl-ACP Reductase, Figure S9A')
