%All units are uM and sec.

%optimization vector
p_vec_k = [0.235521196381529	75.3601070325288];

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(342,1);
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
enz_conc = [0  1    1    1    2.59018648976186	1   10     0.497550254250183	0.0394246097521289	1.52877593179746];

%run the solver file
FAME = CPS_FAME_solv(init_cond,enz_conc,time_range,p_vec_k);

%Plotting and Data Handling
exp = [907.1441169 114.6680163 0 334.0233473 131.5698543 0];
exp_FA_titer = 0.7/256.4*10^6; %uM - reported g/L / MW of C16 * scaing factor
f = 33.85/exp_FA_titer;
exp = f.*exp;

FAME_ratio = FAME./sum(FAME);
figure()
barnames = {'C12','C14','C16','C12:1','C14:1','C16:1'};
bar([exp; FAME]')
set(gca,'xticklabel',barnames)
legend('exp','model')
ylabel('FAME Titer (uM)')
title('FAME, Figure 3C')
