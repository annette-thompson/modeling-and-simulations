% %All units are uM and sec.

%optimization vector
p_vec_k = [0.183850387798547	2.02055368400942];



%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(368,1);
init_cond(1) = 1000; %s1 ATP
init_cond(3) = 500;%s3 (Acetyl-CoA)
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1000;%s7 (NADP)
init_cond(6) = 1000;%s8 (NADPH)
init_cond(8) = 500;%p2 (malonyl-CoA)
init_cond(9) = 0; %p3 (CoA)

compart=1E-3*6.7e-16;
init_cond(309) = 9.61452e-19/compart;%s3 (FtsH)
init_cond(310) = 6.39307e-19/compart;%s6 (LpxC)
init_cond(311) = 2.54062e-19/compart;%s8 (KdtA)

%Specify time range to solve system (seconds)
time_range = [0 720];

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    2.59018648976186	1   p_vec_k(2)     0.497550254250183	0.0394246097521289	1.52877593179746];

%run the solver file
FAEE = CPS_FAEE_solv(init_cond,enz_conc,time_range,p_vec_k);

%Plotting and Data Handling
sum_FAEE = sum(FAEE);
weight_vec = [8,10,12,14,16,18,12,14,16,18]/16;
FAEE_weighted = sum(FAEE.*weight_vec);
exp_FAEE_titer = 484.8326359832636/1000/284.5*10^6; %uM
exp_FA_titer = 0.7/256.4*10^6; %uM
f = 33.85/exp_FA_titer;

exp = [0.08723404255319149 0.2797872340425532 0.36914893617021277 0.251063829787234].*exp_FAEE_titer.*f;
FAEE_comb = [FAEE(3)+FAEE(7) FAEE(4)+FAEE(8) FAEE(5)+FAEE(9) FAEE(6)+FAEE(10)];
figure()
barnames = {'C12','C14','C16','C18'};
bar([exp; FAEE_comb]')
set(gca,'xticklabel',barnames)
legend('exp','model')
ylabel('FAEE Titer (uM)')
title('FAEE, Figure 3D')
