%All units are uM and sec.

%optimization vector when using TesA
p_vec_k = [0.183998491655396	0.680651352240870	132.988017008123 0.146865102220244]; %% TesA

% optimization vector for when using BTE (uncomment as necessary)
% p_vec_k = [0.802490234375000  0.680651352240870	132.988017008123 0.146865102220244]; %% BTE


%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(380,1);
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
enz_conc = [0  1    1    1    2.59018648976186	1   p_vec_k(1)   0.497550254250183	   0.0394246097521289	   1.52877593179746];

%run the solver file
FA_alc = CPS_Alc_solv(init_cond,enz_conc,time_range,p_vec_k);

%Plotting and Data Handling
FA_alc_fig = [FA_alc(1:3) FA_alc(4)+FA_alc(8) FA_alc(5)+FA_alc(9) FA_alc(6)+FA_alc(10) FA_alc(7)+FA_alc(11)];

exp_T = [0 0 0 33.54240849 185.7939824 45.44408747 0];
exp_FA_titer = 0.7/256.4*10^6; %uM - FAEE paper
f = 33.85/exp_FA_titer;
exp_T = f.*exp_T;

figure()
barnames = {'C6','C8','C10','C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
bar([exp_T; FA_alc_fig]')
set(gca,'xticklabel',barnames)
legend('exp','model')
ylabel('Fatty Alcohol Titer - TesA')
title('Fatty Alcohol through FadD/ACR1/AHR, Figure S9D')

%BTE plotting (uncomment as necessary)
% exp_BTE = [0 0 17.93397764 284.9208219 6.865345734 0 0];
% exp_FA_titer = 0.7/256.4*10^6; %uM
% f = 33.85/exp_FA_titer*0.35;
% exp_BTE = f.*exp_BTE;
% 
% figure()
% barnames = {'C6','C8','C10','C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
% bar([exp_BTE; FA_alc_fig]')
% set(gca,'xticklabel',barnames)
% legend('exp','model')
% ylabel('Fatty Alcohol Titer - BTE')
% title('Fatty Alcohol through FadD/ACR1/AHR')
