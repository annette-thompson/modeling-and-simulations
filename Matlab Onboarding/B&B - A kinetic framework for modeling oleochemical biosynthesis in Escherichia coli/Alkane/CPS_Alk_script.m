%All units are uM and sec.

%optimization vector
p_vec_k = [0.0154482766181948	94.5603925895629];

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(346,1);
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
enz_conc1 = [0  1    1    1    2.59018648976186	1   0   0.497550254250183	0.0394246097521289	1.52877593179746];

%Run the solver file
Alk = CPS_Alk_solv(init_cond,enz_conc1,time_range,p_vec_k);

%Plotting and Data Handling
exp = [0.24347826086956526 0.04766505636070858 0 0.7111111111111111];
Alk_ratio = Alk./sum(Alk);
exp_raw = [27.49347488 5.299014145 0 70.7286813];
exp_FA_titer = 0.7/256.4*10^6; %uM - reported g/L / MW of C16 * scaing factor
f = 33.85/exp_FA_titer;
exp_raw = f.*exp_raw;

figure()
barnames = {'C15','C17','C15:1','C17:1'};
bar([Alk./sum(Alk); exp]')
set(gca,'xticklabel',barnames)
legend('model','exp','location','best')
ylabel('Alkane Titer (uM)')
title('Alkanes, Figure 2A')
