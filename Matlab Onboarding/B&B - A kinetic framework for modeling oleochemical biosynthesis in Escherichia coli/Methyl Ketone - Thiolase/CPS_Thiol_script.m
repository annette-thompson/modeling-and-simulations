%optimization Vector
p_vec_k = [11.7210624038738	60.0916342073819	0.165729616633289	0.0578507289945410];

%Optimization value for FadA
p_vec_new = 3.48730468750001;
ket_cond = readtable('Thiol_fitting_cond.xlsx','ReadVariableNames',true);
field1 = 'ket_cond';val1 = ket_cond;
field3 = 'tot_count';val3 = 1;
field4 = 'count';val4=1;
struct_k = struct(field1,val1,field3,val3,field4,val4);

%All units are uM and sec.

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(406,1);
init_cond(1) = 1000;
init_cond(3) = 500;%s3 (Acetyl-CoA)
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1000;%s7 (NADPH)
init_cond(6) = 1000;%s8 (NADH)
init_cond(8) = 500;%p2 (malonyl-CoA)
init_cond(9) = 0;

compart=1E-3*6.7e-16;
init_cond(309) = 9.61452e-19/compart;%s3 (FtsH)
init_cond(310) = 6.39307e-19/compart;%s6 (LpxC)
init_cond(311) = 2.54062e-19/compart;%s8 (KdtA)

%Specify time range to solve system (seconds)
time_range = [0 720];

%Run the model
cond_num = struct_k.tot_count;
ket_amt = zeros(cond_num,13);
FA_ratio = zeros(cond_num,5);
for i = 1:cond_num
    struct_k.count = i;
    enz_conc = [0  1    1    1    2.59018648976186	1   struct_k.ket_cond.TE_conc(i)   0.497550254250183	0.0394246097521289	1.52877593179746];
    [ket_amt(i,:),FA_ratio(i,:)] = CPS_Thiol_solv(init_cond,enz_conc,time_range,struct_k,p_vec_k,p_vec_new);
end

ket_amt_cond = [ket_amt(1:3) ket_amt(4)+ket_amt(8) ket_amt(5)+ket_amt(9) ket_amt(6)+ket_amt(10) ket_amt(7)+ket_amt(11) ket_amt(12)+ket_amt(13)];

A = 3000; %From Katie's paper (uM TesA Experiemntal)
D = 33.85; % from base model (uM TesA model)
f = D/A;
exp_titers = [0	1824.839102	8119.965712	241.6824604	0	0	0 0 ; 0 0 0 0 0 0 0 0];
% exp_titers_scaled = f.*exp_titers;

figure()
subplot(1,2,1)
bar(exp_titers, 'stacked')
ylabel('experimental titers (\muM)')
title('FadA experimental, Figure S8')
legend('C5','C7','C9','C11','C13','C15','C17','C19')

subplot(1,2,2)
bar([ket_amt_cond; 0 0 0 0 0 0 0 0],'stacked')
ylabel('model titers (\muM)')
title('FadA model, Figure S8')
legend('C5','C7','C9','C11','C13','C15','C17','C19')



