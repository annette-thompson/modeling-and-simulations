%optimization vector
p_vec_k = [11.7210624038738	60.0916342073819	0.165729616633289	0.0578507289945410];

ket_cond = readtable('ketone_fitting_cond.xlsx','ReadVariableNames',true);
field1 = 'ket_cond';val1 = ket_cond;
field3 = 'tot_count';val3 = 7;
field4 = 'count';val4=7;
struct_k = struct(field1,val1,field3,val3,field4,val4);


%All units are uM and sec.

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(402,1);
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

%Run all models
cond_num = struct_k.tot_count;
ket_amt = zeros(cond_num,11);
FA_ratio = zeros(cond_num,5);
for i = 1:cond_num
    struct_k.count = i;
    enz_conc = [0  1    1    1    2.59018648976186	1   struct_k.ket_cond.TE_conc(i)   0.497550254250183	0.0394246097521289	1.52877593179746];
    [ket_amt(i,:),FA_ratio(i,:)] = CPS_Ket_solv(init_cond,enz_conc,time_range,struct_k,p_vec_k);
end

%Data Handling and Plotting
figure()
fig2b1exp = [0.5000	0.4000	0.1300	0.9700	0];
fig2b2exp = [0.0780	0.0410	0.0180	0.1700	0.0052];
fig2b1 = [fig2b1exp; FA_ratio(6,:)]'; 
subplot(1,2,1)
barnames = {'C12/C11','C14/C13','C14:1/C13:1','C16/C15','C16:1/C15:1'};
bar(fig2b1)
set(gca,'xticklabel',barnames)
ylabel('Ratio of Fatty Acid to Ketone')
legend('exp','model')
title('Ketones, Figure 2B')

subplot(1,2,2)
barnames = {'C12/C11','C14/C13','C14:1/C13:1','C16/C15','C16:1/C15:1'};
fig2b2 = [fig2b2exp; FA_ratio(7,:)]';
bar(fig2b2)
set(gca,'xticklabel',barnames)
ylabel('Ratio of Fatty Acid to Ketone')
legend('exp','model')
title('Ketones, Figure S6')

exp_titers =    [0 22.77106323	0	0	0 0 0;...
                0 19.79330881	23.83312594	26.54295613	0 0 0;...
                0 192.9409704	23.83312594	36.52592636	0 0 0;...
                0 1875.284638	0	36.52592636	0 0 0;...
                0 170.5202312	516.3843953	1221.269599	289.9998941 0 0]; %uM

ket_amt_cond = zeros(5,7);
for i = 1:5
    ket_amt_cond(i,:) = [ket_amt(i,1:3) ket_amt(i,4)+ket_amt(i,8) ket_amt(i,5)+ket_amt(i,9) ket_amt(i,6)+ket_amt(i,10) ket_amt(i,7)+ket_amt(i,11)];
end

A = 3000; %From Katie paper (um TesA Experiemntal)
D = 33.85; % from base model (uM TesA model)
f = D/A;

exp_titers_scaled = f.*exp_titers;

figure()
subplot(1,2,1)
bar(exp_titers_scaled, 'stacked')
ylabel('experimental titres, scaled (\muM)')
title('Figure 2A experimental ')
legend('C5','C7','C9','C11','C13','C15','C17')


subplot(1,2,2)
bar(ket_amt_cond(1:5,:),'stacked')
ylabel('model titres (\muM)')
title('Figure 2A model')
legend('C5','C7','C9','C11','C13','C15','C17')



