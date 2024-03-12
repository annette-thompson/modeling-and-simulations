%Script for manipulating model inputs easily
%All units are uM and sec.

%   palm_equiv: vector of palmitic acid equivalents at each time point
%   in time vals.
%   dist: vector of the amount of fatty acid of each chain length at each
%   time point in time vals (labels below).
%   unsat: fraction of unsaturated fatty acid at each time point in time
%   vals.
%   time_vals: vector of time points in the numerical solution (seconds)
%   conc_vals: matrix of concentration values for all components of the system
%   at the time points in time vals (time_vals,component). Use the labels
%   below to identify the component.

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(304,1);

init_cond(3) = 500;%s3 (Acetyl-CoA)
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1000;%s7 (NADPH)
init_cond(6) = 1000;%s8 (NADH)
init_cond(8) = 500;%p2 (malonyl-CoA)

%Specify time range to solve system (seconds)
time_range = [0 720];

%% Wild Type FabB = 1 uM

%Specify enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    1    1    10   1    1    1];

%Run the model
[~,dist,~,~,~] = Combined_Pathway_Solver_function(init_cond,enz_conc,time_range);

%labels for distribution
%dist_labels = {'4','6','8','10','12','12:1','14','14:1','16','16:1','18','18:1','20','20:1'};

frac_FA_WT = zeros(1,9);
frac_FA_WT(1:4) = dist(end,1:4);
frac_FA_WT(5) = dist(end,5)+dist(end,6);
frac_FA_WT(6) = dist(end,7)+dist(end,8);
frac_FA_WT(7) = dist(end,9)+dist(end,10);
frac_FA_WT(8) = dist(end,11)+dist(end,12);
frac_FA_WT(9) = dist(end,13)+dist(end,14);
frac_FA_WT = frac_FA_WT./sum(frac_FA_WT);

%% C8 Mutant FabB = 1 uM

%Specify enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    1    1    10   1    1    1];

%Run the model
[~,dist,~,~,~] = Combined_Pathway_Solver_function_mutant(init_cond,enz_conc,time_range);

frac_FA_C8_1 = zeros(1,9);
frac_FA_C8_1(1:4) = dist(end,1:4);
frac_FA_C8_1(5) = dist(end,5)+dist(end,6);
frac_FA_C8_1(6) = dist(end,7)+dist(end,8);
frac_FA_C8_1(7) = dist(end,9)+dist(end,10);
frac_FA_C8_1(8) = dist(end,11)+dist(end,12);
frac_FA_C8_1(9) = dist(end,13)+dist(end,14);
frac_FA_C8_1 = frac_FA_C8_1./sum(frac_FA_C8_1);

%% C8 Mutant FabB = 0.25 uM

%Specify enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    1    1    10   1    1    0.25];

%Run the model
[~,dist,~,~,~] = Combined_Pathway_Solver_function_mutant(init_cond,enz_conc,time_range);

frac_FA_C8_25 = zeros(1,9);
frac_FA_C8_25(1:4) = dist(end,1:4);
frac_FA_C8_25(5) = dist(end,5)+dist(end,6);
frac_FA_C8_25(6) = dist(end,7)+dist(end,8);
frac_FA_C8_25(7) = dist(end,9)+dist(end,10);
frac_FA_C8_25(8) = dist(end,11)+dist(end,12);
frac_FA_C8_25(9) = dist(end,13)+dist(end,14);
frac_FA_C8_25 = frac_FA_C8_25./sum(frac_FA_C8_25);

%% C8 Mutant FabB = 0 uM

%Specify enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    1    1    10   1    1    0];

%Run the model
[~,dist,~,~,~] = Combined_Pathway_Solver_function_mutant(init_cond,enz_conc,time_range);

frac_FA_C8_0 = zeros(1,9);
frac_FA_C8_0(1:4) = dist(end,1:4);
frac_FA_C8_0(5) = dist(end,5)+dist(end,6);
frac_FA_C8_0(6) = dist(end,7)+dist(end,8);
frac_FA_C8_0(7) = dist(end,9)+dist(end,10);
frac_FA_C8_0(8) = dist(end,11)+dist(end,12);
frac_FA_C8_0(9) = dist(end,13)+dist(end,14);
frac_FA_C8_0 = frac_FA_C8_0./sum(frac_FA_C8_0);

%% Plots

%x values for bar chart
x = 4:2:20;

figure;
bar(x,[frac_FA_WT;frac_FA_C8_1])
ylim([0,0.8])
legend('WT','C8 Mutant')
xlabel('Fatty Acid Chain Length')
ylabel('Mole Fraction')
legend('WT','C8 Mutant')


figure;
bar(x,[frac_FA_C8_1;frac_FA_C8_25;frac_FA_C8_0])
ylim([0,0.8])
legend('WT','C8 Mutant')
xlabel('Fatty Acid Chain Length')
ylabel('Mole Fraction')
legend('FabB = 1 um','FabB = 0.25 um','FabB = 0 um')
