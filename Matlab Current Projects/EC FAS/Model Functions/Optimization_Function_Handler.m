%% Optimization_Function_Handler
function [total_obj] = Optimization_Function_Handler(p_vec,ODE_options)
% Calculates objectives for optimizing model output with experimental data
%   Input:
%       p_vec: vector guess with variables to be optimized
%       ODE_options: options for ODE solver
%   Output:
%       total_obj: objective to be minimized

% Assign kinetic parameters based on p_vec guess
S = set_vars_opt(p_vec);

% Limiting specific values needed to be controlled (based on what is being
% parametrized) - set to make sure largest value being scaled doesn't go
% over its diffusion limit
tp_vec = p_vec;
if tp_vec(7) < 0 || tp_vec(19) > 18.575 || tp_vec(19) < 0
%if tp_vec(19) > 18.575 || tp_vec(19) < 0
    total_obj = 1E8;
    return
end

% Initial rate (2.5 mins)
range_rate = [0 150];

% Final production (12 mins)
range_prod = [0 720];

% Initial conditions
% Reference
init_cond_1 = zeros(S.num,1);
init_cond_1(3) =    500; %  Acetyl-CoA
init_cond_1(12) =   10; %   holo ACP
init_cond_1(13) =   1000; % NADPH
init_cond_1(15) =   1000; % NADH
init_cond_1(18) =   500; %  Malonyl-CoA
% No Acetyl-CoA
init_cond_2 = zeros(S.num,1);
init_cond_2(3) =    0; %    Acetyl-CoA
init_cond_2(12) =   10; %   holo ACP
init_cond_2(13) =   1000; % NADPH
init_cond_2(15) =   1000; % NADH
init_cond_2(18) =   500; %  Malonyl-CoA
% Add ACC substrates, no Malonyl-CoA
init_cond_3 = zeros(S.num,1);
init_cond_3(1) =    1000; % ATP
init_cond_3(2) =    500; % Bicarbonate
init_cond_3(3) =    1000; %  Acetyl-CoA
init_cond_3(12) =   10; %   holo ACP
init_cond_3(13) =   1000; % NADPH
init_cond_3(15) =   1000; % NADH
init_cond_3(18) =   0; %    Malonyl-CoA

% Enzyme compositions
% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc_ref =  [0 1 1 1 1 1 10 1 1 1];%    Reference composition (no ACC)
enz_conc_H =    [0 1 0 1 1 1 10 1 1 1]; %   No FabH
enz_conc_HF =   [0 1 0 1 1 1 10 0 1 1]; %   No FabH or FabF
enz_conc_HB =   [0 1 0 1 1 1 10 1 1 0]; %   No FabH or FabB
enz_conc_ACC =  [1 1 1 1 1 1 10 1 1 1];%   Add ACC

%% Objective 1 (Initial Rates)
% Reference
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_ref,range_rate,init_cond_1,ODE_options);
rate(1) = F_weighted(end,1)/T(end)*60;
% No FabH
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_H,range_rate,init_cond_1,ODE_options);
rate(2) = F_weighted(end,1)/T(end)*60;
% No FabH or FabF
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_HF,range_rate,init_cond_1,ODE_options);
rate(3) = F_weighted(end,1)/T(end)*60;
% No FabH or FabB
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_HB,range_rate,init_cond_1,ODE_options);
rate(4) = F_weighted(end,1)/T(end)*60;
% No FabH or Acetyl-CoA
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_H,range_rate,init_cond_2,ODE_options);
rate(5) = F_weighted(end,1)/T(end)*60;
% No FabH or FabF or Acetyl-CoA
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_HF,range_rate,init_cond_2,ODE_options);
rate(6) = F_weighted(end,1)/T(end)*60;
% No FabH or FabB or Acetyl-CoA
[T,F_weighted,~] = Optimization_ODE_Function(S,enz_conc_HB,range_rate,init_cond_2,ODE_options);
rate(7) = F_weighted(end,1)/T(end)*60;

% Experimental values
rate_exp = S.opt_init_rate_data;

% SSE model and experimental inital rates
obj1 = sum((rate_exp - rate).^2);

% if obj1>100
%     total_obj = 1E6;
%     return
% end

%% Objective 2 (Palmitic Acid Equiv. vs. Time)
[T,F_weighted,F_raw_1] = Optimization_ODE_Function(S,enz_conc_ref,range_prod,init_cond_1,ODE_options);

time_val = T/60; % conversion to minutes

% SSE model and experimental total production
obj2 = LeastSquaresCalc(time_val,F_weighted,S.opt_tot_prod_file);

% if obj2>100
%     total_obj = 1E6;
%     return
% end

%% Objective 3 (FFA Profile)
[~,~,F_raw_2] = Optimization_ODE_Function_ACC_MM(S,enz_conc_ACC,range_prod,init_cond_3,ODE_options);

total_FA = sum(F_raw_1); % total fatty acid (concentration) from in vitro for scaling
fit_dist = total_FA.*S.opt_prod_dist_data; % Expected distribution (experimental)

% SSE model and experimental profile
%obj3 = sum((F_raw_1 - fit_dist).^2);
obj3 = sum((F_raw_2 - fit_dist).^2);

% if obj3>100
%     total_obj = 1E6;
%     return
% end

%% Total Objective
total_obj = obj1*obj2*obj3;

