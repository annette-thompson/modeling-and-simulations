%% Optimization_Function_Handler
function [total_obj] = Optimization_Function_Handler(S,p_vec,ODE_options)
% Calculates objectives for optimizing model output with experimental data
%   Input:
%       S: structure containing all kinetic parameters
%       p_vec:
%       ODE_options:
%   Output:
%       total_obj:


%Calculates kinetic parameters
P = Param_Function(S);

% Various checks to make sure that parameters aren't becoming unreasonable
tp_vec = S.p_vec;
tp_vec(10) = abs(tp_vec(10));
tp_vec(11) = abs(tp_vec(11));
if any(tp_vec(tp_vec<=0))
    total_obj = 1E8;
    return
elseif tp_vec(4) > 6.29E4 || tp_vec(5) > 240 || tp_vec(6) > 1 || tp_vec(7) > 15
    total_obj = 1E8;
    return
elseif tp_vec(1) < .001 || tp_vec(2) < .001 || tp_vec(3) < .001 || tp_vec(5) < 0.01 || tp_vec(6) < 7.41E-05 || tp_vec(7) < 0.001 || tp_vec(8) < 2.46E-7 || tp_vec(9) < 2.46E-7
    total_obj = 1E8;
    return
end

P = orderfields(P);
fields = fieldnames(P);
for i = 1:length(fields)
    if ismember(fields(i),{'k2_1f','k2_3f','k3_1f','k3_3f','k4_1f','k4_2f','k5_1f','k6_1f','k6_2f','k7_1f','k8_1f','k8_3f'})
        if ismember(fields(i),{'k2_1f','k2_3f','k3_1f'})
            if P.(cell2mat(fields(i))) > 1650
                total_obj = 1E8;
                return
            elseif P.(cell2mat(fields(i+1)))/P.(cell2mat(fields(i))) < .01
                total_obj = 1E8;
                return 
            end
        elseif ismember(fields(i),{'k3_3f','k4_2f','k5_1f','k6_2f','k8_1f','k8_3f'})
            if P.(cell2mat(fields(i))) > 629
                total_obj = 1E8;
                return
            elseif P.(cell2mat(fields(i+1)))/P.(cell2mat(fields(i))) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(fields(i),{'k4_1f','k6_1f'})
            if P.(cell2mat(fields(i))) > 1650
                total_obj = 1E8;
                return
            elseif P.(cell2mat(fields(i+1)))/P.(cell2mat(fields(i))) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(fields(i),{'k7_1f'})
            if max(P.(cell2mat(fields(i)))) > 629
                total_obj = 1E8;
                return
            end
        end
    end
end

%% Calculation of Obj1 (New data with various combinations of FabH, FabF, FabB, and AcCoA removed)
% Finding initial rate (2.5 mins)
range_init = [0 150];

init_cond_1 = zeros(S.num,1);
init_cond_1(3) = 500;%s3 (Acetyl-CoA, 0.5mM)
init_cond_1(4) = 10;%s6 (holo ACP)
init_cond_1(5) = 1000;%s7 (NADPH)
init_cond_1(6) = 1000;%s8 (NADH)
init_cond_1(8) = 500;%p2 (malonyl-CoA)

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc_ref = [0 1 1 1 1 1 10 1 1 1];

[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_ref,range_init,init_cond_1,ODE_options);
rate(1) = F_weighted(end,1);
enz_conc_H = [0 1 0 1 1 1 10 1 1 1];
[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_H,range_init,init_cond_1,ODE_options);
rate(2) = F_weighted(end,1);
enz_conc_HF = [0 1 0 1 1 1 10 0 1 1];
[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_HF,range_init,init_cond_1,ODE_options);
rate(3) = F_weighted(end,1);
enz_conc_HB = [0 1 0 1 1 1 10 1 1 0];
[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_HB,range_init,init_cond_1,ODE_options);
rate(4) = F_weighted(end,1);

init_cond_2 = zeros(S.num,1);
init_cond_2(3) = 0;%s3 (Acetyl-CoA, 0.5mM)
init_cond_2(4) = 10;%s6 (holo ACP)
init_cond_2(5) = 1000;%s7 (NADPH)
init_cond_2(6) = 1000;%s8 (NADH)
init_cond_2(8) = 500;%p2 (malonyl-CoA)

[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_H,range_init,init_cond_2,ODE_options);
rate(5) = F_weighted(end,1);
[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_HF,range_init,init_cond_2,ODE_options);
rate(6) = F_weighted(end,1);
[~,F_weighted,~] = Optimization_ODE_Function(S,p_vec,enz_conc_HB,range_init,init_cond_2,ODE_options);
rate(7) = F_weighted(end,1);

rate_exp = S.opt_init_rate_data;

% SSE model and experimental inital rates
obj1 = sum((rate_exp - rate).^2);

%% Obj2 (Palmitic Acid Equiv. vs. Time)
% Final production (12 mins)
range_init_2 = [0 720];

[T,F_weighted,F_raw] = Optimization_ODE_Function(S,p_vec,enz_conc_ref,range_init_2,init_cond_1,ODE_options);

time_val = T/60;%conversion to minutes

% SSE model and experimental total production
ob2 = LeastSquaresCalc(time_val,F_weighted,S.opt_tot_prod_file);

%% Obj3 (FFA Profile)
% Final concentration (12 mins)
total_FA = sum(F_raw);%total fatty acid (concentration)
fit_dist = total_FA.*S.opt_prod_dist_data; %Expected distribution (experimental)

% SSE model and experimental profile
obj3 = sum((F_raw - fit_dist).^2);

%%
total_obj = obj1*ob2*obj3;

end