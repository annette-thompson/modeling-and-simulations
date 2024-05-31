%% Optimization_ODE_Function
function [T,F_weighted,F_raw] = Optimization_ODE_Function(S,p_vec,enz_conc,range,init_cond,ODE_options)
% Runs the ODE_Function for optimization and calculates final FA production
% for optimization comparison
%   Input:
%       S: structure containing all kinetic parameters
%       p_vec: 
%       enz_conc:
%       range:
%       init_cond:
%       ODE_options:
%   Output:
%       T:
%       F_weighted:
%       F_raw:


S.p_vec = p_vec;
S.enzyme_conc = enz_conc;
S.range = range;
S.init_cond = init_cond;

% Calculate kinetic parameters
P = Param_Function(S);

% Make ODEs with new params
parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);

% Intial Conditions
init_cond = S.init_cond;

%Numerical solution is found with solver ode15s (ideal for stiff systems)
%T: Time (sec)
%C: Matrix of concentration values of each species (columns) at each time
%point in T (rows)
[T,C] = ode15s(parameterized_ODEs,range,init_cond,ODE_options);

%Calculates the total concentration of each intermediate from the
%concentration matrix (for all chain lengths)

%F: Fatty Acid
%Q: Beta-ketoacyl-ACP
%M: Beta-hydroxyacyl-ACP
%R: Enoyl-Acyl-ACP
%T: Acyl-ACP
F_total = zeros(length(T),1);
Q_total = zeros(length(T),1);
M_total = zeros(length(T),1);
R_total = zeros(length(T),1);
T_total = zeros(length(T),1);
for ind = 1:length(S.labels)
        label_val = char(S.labels(ind));
        first_char = label_val(1);
        if first_char == char('F')
            F_total = C(:,ind) + F_total;
        end
        if first_char == char('Q')
            Q_total = C(:,ind) + Q_total;
        end
        if first_char == char('M')
            M_total = C(:,ind) + M_total;
        end
        if first_char == char('R')
            R_total = C(:,ind) + R_total;
        end
        if first_char == char('T')
            T_total = C(:,ind) + T_total;
        end
end

F_weighted = zeros(length(T),1);
F_saved = zeros(1,length(S.FA_dist));
F_raw = zeros(1,length(S.FA_dist));
weight_vec = S.FA_dist/16; %Palmitic Acid equivalents (C16)
count = 1;
for ind = 1:length(S.labels)
    label_val = char(S.labels(ind));
    first_char = label_val(1);
    if first_char == char('F')
        F_saved(count) = weight_vec(count)*(C(end,ind));
        F_raw(count) = C(end,ind);
        F_weighted = weight_vec(count)*(C(:,ind)) + F_weighted;
        count = count + 1;
    end
end
