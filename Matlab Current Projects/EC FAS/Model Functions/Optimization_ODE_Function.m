%% Optimization_ODE_Function
function [T,F_weighted,F_raw_new] = Optimization_ODE_Function(S,enz_conc,range,init_cond,ODE_options)
% Runs the ODE_Function for optimization and calculates final FA production
% for optimization comparison
%   Input:
%       S: structure containing all kinetic parameters
%       enz_conc: enzyme composition
%       range: time to run
%       init_cond: initial species conditions
%       ODE_options: options for ODE solver
%   Output:
%       T: time points
%       F_weighted: total palmitic acid equivalents at each time point
%       F_raw_new: product profile concentrations

% Set run conditions
S.enzyme_conc = enz_conc;
S.range = range;
S.init_cond = init_cond;

% Calculate kinetic parameters
P = Param_Function(S);

% Make ODEs with new params
parameterized_ODEs = @(t,c) ODE_Function(t,c,P);

% Numerical solution is found with solver ode15s (ideal for stiff systems)
% T: Time (sec)
% C: Matrix of concentration values of each species (columns) at each time
% point in T (rows)
[T,C] = ode15s(parameterized_ODEs,range,init_cond,ODE_options);

% Calculates the total concentration of each FA
F_weighted = zeros(length(T),1);
F_saved = zeros(1,length(S.FA_dist));
F_raw = zeros(1,length(S.FA_dist));
weight_vec = S.FA_dist/16; % Palmitic Acid equivalents (C16)
count = 1;
for ind = 1:length(S.labels)
    label_val = char(S.labels(ind));
    if contains(label_val, '_FA','IgnoreCase',false)
        F_saved(count) = weight_vec(count)*(C(end,ind));
        F_raw(count) = C(end,ind);
        F_weighted = weight_vec(count)*(C(:,ind)) + F_weighted;
        count = count + 1;
    end
end

% Rearrange to match experimental data
F_raw_new(1,1:4) = F_raw(1,1:4);
j=10;
for i=5:2:13
    F_raw_new(1,i) = F_raw(j-5);
    F_raw_new(1,i+1) = F_raw(j);
    j = j+1;
end

end
