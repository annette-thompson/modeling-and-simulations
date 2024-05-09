%% Optimization_Function
function [p_vec_sol] = Optimization_Function(S,p_opt,opt_options,ODE_options)
% Calculates objectives for optimizing model output with experimental data
%   Input:
%       S: structure containing all kinetic parameters
%       p_opt: contains positions of p_vec to be optimized
%       opt_options: contains options for fminsearch
%       ODE_options: contains options for ode15s
%   Output:
%       p_vec_sol: solution of optimized p_vec parameters

p_vec0 = S.p_vec;

% See if adding any new terms (lengthen vector)
if any(p_opt>length(p_vec0))
    sz = max(p_opt);
else
    sz = length(p_vec0);
end

% Positions to change in the vector
new = cell(1,sz);
for i = 1:sz
    if any(i==p_opt)
        new{i} = sprintf('x(%d)', i);
    else 
        new{i} = num2str(p_vec0(i));
    end
end

% Compose new vector with x variables as values to be optimized
new = strjoin(new,' ');
new = @(x) str2func(new);

fitfunc = @(x) Optimization_Function_Handler(S,new,ODE_options);

[p_vec_sol,~,~,~] = fminsearch(fitfunc,p_vec0,opt_options);
