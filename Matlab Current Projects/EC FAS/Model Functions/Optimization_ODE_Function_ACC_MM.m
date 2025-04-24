%% Optimization_ODE_Function
function [F_raw] = Optimization_ODE_Function_ACC_MM(S,enz_conc,range,init_cond,ODE_options)
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
parameterized_ODEs = @(t,c) ODE_Function_ACC_MM(t,c,P);

% Numerical solution is found with solver ode15s (ideal for stiff systems)
% T: Time (sec)
% C: Matrix of concentration values of each species (columns) at each time
% point in T (rows)
[~,C] = ode15s(parameterized_ODEs,range,init_cond,ODE_options);

% Find indices matching the patterns
pattern = '^c_C\d{1,2}_(BHyAcACP|EnAcACP|AcACP)';
indices = find(cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels));

% Extract labels and chain lengths
labels = S.labels(indices);
chain_lengths = cellfun(@(x) str2double(regexp(x, '(?<=c_C)\d{1,2}(?=_(BHyAcACP|EnAcACP|AcACP))', 'match', 'once')), labels, 'UniformOutput', false);
chain_lengths = cell2mat(chain_lengths);
chain_lengths(isnan(chain_lengths)) = 0;

% Sort labels and indices based on chain lengths
[~, sortedOrder] = sort(chain_lengths);
sorted_indices = indices(sortedOrder);
plottedSpecies = [2,3,5,6,8,9,11,12,14,16,19,15,20,22,25,21,28];
plot_indices = sorted_indices(plottedSpecies);
F_raw = C(end,plot_indices); % Model data
total_F_raw = sum(F_raw);
F_raw = F_raw/total_F_raw;