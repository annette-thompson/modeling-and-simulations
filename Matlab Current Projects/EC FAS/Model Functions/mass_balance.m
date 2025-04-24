function [balance_conc, balances, total_conc, carbon] = mass_balance(C, P)
% Mass Balance
% Inputs:
%   C - Concentration matrix (rows: time points, cols: species)
%   P - Struct with species labels
% Outputs:
%   balance_conc - Species grouped by each balance category
%   balances - Net final-initial change for each balance
%   total_conc - Total concentration of each balance over time
%   carbon - Carbon species scaled by length over time

% Define balances to search for
balance_search = {'CoA', 'ACP', 'NAD(?!P)', 'NADP', 'ATP.*|.*ADP', 'C(\d+)_'};
n_balances = numel(balance_search);

% Initialize outputs
balances = zeros(1, n_balances);
total_conc = cell(1, n_balances);
balance_conc = cell(1, n_balances);
carbon = [];

for i = 1:n_balances
    % Find species matching the current balance
    matches = cellfun(@(label) regexp(label, ['.*' balance_search{i} '.*'], 'match'), P.labels, 'UniformOutput', false);
    balance_conc{i} = [matches{:}];
    
    % Handle carbon-specific balances
    if strcmp(balance_search{i}, 'C(\d+)_')
        carbon_length = cellfun(@(label) str2double(regexp(label, 'C(\d+)_', 'tokens', 'once')), balance_conc{i});
    else
        carbon_length = ones(1, numel(balance_conc{i}));
    end

    % Find indices of matched species in the concentration matrix
    conc_index = find(ismember(P.labels, balance_conc{i}));
    
    % Calculate the net change for each balance
    final_change = C(end, conc_index) - C(1, conc_index);
    balances(i) = sum(final_change .* carbon_length);

    % Calculate total concentration over time
    total_conc{i} = sum(C(:, conc_index) .* carbon_length, 2);

    % Store scaled carbon concentrations over time (only for carbon balance)
    if strcmp(balance_search{i}, 'C(\d+)_')
        carbon = C(:, conc_index) .* carbon_length;
    end
end
