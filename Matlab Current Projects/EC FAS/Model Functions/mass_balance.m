% Mass Balance
function [balance_conc, balances, total_conc, carbon] = mass_balance(C,P)
% C: concentration matrix (rows = time points, cols = species)
% P: struct containing list of lables
% balance_conc: lists of species pulled for each balance
% balances: overall final-inital for each balance
% total_conc: total concentration of balance vs time
% carbon: each carbon species (scaled by their length) versus time

% What to search for for each balance
balance_search = {'CoA','ACP','NAD(?!P)','NADP','ATP.*|.*_ADP','_C(\d+)_'};

% Where final difference between beginning and end will be stored
balances = zeros(1,numel(balance_search));

% Where total concentration over time will be stored
total_conc = cell(1,numel(balance_search));

for i=1:numel(balance_search)

    % Find all species that are included in the balance
    balance_conc{i} = {};
    for j = 1:numel(P.labels)
        matches = regexp(P.labels{j}, ['.*' balance_search{i} '.*'], 'match');
        if ~isempty(matches)
            balance_conc{i} = [balance_conc{i}, matches];
        end
    end

    % If doing carbon balance, find length of carbon for each species
    if strcmp(balance_search{i},'_C(\d+)_')
        carbon_length = [];
        for k = 1:numel(balance_conc{i})
            % Find length from name
            match = regexp(balance_conc{i}{k}, 'c_C(\d+)_', 'tokens');
            number = str2double(match{1}{1});
            carbon_length = [carbon_length, number];
        end
    end

    % Find what number in the concentration matrix corresponds to each
    % species
    conc_index = zeros(1,numel(balance_conc{i}));
    for j=1:numel(balance_conc{i})
        conc_index(j) = find(ismember(P.labels, balance_conc{i}{j}), 1);
    end

    % Calculate difference between end and start (should be as close to 0
    % as possible)
    balances(i) = 0;
    for j = 1:numel(conc_index)
            final_change = C(end,conc_index(j)) - C(1,conc_index(j));
            if strcmp(balance_search{i},'_C(\d+)_')
                balances(i) = balances(i) + final_change*carbon_length(j);
            else
                balances(i) = balances(i) + final_change;
            end
    end

    % Calculate concentration at each time point (should stay constant)
    for j=1:numel(C(:,1))
            sum=0;
            for k=1:numel(conc_index)
                if strcmp(balance_search{i},'_C(\d+)_')
                    sum = sum +C(j,conc_index(k))*carbon_length(k);
                    carbon(j,k) = C(j,conc_index(k))*carbon_length(k);
                else
                    sum = sum +C(j,conc_index(k));
                end
            end
            total_conc{i}(j)=sum;
    end
end
