%% Calc_Function
function [F_raw, rel_rate] = Calc_Function(T,C,S)
% Calculates production profile and rate
%   Input:
%       T: time points [sec]
%       C: concentration matrix [uM]
%       S: structure storing all variables
%   Output:
%       F_raw: concentrations of all free fatty acids [uM]
%       rel_rate: rate of fatty acid production [uM/min]
%       F_weighted (not currently outputed): palmitic acid equivalents at
%       each time point [uM]

% Initialize vectors
F_weighted = zeros(length(T),1);
F_saved = zeros(1,length(S.FA_dist));
F_raw = zeros(1,length(S.FA_dist));

% Find final concentration of fatty acids and palmitic acid equivalents at each time point
weight_vec = S.FA_dist/16; % Palmitic acid equivalents (C16)
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

% Calculate rate
rel_rate = F_weighted(end,1)/T(end)*60;

end