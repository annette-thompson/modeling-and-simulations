%% Calc_Function
function [F_raw, rel_rate] = Calc_Function(T,C,S)

F_weighted = zeros(length(T),1);
F_saved = zeros(1,length(S.FA_dist));
F_raw = zeros(1,length(S.FA_dist));
weight_vec = S.FA_dist/16; %Palmitic Acid equivalents (C16)
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

rel_rate = F_weighted(end,1)/150*60;

end