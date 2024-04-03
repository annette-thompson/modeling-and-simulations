%% Calc_Function
function [F_raw, rel_rate] = Calc_Function(T,C,S)

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

rel_rate = F_weighted(end,1)/150*60;

end