% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% Make testing matrix
% Testing starting substrate concentration
bicarbonate = linspace(100, 1000, 10);
acetylcoa = linspace(100, 1000, 10);

% Experimental data (guessed from Allen plots)
exp_pmol_data = [160 78 5 6 25 12 12 7 3 3 1 0.5];
exp_peak_data = [1.5e4 1.4e4 1.6e4 4e4 6e4 1.2e5 1.4e4 2.5e4 1.4e5 2e4 0.5e4 0.5e4 2e5 3e4 0.05e4];
pmol_total = sum(exp_pmol_data);
peak_total = sum(exp_peak_data);
pmol_dist = exp_pmol_data / pmol_total;
peak_dist = exp_peak_data / peak_total;
exp_dist = [pmol_dist peak_dist];
exp_indices = [65 66 67 68 69 70 71 76 72 77 36 50 37 51 38 52 39 53 40 54 74 45 41 55 75 46 56];
%'c_C4_AcACP','c_C6_AcACP','c_C8_AcACP','c_C10_AcACP','c_C12_AcACP','c_C14_AcACP','c_C16_AcACP','c_C16_AcACP_un','c_C18_AcACP','c_C18_AcACP_un'...
% 'c_C6_BHyAcACP','c_C6_EnAcACP','c_C8_BHyAcACP','c_C8_EnAcACP','c_C10_BHyAcACP','c_C10_EnAcACP','c_C12_BHyAcACP','c_C12_EnAcACP','c_C12_AcACP_un','c_C12_BHyAcACP_un','c_C14_BHyAcACP','c_C14_EnAcACP','c_C14_AcACP_un','c_C14_BHyAcACP_un','c_C16_EnAcACP'

% Optimization in vivo data - probably should switch to Katie's
% (Grisewood, M.; Hernández-Lozada, N.; Thoden, J.; Gifford, N.; Mendez-Perez, 
% D.;Lai, R.; Holden, H.; Pfleger, Et. al. ACS Catalysis 2017, 7, 3837-3849.)
vivo_dist = [0 0 0.059771046 0.011042448 0.202613717 0.075980144 0.29378989 0.07192787 0.143855739 0.117515956 0.002228751 0.02127444 0 0];
vivo_indices = [79 80 81 82 83 88 84 89 85 90 86 91 87 92];

[Y1, Y2] = ndgrid(bicarbonate, acetylcoa);
Y = [Y1(:), Y2(:)];
n_sim = size(Y, 1);

results = zeros(n_sim, 8);

%% Run all tests
parfor n=1:n_sim

    % Allen data (no TesA)
    S = set_vars();
    S.init_cond = zeros(S.num, 1);
    S.init_cond(1) = 1000; % ATP
    S.init_cond(12) = 10; % holo ACP
    S.init_cond(13) = 1000; % NADPH
    S.init_cond(15) = 1000; % NADH
    S.init_cond(18) = 0; % Malonyl-CoA
    S.enzyme_conc = [10 1 1 1 1 1 0 1 1 1];
    S.range = [0 720];
    
    S.init_cond(2) = Y(n, 1); % Bicarbonate
    S.init_cond(3) = Y(n, 2); % Acetyl-CoA

    P = Param_Function(S);
    
    % Turn off FabB/F Initiation
    P.kcat8_H = 0;
    P.kcat10_H = 0;
    
    % Set Michaelis Menton parameters for ACC
    P.kcat1_1 = 85.17/60; % s^-1
    P.Km1_1 = 170; % uM
    P.kcat1_2 = 73.8/60; % s^-1
    P.Km1_2 = 370; % uM
    P.kcat1_3 = 1000.8/60*S.scaling_factor_kcat_init; % s^-1
    P.Km1_3 = 160; % uM
    P.kcat1_4 = 2031.8/60*S.scaling_factor_kcat_init; % s^-1
    P.Km1_4 = 450; % uM
    P.kcat1_5 = 30.1; % s^-1
    P.Km1_5 = 48.7; % uM

    parameterized_ODEs = @(t, c) ODE_Function_ACC_MM(t, c, P);
    [T1, C1] = ode15s(parameterized_ODEs, S.range, S.init_cond, ODE_options);

    F_raw_1 = C1(:, exp_indices);
    total_F_raw_1 = sum(F_raw_1, 2);
    F_raw_dist_1 = F_raw_1 ./ total_F_raw_1;
    obj1 = sum((F_raw_dist_1 - exp_dist).^2, 2);

    % In vivo data (with TesA)
    S.enzyme_conc(7) = 10;

    P = Param_Function(S);

    % Turn off FabB/F Initiation
    P.kcat8_H = 0;
    P.kcat10_H = 0;
    
    % Set Michaelis Menton parameters for ACC
    P.kcat1_1 = 85.17/60; % s^-1
    P.Km1_1 = 170; % uM
    P.kcat1_2 = 73.8/60; % s^-1
    P.Km1_2 = 370; % uM
    P.kcat1_3 = 1000.8/60*S.scaling_factor_kcat_init; % s^-1
    P.Km1_3 = 160; % uM
    P.kcat1_4 = 2031.8/60*S.scaling_factor_kcat_init; % s^-1
    P.Km1_4 = 450; % uM
    P.kcat1_5 = 30.1; % s^-1
    P.Km1_5 = 48.7; % uM

    parameterized_ODEs = @(t, c) ODE_Function_ACC_MM(t, c, P);
    [T2, C2] = ode15s(parameterized_ODEs, S.range, S.init_cond, ODE_options);

    F_raw = C2(:, vivo_indices);
    total_F_raw = sum(F_raw, 2);
    F_raw_dist = F_raw ./ total_F_raw;
    obj2 = sum((F_raw_dist - vivo_dist).^2, 2);

    % Find common time points and interpolate
    T_common = linspace(S.range(1), S.range(2), 500);
    obj1_interp = interp1(T1, obj1, T_common, 'linear');
    obj2_interp = interp1(T2, obj2, T_common, 'linear');

    % Compute product of objectives and find minimum
    product_obj = obj2_interp .* obj1_interp;
    min_product = min(product_obj);
    min_indice = find(product_obj == min_product, 1, 'last');
    time_at_min = T_common(min_indice);
    obj1_min = obj1_interp(min_indice);
    count1 = sum(obj1_interp == obj1_min);
    obj2_min = obj2_interp(min_indice);
    count2 = sum(obj2_interp == obj2_min);

    results(n, :) = [Y(n,:) time_at_min min_product obj1_min obj2_min count1 count2];

    fprintf('%.4g\t', results(n, :))
    fprintf('\n')
end

[val,loc]=min(results(:,4));
disp(results(loc,1:6))
results_both=results;

%% Run based on results

% Run variable code
S = set_vars();

% Time
S.range = [0 180];

% Initial conditions
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 500; % Bicarbonate
S.init_cond(3) = 700; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(18) = 0; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc = [10 1 1 1 1 1 10 1 1 1];

P = Param_Function(S);

% Turn off FabB/F Initiation
P.kcat8_H = 0;
P.kcat10_H = 0;

% Set Michaelis Menton parameters for ACC
P.kcat1_1 = 85.17/60; % s^-1
P.Km1_1 = 170; % uM
P.kcat1_2 = 73.8/60; % s^-1
P.Km1_2 = 370; % uM
P.kcat1_3 = 1000.8/60*S.scaling_factor_kcat_init; % s^-1
P.Km1_3 = 160; % uM
P.kcat1_4 = 2031.8/60*S.scaling_factor_kcat_init; % s^-1
P.Km1_4 = 450; % uM
P.kcat1_5 = 30.1; % s^-1
P.Km1_5 = 48.7; % uM

parameterized_ODEs = @(t,c) ODE_Function_ACC_MM(t,c,P);
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

F_raw = C(end, vivo_indices);
total_F_raw = sum(F_raw);
F_raw_dist = F_raw/total_F_raw;

% In vivo opt data plotting
figure()
b = bar([F_raw_dist;vivo_dist]','grouped');
hold on

vivo_labels = {'4','6','8','10','12','12:1','14',...
    '14:1','16','16:1','18','18:1','20','20:1'};

ylim([0 max(max(F_raw_dist,vivo_dist))*1.1]);
ylabel('Mole %');
set(gca, 'XTickLabel', vivo_labels, 'XTick', 1:length(vivo_indices));
xlabel('acyl-ACPs')
legend('Model','Experiment')

%% Allen lab plotting
% Panel A: pmol/Total OD600
% Find indices matching the patterns
pattern = '^c_C\d{1,2}_AcACP|^c_ACP$';
indices = find(cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels));

% Extract labels and chain lengths
labels = S.labels(indices);
chain_lengths = cellfun(@(x) str2double(regexp(x, '(?<=c_C)\d{1,2}(?=_AcACP)', 'match', 'once')), labels, 'UniformOutput', false);
chain_lengths = cell2mat(chain_lengths);
chain_lengths(isnan(chain_lengths)) = 0;

% Sort labels and indices based on chain lengths
[~, sortedOrder] = sort(chain_lengths);
sorted_indices = indices(sortedOrder);
plottedSpecies = [1:7,9,11:14];
plot_indices = sorted_indices(plottedSpecies);
acyl_ACPs_A = {'Holo', 'C2', 'C4', 'C6', 'C8', 'C10', 'C12', 'C14', 'C16', 'C16:1', 'C18', 'C18:1'};
model_A = C(end,plot_indices); % Model data
total_model_A = sum(model_A);
model_A = model_A/total_model_A;
pmol_37C = [160, 78, 5, 6, 25, 12, 12, 7, 3, 3, 1, 0.5];  % 37°C data
total_pmol_37C = sum(pmol_37C);
pmol_37C = pmol_37C/total_pmol_37C;

figure;
subplot(9,1,1:3);
hold on;
b_A = bar([model_A; pmol_37C]', 'grouped');

% Set colors for 30°C and 37°C
b_A(1).FaceColor = [0 0 0];  % Black for model
b_A(2).FaceColor = [0.9, 0.8, 0.6];  % Gold for 37°C

ylim([0 max(max(model_A,pmol_37C))*1.1]);
ylabel('Mole %');
set(gca, 'XTickLabel', acyl_ACPs_A, 'XTick', 1:length(acyl_ACPs_A));
legend('Model', 'Experiment');
set(gca,'FontSize', 18)

hold off;

% Panel B: Peak area/Total OD600
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
sorted_labels = labels(sortedOrder);
sorted_indices = indices(sortedOrder);
plottedSpecies = [2,3,5,6,8,9,11,12,14,16,19,15,20,22,25,21,28];
plot_indices = sorted_indices(plottedSpecies);

acyl_ACPs_B = {'βhyd-C4', 'C4-enoyl', 'βhyd-C6', 'C6-enoyl', 'βhyd-C8', ...
               'C8-enoyl', 'βhyd-C10', 'C10-enoyl', 'C10:1', 'βhyd-C12', 'C12-enoyl', ...
               'C12:1', 'βhyd-C12:1', 'βhyd-C14', 'C14-enoyl', 'C14:1', 'βhyd-C14:1', 'C16-enoyl'};
model_B = [C(end,plot_indices(1:8)),0,C(end,plot_indices(9:end))]; % Model data
total_model_B = sum(model_B);
model_B = model_B/total_model_B;
peak_37C = [3e5, 1.4e4, 1.5e4, 1.4e4, 1.6e4, 4e4, 6e4, 1.2e5, 0e4, 1.4e4, 2.5e4, 1.4e5, 2e4, 0.5e4, 0.5e4, 2e5, 3e4, 0.05e4];  % 37°C data
total_peak_37C = sum(peak_37C);
peak_37C = peak_37C/total_peak_37C;

subplot(9,1,5:7);
hold on;
b_B = bar([model_B; peak_37C]', 'grouped');

% Set colors for 30°C and 37°C
b_B(1).FaceColor = [0 0 0];  % Black for model
b_B(2).FaceColor = [0.9, 0.8, 0.6];  % Gold for 37°C

ylim([0 max(max(model_B,peak_37C))*1.1]);
ylabel('Mole %');
set(gca, 'XTickLabel', acyl_ACPs_B, 'XTick', 1:length(acyl_ACPs_B));
set(gca,'FontSize', 18)
xlabel('acyl-ACPs')

hold off;