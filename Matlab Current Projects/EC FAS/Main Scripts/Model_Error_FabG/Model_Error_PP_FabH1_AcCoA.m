% PP FabH1 Acetyl CoA

% Set mean and standard deviation
% mu= mean, sigma = stdev
% AcCoA OcCoA ACP NADPH NADH MalCoA FabD FabH FabG FabI TesA FabF FabA FabB
mu=[100 0 10 1000 1000 500  1    1    1    1    1    10   1    1    1];
sigma = 0.05.*mu;

% Get matrix of concentrations for given number of runs
samples=300;
r=zeros(samples,length(mu));
for j=1:length(mu)
    r(:,j) = normrnd(mu(j),sigma(j),[samples,1]);
end

% Run variable code
S = set_vars();
S.range = [0 150]; % 2.5 mins (initial rate)

% Scaling for FabH
EC_kcat3_scaling = [1,0,0,0,0,0,0,0,0];
PP_H1_kcat3_scaling = [5.5468,0,0,0,0,0,0,0,0]; 
PP_H2_kcat3_scaling = [0,0,0,3.7729,0,0,0,0,0]; 

EC_FabH_kon_scaling = [1,0,0,0,0,0,0,0,0];
PP_FabH1_kon_scaling = [1,0,0,0,0,0,0,0,0];
PP_FabH2_kon_scaling = [0,0,0,1,0,0,0,0,0];

% Using E. coli FabH
S.kcat_scaling_fabH = PP_H1_kcat3_scaling;  
S.kon_scaling_fabH = PP_FabH1_kon_scaling; 

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Run model at each set of concetrations
rel_rate = zeros(1,samples);
NAD_rate = zeros(1,samples);
for i=1:samples
    % Set initial conditions
    S.init_cond = zeros(S.num,1);
    S.init_cond(3) = r(i,1); % AcCoA
    S.init_cond(6) = r(i,2); % OcCoA
    S.init_cond(12) = r(i,3); % ACP
    S.init_cond(13) = r(i,4); % NADPH
    S.init_cond(15) = r(i,5); % NADH
    S.init_cond(18) = r(i,6); % MalCoA
    S.enzyme_conc =[0 r(i,7:end)]; % ACC FabD FabH FabG FabI TesA FabF FabA FabB
    
    % Calculate parameters
    P = Param_Function(S);

    % Solve system of ODEs
    parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
    [T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);

    % Calculate initial rate
    [~, rel_rate(i), NAD_rate(i)] = Calc_Function(T,C,S);
end

% Calculate mean and stdev
mean_rel_rate = mean(rel_rate);
std_rel_rate = std(rel_rate);
mean_NAD_rate = mean(NAD_rate);
std_NAD_rate = std(NAD_rate);

% Save results to CSV
output_data = table(r(:,1), r(:,2), r(:,3), r(:,4), r(:,5), r(:,6), r(:,7),...
    r(:,8), r(:,9), r(:,10), r(:,11), r(:,12), r(:,13), r(:,14), r(:,15), rel_rate', NAD_rate', ...
    'VariableNames', {'AcCoA', 'OcCoA', 'ACP', 'NADPH', 'NADH', 'MalCoA', 'FabD',...
    'FabH', 'FabG', 'FabZ', 'FabI', 'TesA', 'FabF', 'FabA', 'FabB', 'RelRate', 'NADRate'});

summary_data = table(mean_rel_rate, std_rel_rate, mean_NAD_rate, std_NAD_rate,...
    'VariableNames', {'MeanRelRate', 'StDevRelRate', 'MeanNADRate', 'StDevNADRate'});

writetable(output_data, 'Model_error_output_PP_FabH1_Acetyl_CoA.csv');
writetable(summary_data, 'Model_error_summary_PP_FabH1_Acetyl_CoA.csv');