function [total_AA,average_chain_length,En_BH_ratio,unsaturated_frac]=Intermediate_Solver_func(X)

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Time
%S.range = [0 720]; % 12 mins
%S.range = [0 150]; % 2.5 mins
S.range = [0 30]; % 30 sec

% Initial conditions
S.init_cond = zeros(S.num,1);
% S.init_cond(1) = 0; % ATP
% S.init_cond(2) = 0; % Bicarbonate
% S.init_cond(3) = 500; % Acetyl-CoA
% S.init_cond(12) = 10; % holo ACP
% S.init_cond(13) = 1000; % NADPH
% S.init_cond(15) = 1000; % NADH
% S.init_cond(18) = 500; % Malonyl-CoA
% 
% % (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
% S.enzyme_conc = [0, X];

S.init_cond(1) = X(1); % ATP
S.init_cond(2) = X(2); % Bicarbonate
S.init_cond(3) = X(3); % Acetyl-CoA
S.init_cond(12) = X(4); % holo ACP
S.init_cond(13) = X(5); % NADPH
S.init_cond(15) = X(6); % NADH
S.init_cond(18) = X(7); % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc = [1 1 1 1 1 1 0 1 1 1];

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
% parameterized_ODEs = @(t,c) ODE_Function(t,c,P);

[~,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);

% Find average chain length
pattern = '^c_C\d{1,2}_(BHyAcACP|EnAcACP|AcACP)';
indices = find(cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels));
labels = S.labels(indices);
F_raw = C(end,indices);

chain_lengths = zeros(size(F_raw));
conc = zeros(size(F_raw));
for i = 1:length(labels)
    chain_length = regexp(labels(i), '(?<=c_C)\d+', 'match', 'once');
    if ~isempty(chain_length)
        chain_lengths(i) = str2double(chain_length);
    else
        chain_lengths(i) = NaN;
    end
    conc(i) = F_raw(i); 
end

chain_lengths(isnan(chain_lengths)) = 0;

average_chain_length = sum(chain_lengths .* conc) / sum(conc);

% Find total AcACPs
pattern = '^c_C\d{1,2}_AcACP';
indices = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels);
F_raw = C(end,indices);
total_AA = sum(F_raw);

% Find EnAcACP:BHyAcACP ratio
pattern = '^c_C\d{1,2}_EnAcACP';
indices = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels);
F_raw = C(end,indices);
total_En = sum(F_raw);

pattern = '^c_C\d{1,2}_BHyAcACP';
indices = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels);
F_raw = C(end,indices);
total_BH = sum(F_raw);
En_BH_ratio = total_En/total_BH;

% Find unsaturated fraction
pattern = '^c_C\d{1,2}_(AcACP|EnAcACP|BHyAcACP)$';
indices = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels);
F_raw = C(end,indices);
total_sat = sum(F_raw);

pattern = '^c_C\d{1,2}_(AcACP|EnAcACP|BHyAcACP)_un$';
indices = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels);
F_raw = C(end,indices);
total_unsat = sum(F_raw);

unsaturated_frac = total_unsat/total_sat;