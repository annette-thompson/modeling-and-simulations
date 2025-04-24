% Find kcats for FabH to fit Katie's data

my_dir = '/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
addpath(genpath(my_dir));
%%
x0 = 3.6464;

S = set_vars();
S.range = [0 150]; % 2.5 min

S.init_cond = zeros(S.num, 1);
starting_substrates = [3, 6, 12, 13, 15, 18]; % AcCoA, OcCoA, ACP, NADPH, NADH, MalCoA
% starting_substrate_concentrations = [100, 0, 10, 1000, 1000, 500]; % uM
starting_substrate_concentrations = [0, 100, 10, 1000, 1000, 500]; % uM
S.init_cond(starting_substrates) = starting_substrate_concentrations;

S.enzyme_conc = [0 1 1 1 1 1 10 1 1 1]; % ACC, FabD, FabH, FabG, FabZ, FabI, TesA, FabF, FabA, FabB

% S.kon_scaling_fabH = [1,0,0,0,0,0,0,0,0];
S.kon_scaling_fabH = [0,0,0,1,0,0,0,0,0];

fitfunc = @(kcat_scaling) Fit_kcat(kcat_scaling, S);

opt_options = optimset('Display','iter','OutputFcn',@optout);

[PP_FabH2_B_w_init_kcat, PP_FabH2_B_w_init_obj] = fminsearch(fitfunc, x0,opt_options);

function obj = Fit_kcat(kcat_scaling, S)
% S.kcat_scaling_fabH = [kcat_scaling,0,0,0,0,0,0,0,0];
S.kcat_scaling_fabH = [0,0,0,kcat_scaling,0,0,0,0,0];
P = Param_Function(S);
parameterized_ODEs = @(t, c) ODE_Function(t, c, P);

ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

[T, C] = ode15s(parameterized_ODEs, S.range, S.init_cond, ODE_options);

% FA_indices = contains(S.labels, '_FA', 'IgnoreCase', false);
% FA_raw = C(end, FA_indices);
% FA_weighted = sum(S.FA_dist / 16 .* FA_raw);
% rate = FA_weighted / T(end) * 60;

NAD_converted = sum(C(end, strcmp(S.labels, 'NADP') | strcmp(S.labels, 'NAD')));
rate = NAD_converted / T(end) * 60;

% Katie's data
% exp_rate = 8.857864358; % EC FabH, A
% exp_rate = 4.886394558; % EC FabH, D
% exp_rate = 6.65303832; % PP FabH1, A
% exp_rate = 4.757816258; % PP FabH2, B
% exp_rate = 4.757816258*14/8; % PP FabH2, B fixed

% Katie's NADPH data
% exp_rate = 124.010101; % EC FabH, A
% exp_rate = 121.6416867; % EC FabH, D
% exp_rate = 93.14253648; % PP FabH1, A
exp_rate = 66.60942761; % PP FabH2, B

obj = (exp_rate - rate)^2;
end
%%
x0 = 15;

S = set_vars();
S.range = [0 150]; % 2.5 min

S.init_cond = zeros(S.num, 1);
starting_substrates = [3, 6, 12, 13, 15, 18]; % AcCoA, OcCoA, ACP, NADPH, NADH, MalCoA
starting_substrate_concentrations = [100, 0, 10, 1000, 1000, 500]; % uM
S.init_cond(starting_substrates) = starting_substrate_concentrations;

S.enzyme_conc = [0 1 1 1 1 1 10 1 1 1]; % ACC, FabD, FabH, FabG, FabZ, FabI, TesA, FabF, FabA, FabB

S.kon_scaling_fabH = [1,0,0,0,0,0,0,0,0];
S.kcat_scaling_fabH = [1,0,0,0,0,0,0,0,0];

fitfunc = @(FabZ_conc) Fit_FabZ(FabZ_conc, S);

opt_options = optimset('Display','iter','OutputFcn',@optout);

[FabZ, obj] = fminsearch(fitfunc, x0,opt_options);

function obj = Fit_FabZ(FabZ_conc, S)
S.enzyme_conc(5) = FabZ_conc;

P = Param_Function(S);
parameterized_ODEs = @(t, c) ODE_Function(t, c, P);

ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

[T, C] = ode15s(parameterized_ODEs, S.range, S.init_cond, ODE_options);

FA_indices = contains(S.labels, '_FA', 'IgnoreCase', false);
FA_raw = C(end, FA_indices);
FA_weighted = sum(S.FA_dist / 16 .* FA_raw);
rate = FA_weighted / T(end) * 60;

% Katie's data
exp_rate = 8.857864358; % EC FabH, A
% exp_rate = 4.886394558; % EC FabH, D
% exp_rate = 6.65303832; % PP FabH1, A
% exp_rate = 4.757816258; % PP FabH2, B
% exp_rate = 4.757816258*14/8; % PP FabH2, B fixed?

obj = (exp_rate - rate)^2;
end

%%
x0 = [1, 1];

S = set_vars();
S.range = [0 150]; % 2.5 min

S.init_cond = zeros(S.num, 1);
starting_substrates = [3, 6, 12, 13, 15, 18]; % AcCoA, OcCoA, ACP, NADPH, NADH, MalCoA
starting_substrate_concentrations = [100, 0, 10, 1000, 1000, 500]; % uM
S.init_cond(starting_substrates) = starting_substrate_concentrations;

S.enzyme_conc = [0 1 1 1 1 1 10 1 1 1]; % ACC, FabD, FabH, FabG, FabZ, FabI, TesA, FabF, FabA, FabB

S.kon_scaling_fabH = [1,0,0,0,0,0,0,0,0];
S.kcat_scaling_fabH = [1,0,0,0,0,0,0,0,0];

fitfunc = @(x) Fit_Both(x, S);

opt_options = optimset('Display','iter','OutputFcn',@optout);

[sol, obj] = fminsearch(fitfunc, x0,opt_options);

function obj = Fit_Both(x, S)
S.kcat_scaling_fabH = [x(1),0,0,0,0,0,0,0,0];
S.enzyme_conc(5) = x(2);

P = Param_Function(S);
parameterized_ODEs = @(t, c) ODE_Function(t, c, P);

ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

[T, C] = ode15s(parameterized_ODEs, S.range, S.init_cond, ODE_options);

FA_indices = contains(S.labels, '_FA', 'IgnoreCase', false);
FA_raw = C(end, FA_indices);
FA_weighted = sum(S.FA_dist / 16 .* FA_raw);
rate = FA_weighted / T(end) * 60;

% Katie's data
exp_rate = 8.857864358; % EC FabH, A
% exp_rate = 4.886394558; % EC FabH, D
% exp_rate = 6.65303832; % PP FabH1, A
% exp_rate = 4.757816258; % PP FabH2, B
% exp_rate = 4.757816258*14/8; % PP FabH2, B fixed?

obj = (exp_rate - rate)^2;
end