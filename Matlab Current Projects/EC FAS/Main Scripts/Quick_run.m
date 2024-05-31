% Trying to see if parameters run for ACC

% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% Simulation

% Time
S.range = [0 150]; %2.5 mins (initial rate)

% Initial conditions
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 500; % Bicarbonate
S.init_cond(3) = 600; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(15) = 1300; % NADH
S.init_cond(18) = 0; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc = [1 1 1 1 1 1 10 1 1 1];

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw, rel_rate] = Calc_Function(T,C,S);

[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);
