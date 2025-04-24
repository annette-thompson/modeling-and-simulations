% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/Pennycress FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% Simulation

% Time
S.range = [0 720]; % 12 mins

% Initial conditions
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 0; % ATP
S.init_cond(2) = 0; % Bicarbonate
S.init_cond(3) = 500; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(17) = 1000; % Fd+
S.init_cond(20) = 500; % Malonyl-CoA

% (ACC,MCMT,KASIII,KAR,HAD,ER,FatA,KASI,SAD,KASII)
S.enzyme_conc = [0 1 1 1 1 1 10 1 0 1];

% S.kcat_scaling_KASI = [0.914,0.914,0.901,1,0.9,0.289,0.222,0.222,0.222]; % FabF scaling
% S.kcat_scaling_KASII = [0.855,0.855,0.975,0.967,1,0.125,0.0208,0.0208,0.0208]; % FabB scaling
S.kcat_scaling_KASI = [0.914,0.914,0.901,1,0.9,0.289,0.222,0,0];
S.kcat_scaling_KASII = [0.855,0.855,0.975,0.967,1,0.125,0.0208,0,0];

% S.FatA_fitting_source = 'FatA'; % kcat_scaling_FatA controls kon and kcat instead of Pf
% To get to Pf kcat, scale to 0.0569    0.0509    0.1035    0.0158    0.2526    0.4582    1.0000    1.2210    1.5369
% To get to Pf kon, scale to 0.0172    0.0280    0.1573    0.5625    1.0830    1.9170    3.3955    6.0119   10.6455
% S.kcat_scaling_FatA = [0.1, 0.1, 0.1, 10, 0.01, 0.01, 0.01, 0.01, 0.01]; % "Cuphea" TE

P = Param_Function(S);

% Solve ODEs
parameterized_ODEs = @(t,c) ODE_Function_ACC_MM(t,c,P);
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[F_raw, rel_rate] = Calc_Function(T,C,S);
[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

%% Plotting
% Profile
figure()
F_raw_new(1,1:4) = F_raw(1,1:4);
j=10;
for i=5:2:13
    F_raw_new(1,i) = F_raw(j-5);
    F_raw_new(1,i+1) = F_raw(j);
    j = j+1;
end
total = sum(F_raw_new);
bar(F_raw_new)
xticklabels(['  4 ';'  6 ';'  8 ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1';])
ylabel('Production (uM)')
xlabel('Chain Length')
ylim([0 ceil(max(F_raw_new))])
set(gca,'FontSize', 18)
