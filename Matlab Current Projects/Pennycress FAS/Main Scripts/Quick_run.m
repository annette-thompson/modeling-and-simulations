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
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 1000; % Bicarbonate
S.init_cond(3) = 500; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(17) = 1000; % Fd+
S.init_cond(20) = 500; % Malonyl-CoA

% (ACC,MCMT,KASIII,KAR,HAD,ER,FatA,KASI,SAD,KASII)
S.enzyme_conc = [1 1 1 1 1 1 10 1 1 1];

P = Param_Function(S);

% Turn off KASI/II Initiation
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

% Solve ODEs
parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[F_raw_12, rel_rate] = Calc_Function(T,C,S);
[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

%% Plotting
% Profile
figure()
F_raw_new(1,1:4) = F_raw_12(1,1:4);
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
set(gca,'FontSize', 24)
