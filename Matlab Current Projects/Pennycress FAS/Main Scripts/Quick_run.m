% Testing pennycress model (no ACC)

% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/Pennycress FAS';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% Simulation

% Time
S.range = [0 7200]; % 2 hrs (total production)

% Initial conditions
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 0; % ATP
S.init_cond(2) = 0; % Bicarbonate
S.init_cond(3) = 300; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 2600; % NADPH
S.init_cond(15) = 2600; % NADH
S.init_cond(17) = 2600; % Fd+
S.init_cond(20) = 1500; % Malonyl-CoA

% (ACC,MCMT,KASIII,KAR,HAD,ER,FatA,KASI,SAD,KASII)
S.enzyme_conc = [0 1 1 1 1 1 10 1 1 1];

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);

tic
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc

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
xticklabels([' 4  ';' 6  ';' 8  ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1';])
ylabel('Production (uM)')
title("2 hrs, No ACC")
xlabel('Chain Length')
ylim([0 90])

% Ac/MalCoA
% figure()
% plot(T/60,C(:,18))
% hold on
% plot(T/60,C(:,3))
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% legend('Malonyl-CoA','Acetyl-CoA')
% title("With ACC 2 hrs")