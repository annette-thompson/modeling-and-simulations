% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

% Variables

% Run variable code
S = set_vars();

% Set FabG/FabH scaling
EC_kcat3_scaling = [1,0,0,0,0,0,0,0,0];
PP_H1_kcat3_scaling = [0.5,0,0,0,0,0,0,0,0]; % changed
PP_H2_kcat3_scaling = [0,0,0,0.4,0,0,0,0,0]; % changed

EC_kcat4_scaling = [1,1,1,1,1,1,1,1,1];
PP_1914_kcat4_scaling = [.1,.1,.1,.1,.1,.1,.1,.1,.1]; % changed
PP_2783_kcat4_scaling = [0,0,0,.025,.025,.025,.025,.025,.025]; % changed

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Figure A AcCoA
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_A = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 1000; % Bicarbonate
S.init_cond(3) = 600; % Acetyl-CoA
S.init_cond(6) = 0; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(15) = 1300; % NADH
S.init_cond(18) = 0; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [1 1 0 1 1 1 10 1 1 1;
                   1 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta1,Ca1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(1)] = Calc_Function(Ta1,Ca1,S);

[balance_conc_a1, balances_a1, total_conc_a1, carbon_a1] = mass_balance(Ca1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling;  % Using E. coli FabH

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta2,Ca2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(2)] = Calc_Function(Ta2,Ca2,S);

[balance_conc_a2, balances_a2, total_conc_a2, carbon_a2] = mass_balance(Ca2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling; % Using PP FabH1

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta3,Ca3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(3)] = Calc_Function(Ta3,Ca3,S);

[balance_conc_a3, balances_a3, total_conc_a3, carbon_a3] = mass_balance(Ca3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta4,Ca4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(4)] = Calc_Function(Ta4,Ca4,S);

[balance_conc_a4, balances_a4, total_conc_a4, carbon_a4] = mass_balance(Ca4,P);

% Plot
figure('Position',[500 600 250 175])
bar(rel_rate_A,'magenta')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
ylim([0 15])
ax = gca;
ax.FontSize = 10; 
text(0.1, 14, 'Acetyl-CoA','FontSize',10)


%% Figure B OcCoA
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_B = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 1000; % Bicarbonate
S.init_cond(3) = 500; % Acetyl-CoA
S.init_cond(6) = 100; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(15) = 1300; % NADH
S.init_cond(18) = 0; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [1 1 0 1 1 1 10 1 1 1;
                   1 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb1,Cb1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(1)] = Calc_Function(Tb1,Cb1,S);

[balance_conc_b1, balances_b1, total_conc_b1, carbon_b1] = mass_balance(Cb1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling; % Using E. coli FabH

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb2,Cb2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(2)] = Calc_Function(Tb2,Cb2,S);

[balance_conc_b2, balances_b2, total_conc_b2, carbon_b2] = mass_balance(Cb2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling; % Using PP FabH1

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb3,Cb3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(3)] = Calc_Function(Tb3,Cb3,S);

[balance_conc_b3, balances_b3, total_conc_b3, carbon_b3] = mass_balance(Cb3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb4,Cb4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(4)] = Calc_Function(Tb4,Cb4,S);

[balance_conc_b4, balances_b4, total_conc_b4, carbon_b4] = mass_balance(Cb4,P);

% Plot
figure('Position',[500 350 250 175])
bar(rel_rate_B,'magenta')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
ylim([0 15])
ax = gca;
ax.FontSize = 10; 
text(0.1, 14, 'Octanoyl-CoA','FontSize',10)


%% Figure C No Acyl-CoA
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_C = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 1000; % Bicarbonate
S.init_cond(3) = 500; % Acetyl-CoA
S.init_cond(6) = 0; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(15) = 1300; % NADH
S.init_cond(18) = 0; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [1 1 0 1 1 1 10 1 1 1;
                   1 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc1,Cc1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(1)] = Calc_Function(Tc1,Cc1,S);

[balance_conc_c1, balances_c1, total_conc_c1, carbon_c1] = mass_balance(Cc1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling; % Using E. coli FabH

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc2,Cc2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(2)] = Calc_Function(Tc2,Cc2,S);

[balance_conc_c2, balances_c2, total_conc_c2, carbon_c2] = mass_balance(Cc2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling; % Using PP FabH1

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc3,Cc3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(3)] = Calc_Function(Tc3,Cc3,S);

[balance_conc_c3, balances_c3, total_conc_c3, carbon_c3] = mass_balance(Cc3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc4,Cc4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(4)] = Calc_Function(Tc4,Cc4,S);

[balance_conc_c4, balances_c4, total_conc_c4, carbon_c4] = mass_balance(Cc4,P);

% Plot
figure('Position',[500 100 250 175])
bar(rel_rate_C,'magenta')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
ylim([0 15])
ax = gca;
ax.FontSize = 10; 
text(0.1, 14, 'No acyl-CoA','FontSize',10)
