% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
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

% Figure D EC_FabH AcCoA
S.kcat_scaling_fabH = EC_kcat3_scaling;  % Using E. coli FabH

S.range = [0 150]; % 2.5 mins (initial rate)

rel_rate_D = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 1000; % Bicarbonate
S.init_cond(3) = 100; % Acetyl-CoA
S.init_cond(6) = 0; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(15) = 1300; % NADH
S.init_cond(18) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 1 0 1 1 10 1 1 1;
                   0 1 1 1 1 1 10 1 1 1]; 

% No FabG
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td1,Cd1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(1)] = Calc_Function(Td1,Cd1,S);

[balance_conc_d1, balances_d1, total_conc_d1, carbon_d1] = mass_balance(Cd1,P);

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling;  % Using E. coli FabG

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td2,Cd2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_D(2,:), rel_rate_D(2)] = Calc_Function(Td2,Cd2,S);

[balance_conc_d2, balances_d2, total_conc_d2, carbon_d2] = mass_balance(Cd2,P);

% Profile plot - normal EC setup
% figure()
% F_raw_D_new(1,1:9) = F_raw_D(2,1:9);
% for i=10:14
%     F_raw_D_new(1,i-5) = F_raw_D(2,i-5)+F_raw_D(2,i);
% end
% total = sum(F_raw_D_new);
% bar(F_raw_D_new/total)
% xticklabels([' 4  ';' 6  ';' 8  ';' 10 ';' 12 ';' 14 ';' 16 ';' 18 ';' 20 ';])
% ylabel('Production (mole frac)')
% title("No ACC 2.5 mins")
% xlabel('Chain Length')
% ylim([0 1])

% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling; % Using PP 1914 FabG

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td3,Cd3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(3)] = Calc_Function(Td3,Cd3,S);

[balance_conc_d3, balances_d3, total_conc_d3, carbon_d3] = mass_balance(Cd3,P);

% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling; % Using PP 2783 FabG

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td4,Cd4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(4)] = Calc_Function(Td4,Cd4,S);

[balance_conc_d4, balances_d4, total_conc_d4, carbon_d4] = mass_balance(Cd4,P);

% Plot
figure('Position',[500 600 250 175])
bar(rel_rate_D,'cyan')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabG';'EC FabG';'PP 1914';'PP 2783'])
ylim([0 15])
ax = gca;
ax.FontSize = 8; 
text(0.1, 14, 'Acetyl-CoA','FontSize',10)
text(0.1, 12.5, '1 uM EC FabH','FontSize',10)

%% Figure E PP_FabH2 OcCoA
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;  % Using PP FabH2

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_E = zeros(1,4);

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
enz_conc = [1 1 10 0 1 1 10 1 1 1;
                   1 1 10 1 1 1 10 1 1 1]; 

% No FabG
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te1,Ce1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(1)] = Calc_Function(Te1,Ce1,S);

[balance_conc_e1, balances_e1, total_conc_e1, carbon_e1] = mass_balance(Ce1,P);

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te2,Ce2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(2)] = Calc_Function(Te2,Ce2,S);

[balance_conc_e2, balances_e2, total_conc_e2, carbon_e2] = mass_balance(Ce2,P);

% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling; % Using PP 1914 FabG

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te3,Ce3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(3)] = Calc_Function(Te3,Ce3,S);

[balance_conc_e3, balances_e3, total_conc_e3, carbon_e3] = mass_balance(Ce3,P);

% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling; % Using PP 2783 FabG

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te4,Ce4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(4)] = Calc_Function(Te4,Ce4,S);

[balance_conc_e4, balances_e4, total_conc_e4, carbon_e4] = mass_balance(Ce4,P);

% Plot
figure('Position',[500 350 250 175])
bar(rel_rate_E,'cyan')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabG';'EC FabG';'PP 1914';'PP 2783'])
ylim([0 15])
ax = gca;
ax.FontSize = 8; 
text(0.1, 14, 'Octanoyl-CoA','FontSize',10)
text(0.1, 12.5, '10 uM PP FabH2','FontSize',10)

%% Figure F PP_FabH2 OcCoA
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

% changed
S.range = [0 720]; % 12 mins (total production)

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
enz_conc = [1 1 10 1 1 1 10 1 1 1]; 

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tf1,Cf1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(1,:),~] = Calc_Function(Tf1,Cf1,S);

[balance_conc_f1, balances_f1, total_conc_f1, carbon_f1] = mass_balance(Cf1,P);

% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling; % Using PP 1914 FabG

S.enzyme_conc = enz_conc;
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tf2,Cf2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(2,:),~] = Calc_Function(Tf2,Cf2,S);

[balance_conc_f2, balances_f2, total_conc_f2, carbon_f2] = mass_balance(Cf2,P);

% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling; % Using PP 2783 FabG

S.enzyme_conc = enz_conc;
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tf3,Cf3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc

[F_raw_F(3,:),~] = Calc_Function(Tf3,Cf3,S);


[balance_conc_f3, balances_f3, total_conc_f3, carbon_f3] = mass_balance(Cf3,P);

F_raw_F_new(:,1:9) = F_raw_F(:,1:9);
for i=10:14
    F_raw_F_new(:,i-5) = F_raw_F(:,i-5)+F_raw_F(:,i);
end
F_raw_F_plot = F_raw_F_new(:,4:8);

% Plot
figure('Position',[500 100 250 175])
stack_labels = {'C10','C12','C14','C16','C18'};
bh = bar(F_raw_F_plot, .9, 'stacked');
xticklabels(['EC FabG';'PP 1914';'PP 2783'])
xtickangle(30)
legend(stack_labels)
ylabel('Production (uM)')
ylim([0 150])
xlim([0.4, 3.6])
set(bh, 'FaceColor', 'Flat')
colors =  [106/255, 173/255, 138/255;...
               238/255, 210/255, 148/255;...
               198/255, 96/255, 93/255;...
               145/255, 145/255, 145/255;...
               5/255, 84/255, 117/255];
colors = mat2cell(colors,ones(5,1),3);
set(bh, {'CData'}, colors)
ax = gca;
ax.FontSize = 8; 

% Plot production vs time
% figure()
% x = find(endsWith(S.labels, '_FA') | endsWith(S.labels, '_FA_un'));
% colors = distinguishable_colors(length(x));
% plot(Tf1/60,Cf1(:,x),'LineWidth',1);
% lineHandles = findobj(gca, 'Type', 'Line');
% set(lineHandles, {'Color'}, num2cell(colors, 2));
% legend(cellfun(@(str) strrep(str(3:end), '_', ' '), {S.labels{x}}, 'UniformOutput', false),'Location','bestoutside')
% ylabel('Production (uM)')
% xlabel('Time (min)')
% axis('tight')
% ax = gca;
% ax.FontSize = 18; 
