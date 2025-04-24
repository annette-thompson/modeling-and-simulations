% Give access to all necessary folders

my_dir = '/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Set FabG/FabH scaling
% 10 uM PP FabH fit
% EC_kcat3_scaling_fit_A = [5.6325,0,0,0,0,0,0,0,0];
% EC_kcat3_scaling_fit_D = [0.7877,0,0,0,0,0,0,0,0]; 
% PP_H1_kcat3_scaling_fit_A = [0.5746,0,0,0,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fit_B = [0,0,0,0.0284,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fixed_B = [0,0,0,0.3981,0,0,0,0,0]; 

% % 1 uM PP FabH fit
% EC_kcat3_scaling_fit_A = [5.6325,0,0,0,0,0,0,0,0];
% EC_kcat3_scaling_fit_D = [0.7877,0,0,0,0,0,0,0,0]; 
% PP_H1_kcat3_scaling_fit_A = [5.6367,0,0,0,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fit_B = [0,0,0,0.2520,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fixed_B = [0,0,0,3.7639,0,0,0,0,0]; 

% % 1 uM PP FabH fit with initiation
% EC_kcat3_scaling_fit_A_init = [5.5516,0,0,0,0,0,0,0,0];
% EC_kcat3_scaling_fit_D_init = [0.6843,0,0,0,0,0,0,0,0]; 
% PP_H1_kcat3_scaling_fit_A_init = [5.5468,0,0,0,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fit_B_init = [0,0,0,0.2339,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fixed_B_init = [0,0,0,3.7729,0,0,0,0,0]; 

% % NADPH Fit
% EC_kcat3_scaling_fit_A = [6.9576,0,0,0,0,0,0,0,0];
% EC_kcat3_scaling_fit_new_D = [6.9576,0,0,0,0,0,0,0,0]; % same
% PP_H1_kcat3_scaling_fit_A = [1.7043,0,0,0,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fit_B = [0,0,0,3.6464,0,0,0,0,0]; 

% % NADPH Fit with initiation
% EC_kcat3_scaling_fit_A = [6.8708,0,0,0,0,0,0,0,0];
% EC_kcat3_scaling_fit_new_D = [6.8708,0,0,0,0,0,0,0,0]; % same
% PP_H1_kcat3_scaling_fit_A = [1.6041,0,0,0,0,0,0,0,0]; 
% PP_H2_kcat3_scaling_fit_B = [0,0,0,3.5350,0,0,0,0,0]; 

% No fit
EC_kcat3_scaling = [1,0,0,0,0,0,0,0,0];
PP_H1_kcat3_scaling = [1,0,0,0,0,0,0,0,0]; 
PP_H2_kcat3_scaling = [0,0,0,1,0,0,0,0,0]; 

EC_FabH_kon_scaling = [1,0,0,0,0,0,0,0,0];
PP_FabH1_kon_scaling = [1,0,0,0,0,0,0,0,0];
PP_FabH2_kon_scaling = [0,0,0,1,0,0,0,0,0];

EC_kcat4_scaling = [1,1,1,1,1,1,1,1,1];
PP_1914_kcat4_scaling = [.1,.1,.1,.1,.1,.1,.1,.1,.1]; % changed
PP_2783_kcat4_scaling = [0,0,0,.025,.025,.025,.025,.025,.025]; % changed

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');


% Figure A AcCoA

% Run variable code
S = set_vars();

S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; % 2.5 mins (initial rate)

rel_rate_A_fit_A_A_B = zeros(1,4);
NAD_rate_A_fit_A_A_B = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(1) = 0; % ATP
S.init_cond(2) = 0; % Bicarbonate
S.init_cond(3) = 100; % Acetyl-CoA
S.init_cond(6) = 0; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(18) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
               0 1 1 1 1 1 10 1 1 1;
               0 1 1 1 1 1 10 1 1 1];

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Ta1,Ca1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_A_fit_A_A_B(1), NAD_rate_A_fit_A_A_B(1)] = Calc_Function(Ta1,Ca1,S);
[balance_conc_a1, balances_a1, total_conc_a1, carbon_a1] = mass_balance(Ca1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling_fit_A;  % Using E. coli FabH
S.kon_scaling_fabH = EC_FabH_kon_scaling; 

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Ta2,Ca2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_A_fit_A_A_B(2), NAD_rate_A_fit_A_A_B(2)] = Calc_Function(Ta2,Ca2,S);
[balance_conc_a2, balances_a2, total_conc_a2, carbon_a2] = mass_balance(Ca2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling_fit_A; % Using PP FabH1
S.kon_scaling_fabH = PP_FabH1_kon_scaling;

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Ta3,Ca3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_A_fit_A_A_B(3), NAD_rate_A_fit_A_A_B(3)] = Calc_Function(Ta3,Ca3,S);
[balance_conc_a3, balances_a3, total_conc_a3, carbon_a3] = mass_balance(Ca3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling_fit_B; % Using PP FabH2
S.kon_scaling_fabH = PP_FabH2_kon_scaling;

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Ta4,Ca4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_A_fit_A_A_B(4), NAD_rate_A_fit_A_A_B(4)] = Calc_Function(Ta4,Ca4,S);
[balance_conc_a4, balances_a4, total_conc_a4, carbon_a4] = mass_balance(Ca4,P);


% Figure B OcCoA

% Run variable code
S = set_vars();

S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_B_fit_A_A_B = zeros(1,4);
NAD_rate_B_fit_A_A_B = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(1) = 0; % ATP
S.init_cond(2) = 0; % Bicarbonate
S.init_cond(3) = 0; % Acetyl-CoA
S.init_cond(6) = 100; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(18) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
               0 1 1 1 1 1 10 1 1 1;
               0 1 10 1 1 1 10 1 1 1];

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tb1,Cb1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_B_fit_A_A_B(1), NAD_rate_B_fit_A_A_B(1)] = Calc_Function(Tb1,Cb1,S);
[balance_conc_b1, balances_b1, total_conc_b1, carbon_b1] = mass_balance(Cb1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling_fit_A; % Using E. coli FabH
S.kon_scaling_fabH = EC_FabH_kon_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tb2,Cb2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_B_fit_A_A_B(2), NAD_rate_B_fit_A_A_B(2)] = Calc_Function(Tb2,Cb2,S);
[balance_conc_b2, balances_b2, total_conc_b2, carbon_b2] = mass_balance(Cb2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling_fit_A; % Using PP FabH1
S.kon_scaling_fabH = PP_FabH1_kon_scaling;

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tb3,Cb3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_B_fit_A_A_B(3), NAD_rate_B_fit_A_A_B(3)] = Calc_Function(Tb3,Cb3,S);
[balance_conc_b3, balances_b3, total_conc_b3, carbon_b3] = mass_balance(Cb3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling_fit_B; % Using PP FabH2
S.kon_scaling_fabH = PP_FabH2_kon_scaling;

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tb4,Cb4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_B_fit_A_A_B(4), NAD_rate_B_fit_A_A_B(4)] = Calc_Function(Tb4,Cb4,S);
[balance_conc_b4, balances_b4, total_conc_b4, carbon_b4] = mass_balance(Cb4,P);


% Figure C No Acyl-CoA

% Run variable code
S = set_vars();

S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_C_fit_A_A_B = zeros(1,4);
NAD_rate_C_fit_A_A_B = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(1) = 0; % ATP
S.init_cond(2) = 0; % Bicarbonate
S.init_cond(3) = 0; % Acetyl-CoA
S.init_cond(6) = 0; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(18) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
               0 1 1 1 1 1 10 1 1 1;
               0 1 10 1 1 1 10 1 1 1];

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tc1,Cc1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_C_fit_A_A_B(1), NAD_rate_C_fit_A_A_B(1)] = Calc_Function(Tc1,Cc1,S);
[balance_conc_c1, balances_c1, total_conc_c1, carbon_c1] = mass_balance(Cc1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling_fit_A; % Using E. coli FabH
S.kon_scaling_fabH = EC_FabH_kon_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tc2,Cc2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_C_fit_A_A_B(2), NAD_rate_C_fit_A_A_B(2)] = Calc_Function(Tc2,Cc2,S);
[balance_conc_c2, balances_c2, total_conc_c2, carbon_c2] = mass_balance(Cc2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling_fit_A; % Using PP FabH1
S.kon_scaling_fabH = PP_FabH1_kon_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tc3,Cc3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_C_fit_A_A_B(3), NAD_rate_C_fit_A_A_B(3)] = Calc_Function(Tc3,Cc3,S);
[balance_conc_c3, balances_c3, total_conc_c3, carbon_c3] = mass_balance(Cc3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling_fit_B; % Using PP FabH2
S.kon_scaling_fabH = PP_FabH2_kon_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
[Tc4,Cc4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[~, rel_rate_C_fit_A_A_B(4), NAD_rate_C_fit_A_A_B(4)] = Calc_Function(Tc4,Cc4,S);
[balance_conc_c4, balances_c4, total_conc_c4, carbon_c4] = mass_balance(Cc4,P);

%% Plotting
% 
% % Figure 2A
% %figure()
% figure('Position',[500 600 250 175])
% b = bar(rel_rate_A);
% set(b, 'FaceColor', 'Flat')
% color = mat2cell([124/255,28/255,108/255],ones(1,1),3);
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM C16/m)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 15])
% ax = gca;
% ax.FontSize = 10; 
% text(0.1, 14, 'Acetyl-CoA','FontSize',10)
% %title("Acetyl-CoA, FabB Scaling Init = 0.1")
% 
% % A - Acetyl-CoA and Malonyl-CoA
% % figure()
% % colors =  [106/255, 173/255, 138/255;...
% %                238/255, 210/255, 148/255;...
% %                198/255, 96/255, 93/255;...
% %                5/255, 84/255, 117/255];
% % La(1) = length(Ta1);
% % La(2) = length(Ta2);
% % La(3) = length(Ta3);
% % La(4) = length(Ta4);
% % clear Ta
% % Ta = zeros(max(La),4);
% % Ta(1:La(1),1)=Ta1;
% % Ta(1:La(2),2)=Ta2;
% % Ta(1:La(3),3)=Ta3;
% % Ta(1:La(4),4)=Ta4;
% % clear CA
% % CA = zeros(max(La),4);
% % CA(1:La(1),1)=Ca1(:,3);
% % CA(1:La(2),2)=Ca2(:,3);
% % CA(1:La(3),3)=Ca3(:,3);
% % CA(1:La(4),4)=Ca4(:,3);
% % clear CM
% % CM = zeros(max(La),4);
% % CM(1:La(1),1)=Ca1(:,18);
% % CM(1:La(2),2)=Ca2(:,18);
% % CM(1:La(3),3)=Ca3(:,18);
% % CM(1:La(4),4)=Ca4(:,18);
% % for i=1:4
% %     plot(Ta(1:La(i),i)/60,CA(1:La(i),i),'Color',colors(i,:),'LineWidth',2)
% %     hold on
% %     plot(Ta(1:La(i),i)/60,CM(1:La(i),i),'Color',colors(i,:),'LineWidth',2)
% % end
% % ylabel('Concentration (uM)')
% % xlabel('Time (min)')
% % legend('No FabH Acetyl-CoA', 'No FabH Malonyl-CoA',...
% %     'EC FabH Acetyl-CoA', 'EC FabH Malonyl-CoA',...
% %     'PP FabH1 Acetyl-CoA', 'PP FabH1 Malonyl-CoA',...
% %     'PP FabH2 Acetyl-CoA', 'PP FabH2 Malonyl-CoA')
% % title("Acetyl-CoA, FabB Scaling Init = 10")
% 
% % Figure 2B
% %figure()
% figure('Position',[500 350 250 175])
% b = bar(rel_rate_B);
% set(b, 'FaceColor', 'Flat')
% color = mat2cell([124/255,28/255,108/255],ones(1,1),3);
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM C16/m)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 15])
% ax = gca;
% ax.FontSize = 10; 
% text(0.1, 14, 'Octanoyl-CoA','FontSize',10)
% %title("Octonoyl-CoA, FabB Scaling Init = 0.1")
% 
% % B - Acetyl-CoA and Malonyl-CoA
% % figure()
% % colors =  [106/255, 173/255, 138/255;...
% %                238/255, 210/255, 148/255;...
% %                198/255, 96/255, 93/255;...
% %                5/255, 84/255, 117/255];
% % Lb(1) = length(Tb1);
% % Lb(2) = length(Tb2);
% % Lb(3) = length(Tb3);
% % Lb(4) = length(Tb4);
% % clear Tb
% % Tb = zeros(max(Lb),4);
% % Tb(1:Lb(1),1)=Tb1;
% % Tb(1:Lb(2),2)=Tb2;
% % Tb(1:Lb(3),3)=Tb3;
% % Tb(1:Lb(4),4)=Tb4;
% % clear CA
% % CA = zeros(max(Lb),4);
% % CA(1:Lb(1),1)=Cb1(:,3);
% % CA(1:Lb(2),2)=Cb2(:,3);
% % CA(1:Lb(3),3)=Cb3(:,3);
% % CA(1:Lb(4),4)=Cb4(:,3);
% % clear CM
% % CM = zeros(max(Lb),4);
% % CM(1:Lb(1),1)=Cb1(:,18);
% % CM(1:Lb(2),2)=Cb2(:,18);
% % CM(1:Lb(3),3)=Cb3(:,18);
% % CM(1:Lb(4),4)=Cb4(:,18);
% % for i=1:4
% %     plot(Tb(1:Lb(i),i)/60,CA(1:Lb(i),i),'Color',colors(i,:),'LineWidth',2)
% %     hold on
% %     plot(Tb(1:Lb(i),i)/60,CM(1:Lb(i),i),'Color',colors(i,:),'LineWidth',2)
% % end
% % ylabel('Concentration (uM)')
% % xlabel('Time (min)')
% % legend('No FabH Acetyl-CoA', 'No FabH Malonyl-CoA',...
% %     'EC FabH Acetyl-CoA', 'EC FabH Malonyl-CoA',...
% %     'PP FabH1 Acetyl-CoA', 'PP FabH1 Malonyl-CoA',...
% %     'PP FabH2 Acetyl-CoA', 'PP FabH2 Malonyl-CoA')
% % title("Octonoyl-CoA, FabB Scaling Init = 10")
% 
% % Figure 2C
% %figure()
% figure('Position',[500 100 250 175])
% b = bar(rel_rate_C);
% set(b, 'FaceColor', 'Flat')
% color = mat2cell([124/255,28/255,108/255],ones(1,1),3);
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM C16/m)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 15])
% ax = gca;
% ax.FontSize = 10;  
% text(0.1, 14, 'No acyl-CoA','FontSize',10)
% %title("No Acyl-CoA, FabB Scaling Init = 0.1")
% 
% % C - Acetyl-CoA and Malonyl-CoA
% % figure()
% % colors =  [106/255, 173/255, 138/255;...
% %                238/255, 210/255, 148/255;...
% %                198/255, 96/255, 93/255;...
% %                5/255, 84/255, 117/255];
% % Lc(1) = length(Tc1);
% % Lc(2) = length(Tc2);
% % Lc(3) = length(Tc3);
% % Lc(4) = length(Tc4);
% % clear Tc
% % Tc = zeros(max(Lc),4);
% % Tc(1:Lc(1),1)=Tc1;
% % Tc(1:Lc(2),2)=Tc2;
% % Tc(1:Lc(3),3)=Tc3;
% % Tc(1:Lc(4),4)=Tc4;
% % clear CA
% % CA = zeros(max(Lc),4);
% % CA(1:Lc(1),1)=Cc1(:,3);
% % CA(1:Lc(2),2)=Cc2(:,3);
% % CA(1:Lc(3),3)=Cc3(:,3);
% % CA(1:Lc(4),4)=Cc4(:,3);
% % clear CM
% % CM = zeros(max(Lc),4);
% % CM(1:Lc(1),1)=Cc1(:,18);
% % CM(1:Lc(2),2)=Cc2(:,18);
% % CM(1:Lc(3),3)=Cc3(:,18);
% % CM(1:Lc(4),4)=Cc4(:,18);
% % for i=1:4
% %     plot(Tc(1:Lc(i),i)/60,CA(1:Lc(i),i),'Color',colors(i,:),'LineWidth',2)
% %     hold on
% %     plot(Tc(1:Lc(i),i)/60,CM(1:Lc(i),i),'Color',colors(i,:),'LineWidth',2)
% % end
% % ylabel('Concentration (uM)')
% % xlabel('Time (min)')
% % legend('No FabH Acetyl-CoA', 'No FabH Malonyl-CoA',...
% %     'EC FabH Acetyl-CoA', 'EC FabH Malonyl-CoA',...
% %     'PP FabH1 Acetyl-CoA', 'PP FabH1 Malonyl-CoA',...
% %     'PP FabH2 Acetyl-CoA', 'PP FabH2 Malonyl-CoA')
% % title("No Acyl-CoA, FabB Scaling Init = 10")

%% Plotting Palmitic Acid Equivalents
% 
% color = mat2cell([124/255,28/255,108/255],ones(1,1),3);
% 
% % Figure 2A
% figure()
% b = bar(rel_rate_A);
% set(b, 'FaceColor', 'Flat')
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM C16/min)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 15])
% 
% % Figure 2B
% figure()
% b = bar(rel_rate_B);
% set(b, 'FaceColor', 'Flat')
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM C16/min)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 15])
% 
% % Figure 2C
% figure()
% b = bar(rel_rate_C);
% set(b, 'FaceColor', 'Flat')
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM C16/min)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 15])

%% Plotting NAD(P)H conversion
% 
% color = mat2cell([124/255,28/255,108/255],ones(1,1),3);
% 
% % Figure 2A
% figure('Position',[500 600 250 175])
% b = bar(NAD_rate_A);
% set(b, 'FaceColor', 'Flat')
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM NAD(P)/min)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 120])
% ax = gca;
% ax.FontSize = 10; 
% 
% % Figure 2B
% figure('Position',[500 350 250 175])
% b = bar(NAD_rate_B);
% set(b, 'FaceColor', 'Flat')
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM NAD(P)/min)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 120])
% ax = gca;
% ax.FontSize = 10; 
% 
% % Figure 2C
% figure('Position',[500 100 250 175])
% b = bar(NAD_rate_C);
% set(b, 'FaceColor', 'Flat')
% set(b, {'CData'},color)
% ylabel('Initial Rate (uM NAD(P)/min)')
% xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
% ylim([0 120])
% ax = gca;
% ax.FontSize = 10; 
