% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
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
S.init_cond(18) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc = [0 1 1 1 1 1 10 1 1 1];

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

% Solve ODEs
parameterized_ODEs = @(t,c) ODE_Function_ACC_MM(t,c,P);
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[F_raw, ~] = Calc_Function(T,C,S);
[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

%% Plotting
% % Substrates vs Time
% figure()
% plot(T/60,C(:,1))
% hold on
% plot(T/60,C(:,2))
% plot(T/60,C(:,3))
% plot(T/60,C(:,13))
% plot(T/60,C(:,15))
% xlabel('Time (min)')
% ylabel('Concentration (uM)')
% legend('ATP','Bicarbonate','Acetyl-CoA','NADPH','NADH')
% 
% % Malonyl-CoA vs Time
% figure()
% plot(T/60,C(:,18))
% xlabel('Time (min)')
% ylabel('Concentration (uM)')
% legend('Malonyl-CoA')
% 
% % Fatty Acids vs Time
% figure()
% indices = find(endsWith(S.labels,'_FA'));
% indices_un = find(endsWith(S.labels,'_FA_un'));
% C_new = zeros(length(C),9);
% C_new(:,1:4) = C(:,indices(1:4));
% C_new(:,5:9) = C(:,indices(5:9))+C(:,indices_un);
% Lables = strrep(S.labels(indices),'c_','');
% Lables = strrep(Lables,'_FA','');
% colors = distinguishable_colors(9);
% plot(T/60,C_new,'LineWidth',1)
% lineHandles = findobj(gca,'Type','Line');
% set(lineHandles, {'Color'}, num2cell(colors, 2));
% xlabel('Time (min)')
% ylabel('Concentration (uM)')
% legend(Lables)
% % Print relevant values
% avg_lgth = sum(F_raw.*S.FA_dist)/sum(F_raw)
% total_prod = sum(F_raw)
% C16_equiv = sum(F_raw.*S.FA_dist)/16
% 
% % Initial rate bar graph
% figure()
% bar([rel_rate_05G;rel_rate_05G_m])
% ylabel('Initial Rate (uM C16/m)')
% xticklabels('0.05 G WT';'0.05 G MT')
% xtickangle(30)
% 
% % Production profile bar plot (unsat separated) vs experiment
% figure()
% bar_labels = {'  4 ','  6 ','  8 ',' 10',' 12 ','12:1',' 14 ','14:1',' 16 ','16:1',' 18 ','18:1',' 20 ','20:1'};
% F_raw_new(1,1:4) = F_raw(1,1:4);
% j=10;
% for i=5:2:13
%     F_raw_new(1,i) = F_raw(j-5);
%     F_raw_new(1,i+1) = F_raw(j);
%     j = j+1;
% end
% dist_total = F_raw_new;
% total_FA = sum(F_raw_new);
% fit_dist = total_FA.*[0,0,0.059771046,0.011042448,0.202613717,0.075980144,0.29378989,0.07192787,0.143855739,0.117515956,0.002228751,0.02127444,0,0];
% y_bar = cat(1,fit_dist,dist_total);
% b = bar(y_bar');
% b(1).FaceColor ='k';
% b(2).FaceColor ='w';
% set(gca,'XTick',1:length(bar_labels),'XTickLabel',bar_labels);
% ylim([0 25])
% xlabel('Fatty Acid Chain Length','FontSize',18)
% ylabel('Concentration of Fatty Acid \muM','FontSize',18)
% legend({'Experimental Data','Model Result'},'FontSize',18)
% title('Fatty Acid Distributions','FontSize',18)
%
% % Production profile bar plot (unsat separated)
% figure()
% F_raw_new(1,1:4) = F_raw(1,1:4);
% j=10;
% for i=5:2:13
%     F_raw_new(1,i) = F_raw(j-5);
%     F_raw_new(1,i+1) = F_raw(j);
%     j = j+1;
% end
% total = sum(F_raw_new);
% bar(F_raw_new/total)
% xticklabels(['  4 ';'  6 ';'  8 ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1'])
% ylabel('Mole Fraction')
% xlabel('Fatty Acid Chain Length')
% ylim([0 1])
% 
% % Production profile bar plot (unsat + sat)
% figure('Position',[440 278 548 300])
% %figure()
% F_raw_new(1,1:9) = F_raw(1,1:9);
% for i=10:14
%     F_raw_new(1,i-5) = F_raw(1,i-5)+F_raw(1,i);
% end
% total = sum(F_raw_new);
% bar(F_raw_new(1:8)/total,'black')
% xticklabels({'4';'6';'8';'10';'12';'14';'16';'18'})
% ylim([0 1])
% yticks(0:0.2:1)
% ylabel('Mole Fraction')
% xlabel('Fatty Acid Chain Length')

%% Optimization
clear all;close all;clc

% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Set optimization options
format shortG % show significant figures in shortest way
max_fun_evals = 1000; %max number of function evaluations
max_iter = 500; %max number of optimiztion iterations
opt_options = optimset('Display','iter','OutputFcn',@optout);

% p_vec = [a1 a2 a3 b1 c2 c3...
%           c1 b2 b3 d1 d2 e...
%           f c4 x1 x2 x3 x4];
% p_vec = [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
%           4.645634978 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673...
%           132.8499358 2180.050007 0.539756276 0.053673263 34.49718991 11.15058888];

% Initial guess for parameters being optimized
% Just 19: 1.2044, 1544.53
% Just 7: 6.3986, 11772 (with ACC now)
% Both with ACC: 6.9393 2.0016, 10187
x0 = [6.3986 1.2044]; 

% p_vec with x(1),x(2),... for parameters being optimized
fitfunc = @(x) Optimization_Function_Handler(...
    [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
    x(1) 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673...
    132.8499358 2180.050007 0.539756276 0.053673263 34.49718991 11.15058888 x(2)],...
    ODE_options);

% Search for parameter values that make the simulation closest to experiments
[x_sol, fval] = fminsearch(fitfunc,x0,opt_options);

% Play sound when done
% load handel
% sound(y,Fs)

%%
total_obj=zeros(1,1000);

parfor i=1:1000
        total_obj(i) = Optimization_Function_Handler([142473.7238 7597.676912...
            4.276689943 40213.92919 88.88525384 0.005388274 4.645634978 0.006677519...
            0.284982219 -0.285700283 3.348915642 2.886607673 132.8499358 2180.050007...
            0.539756276 0.053673263 34.49718991 11.15058888 i/100],ODE_options);
        pause(30)
end
