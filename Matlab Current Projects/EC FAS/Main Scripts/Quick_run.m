% Trying to see if parameters run for ACC
clear all

% Give access to all necessary folders
my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
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
S.init_cond(18) = 1500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc = [0 1 1 1 1 1 10 2 1 1];

P = Param_Function(S);

% P.kcat1_1 = 85.17/60*100; % s^-1
% P.Km1_1 = 170; % uM
% P.kcat1_2 = 73.8/60*100; % s^-1
% P.Km1_2 = 370; % uM
% P.kcat1_3 = 1000.8/60*100; % s^-1
% P.Km1_3 = 160; % uM
% P.kcat1_4 = 2031.8/60*100; % s^-1
% P.Km1_4 = 450; % uM
% P.kcat1_5 = 30.1; % s^-1
% P.Km1_5 = 48.7; % uM

x = 1000;
y = 1000;
z = 0.001;
P.k1_1r = x;
P.k1_1f = x/170;
P.k1_2r = x;
P.k1_2f = x/370;
P.kcat1_1 = 73.8/60*y;
P.k1_3r = x;
P.k1_3f = x/(160*z);
P.kcat1_2 = 1000.8/60*y;
P.k1_4r = x;
P.k1_4f = x/(450*z);
P.kcat1_3 = 2031.8/60*y;
P.k1_5r = x;
P.k1_5f = x/48.7;
P.kcat1_4 = 30.1;

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw, rel_rate] = Calc_Function(T,C,S);

[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

%% Plotting
% figure()
% bar_labels = {'4','6','8','10','12','12:1','14','14:1','16','16:1','18','18:1','20','20:1'};
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
% b(1).FaceColor = 'k';
% b(2).FaceColor = 'w';
% set(gca,'XTick',1:length(bar_labels),'XTickLabel',bar_labels);
% ylim([0 25])
% xlabel('Fatty Acid Chain Length','FontSize',18)
% ylabel('Concentration of Fatty Acid \muM','FontSize',18)
% legend({'Experimental Data','Model Result'},'FontSize',18)
% title('Fatty Acid Distributions','FontSize',18)

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
% xticklabels([' 4  ';' 6  ';' 8  ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1';])
% ylabel('Mole Fraction')
% title("My Model: No ACC")
% xlabel('Fatty Acid Chain Length')
% ax = gca;
% ax.FontSize = 18; 
% ylim([0 1])

figure('Position',[440 278 548 300])
%figure()
F_raw_new(1,1:9) = F_raw(1,1:9);
for i=10:14
    F_raw_new(1,i-5) = F_raw(1,i-5)+F_raw(1,i);
end
total = sum(F_raw_new);
bar(F_raw_new(1:8)/total)
xticklabels([' 4  ';' 6  ';' 8  ';' 10 ';' 12 ';' 14 ';' 16 ';' 18 ';])
ylim([0 1])
yticks(0:0.2:1)
ylabel('Mole Fraction')
xlabel('Fatty Acid Chain Length')
text(2.5, 0.85, '1 uM FabB, 2 uM FabF','FontSize',27,'Color','red')
%title("My Model: No ACC")
 ax = gca;
 ax.FontSize = 27; 

% Ac/MalCoA
% figure()
% plot(T/60,C(:,18))
% hold on
% plot(T/60,C(:,3))
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% legend('Malonyl-CoA','Acetyl-CoA')
% title("With ACC 2 hrs")