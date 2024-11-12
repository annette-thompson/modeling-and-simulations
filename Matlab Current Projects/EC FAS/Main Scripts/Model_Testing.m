% Runs for Model Run Log

%% PNAS Model
% Move to correct folder(s)
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/Previous Models/PNAS';
cd(my_dir)
addpath(genpath(my_dir))
%%
% Time
time_range = [0 720]; % 12 mins

% Initial conditions
init_cond = zeros(304,1);
init_cond(3) = 100; % Acetyl-CoA
init_cond(4) = 10; % holo ACP
init_cond(5) = 1300; % NADPH
init_cond(6) = 1300; % NADH
init_cond(8) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enzyme_conc = [0 1 0 1 1 1 10 1 1 1];

[Total_Production, Product_Profile] = PNAS_Model(time_range,init_cond,enzyme_conc);

% Results
disp('PNAS:')
fprintf('Total Production = %f\n', Total_Production)
fprintf('Product Profile = [%s]\n', num2str(Product_Profile, '%f '))

%% Plotting
% C16 Equivalents vs Time
figure()
plot(T/60,F_weighted)
xlabel('Time [min]')
ylabel('Palmitic Equivalents [uM]')
ylim([0 ceil(Total_Production)])
title('PNAS Model')

% Product Profile
figure()
bar(Product_Profile)
xticklabels(['  4 ';'  6 ';'  8 ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1'])
ylabel('Concentration [uM]')
xlabel('Fatty Acid Chain Length')
title('PNAS Model')

%% ME1 Model
% Move to correct folder(s)
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/Previous Models/ME1';
cd(my_dir)
addpath(genpath(my_dir))

% Time
time_range = [0 720]; % 12 mins

% Initial conditions
init_cond = zeros(315,1);
init_cond(3) = 100; % Acetyl-CoA
init_cond(4) = 10; % holo ACP
init_cond(5) = 1300; % NADPH
init_cond(6) = 1300; % NADH
init_cond(8) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enzyme_conc = [0 1 0 1 1 1 10 1 1 1];

[Total_Production, Product_Profile] = ME1_Model(time_range,init_cond,enzyme_conc);

% Results
disp('ME1:')
fprintf('Total Production = %f\n', Total_Production)
fprintf('Product Profile = [%s]\n', num2str(Product_Profile, '%f '))

%% Plotting
% C16 Equivalents vs Time
figure()
plot(T/60,F_weighted)
xlabel('Time [min]')
ylabel('Palmitic Equivalents [uM]')
ylim([0 ceil(Total_Production)])
title('ME1 Model')

% Product Profile
figure()
bar(Product_Profile)
xticklabels(['  4 ';'  6 ';'  8 ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1'])
ylabel('Concentration [uM]')
xlabel('Fatty Acid Chain Length')
title('ME1 Model')

%% EC FAS Model
% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

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

% Calculate kinetic parameters
P = Param_Function(S);

% Turn off FabB/F Initiation
P.kcat8_H = 0;
P.kcat10_H = 0;

% Make ODEs with new params
parameterized_ODEs = @(t,c) ODE_Function(t,c,P);

% Solve numerically
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);

% Calculates the total concentration of each FA
F_weighted = zeros(length(T),1);
F_raw = zeros(1,length(S.FA_dist));
weight_vec = S.FA_dist/16; % Palmitic Acid equivalents (C16)
count = 1;
for ind = 1:length(S.labels)
    label_val = char(S.labels(ind));
    if contains(label_val, '_FA','IgnoreCase',false)
        F_raw(count) = C(end,ind);
        F_weighted = weight_vec(count)*(C(:,ind)) + F_weighted;
        count = count + 1;
    end
end

% Final production
Total_Production = F_weighted(end);

% Rearrange to match experimental data
% (4,6,8,10,12,12:1,14,14:1,16,16:1,18,18:1,20,20:1)
Product_Profile(1,1:4) = F_raw(1,1:4);
j=10;
for i=5:2:13
    Product_Profile(1,i) = F_raw(j-5);
    Product_Profile(1,i+1) = F_raw(j);
    j = j+1;
end

% Results
disp('EC FAS:')
fprintf('Total Production = %f\n', Total_Production)
fprintf('Product Profile = [%s]\n', num2str(Product_Profile, '%f '))

%% Plotting
% C16 Equivalents vs Time
figure()
plot(T/60,F_weighted,'LineWidth',1)
hold on
plot(save_T/60,save_weight,'LineWidth',1)
warning('off','all')
validation_data = readtable(S.opt_tot_prod_file,'ReadRowNames',false);
warning('on','all')
time_data = validation_data.Time_min_;
conc_data = validation_data.Concentration_um_;
plot(time_data,conc_data,'k.','LineWidth',1)
xlabel('Time [min]')
ylabel('Palmitic Equivalents [uM]')
ylim([0 ceil(Total_Production)])
%title('EC FAS Model')
legend('Model No Change','Model Change 1','Experiment')

% Product Profile
figure()
total_FA = sum(F_raw); % total fatty acid (concentration) from in vitro for scaling
fit_dist = total_FA.*S.opt_prod_dist_data;
profile = [Product_Profile;save_product;fit_dist];
bar(transpose(profile))
xticklabels(['  4 ';'  6 ';'  8 ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1'])
ylabel('Concentration [uM]')
xlabel('Fatty Acid Chain Length')
%title('EC FAS Model')
legend('Model No Change','Model Change 1','Experiment')

%%
y = [0	5.1227864	0	4.1256764;...
    2.6891828	5.2055916	0	0;...
    0	5.1150236	0	0;...
    0.32531666	4.886394558	1.388407888	5.179493346];

figure()
bar(transpose(y))
xticklabels(['   AcCoA No FabH';'   AcCoA    FabH';'No AcCoA No FabH';'No AcCoA    FabH'])
xtickangle(30)
ylabel('uM C16/min')
legend('Change 1','Change 2','Change 3',"Katie's Data")
ylim([0 8])

profile = [0.1500	0.2244	1.9213	1.1289	9.4517	3.1988	12.8771	5.1225	7.3644	2.5405	0.3513	0.7426	0.0108	0.1501;...
    0.2052	0.3053	2.5685	1.5335	10.8556	3.2841	11.1236	4.3659	5.5588	1.9960	0.2165	0.5892	0.0054	0.1171;...
    0.2052	0.3053	2.5684	1.5335	10.8554	3.2840	11.1234	4.3659	5.5588	1.9960	0.2165	0.5892	0.0054	0.1171;...
    0	0	2.553664433	0.471778706	8.656489678	3.246183635	12.55191005	3.073053856	6.146107668	5.020764019	0.095221391	0.908931404	0	0];
figure()
bar(transpose(profile))
xticklabels(['  4 ';'  6 ';'  8 ';' 10 ';' 12 ';'12:1';' 14 ';'14:1';' 16 ';'16:1';' 18 ';'18:1';' 20 ';'20:1'])
ylabel('Concentration [uM]')
xlabel('Fatty Acid Chain Length')
%title('EC FAS Model')
legend('PNAS','ME 2023','Current Model','Experiment scaled to Current Model')