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
% S.enzyme_conc = [0 1 1 1 1 1 10 1 1 0]; % No FabB
% S.enzyme_conc = [0 1 1 1 1 1 10 1 1 0.1]; % 0.1 uM FabB
S.enzyme_conc = [0 1 1 1 1 1 10 1 1 1]; % 1 uM FabB

% Change specificity for FabF mutant
% S.kcat_scaling_fabF = [0.914,0.914,0.901,1,0.9,0.289,0.0222,0.0222,0.0222];
S.kcat_scaling_fabF = [0.914,0.914,0,0,0,0,0,0,0]; % Only work on C4-C6
S.kcat_scaling_fabF_unsat = [1,1,1,1,0.9,0.289,0.34,0.34,0.34]; % Without FabB, no unsat

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
% [F_raw_FabF_WT, ~] = Calc_Function(T,C,S);
% [F_raw_FabF_C8_Mutant, ~] = Calc_Function(T,C,S);
% [F_raw_FabF_C8_Mutant_01, ~] = Calc_Function(T,C,S);
[F_raw_FabF_C8_Mutant_1, ~] = Calc_Function(T,C,S);
[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

%% Plot 1
bgColor = [255, 233, 213] / 255;
lightpink = [249, 149, 217] / 255;
darkpink = [181, 48, 111] / 255;
goldyellow = [246, 180, 9] / 255;
burntorange = [225, 91, 40] / 255;
maroon = [112, 48, 160] / 255;
darkpurple = [48, 57, 181] / 255;
darkblue = [48, 154, 181] / 255;

% Data
categories = {'4', '6', '8', '10', '12', '14', '16', '18'};
data1 = F_raw_FabF_WT(1:8)/sum(F_raw_FabF_WT(1:9));
data2 = F_raw_FabF_C8_Mutant(1:8)/sum(F_raw_FabF_C8_Mutant(1:9));

% Create a figure
figure;

% Set the size of the figure
set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPosition', [1, 1, 6, 4]);

% Creating bar chart
b = bar(categories, [data1; data2]', 'grouped');

% Adjust the bar width
b(1).BarWidth = 1;
b(2).BarWidth = 1;

% Setting bar colors
b(1).FaceColor = lightpink;
b(2).FaceColor = darkpink;

% Customizing the y-axis ticks
yticks(0:0.2:1);  % Major ticks every 0.2
set(gca, 'YTickLabel', yticks); % Label the major ticks
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off'); % Turn on minor ticks but not minor grid
ax = gca;
ax.YAxis.MinorTickValues = 0.1:0.2:0.9;  % Set minor ticks between major ticks
set(gca, 'TickLength', [0.02, 0.02]);  % Adjust tick length
set(gca, 'TickDir', 'out', 'TickDirMode', 'manual');  % Ticks outside

% Customizing the appearance to match the image
set(gca, 'FontSize', 18, 'FontName', 'Open Sans Medium', 'LineWidth', 1.5);
set(gca, 'YColor', 'k', 'XColor', 'k');

% Legend
legend({'FabF WT', 'C8 Mutant'}, 'Location', 'northeast', 'FontSize', 18, 'Box', 'off');

% Axis labels and title
ylabel('Mole Fraction', 'FontSize', 18, 'FontName', 'Open Sans Medium');
xlabel('Fatty Acid Chain Length', 'FontSize', 18, 'FontName', 'Open Sans Medium');
ylim([0 1])

% Adjusting the figure background
set(gcf, 'Color', 'none');  % Transparent background for the figure
set(gca, 'Color', 'none');  % Transparent background for the axes

% Saving the figure as an SVG file
% export_fig Model_F_C8_Mutant.png -png -transparent -r300

%% Plot 2
% Data
categories = {'4', '6', '8', '10', '12', '14', '16', '18'};
data1 = F_raw_FabF_C8_Mutant(1:8)/sum(F_raw_FabF_C8_Mutant(1:9));
data2 = F_raw_FabF_C8_Mutant_01(1:8)/sum(F_raw_FabF_C8_Mutant_01(1:9));
data3 = F_raw_FabF_C8_Mutant_1(1:8)/sum(F_raw_FabF_C8_Mutant_1(1:9));

% Create a figure
figure;

% Set the size of the figure
set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPosition', [1, 1, 6, 4]);

% Creating bar chart
b = bar(categories, [data1; data2; data3]', 'grouped');

% Adjust the bar width
b(1).BarWidth = 1;
b(2).BarWidth = 1;
b(3).BarWidth = 1;

% Setting bar colors
b(1).FaceColor = darkpink;
b(2).FaceColor = maroon;
b(3).FaceColor = darkpurple;

% Customizing the y-axis ticks
yticks(0:0.2:1);  % Major ticks every 0.2
set(gca, 'YTickLabel', yticks); % Label the major ticks
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off'); % Turn on minor ticks but not minor grid
ax = gca;
ax.YAxis.MinorTickValues = 0.1:0.2:0.9;  % Set minor ticks between major ticks
set(gca, 'TickLength', [0.02, 0.02]);  % Adjust tick length
set(gca, 'TickDir', 'out', 'TickDirMode', 'manual');  % Ticks outside

% Customizing the appearance to match the image
set(gca, 'FontSize', 18, 'FontName', 'Open Sans Medium', 'LineWidth', 1.5);
set(gca, 'YColor', 'k', 'XColor', 'k');

% Legend
legend({'0 uM', '0.1 uM', '1 uM'}, 'Location', 'northeast', 'FontSize', 18, 'Box', 'off');

% Axis labels and title
ylabel('Mole Fraction', 'FontSize', 18, 'FontName', 'Open Sans Medium');
xlabel('Fatty Acid Chain Length', 'FontSize', 18, 'FontName', 'Open Sans Medium');
ylim([0 1])

% Adjusting the figure background
set(gcf, 'Color', 'none');  % Transparent background for the figure
set(gca, 'Color', 'none');  % Transparent background for the axes

% Saving the figure as an SVG file
% export_fig Model_F_B_C8_Mutant.png -png -transparent -r300
