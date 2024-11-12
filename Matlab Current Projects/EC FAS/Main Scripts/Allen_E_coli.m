% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Time
S.range = [0 30];

% Initial conditions
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 100; % Bicarbonate
S.init_cond(3) = 200; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1000; % NADPH
S.init_cond(15) = 1000; % NADH
S.init_cond(18) = 0; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc = [10 1 1 1 1 1 0 1 1 1];

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

parameterized_ODEs = @(t,c) ODE_Function_ACC_MM(t,c,P);
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);

%% Acyl-ACPs over Time
% Find indices matching the patterns
pattern = '^c_C\d{1,2}_AcACP|^c_ACP$';
indices = find(cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels));

% Extract labels and chain lengths
labels = S.labels(indices);
chain_lengths = cellfun(@(x) str2double(regexp(x, '(?<=c_C)\d{1,2}(?=_AcACP)', 'match', 'once')), labels, 'UniformOutput', false);
chain_lengths = cell2mat(chain_lengths);
chain_lengths(isnan(chain_lengths)) = 0;

% Sort labels and indices based on chain lengths
[~, sortedOrder] = sort(chain_lengths);
sorted_labels = labels(sortedOrder);
sorted_indices = indices(sortedOrder);

% Animation setup
numFrames = 1441;
frameTimes = linspace(T(2), T(end), numFrames);
figure("Position", [662,379,779,418]);
[~, timeIndice] = min(abs(T - frameTimes(61)));
F_raw = C(timeIndice, sorted_indices);
h = bar(F_raw, 'FaceColor', [0.9, 0.8, 0.6]);
bar_labels = {'holo', '  2 ', '  4 ', '  6 ', '  8 ', ' 10 ', ' 12 ', '12:1', ' 14 ', '14:1', ' 16 ', '16:1', ' 18 ', '18:1', ' 20 ', '20:1'};
set(gca, 'XTick', 1:length(sorted_indices), 'XTickLabel', bar_labels);
ylabel('Concentration [uM]');
xlabel('Acyl-ACPs');
timeText = text(0.9, 0.9, '', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
set(timeText, 'String', sprintf('Time: %.2f', frameTimes(61)));
%ratioText = text(0.9, 0.85, '', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
%set(ratioText, 'String', sprintf('holo/C8: %.2f', F_raw(1) / F_raw(5)));
frameText = text(0.9, 0.95, '', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

% Calculate the ratio and find the longest continuous section
%ratios = C(:, sorted_indices(1)) ./ C(:, sorted_indices(5));
%within_range = (ratios >= 2) & (ratios <= 10) & C(:, sorted_indices(5))>C(:, sorted_indices(4));
%times = find(within_range==1);

% Video setup
v = VideoWriter('Acyl_ACPs_WithoutTesA.mp4','MPEG-4');
open(v);

% Update the bar chart in a loop
for k = 1:182
    % Capture the plot as an image and write it to the video
    frame = getframe(gcf);
    writeVideo(v, frame);

    %F_raw = C(times(k), sorted_indices);
    [~, timeIndice] = min(abs(T - frameTimes(k)));
    F_raw = C(timeIndice,sorted_indices);

    set(h, 'YData', F_raw);
    %set(timeText, 'String', sprintf('Time: %.2f', T(times(k))));
    set(timeText, 'String', sprintf('Time: %.2f', frameTimes(k)));
    %set(ratioText, 'String', sprintf('holo/C8: %.2f', F_raw(1) / F_raw(5)));
    set(gca, 'XTick', 1:length(sorted_indices), 'XTickLabel', bar_labels);
    set(gca,'YLim',[0 max(F_raw)])
    set(frameText, 'String', sprintf('Frame: %i', k));

    %pause(.1);
end

% Close the video file
close(v);

%% Intermediate ACPs over Time
% Find indices matching the patterns
pattern = '^c_C\d{1,2}_(BHyAcACP|EnAcACP|AcACP)';
indices = find(cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), S.labels));

% Extract labels and chain lengths
labels = S.labels(indices);
chain_lengths = cellfun(@(x) str2double(regexp(x, '(?<=c_C)\d{1,2}(?=_(BHyAcACP|EnAcACP|AcACP))', 'match', 'once')), labels, 'UniformOutput', false);
chain_lengths = cell2mat(chain_lengths);
chain_lengths(isnan(chain_lengths)) = 0;

% Sort labels and indices based on chain lengths
[~, sortedOrder] = sort(chain_lengths);
sorted_labels = labels(sortedOrder);
sorted_indices = indices(sortedOrder);
plottedSpecies = [2,3,5,6,8,9,11,12,14,16,19,15,20,22,25,21,28];
plot_indices = sorted_indices(plottedSpecies);
plot_labels = sorted_labels(plottedSpecies);

% Animation setup
numFrames = 1441;
frameTimes = linspace(T(2), T(end), numFrames);
figure("Position", [662,73,779,418]);
[~, timeIndice] = min(abs(T - frameTimes(61)));
F_raw = C(timeIndice, plot_indices);
h = bar(F_raw, 'FaceColor', [0.9, 0.8, 0.6]);
set(gca,'YLim',[0 max(F_raw)])
bar_labels = strrep(plot_labels,'c_','');
bar_labels = strrep(bar_labels,'_',' ');
set(gca,'XTick',1:length(plot_indices),'XTickLabel',bar_labels);
ylabel('Concentration [uM]')
xlabel('Acyl-ACPs')
timeText = text(0.9, 0.9, '', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
set(timeText, 'String', sprintf('Time: %.2f', frameTimes(61)));
frameText = text(0.9, 0.95, '', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

% Video setup
v = VideoWriter('Intermediates_WithoutTesA.mp4','MPEG-4');
open(v);

% Update the bar chart in a loop
for k = 1:182

    % Capture the plot as an image and write it to the video
    frame = getframe(gcf);
    writeVideo(v, frame);

    % Update data
    [~, timeIndice] = min(abs(T - frameTimes(k)));
    F_raw = C(timeIndice,plot_indices);

    % Update bar chart
    set(h, 'YData', F_raw);
    set(gca,'XTick',1:length(plot_indices),'XTickLabel',bar_labels);
    if max(F_raw)>0
        set(gca,'YLim',[0 max(F_raw)])
    end
    set(timeText, 'String', sprintf('Time: %.2f', frameTimes(k)));
    set(frameText, 'String', sprintf('Frame: %i', k));


    %pause(.1);
end

% Close the video file
close(v);


%% Maneesh's Data
% Panel A: pmol/Total OD600
acyl_ACPs_A = {'Holo', 'C2', 'C4', 'C6', 'C8', 'C10', 'C12', 'C14', 'C16', 'C16:1', 'C18', 'C18:1'};
pmol_30C = [220, 63, 8, 4, 18, 10, 8, 2, 2, 10, 0.5, 2];  % 30°C data
pmol_37C = [160, 78, 5, 6, 25, 12, 12, 7, 3, 3, 1, 0.5];  % 37°C data
err_30C_A = [112, 5, 1, 1, 1, 1, 0.5, 0.1, 0.1, 1, 0.1, 0.5];  % Error bars for 30°C
err_37C_A = [85, 15, .5, 1, 3, 2, 2, 2, 1, 1, 0.5, 0.1];  % Error bars for 37°C

% Split y-axis ranges
ylow = [0 35];  % Lower y-axis range
yhigh = [35 375];  % Higher y-axis range

figure;
% First subplot for lower y-axis range
subplot(9,1,1:2); % Top part of the plot
hold on;
b_A_high = bar([pmol_30C; pmol_37C]', 'grouped');

% Set colors for 30°C and 37°C
b_A_high(1).FaceColor = [0 0.4 0.7];  % Blue for 30°C
b_A_high(2).FaceColor = [0.9, 0.8, 0.6];  % Gold for 37°C

% Error bars for the high part
ngroups_A = length(acyl_ACPs_A);
nbars_A = size([pmol_30C; pmol_37C], 1);
groupwidth_A = min(0.8, nbars_A/(nbars_A + 1.5));
for i = 1:nbars_A
    x_A = (1:ngroups_A) - groupwidth_A/2 + (2*i-1) * groupwidth_A / (2*nbars_A);
    if i == 1
        errorbar(x_A, pmol_30C, err_30C_A, 'k', 'linestyle', 'none');
    else
        errorbar(x_A, pmol_37C, err_37C_A, 'k', 'linestyle', 'none');
    end
end

% Statistical significance stars (high part)
sig_indices_A_high = [2];
for i = sig_indices_A_high
    text(i, max(pmol_30C(i) + err_30C_A(i), pmol_37C(i) + err_37C_A(i)), '*', 'FontSize', 18, 'HorizontalAlignment', 'center','Color','red');
end

% Customize upper y-axis range
ylim(yhigh);  % Set the y-axis limit for the high plot
set(gca, 'XTickLabel', []);  % Remove x-tick labels to simulate continuity
ylabel('pmol/Total OD_{600}');
legend('30°C', '37°C');

hold off;

% Second subplot for lower y-axis range
subplot(9,1,3:4); % Bottom part of the plot
hold on;
b_A_low = bar([pmol_30C; pmol_37C]', 'grouped');

% Set colors for 30°C and 37°C
b_A_low(1).FaceColor = [0 0.4 0.7];  % Blue for 30°C
b_A_low(2).FaceColor = [0.9, 0.8, 0.6];  % Gold for 37°C

% Error bars for the low part
for i = 1:nbars_A
    x_A = (1:ngroups_A) - groupwidth_A/2 + (2*i-1) * groupwidth_A / (2*nbars_A);
    if i == 1
        errorbar(x_A, pmol_30C, err_30C_A, 'k', 'linestyle', 'none');
    else
        errorbar(x_A, pmol_37C, err_37C_A, 'k', 'linestyle', 'none');
    end
end

% Statistical significance stars (high part)
sig_indices_A_low = [3,5,7,8,9,10,12];
for i = sig_indices_A_low
    text(i, max(pmol_30C(i) + err_30C_A(i), pmol_37C(i) + err_37C_A(i)), '*', 'FontSize', 18, 'HorizontalAlignment', 'center','Color','red');
end

% Customize lower y-axis range
ylim(ylow);  % Set the y-axis limit for the low plot
set(gca, 'XTickLabel', acyl_ACPs_A, 'XTick', 1:ngroups_A);
ylabel('pmol/Total OD_{600}');
xlabel('acyl-ACPs')
hold off;

% Panel B: Peak area/Total OD600
acyl_ACPs_B = {'β-hydroxy-C4', 'C4-enoyl', 'β-hydroxy-C6', 'C6-enoyl', 'β-hydroxy-C8', ...
               'C8-enoyl', 'β-hydroxy-C10', 'C10-enoyl', 'C10:1', 'β-hydroxy-C12', 'C12-enoyl', ...
               'C12:1', 'β-hydroxy-C12:1', 'β-hydroxy-C14', 'C14-enoyl', 'C14:1', 'β-hydroxy-C14:1', 'C16-enoyl'};
peak_30C = [3.5e5, 1.5e4, 1.3e4, 1.5e4, 1.7e4, 5e4, 5.2e4, 1.1e5, 0.5e4, 1.3e4, 2e4, 1.5e5, 4e4, 0.5e4, 0.5e4, 3e5, 6e4, 0.05e4];  % 30°C data
peak_37C = [3e5, 1.4e4, 1.5e4, 1.4e4, 1.6e4, 4e4, 6e4, 1.2e5, 0e4, 1.4e4, 2.5e4, 1.4e5, 2e4, 0.5e4, 0.5e4, 2e5, 3e4, 0.05e4];  % 37°C data
err_30C_B = [5e4, 1e3, 1e3, 1e3, 3e3, 9e3, 12e3, 20e3, 1e3, 2e3, 3e3, 20e3, 3e3, 0.5e3, 1e3, 20e3, 3e3, 0.1e3];  % Error bars for 30°C
err_37C_B = [10e4, 5e3, 8e3, 8e3, 10e3, 15e3, 15e3, 40e3, 0e3, 3e3, 6e3, 20e3, 4e3, 0.5e3, 1e3, 30e3, 6e3, 0.1e3];  % Error bars for 37°C

% Split y-axis ranges
ylow = [0 1e5];  % Lower y-axis range
yhigh = [1e5 5e5];  % Higher y-axis range

% First subplot for lower y-axis range
subplot(9,1,6:7); % Top part of the plot
hold on;
b_B_high = bar([peak_30C; peak_37C]', 'grouped');

% Set colors for 30°C and 37°C
b_B_high(1).FaceColor = [0 0.4 0.7];  % Blue for 30°C
b_B_high(2).FaceColor = [0.9, 0.8, 0.6];  % Gold for 37°C

% Error bars for the high part
ngroups_B = length(acyl_ACPs_B);
nbars_B = size([peak_30C; peak_37C], 1);
groupwidth_B = min(0.8, nbars_B/(nbars_B + 1.5));
for i = 1:nbars_B
    x_B = (1:ngroups_B) - groupwidth_B/2 + (2*i-1) * groupwidth_B / (2*nbars_B);
    if i == 1
        errorbar(x_B, peak_30C, err_30C_B, 'k', 'linestyle', 'none');
    else
        errorbar(x_B, peak_37C, err_37C_B, 'k', 'linestyle', 'none');
    end
end

% Statistical significance stars (high part)
sig_indices_B_high = [16];
for i = sig_indices_B_high
    text(i, max(peak_30C(i) + err_30C_B(i), peak_37C(i) + err_37C_B(i)), '*', 'FontSize', 18, 'HorizontalAlignment', 'center','Color','red');
end

% Customize upper y-axis range
ylim(yhigh);  % Set the y-axis limit for the high plot
set(gca, 'XTickLabel', []);  % Remove x-tick labels to simulate continuity
ylabel('Peak area/Total OD_{600}');
legend('30°C', '37°C');

hold off;

% Second subplot for lower y-axis range
subplot(9,1,8:9); % Bottom part of the plot
hold on;
b_B_low = bar([peak_30C; peak_37C]', 'grouped');

% Set colors for 30°C and 37°C
b_B_low(1).FaceColor = [0 0.4 0.7];  % Blue for 30°C
b_B_low(2).FaceColor = [0.9, 0.8, 0.6];  % Gold for 37°C

% Error bars for the low part
for i = 1:nbars_B
    x_B = (1:ngroups_B) - groupwidth_B/2 + (2*i-1) * groupwidth_B / (2*nbars_B);
    if i == 1
        errorbar(x_B, peak_30C, err_30C_B, 'k', 'linestyle', 'none');
    else
        errorbar(x_B, peak_37C, err_37C_B, 'k', 'linestyle', 'none');
    end
end

% Statistical significance stars (low part)
sig_indices_B_low = [9,13,17];
for i = sig_indices_B_low
    text(i, max(peak_30C(i) + err_30C_B(i), peak_37C(i) + err_37C_B(i)), '*', 'FontSize', 18, 'HorizontalAlignment', 'center','Color','red');
end

% Customize lower y-axis range
ylim(ylow);  % Set the y-axis limit for the low plot
set(gca, 'XTickLabel', acyl_ACPs_B, 'XTick', 1:ngroups_B);
ylabel('Peak area/Total OD_{600}');
xlabel('acyl-ACPs')
hold off;