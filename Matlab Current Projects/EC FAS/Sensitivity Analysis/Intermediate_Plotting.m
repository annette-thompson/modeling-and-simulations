% Create a matrix to store your y values
table = transpose([miA,miC,miU,miR]);

miA_frac = miA/max(miA);
miC_frac = miC/max(miC);
miU_frac = miU/max(miU);
miR_frac = miR/max(miR);
yValues = [miA_frac;miC_frac;miU_frac;miR_frac];
yValues = transpose(yValues);

% Define the categories
categories = {'Total AcACP Production', 'Average Intermediate Chain Length', 'Unsaturated Fraction', 'EnAcACP to BHyAcACP Ratio'};

% Create a bar chart
figure("WindowState","fullscreen")
h=bar(yValues);
%xlabel('Enzyme');
xlabel('Substrate');
ylabel('Norm. Elementary Effect');
legend(categories);
set(gca,'FontSize',18)
xticks(1:9);
%xticklabels({'D','H','G','Z','I','T','F','A','B'});
xticklabels({'ATP','Bicarbonate','Acetyl-CoA','ACP','NADPH','NADH','Malonyl-CoA'});
