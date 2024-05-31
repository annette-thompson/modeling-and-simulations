% Create a matrix to store your y values
miS2 = miS;
miS2(5) = miS(8);
miS2(6) = miS(7);
miS2(7) = miS(9);
miS2(8) = miS(5);
miS2(9) = miS(6);

miP2 = miP;
miP2(5) = miP(8);
miP2(6) = miP(7);
miP2(7) = miP(9);
miP2(8) = miP(5);
miP2(9) = miP(6);

miA2 = miA;
miA2(5) = miA(8);
miA2(6) = miA(7);
miA2(7) = miA(9);
miA2(8) = miA(5);
miA2(9) = miA(6);

miS3 = miS2;
miS3(6) = miS2(8);
miS3(8) = miS2(6);

miP3 = miP2;
miP3(6) = miP2(8);
miP3(8) = miP2(6);

miA3 = miA2;
miA3(6) = miA2(8);
miA3(8) = miA2(6);

table = transpose([miS3,miA3,miP3]);

miS_new = miS2/max(miS);
miA_new = miA2/max(miA);
miP_new = miP2/max(miP);
yValues = [miP_new;miS_new;miA_new];
yValues = transpose(yValues);

% Define the categories
categories = {'Unsaturated Fraction', 'Total Production', 'Average Chain Length'};

% Define custom colors for each category
colors = [0,0,0;1,0,0;0,0,1];

% Create a bar chart
h=bar(yValues);

% Customize the plot
xlabel('Enzyme');
ylabel('Norm. Elementary Effect');
legend(categories);

% Assign custom colors to each bar
for i = 1:3
    h(i).FaceColor = colors(i,:);
end

% Set the x-axis tick labels
xticks(1:9);  % Assuming you have 9 variables
% Function  {'D','H','G','Z','I','T','F','A','B'} ;
xticklabels({'D','H','G','Z','A','F','B','I','T'}); %plot
% Table     {'D','H','G','Z','A','I','B','F','T'} ; 

