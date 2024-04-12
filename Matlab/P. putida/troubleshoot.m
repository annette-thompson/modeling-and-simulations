% List of variable names
var_names = {'ATP', 'C1_Bicarbonate', 'C2_AcCoA', 'C4_SucCoA', 'C6_HexCoA', ...
    'C8_OcCoA', 'C10_DecCoA', 'C12_LauCoA', 'C14_EthCoA', 'C16_PalCoA', 'C18_OcDecCoA'};

% Assuming 'c' is your data array
c = rand(11, 10);  % Example data, replace with your actual data

% Loop to assign variables
for i = 1:numel(var_names)
    var_name = ['c_' var_names{i}];
    eval([var_name ' = c(i, :);']);
end

% Testing the assignment
disp(c_ATP);
disp(c_C1_Bicarbonate);