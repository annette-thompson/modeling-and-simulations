function dcdt = ODE_Function_troubleshoot(t,c,P,num)
% Contains all the differential equations and enzyme balances that define
% the FAS model
% Input:
% t: time values (required as input for MATLAB ODE solver, sec)
% c: concentration values (all components and intermediates, uM)
% P: structure containing all kinetic parameters
% num: number of variables (y)
% Output:
% dcdt: values of differential equations for given 
% concentrations and kinetic parameters

%% ODEs

conc = {'c_ATP', 'c_Bicarbonate'};


% Assign values to each variable
for i = 1:numel(conc)
    % Get the variable name
    var_name = conc{i};
    
    % Get the corresponding index for the value from vector c
    index = c(i);
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = ' num2str(index) ';']);
end

% Set of differential equations
% ATP
d_ATP = -P.k1_2r.*c_ATP + P.k1_2f.*c_Bicarbonate;

% Bicarbonate
d_Bicarbonate = P.k1_2r.*c_ATP - P.k1_2f.*c_Bicarbonate;


dcdt = [d_ATP; d_Bicarbonate];

