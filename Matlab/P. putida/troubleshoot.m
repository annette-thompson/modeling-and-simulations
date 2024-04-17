%% Balance
% Need to have ran the beginning of the real code

S.labels = {};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 100;

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 0 0 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_troubleshoot(t,c,P);
tic
[T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate] = Calc_Function(T,C,S);

[balance_conc, balances, total_conc, carbon] = mass_balance(C,P);


function dcdt = ODE_Function_troubleshoot(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end


% Set of differential equations


dcdt = [];

end
