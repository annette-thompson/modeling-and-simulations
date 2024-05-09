%% FabD Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ACP', 'c_C3_MalCoA', 'c_CoA', 'c_C3_MalACP', 'c_C3_FabD_MalCoA', 'c_C3_FabD_Act', 'c_C3_FabD_Act_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 10; % ACP
S.init_cond(2) = 100; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 0 0 0 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabD(t,c,P);
tic
[T_FabD,C_FabD] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabD] = Calc_Function(T_FabD,C_FabD,S);

[balance_conc_FabD, balances_FabD, total_conc_FabD, carbon_FabD] = mass_balance(C_FabD,P);


function dcdt = ODE_Function_FabD(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabD
c_FabD = P.FabDtot - c_C3_FabD_MalCoA - c_C3_FabD_Act - c_C3_FabD_Act_ACP;

% Set of differential equations
% ACP
d_ACP = P.k2_3r.*c_C3_FabD_Act_ACP - P.k2_3f.*c_C3_FabD_Act.*c_ACP;

% Malonyl-CoA
d_C3_MalCoA = P.k2_1r.*c_C3_FabD_MalCoA - P.k2_1f.*c_FabD.*c_C3_MalCoA; 

% CoA % changed
d_CoA = P.k2_2f.*c_C3_FabD_MalCoA - P.k2_2r.*c_C3_FabD_Act.*c_CoA;

% Malonyl-ACP % changed
d_C3_MalACP = P.k2_4f.*c_C3_FabD_Act_ACP - P.k2_4r.*c_FabD.*c_C3_MalACP;

% FabD-Malonyl-CoA
d_C3_FabD_MalCoA = P.k2_1f.*c_FabD.*c_C3_MalCoA - P.k2_1r.*c_C3_FabD_MalCoA + P.k2_2r.*c_C3_FabD_Act.*c_CoA - P.k2_2f.*c_C3_FabD_MalCoA; 

% FabD*
d_C3_FabD_Act = P.k2_2f.*c_C3_FabD_MalCoA - P.k2_2r.*c_C3_FabD_Act.*c_CoA + P.k2_3r.*c_C3_FabD_Act_ACP - P.k2_3f.*c_C3_FabD_Act.*c_ACP;

% FabD*-ACP
d_C3_FabD_Act_ACP = P.k2_3f.*c_C3_FabD_Act.*c_ACP - P.k2_3r.*c_C3_FabD_Act_ACP + P.k2_4r.*c_FabD.*c_C3_MalACP - P.k2_4f.*c_C3_FabD_Act_ACP;


dcdt = [d_ACP;d_C3_MalCoA;d_CoA;d_C3_MalACP;d_C3_FabD_MalCoA;d_C3_FabD_Act;d_C3_FabD_Act_ACP];

end
