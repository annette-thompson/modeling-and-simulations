%% ACC Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ATP', 'c_C1_Bicarbonate', 'c_C2_AcCoA', 'c_ADP', 'c_C3_MalCoA', 'c_ACC_s1', 'c_C1_ACC_s2', 'c_C1_ACC_s3', 'c_C3_ACC_s4'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 100; % ATP
S.init_cond(2) = 100; % Bicarbonate
S.init_cond(3) = 100; % Acetyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [1 0 0 0 0 0 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_ACC(t,c,P);
tic
[T_ACC,C_ACC] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_ACC] = Calc_Function(T_ACC,C_ACC,S);

[balance_conc_ACC, balances_ACC, total_conc_ACC, carbon_ACC] = mass_balance(C_ACC,P);


function dcdt = ODE_Function_ACC(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

c_ACC = P.ACCtot - c_ACC_s1 - c_C1_ACC_s2 - c_C1_ACC_s3 - c_C3_ACC_s4;

% Set of differential equations
% ATP
d_ATP = P.k1_1r.*c_ACC_s1 - P.k1_1f.*c_ACC.*c_ATP;

% Bicarbonate
d_C1_Bicarbonate = P.k1_2r.*c_C1_ACC_s2 - P.k1_2f.*c_ACC_s1.*c_C1_Bicarbonate;

% C2n (n=1:9)-CoA % changed 
d_C2_AcCoA = P.k1_3r.*c_C3_ACC_s4 - P.k1_3f.*c_C1_ACC_s3.*c_C2_AcCoA;

% ADP
d_ADP = P.kcat1_2.*c_C3_ACC_s4;

% Malonyl-CoA
d_C3_MalCoA = P.kcat1_2.*c_C3_ACC_s4;

% ACC Step 1
d_ACC_s1 = P.k1_1f.*c_ACC.*c_ATP - P.k1_1r.*c_ACC_s1 + P.k1_2r.*c_C1_ACC_s2 - P.k1_2f.*c_ACC_s1.*c_C1_Bicarbonate;

% ACC Step 2
d_C1_ACC_s2 = P.k1_2f.*c_ACC_s1.*c_C1_Bicarbonate - P.k1_2r.*c_C1_ACC_s2 - P.kcat1_1.*c_C1_ACC_s2;

% ACC Step 3
d_C1_ACC_s3 = P.kcat1_1.*c_C1_ACC_s2 + P.k1_3r.*c_C3_ACC_s4 - P.k1_3f.*c_C1_ACC_s3.*c_C2_AcCoA;

% ACC Step 4
d_C3_ACC_s4 = P.k1_3f.*c_C1_ACC_s3.*c_C2_AcCoA - P.k1_3r.*c_C3_ACC_s4 - P.kcat1_2.*c_C3_ACC_s4;

dcdt = [d_ATP;d_C1_Bicarbonate;d_C2_AcCoA;d_ADP;d_C3_MalCoA;d_ACC_s1;d_C1_ACC_s2;d_C1_ACC_s3;d_C3_ACC_s4];

end
