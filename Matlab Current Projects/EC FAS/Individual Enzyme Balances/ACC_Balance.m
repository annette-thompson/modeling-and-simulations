% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% ACC Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ATP','c_C1_Bicarbonate','c_C2_AcCoA','c_ADP','c_C3_MalCoA','c_BC_ATP','c_C1_BC_ATP_HCO3','c_C1_BC_Pi_HCO3',...
    'c_C1_BC_Pi_HCO3_BCCP_Biotin','c_C1_BCCP_Biotin_CO2','c_C1_CT_BCCP_Biotin_CO2','c_C1_CT_Act','c_C3_CT_Act_AcCoA'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1000; % ATP
S.init_cond(2) = 1000; % Bicarbonate
S.init_cond(3) = 600; % Acetyl-CoA

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

figure()
plot(T_ACC,C_ACC)
legend(S.labels)

function dcdt = ODE_Function_ACC(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% ACC % changed for ACC
% BC (AccC)
c_ACC_C = P.ACC_Ctot - c_BC_ATP - c_C1_BC_ATP_HCO3 - c_C1_BC_Pi_HCO3 - c_C1_BC_Pi_HCO3_BCCP_Biotin;
% BCCP-Biotin (AccB-Biotin)
c_ACC_B = P.ACC_Btot - c_C1_BC_Pi_HCO3_BCCP_Biotin - c_C1_BCCP_Biotin_CO2 - c_C1_CT_BCCP_Biotin_CO2;
% CT (AccAD)
c_ACC_AD = P.ACC_ADtot - c_C1_CT_BCCP_Biotin_CO2 - c_C1_CT_Act - c_C3_CT_Act_AcCoA;

% Set of differential equations
% ATP % changed for ACC
d_ATP = P.k1_1r.*c_BC_ATP - P.k1_1f.*c_ACC_C.*c_ATP;

% Bicarbonate % changed for ACC
d_C1_Bicarbonate = P.k1_2r.*c_C1_BC_ATP_HCO3 - P.k1_2f.*c_BC_ATP.*c_C1_Bicarbonate;

% C2n (n=1:9)-CoA % changed for ACC
d_C2_AcCoA = P.k1_5r.*c_C3_CT_Act_AcCoA - P.k1_5f.*c_C1_CT_Act.*c_C2_AcCoA;

% ADP % changed for ACC
d_ADP = P.kcat1_1.*c_C1_BC_ATP_HCO3;

% Malonyl-CoA kcat1_4 % changed for ACC
d_C3_MalCoA = P.kcat1_4.*c_C3_CT_Act_AcCoA; 

% BC-ATP k1_1 % changed for ACC
d_BC_ATP = P.k1_1f.*c_ACC_C.*c_ATP - P.k1_1r.*c_BC_ATP + P.k1_2r.*c_C1_BC_ATP_HCO3 - P.k1_2f.*c_BC_ATP.*c_C1_Bicarbonate;

% BC-ATP-HCO3 k1_2 % changed for ACC
d_C1_BC_ATP_HCO3 = P.k1_2f.*c_BC_ATP.*c_C1_Bicarbonate - P.k1_2r.*c_C1_BC_ATP_HCO3 - P.kcat1_1.*c_C1_BC_ATP_HCO3;

% BC-Pi-HCO3 kcat1_1 % changed for ACC
d_C1_BC_Pi_HCO3 = P.kcat1_1.*c_C1_BC_ATP_HCO3 + P.k1_3r.*c_C1_BC_Pi_HCO3_BCCP_Biotin - P.k1_3f.*c_ACC_B.*c_C1_BC_Pi_HCO3;

% BC-Pi-HCO3-BCCP-Biotin k1_3 % changed for ACC
d_C1_BC_Pi_HCO3_BCCP_Biotin = P.k1_3f.*c_ACC_B.*c_C1_BC_Pi_HCO3 - P.k1_3r.*c_C1_BC_Pi_HCO3_BCCP_Biotin - P.kcat1_2.*c_C1_BC_Pi_HCO3_BCCP_Biotin;

% BCCP-Biotin-CO2 kcat1_2 % changed for ACC
d_C1_BCCP_Biotin_CO2 = P.kcat1_2.*c_C1_BC_Pi_HCO3_BCCP_Biotin + P.k1_4r.*c_C1_CT_BCCP_Biotin_CO2 - P.k1_4f.*c_ACC_AD.*c_C1_BCCP_Biotin_CO2;

% CT-BCCP-Biotin-CO2 k1_4 % changed for ACC
d_C1_CT_BCCP_Biotin_CO2 = P.k1_4f.*c_ACC_AD.*c_C1_BCCP_Biotin_CO2 - P.k1_4r.*c_C1_CT_BCCP_Biotin_CO2 - P.kcat1_3.*c_C1_CT_BCCP_Biotin_CO2;

% CT* kcat1_3 % changed for ACC
d_C1_CT_Act = P.kcat1_3.*c_C1_CT_BCCP_Biotin_CO2 + P.k1_5r.*c_C3_CT_Act_AcCoA - P.k1_5f.*c_C1_CT_Act.*c_C2_AcCoA;

% CT*-AcCoA k1_5 % changed for ACC
d_C3_CT_Act_AcCoA = P.k1_5f.*c_C1_CT_Act.*c_C2_AcCoA - P.k1_5r.*c_C3_CT_Act_AcCoA - P.kcat1_4.*c_C3_CT_Act_AcCoA;

dcdt = [d_ATP;d_C1_Bicarbonate;d_C2_AcCoA;d_ADP;d_C3_MalCoA;d_BC_ATP;d_C1_BC_ATP_HCO3;d_C1_BC_Pi_HCO3;...
    d_C1_BC_Pi_HCO3_BCCP_Biotin;d_C1_BCCP_Biotin_CO2;d_C1_CT_BCCP_Biotin_CO2;d_C1_CT_Act;d_C3_CT_Act_AcCoA];

end
