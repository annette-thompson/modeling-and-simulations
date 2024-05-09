% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% FabH Balance
% Need to have ran the beginning of the real code

S.labels = {'c_C2_AcCoA', 'c_C4_SucCoA', 'c_C6_HexCoA', 'c_C8_OcCoA', 'c_C10_DecCoA',...
    'c_C12_LauCoA', 'c_C14_EthCoA', 'c_C16_PalCoA', 'c_C18_OcDecCoA', 'c_ACP', 'c_CoA',...
    'c_C3_MalACP', 'c_C1_CO2', 'c_C4_BKeACP', 'c_C6_BKeACP', 'c_C8_BKeACP', 'c_C10_BKeACP',...
    'c_C12_BKeACP', 'c_C14_BKeACP', 'c_C16_BKeACP', 'c_C18_BKeACP', 'c_C20_BKeACP',...
    'c_C4_AcACP', 'c_C6_AcACP', 'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP', 'c_C14_AcACP',...
    'c_C16_AcACP', 'c_C18_AcACP', 'c_C20_AcACP', 'c_C12_AcACP_un', 'c_C14_AcACP_un',...
    'c_C16_AcACP_un', 'c_C18_AcACP_un', 'c_C20_AcACP_un', 'c_C2_FabH_CoA', 'c_C4_FabH_CoA',...
    'c_C6_FabH_CoA', 'c_C8_FabH_CoA', 'c_C10_FabH_CoA', 'c_C12_FabH_CoA', 'c_C14_FabH_CoA',...
    'c_C16_FabH_CoA', 'c_C18_FabH_CoA', 'c_C2_FabH_Act', 'c_C4_FabH_Act', 'c_C6_FabH_Act',...
    'c_C8_FabH_Act', 'c_C10_FabH_Act', 'c_C12_FabH_Act', 'c_C14_FabH_Act', 'c_C16_FabH_Act',...
    'c_C18_FabH_Act', 'c_C5_FabH_Act_MalACP', 'c_C7_FabH_Act_MalACP', 'c_C9_FabH_Act_MalACP',...
    'c_C11_FabH_Act_MalACP', 'c_C13_FabH_Act_MalACP', 'c_C15_FabH_Act_MalACP',...
    'c_C17_FabH_Act_MalACP', 'c_C19_FabH_Act_MalACP', 'c_C21_FabH_Act_MalACP',...
    'c_C4_FabH_AcACP', 'c_C6_FabH_AcACP', 'c_C8_FabH_AcACP', 'c_C10_FabH_AcACP',...
    'c_C12_FabH_AcACP', 'c_C14_FabH_AcACP', 'c_C16_FabH_AcACP', 'c_C18_FabH_AcACP',...
    'c_C20_FabH_AcACP', 'c_C12_FabH_AcACP_un', 'c_C14_FabH_AcACP_un',...
    'c_C16_FabH_AcACP_un', 'c_C18_FabH_AcACP_un', 'c_C20_FabH_AcACP_un',...
    'c_C6_FabH_Act_AcACP', 'c_C8_FabH_Act_AcACP', 'c_C10_FabH_Act_AcACP',...
    'c_C12_FabH_Act_AcACP', 'c_C14_FabH_Act_AcACP', 'c_C16_FabH_Act_AcACP',...
    'c_C18_FabH_Act_AcACP', 'c_C20_FabH_Act_AcACP', 'c_C22_FabH_Act_AcACP',...
    'c_C14_FabH_Act_AcACP_un', 'c_C16_FabH_Act_AcACP_un', 'c_C18_FabH_Act_AcACP_un',...
    'c_C20_FabH_Act_AcACP_un', 'c_C22_FabH_Act_AcACP_un', 'c_FabH_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 100; % Acetyl-CoA
S.init_cond(12) = 10; % Malonyl-ACP
%S.init_cond(4) = 100; % Octanoyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 1 0 0 0 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabH(t,c,P);
tic
[T_FabH,C_FabH] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabH] = Calc_Function(T_FabH,C_FabH,S);

[balance_conc_FabH, balances_FabH, total_conc_FabH, carbon_FabH] = mass_balance(C_FabH,P);


function dcdt = ODE_Function_FabH(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabH % changed
c_FabH = P.FabHtot - c_FabH_ACP...
    - c_C2_FabH_CoA - c_C4_FabH_CoA - c_C6_FabH_CoA - c_C8_FabH_CoA - c_C10_FabH_CoA - c_C12_FabH_CoA - c_C14_FabH_CoA - c_C16_FabH_CoA - c_C18_FabH_CoA...
    - c_C2_FabH_Act - c_C4_FabH_Act - c_C6_FabH_Act - c_C8_FabH_Act - c_C10_FabH_Act - c_C12_FabH_Act - c_C14_FabH_Act - c_C16_FabH_Act - c_C18_FabH_Act...
    - c_C5_FabH_Act_MalACP - c_C7_FabH_Act_MalACP - c_C9_FabH_Act_MalACP - c_C11_FabH_Act_MalACP - c_C13_FabH_Act_MalACP - c_C15_FabH_Act_MalACP - c_C17_FabH_Act_MalACP - c_C19_FabH_Act_MalACP - c_C21_FabH_Act_MalACP...
    - c_C4_FabH_AcACP - c_C6_FabH_AcACP - c_C8_FabH_AcACP - c_C10_FabH_AcACP - c_C12_FabH_AcACP - c_C14_FabH_AcACP - c_C16_FabH_AcACP - c_C18_FabH_AcACP - c_C20_FabH_AcACP...
    - c_C12_FabH_AcACP_un - c_C14_FabH_AcACP_un - c_C16_FabH_AcACP_un - c_C18_FabH_AcACP_un - c_C20_FabH_AcACP_un...
    - c_C6_FabH_Act_AcACP - c_C8_FabH_Act_AcACP - c_C10_FabH_Act_AcACP - c_C12_FabH_Act_AcACP - c_C14_FabH_Act_AcACP - c_C16_FabH_Act_AcACP - c_C18_FabH_Act_AcACP - c_C20_FabH_Act_AcACP - c_C22_FabH_Act_AcACP...
    - c_C14_FabH_Act_AcACP_un - c_C16_FabH_Act_AcACP_un - c_C18_FabH_Act_AcACP_un - c_C20_FabH_Act_AcACP_un - c_C22_FabH_Act_AcACP_un;

% Set of differential equations
% C2n (n=1:9)-CoA % changed 
d_C2_AcCoA        = P.k3_1r(1).*c_C2_FabH_CoA   - P.k3_1f(1).*c_FabH.*c_C2_AcCoA;
d_C4_SucCoA      = P.k3_1r(2).*c_C4_FabH_CoA   - P.k3_1f(2).*c_FabH.*c_C4_SucCoA;
d_C6_HexCoA      = P.k3_1r(3).*c_C6_FabH_CoA   - P.k3_1f(3).*c_FabH.*c_C6_HexCoA;
d_C8_OcCoA        = P.k3_1r(4).*c_C8_FabH_CoA   - P.k3_1f(4).*c_FabH.*c_C8_OcCoA;
d_C10_DecCoA    = P.k3_1r(5).*c_C10_FabH_CoA  - P.k3_1f(5).*c_FabH.*c_C10_DecCoA;
d_C12_LauCoA     = P.k3_1r(6).*c_C12_FabH_CoA - P.k3_1f(6).*c_FabH.*c_C12_LauCoA;
d_C14_EthCoA     = P.k3_1r(7).*c_C14_FabH_CoA  - P.k3_1f(7).*c_FabH.*c_C14_EthCoA;
d_C16_PalCoA      = P.k3_1r(8).*c_C16_FabH_CoA - P.k3_1f(8).*c_FabH.*c_C16_PalCoA;
d_C18_OcDecCoA = P.k3_1r(9).*c_C18_FabH_CoA - P.k3_1f(9).*c_FabH.*c_C18_OcDecCoA;

% ACP
d_ACP = P.k3_inh_r.*c_FabH_ACP   - P.k3_inh_f.*c_FabH.*c_ACP; 

% CoA % changed
d_CoA = P.k3_2f(1).*c_C2_FabH_CoA - P.k3_2r(1).*c_C2_FabH_Act.*c_CoA...
    + P.k3_2f(2).*c_C4_FabH_CoA   - P.k3_2r(2).*c_C4_FabH_Act.*c_CoA...
    + P.k3_2f(3).*c_C6_FabH_CoA   - P.k3_2r(3).*c_C6_FabH_Act.*c_CoA...
    + P.k3_2f(4).*c_C8_FabH_CoA   - P.k3_2r(4).*c_C8_FabH_Act.*c_CoA...
    + P.k3_2f(5).*c_C10_FabH_CoA - P.k3_2r(5).*c_C10_FabH_Act.*c_CoA...
    + P.k3_2f(6).*c_C12_FabH_CoA - P.k3_2r(6).*c_C12_FabH_Act.*c_CoA...
    + P.k3_2f(7).*c_C14_FabH_CoA - P.k3_2r(7).*c_C14_FabH_Act.*c_CoA...
    + P.k3_2f(8).*c_C16_FabH_CoA - P.k3_2r(8).*c_C16_FabH_Act.*c_CoA...
    + P.k3_2f(9).*c_C18_FabH_CoA - P.k3_2r(9).*c_C18_FabH_Act.*c_CoA;

% Malonyl-ACP % changed
d_C3_MalACP = P.k3_3r(1).*c_C5_FabH_Act_MalACP   - P.k3_3f(1).*c_C2_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(2).*c_C7_FabH_Act_MalACP   - P.k3_3f(2).*c_C4_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(3).*c_C9_FabH_Act_MalACP   - P.k3_3f(3).*c_C6_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(4).*c_C11_FabH_Act_MalACP - P.k3_3f(4).*c_C8_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(5).*c_C13_FabH_Act_MalACP - P.k3_3f(5).*c_C10_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(6).*c_C15_FabH_Act_MalACP - P.k3_3f(6).*c_C12_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(7).*c_C17_FabH_Act_MalACP - P.k3_3f(7).*c_C14_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(8).*c_C19_FabH_Act_MalACP - P.k3_3f(8).*c_C16_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(9).*c_C21_FabH_Act_MalACP - P.k3_3f(9).*c_C18_FabH_Act.*c_C3_MalACP;

% CO2 % changed
d_C1_CO2 = P.kcat3(1).*c_C5_FabH_Act_MalACP + P.kcat3(2).*c_C7_FabH_Act_MalACP + P.kcat3(3).*c_C9_FabH_Act_MalACP + P.kcat3(4).*c_C11_FabH_Act_MalACP...
    + P.kcat3(5).*c_C13_FabH_Act_MalACP + P.kcat3(6).*c_C15_FabH_Act_MalACP + P.kcat3(7).*c_C17_FabH_Act_MalACP + P.kcat3(8).*c_C19_FabH_Act_MalACP + P.kcat3(9).*c_C21_FabH_Act_MalACP;

% C2n (n=2:10) B-ketoacyl-ACPs (FabH + FabF + FabB - FabG) % changed
d_C4_BKeACP    = P.kcat3(1).*c_C5_FabH_Act_MalACP;
d_C6_BKeACP    = P.kcat3(2).*c_C7_FabH_Act_MalACP;
d_C8_BKeACP    = P.kcat3(3).*c_C9_FabH_Act_MalACP;
d_C10_BKeACP  = P.kcat3(4).*c_C11_FabH_Act_MalACP;
d_C12_BKeACP  = P.kcat3(5).*c_C13_FabH_Act_MalACP;
d_C14_BKeACP  = P.kcat3(6).*c_C15_FabH_Act_MalACP;
d_C16_BKeACP  = P.kcat3(7).*c_C17_FabH_Act_MalACP;
d_C18_BKeACP  = P.kcat3(8).*c_C19_FabH_Act_MalACP;
d_C20_BKeACP  = P.kcat3(9).*c_C21_FabH_Act_MalACP;

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C4_AcACP   = P.k3_4r(1).*c_C4_FabH_AcACP  - P.k3_4f(1).*c_FabH.*c_C4_AcACP   + P.k3_5r(1).*c_C6_FabH_Act_AcACP   - P.k3_5f(1).*c_C2_FabH_Act.*c_C4_AcACP;
d_C6_AcACP   = P.k3_4r(2).*c_C6_FabH_AcACP  - P.k3_4f(2).*c_FabH.*c_C6_AcACP   + P.k3_5r(2).*c_C8_FabH_Act_AcACP   - P.k3_5f(2).*c_C2_FabH_Act.*c_C6_AcACP;
d_C8_AcACP   = P.k3_4r(3).*c_C8_FabH_AcACP  - P.k3_4f(3).*c_FabH.*c_C8_AcACP   + P.k3_5r(3).*c_C10_FabH_Act_AcACP - P.k3_5f(3).*c_C2_FabH_Act.*c_C8_AcACP;
d_C10_AcACP = P.k3_4r(4).*c_C10_FabH_AcACP - P.k3_4f(4).*c_FabH.*c_C10_AcACP + P.k3_5r(4).*c_C12_FabH_Act_AcACP - P.k3_5f(4).*c_C2_FabH_Act.*c_C10_AcACP;
d_C12_AcACP = P.k3_4r(5).*c_C12_FabH_AcACP - P.k3_4f(5).*c_FabH.*c_C12_AcACP + P.k3_5r(5).*c_C14_FabH_Act_AcACP - P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP;
d_C14_AcACP = P.k3_4r(6).*c_C14_FabH_AcACP - P.k3_4f(6).*c_FabH.*c_C14_AcACP + P.k3_5r(6).*c_C16_FabH_Act_AcACP - P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP;
d_C16_AcACP = P.k3_4r(7).*c_C16_FabH_AcACP - P.k3_4f(7).*c_FabH.*c_C16_AcACP + P.k3_5r(7).*c_C18_FabH_Act_AcACP - P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP;
d_C18_AcACP = P.k3_4r(8).*c_C18_FabH_AcACP - P.k3_4f(8).*c_FabH.*c_C18_AcACP + P.k3_5r(8).*c_C20_FabH_Act_AcACP - P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP;

% C2n (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP = P.k3_4r(9).*c_C20_FabH_AcACP - P.k3_4f(9).*c_FabH.*c_C20_AcACP + P.k3_5r(9).*c_C22_FabH_Act_AcACP - P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C12_AcACP_un = P.k3_4r(5).*c_C12_FabH_AcACP_un - P.k3_4f(5).*c_FabH.*c_C12_AcACP_un + P.k3_5r(5).*c_C14_FabH_Act_AcACP_un - P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP_un;
d_C14_AcACP_un = P.k3_4r(6).*c_C14_FabH_AcACP_un - P.k3_4f(6).*c_FabH.*c_C14_AcACP_un + P.k3_5r(6).*c_C16_FabH_Act_AcACP_un - P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP_un;
d_C16_AcACP_un = P.k3_4r(7).*c_C16_FabH_AcACP_un - P.k3_4f(7).*c_FabH.*c_C16_AcACP_un + P.k3_5r(7).*c_C18_FabH_Act_AcACP_un - P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP_un;
d_C18_AcACP_un = P.k3_4r(8).*c_C18_FabH_AcACP_un - P.k3_4f(8).*c_FabH.*c_C18_AcACP_un + P.k3_5r(8).*c_C20_FabH_Act_AcACP_un - P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP_un;

% C2n:1 (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP_un = P.k3_4r(9).*c_C20_FabH_AcACP_un - P.k3_4f(9).*c_FabH.*c_C20_AcACP_un + P.k3_5r(9).*c_C22_FabH_Act_AcACP_un - P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP_un;

% C2n (n=1:9) FabH-CoA % changed
d_C2_FabH_CoA   = P.k3_1f(1).*c_FabH.*c_C2_AcCoA        - P.k3_1r(1).*c_C2_FabH_CoA   + P.k3_2r(1).*c_C2_FabH_Act.*c_CoA   - P.k3_2f(1).*c_C2_FabH_CoA; 
d_C4_FabH_CoA   = P.k3_1f(2).*c_FabH.*c_C4_SucCoA      - P.k3_1r(2).*c_C4_FabH_CoA   + P.k3_2r(2).*c_C4_FabH_Act.*c_CoA   - P.k3_2f(2).*c_C4_FabH_CoA; 
d_C6_FabH_CoA   = P.k3_1f(3).*c_FabH.*c_C6_HexCoA      - P.k3_1r(3).*c_C6_FabH_CoA   + P.k3_2r(3).*c_C6_FabH_Act.*c_CoA   - P.k3_2f(3).*c_C6_FabH_CoA; 
d_C8_FabH_CoA   = P.k3_1f(4).*c_FabH.*c_C8_OcCoA        - P.k3_1r(4).*c_C8_FabH_CoA   + P.k3_2r(4).*c_C8_FabH_Act.*c_CoA   - P.k3_2f(4).*c_C8_FabH_CoA; 
d_C10_FabH_CoA = P.k3_1f(5).*c_FabH.*c_C10_DecCoA     - P.k3_1r(5).*c_C10_FabH_CoA + P.k3_2r(5).*c_C10_FabH_Act.*c_CoA - P.k3_2f(5).*c_C10_FabH_CoA; 
d_C12_FabH_CoA = P.k3_1f(6).*c_FabH.*c_C12_LauCoA     - P.k3_1r(6).*c_C12_FabH_CoA + P.k3_2r(6).*c_C12_FabH_Act.*c_CoA - P.k3_2f(6).*c_C12_FabH_CoA;
d_C14_FabH_CoA = P.k3_1f(7).*c_FabH.*c_C14_EthCoA     - P.k3_1r(7).*c_C14_FabH_CoA + P.k3_2r(7).*c_C14_FabH_Act.*c_CoA - P.k3_2f(7).*c_C14_FabH_CoA; 
d_C16_FabH_CoA = P.k3_1f(8).*c_FabH.*c_C16_PalCoA      - P.k3_1r(8).*c_C16_FabH_CoA + P.k3_2r(8).*c_C16_FabH_Act.*c_CoA - P.k3_2f(8).*c_C16_FabH_CoA; 
d_C18_FabH_CoA = P.k3_1f(9).*c_FabH.*c_C18_OcDecCoA - P.k3_1r(9).*c_C18_FabH_CoA + P.k3_2r(9).*c_C18_FabH_Act.*c_CoA - P.k3_2f(9).*c_C18_FabH_CoA; 

% C2n (n=1:9) FabH* % changed
% making FabH* - using FabH* - inhibition from Acyl ACPs (only Acetyl-CoA derived FabH*)
d_C2_FabH_Act   = P.k3_2f(1).*c_C2_FabH_CoA - P.k3_2r(1).*c_C2_FabH_Act.*c_CoA... 
 + P.k3_3r(1).*c_C5_FabH_Act_MalACP       - P.k3_3f(1).*c_C2_FabH_Act.*c_C3_MalACP...
 + P.k3_5r(1).*c_C6_FabH_Act_AcACP        - P.k3_5f(1).*c_C2_FabH_Act.*c_C4_AcACP...
 + P.k3_5r(2).*c_C8_FabH_Act_AcACP        - P.k3_5f(2).*c_C2_FabH_Act.*c_C6_AcACP...
 + P.k3_5r(3).*c_C10_FabH_Act_AcACP        - P.k3_5f(3).*c_C2_FabH_Act.*c_C8_AcACP...
 + P.k3_5r(4).*c_C12_FabH_Act_AcACP      - P.k3_5f(4).*c_C2_FabH_Act.*c_C10_AcACP...
 + P.k3_5r(5).*c_C14_FabH_Act_AcACP      - P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP...
 + P.k3_5r(6).*c_C16_FabH_Act_AcACP      - P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP...
 + P.k3_5r(7).*c_C18_FabH_Act_AcACP      - P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP...
 + P.k3_5r(8).*c_C20_FabH_Act_AcACP      - P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP...
 + P.k3_5r(9).*c_C22_FabH_Act_AcACP      - P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP...
 + P.k3_5r(5).*c_C14_FabH_Act_AcACP_un - P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP_un...
 + P.k3_5r(6).*c_C16_FabH_Act_AcACP_un - P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP_un...
 + P.k3_5r(7).*c_C18_FabH_Act_AcACP_un - P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP_un...
 + P.k3_5r(8).*c_C20_FabH_Act_AcACP_un - P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP_un...
 + P.k3_5r(9).*c_C22_FabH_Act_AcACP_un - P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP_un;
d_C4_FabH_Act   = P.k3_2f(2).*c_C4_FabH_CoA  - P.k3_2r(2).*c_C4_FabH_Act.*c_CoA   + P.k3_3r(2).*c_C7_FabH_Act_MalACP   - P.k3_3f(2).*c_C4_FabH_Act.*c_C3_MalACP;
d_C6_FabH_Act   = P.k3_2f(3).*c_C6_FabH_CoA  - P.k3_2r(3).*c_C6_FabH_Act.*c_CoA   + P.k3_3r(3).*c_C9_FabH_Act_MalACP   - P.k3_3f(3).*c_C6_FabH_Act.*c_C3_MalACP;
d_C8_FabH_Act   = P.k3_2f(4).*c_C8_FabH_CoA  - P.k3_2r(4).*c_C8_FabH_Act.*c_CoA   + P.k3_3r(4).*c_C11_FabH_Act_MalACP   - P.k3_3f(4).*c_C8_FabH_Act.*c_C3_MalACP;
d_C10_FabH_Act = P.k3_2f(5).*c_C10_FabH_CoA - P.k3_2r(5).*c_C10_FabH_Act.*c_CoA + P.k3_3r(5).*c_C13_FabH_Act_MalACP - P.k3_3f(5).*c_C10_FabH_Act.*c_C3_MalACP;
d_C12_FabH_Act = P.k3_2f(6).*c_C12_FabH_CoA - P.k3_2r(6).*c_C12_FabH_Act.*c_CoA + P.k3_3r(6).*c_C15_FabH_Act_MalACP - P.k3_3f(6).*c_C12_FabH_Act.*c_C3_MalACP;
d_C14_FabH_Act = P.k3_2f(7).*c_C14_FabH_CoA - P.k3_2r(7).*c_C14_FabH_Act.*c_CoA + P.k3_3r(7).*c_C17_FabH_Act_MalACP - P.k3_3f(7).*c_C14_FabH_Act.*c_C3_MalACP;
d_C16_FabH_Act = P.k3_2f(8).*c_C16_FabH_CoA - P.k3_2r(8).*c_C16_FabH_Act.*c_CoA + P.k3_3r(8).*c_C19_FabH_Act_MalACP - P.k3_3f(8).*c_C16_FabH_Act.*c_C3_MalACP;
d_C18_FabH_Act = P.k3_2f(9).*c_C18_FabH_CoA - P.k3_2r(9).*c_C18_FabH_Act.*c_CoA + P.k3_3r(9).*c_C21_FabH_Act_MalACP - P.k3_3f(9).*c_C18_FabH_Act.*c_C3_MalACP;

% C2n (n=1:9) FabH*-Malonyl-ACP % changed
d_C5_FabH_Act_MalACP   = P.k3_3f(1).*c_C2_FabH_Act.*c_C3_MalACP   - P.k3_3r(1).*c_C5_FabH_Act_MalACP  - P.kcat3(1).*c_C5_FabH_Act_MalACP; 
d_C7_FabH_Act_MalACP   = P.k3_3f(2).*c_C4_FabH_Act.*c_C3_MalACP   - P.k3_3r(2).*c_C7_FabH_Act_MalACP  - P.kcat3(2).*c_C7_FabH_Act_MalACP; 
d_C9_FabH_Act_MalACP   = P.k3_3f(3).*c_C6_FabH_Act.*c_C3_MalACP   - P.k3_3r(3).*c_C9_FabH_Act_MalACP  - P.kcat3(3).*c_C9_FabH_Act_MalACP; 
d_C11_FabH_Act_MalACP = P.k3_3f(4).*c_C8_FabH_Act.*c_C3_MalACP   - P.k3_3r(4).*c_C11_FabH_Act_MalACP  - P.kcat3(4).*c_C11_FabH_Act_MalACP; 
d_C13_FabH_Act_MalACP = P.k3_3f(5).*c_C10_FabH_Act.*c_C3_MalACP - P.k3_3r(5).*c_C13_FabH_Act_MalACP - P.kcat3(5).*c_C13_FabH_Act_MalACP; 
d_C15_FabH_Act_MalACP = P.k3_3f(6).*c_C12_FabH_Act.*c_C3_MalACP - P.k3_3r(6).*c_C15_FabH_Act_MalACP - P.kcat3(6).*c_C15_FabH_Act_MalACP; 
d_C17_FabH_Act_MalACP = P.k3_3f(7).*c_C14_FabH_Act.*c_C3_MalACP - P.k3_3r(7).*c_C17_FabH_Act_MalACP - P.kcat3(7).*c_C17_FabH_Act_MalACP; 
d_C19_FabH_Act_MalACP = P.k3_3f(8).*c_C16_FabH_Act.*c_C3_MalACP - P.k3_3r(8).*c_C19_FabH_Act_MalACP - P.kcat3(8).*c_C19_FabH_Act_MalACP; 
d_C21_FabH_Act_MalACP = P.k3_3f(9).*c_C18_FabH_Act.*c_C3_MalACP - P.k3_3r(9).*c_C21_FabH_Act_MalACP - P.kcat3(9).*c_C21_FabH_Act_MalACP; 

% C2n (n=2:10) FabH-Acyl-ACPs %no change
d_C4_FabH_AcACP   = P.k3_4f(1).*c_FabH.*c_C4_AcACP   - P.k3_4r(1).*c_C4_FabH_AcACP;
d_C6_FabH_AcACP   = P.k3_4f(2).*c_FabH.*c_C6_AcACP   - P.k3_4r(2).*c_C6_FabH_AcACP;
d_C8_FabH_AcACP   = P.k3_4f(3).*c_FabH.*c_C8_AcACP   - P.k3_4r(3).*c_C8_FabH_AcACP;
d_C10_FabH_AcACP = P.k3_4f(4).*c_FabH.*c_C10_AcACP - P.k3_4r(4).*c_C10_FabH_AcACP;
d_C12_FabH_AcACP = P.k3_4f(5).*c_FabH.*c_C12_AcACP - P.k3_4r(5).*c_C12_FabH_AcACP;
d_C14_FabH_AcACP = P.k3_4f(6).*c_FabH.*c_C14_AcACP - P.k3_4r(6).*c_C14_FabH_AcACP;
d_C16_FabH_AcACP = P.k3_4f(7).*c_FabH.*c_C16_AcACP - P.k3_4r(7).*c_C16_FabH_AcACP;
d_C18_FabH_AcACP = P.k3_4f(8).*c_FabH.*c_C18_AcACP - P.k3_4r(8).*c_C18_FabH_AcACP;
d_C20_FabH_AcACP = P.k3_4f(9).*c_FabH.*c_C20_AcACP - P.k3_4r(9).*c_C20_FabH_AcACP;

% C2n:1 (n=6:10) FabH-Acyl-ACPs %no change
d_C12_FabH_AcACP_un = P.k3_4f(5).*c_FabH.*c_C12_AcACP_un - P.k3_4r(5).*c_C12_FabH_AcACP_un;
d_C14_FabH_AcACP_un = P.k3_4f(6).*c_FabH.*c_C14_AcACP_un - P.k3_4r(6).*c_C14_FabH_AcACP_un;
d_C16_FabH_AcACP_un = P.k3_4f(7).*c_FabH.*c_C16_AcACP_un - P.k3_4r(7).*c_C16_FabH_AcACP_un;
d_C18_FabH_AcACP_un = P.k3_4f(8).*c_FabH.*c_C18_AcACP_un - P.k3_4r(8).*c_C18_FabH_AcACP_un;
d_C20_FabH_AcACP_un = P.k3_4f(9).*c_FabH.*c_C20_AcACP_un - P.k3_4r(9).*c_C20_FabH_AcACP_un;

% C2n (n=2:10) FabH*-Acyl-ACPs %no change
d_C6_FabH_Act_AcACP   = P.k3_5f(1).*c_C2_FabH_Act.*c_C4_AcACP   - P.k3_5r(1).*c_C6_FabH_Act_AcACP;
d_C8_FabH_Act_AcACP   = P.k3_5f(2).*c_C2_FabH_Act.*c_C6_AcACP   - P.k3_5r(2).*c_C8_FabH_Act_AcACP;
d_C10_FabH_Act_AcACP = P.k3_5f(3).*c_C2_FabH_Act.*c_C8_AcACP   - P.k3_5r(3).*c_C10_FabH_Act_AcACP;
d_C12_FabH_Act_AcACP = P.k3_5f(4).*c_C2_FabH_Act.*c_C10_AcACP - P.k3_5r(4).*c_C12_FabH_Act_AcACP;
d_C14_FabH_Act_AcACP = P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP - P.k3_5r(5).*c_C14_FabH_Act_AcACP;
d_C16_FabH_Act_AcACP = P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP - P.k3_5r(6).*c_C16_FabH_Act_AcACP;
d_C18_FabH_Act_AcACP = P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP - P.k3_5r(7).*c_C18_FabH_Act_AcACP;
d_C20_FabH_Act_AcACP = P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP - P.k3_5r(8).*c_C20_FabH_Act_AcACP;
d_C22_FabH_Act_AcACP = P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP - P.k3_5r(9).*c_C22_FabH_Act_AcACP;

% C2n:1 (n=6:10) FabH*-Acyl-ACPs %no change
d_C14_FabH_Act_AcACP_un = P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP_un - P.k3_5r(5).*c_C14_FabH_Act_AcACP_un;
d_C16_FabH_Act_AcACP_un = P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP_un - P.k3_5r(6).*c_C16_FabH_Act_AcACP_un;
d_C18_FabH_Act_AcACP_un = P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP_un - P.k3_5r(7).*c_C18_FabH_Act_AcACP_un;
d_C20_FabH_Act_AcACP_un = P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP_un - P.k3_5r(8).*c_C20_FabH_Act_AcACP_un;
d_C22_FabH_Act_AcACP_un = P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP_un - P.k3_5r(9).*c_C22_FabH_Act_AcACP_un;

% FabH-ACP
d_FabH_ACP = P.k3_inh_f.*c_FabH.*c_ACP - P.k3_inh_r.*c_FabH_ACP;


dcdt = [d_C2_AcCoA;d_C4_SucCoA;d_C6_HexCoA;d_C8_OcCoA;d_C10_DecCoA;...
    d_C12_LauCoA;d_C14_EthCoA;d_C16_PalCoA;d_C18_OcDecCoA;d_ACP;d_CoA;...
    d_C3_MalACP;d_C1_CO2;d_C4_BKeACP;d_C6_BKeACP;d_C8_BKeACP;d_C10_BKeACP;...
    d_C12_BKeACP;d_C14_BKeACP;d_C16_BKeACP;d_C18_BKeACP;d_C20_BKeACP;...
    d_C4_AcACP;d_C6_AcACP;d_C8_AcACP;d_C10_AcACP;d_C12_AcACP;d_C14_AcACP;...
    d_C16_AcACP;d_C18_AcACP;d_C20_AcACP;d_C12_AcACP_un;d_C14_AcACP_un;...
    d_C16_AcACP_un;d_C18_AcACP_un;d_C20_AcACP_un;d_C2_FabH_CoA;d_C4_FabH_CoA;...
    d_C6_FabH_CoA;d_C8_FabH_CoA;d_C10_FabH_CoA;d_C12_FabH_CoA;d_C14_FabH_CoA;...
    d_C16_FabH_CoA;d_C18_FabH_CoA;d_C2_FabH_Act;d_C4_FabH_Act;d_C6_FabH_Act;...
    d_C8_FabH_Act;d_C10_FabH_Act;d_C12_FabH_Act;d_C14_FabH_Act;d_C16_FabH_Act;...
    d_C18_FabH_Act;d_C5_FabH_Act_MalACP;d_C7_FabH_Act_MalACP;d_C9_FabH_Act_MalACP;...
    d_C11_FabH_Act_MalACP;d_C13_FabH_Act_MalACP;d_C15_FabH_Act_MalACP;...
    d_C17_FabH_Act_MalACP;d_C19_FabH_Act_MalACP;d_C21_FabH_Act_MalACP;...
    d_C4_FabH_AcACP;d_C6_FabH_AcACP;d_C8_FabH_AcACP;d_C10_FabH_AcACP;...
    d_C12_FabH_AcACP;d_C14_FabH_AcACP;d_C16_FabH_AcACP;d_C18_FabH_AcACP;...
    d_C20_FabH_AcACP;d_C12_FabH_AcACP_un;d_C14_FabH_AcACP_un;...
    d_C16_FabH_AcACP_un;d_C18_FabH_AcACP_un;d_C20_FabH_AcACP_un;...
    d_C6_FabH_Act_AcACP;d_C8_FabH_Act_AcACP;d_C10_FabH_Act_AcACP;...
    d_C12_FabH_Act_AcACP;d_C14_FabH_Act_AcACP;d_C16_FabH_Act_AcACP;...
    d_C18_FabH_Act_AcACP;d_C20_FabH_Act_AcACP;d_C22_FabH_Act_AcACP;...
    d_C14_FabH_Act_AcACP_un;d_C16_FabH_Act_AcACP_un;d_C18_FabH_Act_AcACP_un;...
    d_C20_FabH_Act_AcACP_un;d_C22_FabH_Act_AcACP_un;d_FabH_ACP];

end
