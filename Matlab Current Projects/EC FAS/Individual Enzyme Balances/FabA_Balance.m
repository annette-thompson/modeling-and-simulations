% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% FabA Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ACP', 'c_C4_BHyAcACP', 'c_C6_BHyAcACP', 'c_C8_BHyAcACP', 'c_C10_BHyAcACP',...
    'c_C12_BHyAcACP', 'c_C14_BHyAcACP', 'c_C16_BHyAcACP', 'c_C18_BHyAcACP', 'c_C20_BHyAcACP',...
    'c_C12_BHyAcACP_un', 'c_C14_BHyAcACP_un', 'c_C16_BHyAcACP_un', 'c_C18_BHyAcACP_un',...
    'c_C20_BHyAcACP_un', 'c_C4_EnAcACP', 'c_C6_EnAcACP', 'c_C8_EnAcACP', 'c_C10_EnAcACP',...
    'c_C12_EnAcACP', 'c_C14_EnAcACP', 'c_C16_EnAcACP', 'c_C18_EnAcACP', 'c_C20_EnAcACP',...
    'c_C10_cis3EnAcACP', 'c_C12_EnAcACP_un', 'c_C14_EnAcACP_un', 'c_C16_EnAcACP_un',...
    'c_C18_EnAcACP_un', 'c_C20_EnAcACP_un', 'c_C4_FabA_BHyAcACP', 'c_C6_FabA_BHyAcACP',...
    'c_C8_FabA_BHyAcACP', 'c_C10_FabA_BHyAcACP', 'c_C12_FabA_BHyAcACP', 'c_C14_FabA_BHyAcACP',...
    'c_C16_FabA_BHyAcACP', 'c_C18_FabA_BHyAcACP', 'c_C20_FabA_BHyAcACP',...
    'c_C12_FabA_BHyAcACP_un', 'c_C14_FabA_BHyAcACP_un', 'c_C16_FabA_BHyAcACP_un',...
    'c_C18_FabA_BHyAcACP_un', 'c_C20_FabA_BHyAcACP_un', 'c_C4_FabA_EnAcACP',...
    'c_C6_FabA_EnAcACP', 'c_C8_FabA_EnAcACP', 'c_C10_FabA_EnAcACP', 'c_C12_FabA_EnAcACP',...
    'c_C14_FabA_EnAcACP', 'c_C16_FabA_EnAcACP', 'c_C18_FabA_EnAcACP', 'c_C20_FabA_EnAcACP',...
    'c_C10_FabA_cis3EnAcACP', 'c_C12_FabA_EnAcACP_un', 'c_C14_FabA_EnAcACP_un',...
    'c_C16_FabA_EnAcACP_un', 'c_C18_FabA_EnAcACP_un', 'c_C20_FabA_EnAcACP_un','c_FabA_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(2:15) = 0.7143; % Beta-hydroxy-acyl-ACPs - evenly distribute 10 ACP to all of them

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 0 0 0 0 1 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabA(t,c,P);
tic
[T_FabA,C_FabA] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabA] = Calc_Function(T_FabA,C_FabA,S);

[balance_conc_FabA, balances_FabA, total_conc_FabA, carbon_FabA] = mass_balance(C_FabA,P);


function dcdt = ODE_Function_FabA(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabA
c_FabA = P.FabAtot - c_FabA_ACP - c_C10_FabA_cis3EnAcACP...
    - c_C4_FabA_BHyAcACP - c_C6_FabA_BHyAcACP - c_C8_FabA_BHyAcACP - c_C10_FabA_BHyAcACP - c_C12_FabA_BHyAcACP - c_C14_FabA_BHyAcACP - c_C16_FabA_BHyAcACP - c_C18_FabA_BHyAcACP - c_C20_FabA_BHyAcACP...
    - c_C12_FabA_BHyAcACP_un - c_C14_FabA_BHyAcACP_un - c_C16_FabA_BHyAcACP_un - c_C18_FabA_BHyAcACP_un - c_C20_FabA_BHyAcACP_un...
    - c_C4_FabA_EnAcACP - c_C6_FabA_EnAcACP - c_C8_FabA_EnAcACP - c_C10_FabA_EnAcACP - c_C12_FabA_EnAcACP - c_C14_FabA_EnAcACP - c_C16_FabA_EnAcACP - c_C18_FabA_EnAcACP - c_C20_FabA_EnAcACP...
    - c_C12_FabA_EnAcACP_un - c_C14_FabA_EnAcACP_un - c_C16_FabA_EnAcACP_un - c_C18_FabA_EnAcACP_un - c_C20_FabA_EnAcACP_un;

% Set of differential equations
% ACP
d_ACP = P.k9_inh_r.*c_FabA_ACP   - P.k9_inh_f.*c_FabA.*c_ACP; 

% C2n (n=2:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C4_BHyAcACP   = P.k9_1r(1).*c_C4_FabA_BHyAcACP   - P.k9_1f(1).*c_FabA.*c_C4_BHyAcACP;
d_C6_BHyAcACP   = P.k9_1r(2).*c_C6_FabA_BHyAcACP   - P.k9_1f(2).*c_FabA.*c_C6_BHyAcACP;
d_C8_BHyAcACP   = P.k9_1r(3).*c_C8_FabA_BHyAcACP   - P.k9_1f(3).*c_FabA.*c_C8_BHyAcACP;
d_C10_BHyAcACP = P.k9_1r(4).*c_C10_FabA_BHyAcACP - P.k9_1f(4).*c_FabA.*c_C10_BHyAcACP;
d_C12_BHyAcACP = P.k9_1r(5).*c_C12_FabA_BHyAcACP - P.k9_1f(5).*c_FabA.*c_C12_BHyAcACP;
d_C14_BHyAcACP = P.k9_1r(6).*c_C14_FabA_BHyAcACP - P.k9_1f(6).*c_FabA.*c_C14_BHyAcACP;
d_C16_BHyAcACP = P.k9_1r(7).*c_C16_FabA_BHyAcACP - P.k9_1f(7).*c_FabA.*c_C16_BHyAcACP;
d_C18_BHyAcACP = P.k9_1r(8).*c_C18_FabA_BHyAcACP - P.k9_1f(8).*c_FabA.*c_C18_BHyAcACP;
d_C20_BHyAcACP = P.k9_1r(9).*c_C20_FabA_BHyAcACP - P.k9_1f(9).*c_FabA.*c_C20_BHyAcACP;

% C2n:1 (n=6:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C12_BHyAcACP_un = P.k9_1r_un(5).*c_C12_FabA_BHyAcACP_un - P.k9_1f_un(5).*c_FabA.*c_C12_BHyAcACP_un;
d_C14_BHyAcACP_un = P.k9_1r_un(6).*c_C14_FabA_BHyAcACP_un - P.k9_1f_un(6).*c_FabA.*c_C14_BHyAcACP_un;
d_C16_BHyAcACP_un = P.k9_1r_un(7).*c_C16_FabA_BHyAcACP_un - P.k9_1f_un(7).*c_FabA.*c_C16_BHyAcACP_un;
d_C18_BHyAcACP_un = P.k9_1r_un(8).*c_C18_FabA_BHyAcACP_un - P.k9_1f_un(8).*c_FabA.*c_C18_BHyAcACP_un;
d_C20_BHyAcACP_un = P.k9_1r_un(9).*c_C20_FabA_BHyAcACP_un - P.k9_1f_un(9).*c_FabA.*c_C20_BHyAcACP_un;

% C2n (n=2:10) Enoyl-Acyl-ACPs (FabZ + FabA - FabI) 
d_C4_EnAcACP   = P.k9_3f(1).*c_C4_FabA_EnAcACP   - P.k9_3r(1).*c_FabA.*c_C4_EnAcACP;
d_C6_EnAcACP   = P.k9_3f(2).*c_C6_FabA_EnAcACP   - P.k9_3r(2).*c_FabA.*c_C6_EnAcACP;
d_C8_EnAcACP   = P.k9_3f(3).*c_C8_FabA_EnAcACP   - P.k9_3r(3).*c_FabA.*c_C8_EnAcACP;
d_C10_EnAcACP = P.k9_3f(4).*c_C10_FabA_EnAcACP - P.k9_3r(4).*c_FabA.*c_C10_EnAcACP;
d_C12_EnAcACP = P.k9_3f(5).*c_C12_FabA_EnAcACP - P.k9_3r(5).*c_FabA.*c_C12_EnAcACP;
d_C14_EnAcACP = P.k9_3f(6).*c_C14_FabA_EnAcACP - P.k9_3r(6).*c_FabA.*c_C14_EnAcACP;
d_C16_EnAcACP = P.k9_3f(7).*c_C16_FabA_EnAcACP - P.k9_3r(7).*c_FabA.*c_C16_EnAcACP;
d_C18_EnAcACP = P.k9_3f(8).*c_C18_FabA_EnAcACP - P.k9_3r(8).*c_FabA.*c_C18_EnAcACP;
d_C20_EnAcACP = P.k9_3f(9).*c_C20_FabA_EnAcACP - P.k9_3r(9).*c_FabA.*c_C20_EnAcACP;

% C10 cis-3-Enoyl-Acyl-ACP (FabA - FabB)
d_C10_cis3EnAcACP = P.k9_3f_un(4).*c_C10_FabA_cis3EnAcACP - P.k9_3r_un(4).*c_FabA.*c_C10_cis3EnAcACP;

% C2n:1 (n=6:10) Enoyl-Acyl-ACPs  (FabZ + FabA - FabI)
d_C12_EnAcACP_un = P.k9_3f(5).*c_C12_FabA_EnAcACP_un - P.k9_3r(5).*c_FabA.*c_C12_EnAcACP_un;
d_C14_EnAcACP_un = P.k9_3f(6).*c_C14_FabA_EnAcACP_un - P.k9_3r(6).*c_FabA.*c_C14_EnAcACP_un;
d_C16_EnAcACP_un = P.k9_3f(7).*c_C16_FabA_EnAcACP_un - P.k9_3r(7).*c_FabA.*c_C16_EnAcACP_un;
d_C18_EnAcACP_un = P.k9_3f(8).*c_C18_FabA_EnAcACP_un - P.k9_3r(8).*c_FabA.*c_C18_EnAcACP_un;
d_C20_EnAcACP_un = P.k9_3f(9).*c_C20_FabA_EnAcACP_un - P.k9_3r(9).*c_FabA.*c_C20_EnAcACP_un;

% C2n (n=2:10) FabA-B-hydroxy-acyl-ACPs
d_C4_FabA_BHyAcACP   = P.k9_1f(1).*c_FabA.*c_C4_BHyAcACP   - P.k9_1r(1).*c_C4_FabA_BHyAcACP   + P.k9_2r(1).*c_C4_FabA_EnAcACP  - P.kcat9(1).*c_C4_FabA_BHyAcACP;
d_C6_FabA_BHyAcACP   = P.k9_1f(2).*c_FabA.*c_C6_BHyAcACP   - P.k9_1r(2).*c_C6_FabA_BHyAcACP   + P.k9_2r(2).*c_C6_FabA_EnAcACP  - P.kcat9(2).*c_C6_FabA_BHyAcACP;
d_C8_FabA_BHyAcACP   = P.k9_1f(3).*c_FabA.*c_C8_BHyAcACP   - P.k9_1r(3).*c_C8_FabA_BHyAcACP   + P.k9_2r(3).*c_C8_FabA_EnAcACP  - P.kcat9(3).*c_C8_FabA_BHyAcACP;
d_C10_FabA_BHyAcACP = P.k9_1f(4).*c_FabA.*c_C10_BHyAcACP - P.k9_1r(4).*c_C10_FabA_BHyAcACP + P.k9_2r(4).*c_C10_FabA_EnAcACP - P.kcat9(4).*c_C10_FabA_BHyAcACP;
d_C12_FabA_BHyAcACP = P.k9_1f(5).*c_FabA.*c_C12_BHyAcACP - P.k9_1r(5).*c_C12_FabA_BHyAcACP + P.k9_2r(5).*c_C12_FabA_EnAcACP - P.kcat9(5).*c_C12_FabA_BHyAcACP;
d_C14_FabA_BHyAcACP = P.k9_1f(6).*c_FabA.*c_C14_BHyAcACP - P.k9_1r(6).*c_C14_FabA_BHyAcACP + P.k9_2r(6).*c_C14_FabA_EnAcACP - P.kcat9(6).*c_C14_FabA_BHyAcACP;
d_C16_FabA_BHyAcACP = P.k9_1f(7).*c_FabA.*c_C16_BHyAcACP - P.k9_1r(7).*c_C16_FabA_BHyAcACP + P.k9_2r(7).*c_C16_FabA_EnAcACP - P.kcat9(7).*c_C16_FabA_BHyAcACP;
d_C18_FabA_BHyAcACP = P.k9_1f(8).*c_FabA.*c_C18_BHyAcACP - P.k9_1r(8).*c_C18_FabA_BHyAcACP + P.k9_2r(8).*c_C18_FabA_EnAcACP - P.kcat9(8).*c_C18_FabA_BHyAcACP;
d_C20_FabA_BHyAcACP = P.k9_1f(9).*c_FabA.*c_C20_BHyAcACP - P.k9_1r(9).*c_C20_FabA_BHyAcACP + P.k9_2r(9).*c_C20_FabA_EnAcACP - P.kcat9(9).*c_C20_FabA_BHyAcACP;

% C2n:1 (n=6:10) FabA-B-hydroxy-acyl-ACPs
d_C12_FabA_BHyAcACP_un = P.k9_1f_un(5).*c_FabA.*c_C12_BHyAcACP_un - P.k9_1r_un(5).*c_C12_FabA_BHyAcACP_un + P.k9_2r(5).*c_C12_FabA_EnAcACP_un - P.kcat9(5).*c_C12_FabA_BHyAcACP_un;
d_C14_FabA_BHyAcACP_un = P.k9_1f_un(6).*c_FabA.*c_C14_BHyAcACP_un - P.k9_1r_un(6).*c_C14_FabA_BHyAcACP_un + P.k9_2r(6).*c_C14_FabA_EnAcACP_un - P.kcat9(6).*c_C14_FabA_BHyAcACP_un;
d_C16_FabA_BHyAcACP_un = P.k9_1f_un(7).*c_FabA.*c_C16_BHyAcACP_un - P.k9_1r_un(7).*c_C16_FabA_BHyAcACP_un + P.k9_2r(7).*c_C16_FabA_EnAcACP_un - P.kcat9(7).*c_C16_FabA_BHyAcACP_un;
d_C18_FabA_BHyAcACP_un = P.k9_1f_un(8).*c_FabA.*c_C18_BHyAcACP_un - P.k9_1r_un(8).*c_C18_FabA_BHyAcACP_un + P.k9_2r(8).*c_C18_FabA_EnAcACP_un - P.kcat9(8).*c_C18_FabA_BHyAcACP_un;
d_C20_FabA_BHyAcACP_un = P.k9_1f_un(9).*c_FabA.*c_C20_BHyAcACP_un - P.k9_1r_un(9).*c_C20_FabA_BHyAcACP_un + P.k9_2r(9).*c_C20_FabA_EnAcACP_un - P.kcat9(9).*c_C20_FabA_BHyAcACP_un;

% C2n (n=2:10) FabA-Enoyl-Acyl-ACPs
d_C4_FabA_EnAcACP   = P.kcat9(1).*c_C4_FabA_BHyAcACP   - P.k9_2r(1).*c_C4_FabA_EnAcACP   + P.k9_3r(1).*c_FabA.*c_C4_EnAcACP  - P.k9_3f(1).*c_C4_FabA_EnAcACP;
d_C6_FabA_EnAcACP   = P.kcat9(2).*c_C6_FabA_BHyAcACP   - P.k9_2r(2).*c_C6_FabA_EnAcACP   + P.k9_3r(2).*c_FabA.*c_C6_EnAcACP  - P.k9_3f(2).*c_C6_FabA_EnAcACP;
d_C8_FabA_EnAcACP   = P.kcat9(3).*c_C8_FabA_BHyAcACP   - P.k9_2r(3).*c_C8_FabA_EnAcACP   + P.k9_3r(3).*c_FabA.*c_C8_EnAcACP  - P.k9_3f(3).*c_C8_FabA_EnAcACP;
d_C10_FabA_EnAcACP = P.kcat9(4).*c_C10_FabA_BHyAcACP - P.k9_2r(4).*c_C10_FabA_EnAcACP + P.k9_3r(4).*c_FabA.*c_C10_EnAcACP - P.k9_3f(4).*c_C10_FabA_EnAcACP + P.k9_2r_un(4).*c_C10_FabA_cis3EnAcACP - P.kcat9_un(4).*c_C10_FabA_EnAcACP;
d_C12_FabA_EnAcACP = P.kcat9(5).*c_C12_FabA_BHyAcACP - P.k9_2r(5).*c_C12_FabA_EnAcACP + P.k9_3r(5).*c_FabA.*c_C12_EnAcACP - P.k9_3f(5).*c_C12_FabA_EnAcACP;
d_C14_FabA_EnAcACP = P.kcat9(6).*c_C14_FabA_BHyAcACP - P.k9_2r(6).*c_C14_FabA_EnAcACP + P.k9_3r(6).*c_FabA.*c_C14_EnAcACP - P.k9_3f(6).*c_C14_FabA_EnAcACP;
d_C16_FabA_EnAcACP = P.kcat9(7).*c_C16_FabA_BHyAcACP - P.k9_2r(7).*c_C16_FabA_EnAcACP + P.k9_3r(7).*c_FabA.*c_C16_EnAcACP - P.k9_3f(7).*c_C16_FabA_EnAcACP;
d_C18_FabA_EnAcACP = P.kcat9(8).*c_C18_FabA_BHyAcACP - P.k9_2r(8).*c_C18_FabA_EnAcACP + P.k9_3r(8).*c_FabA.*c_C18_EnAcACP - P.k9_3f(8).*c_C18_FabA_EnAcACP;
d_C20_FabA_EnAcACP = P.kcat9(9).*c_C20_FabA_BHyAcACP - P.k9_2r(9).*c_C20_FabA_EnAcACP + P.k9_3r(9).*c_FabA.*c_C20_EnAcACP - P.k9_3f(9).*c_C20_FabA_EnAcACP;

% FabA-C10 cis-3-Enoyl-Acyl-ACP
d_C10_FabA_cis3EnAcACP = P.kcat9_un(4).*c_C10_FabA_EnAcACP - P.k9_2r_un(4).*c_C10_FabA_cis3EnAcACP + P.k9_3r_un(4).*c_FabA.*c_C10_cis3EnAcACP - P.k9_3f_un(4).*c_C10_FabA_cis3EnAcACP;

% C2n:1 (n=6:10) FabA-Enoyl-Acyl-ACPs
d_C12_FabA_EnAcACP_un = P.kcat9(5).*c_C12_FabA_BHyAcACP_un - P.k9_2r(5).*c_C12_FabA_EnAcACP_un + P.k9_3r(5).*c_FabA.*c_C12_EnAcACP_un - P.k9_3f(5).*c_C12_FabA_EnAcACP_un;
d_C14_FabA_EnAcACP_un = P.kcat9(6).*c_C14_FabA_BHyAcACP_un - P.k9_2r(6).*c_C14_FabA_EnAcACP_un + P.k9_3r(6).*c_FabA.*c_C14_EnAcACP_un - P.k9_3f(6).*c_C14_FabA_EnAcACP_un;
d_C16_FabA_EnAcACP_un = P.kcat9(7).*c_C16_FabA_BHyAcACP_un - P.k9_2r(7).*c_C16_FabA_EnAcACP_un + P.k9_3r(7).*c_FabA.*c_C16_EnAcACP_un - P.k9_3f(7).*c_C16_FabA_EnAcACP_un;
d_C18_FabA_EnAcACP_un = P.kcat9(8).*c_C18_FabA_BHyAcACP_un - P.k9_2r(8).*c_C18_FabA_EnAcACP_un + P.k9_3r(8).*c_FabA.*c_C18_EnAcACP_un - P.k9_3f(8).*c_C18_FabA_EnAcACP_un;
d_C20_FabA_EnAcACP_un = P.kcat9(9).*c_C20_FabA_BHyAcACP_un - P.k9_2r(9).*c_C20_FabA_EnAcACP_un + P.k9_3r(9).*c_FabA.*c_C20_EnAcACP_un - P.k9_3f(9).*c_C20_FabA_EnAcACP_un;

% FabA-ACP
d_FabA_ACP = P.k9_inh_f.*c_FabA.*c_ACP - P.k9_inh_r.*c_FabA_ACP;


dcdt = [d_ACP;d_C4_BHyAcACP;d_C6_BHyAcACP;d_C8_BHyAcACP;d_C10_BHyAcACP;...
    d_C12_BHyAcACP;d_C14_BHyAcACP;d_C16_BHyAcACP;d_C18_BHyAcACP;d_C20_BHyAcACP;...
    d_C12_BHyAcACP_un;d_C14_BHyAcACP_un;d_C16_BHyAcACP_un;d_C18_BHyAcACP_un;...
    d_C20_BHyAcACP_un;d_C4_EnAcACP;d_C6_EnAcACP;d_C8_EnAcACP;d_C10_EnAcACP;...
    d_C12_EnAcACP;d_C14_EnAcACP;d_C16_EnAcACP;d_C18_EnAcACP;d_C20_EnAcACP;...
    d_C10_cis3EnAcACP;d_C12_EnAcACP_un;d_C14_EnAcACP_un;d_C16_EnAcACP_un;...
    d_C18_EnAcACP_un;d_C20_EnAcACP_un;d_C4_FabA_BHyAcACP;d_C6_FabA_BHyAcACP;...
    d_C8_FabA_BHyAcACP;d_C10_FabA_BHyAcACP;d_C12_FabA_BHyAcACP;d_C14_FabA_BHyAcACP;...
    d_C16_FabA_BHyAcACP;d_C18_FabA_BHyAcACP;d_C20_FabA_BHyAcACP;...
    d_C12_FabA_BHyAcACP_un;d_C14_FabA_BHyAcACP_un;d_C16_FabA_BHyAcACP_un;...
    d_C18_FabA_BHyAcACP_un;d_C20_FabA_BHyAcACP_un;d_C4_FabA_EnAcACP;...
    d_C6_FabA_EnAcACP;d_C8_FabA_EnAcACP;d_C10_FabA_EnAcACP;d_C12_FabA_EnAcACP;...
    d_C14_FabA_EnAcACP;d_C16_FabA_EnAcACP;d_C18_FabA_EnAcACP;d_C20_FabA_EnAcACP;...
    d_C10_FabA_cis3EnAcACP;d_C12_FabA_EnAcACP_un;d_C14_FabA_EnAcACP_un;...
    d_C16_FabA_EnAcACP_un;d_C18_FabA_EnAcACP_un;d_C20_FabA_EnAcACP_un;d_FabA_ACP];

end
