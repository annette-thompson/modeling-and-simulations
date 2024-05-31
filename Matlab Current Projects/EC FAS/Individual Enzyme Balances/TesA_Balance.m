% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% TesA Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ACP', 'c_C4_AcACP', 'c_C6_AcACP', 'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP',...
    'c_C14_AcACP', 'c_C16_AcACP', 'c_C18_AcACP', 'c_C20_AcACP', 'c_C12_AcACP_un',...
    'c_C14_AcACP_un', 'c_C16_AcACP_un', 'c_C18_AcACP_un', 'c_C20_AcACP_un', 'c_C4_FA',...
    'c_C6_FA', 'c_C8_FA', 'c_C10_FA', 'c_C12_FA', 'c_C14_FA', 'c_C16_FA', 'c_C18_FA', 'c_C20_FA',...
    'c_C12_FA_un', 'c_C14_FA_un', 'c_C16_FA_un', 'c_C18_FA_un', 'c_C20_FA_un', 'c_C4_TesA_AcACP',...
    'c_C6_TesA_AcACP', 'c_C8_TesA_AcACP', 'c_C10_TesA_AcACP', 'c_C12_TesA_AcACP',...
    'c_C14_TesA_AcACP', 'c_C16_TesA_AcACP', 'c_C18_TesA_AcACP', 'c_C20_TesA_AcACP',...
    'c_C12_TesA_AcACP_un', 'c_C14_TesA_AcACP_un', 'c_C16_TesA_AcACP_un',...
    'c_C18_TesA_AcACP_un', 'c_C20_TesA_AcACP_un', 'c_TesA_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(2:15) = 0.7143; % Acyl-ACPs - evenly distribute 10 ACP to all of them

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 0 0 10 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_TesA(t,c,P);
tic
[T_TesA,C_TesA] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_TesA] = Calc_Function(T_TesA,C_TesA,S);

[balance_conc_TesA, balances_TesA, total_conc_TesA, carbon_TesA] = mass_balance(C_TesA,P);


function dcdt = ODE_Function_TesA(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

c_TesA = P.TesAtot- c_TesA_ACP...
    - c_C4_TesA_AcACP - c_C6_TesA_AcACP - c_C8_TesA_AcACP - c_C10_TesA_AcACP - c_C12_TesA_AcACP - c_C14_TesA_AcACP - c_C16_TesA_AcACP - c_C18_TesA_AcACP - c_C20_TesA_AcACP...
    - c_C12_TesA_AcACP_un - c_C14_TesA_AcACP_un - c_C16_TesA_AcACP_un - c_C18_TesA_AcACP_un - c_C20_TesA_AcACP_un;

% Set of differential equations
% ACP
d_ACP = P.kcat7(1).*c_C4_TesA_AcACP + P.kcat7(2).*c_C6_TesA_AcACP + P.kcat7(3).*c_C8_TesA_AcACP + P.kcat7(4).*c_C10_TesA_AcACP...
    + P.kcat7(5).*c_C12_TesA_AcACP + P.kcat7(6).*c_C14_TesA_AcACP + P.kcat7(7).*c_C16_TesA_AcACP + P.kcat7(8).*c_C18_TesA_AcACP + P.kcat7(9).*c_C20_TesA_AcACP...
    + P.kcat7(5).*c_C12_TesA_AcACP_un + P.kcat7(6).*c_C14_TesA_AcACP_un + P.kcat7(7).*c_C16_TesA_AcACP_un + P.kcat7(8).*c_C18_TesA_AcACP_un + P.kcat7(9).*c_C20_TesA_AcACP_un...
    + P.k7_inh_r.*c_TesA_ACP   - P.k7_inh_f.*c_TesA.*c_ACP; 

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C4_AcACP   = P.k7_1r(1).*c_C4_TesA_AcACP  - P.k7_1f(1).*c_TesA.*c_C4_AcACP;
d_C6_AcACP   = P.k7_1r(2).*c_C6_TesA_AcACP  - P.k7_1f(2).*c_TesA.*c_C6_AcACP;
d_C8_AcACP   = P.k7_1r(3).*c_C8_TesA_AcACP  - P.k7_1f(3).*c_TesA.*c_C8_AcACP;
d_C10_AcACP = P.k7_1r(4).*c_C10_TesA_AcACP - P.k7_1f(4).*c_TesA.*c_C10_AcACP;
d_C12_AcACP = P.k7_1r(5).*c_C12_TesA_AcACP - P.k7_1f(5).*c_TesA.*c_C12_AcACP;
d_C14_AcACP = P.k7_1r(6).*c_C14_TesA_AcACP - P.k7_1f(6).*c_TesA.*c_C14_AcACP;
d_C16_AcACP = P.k7_1r(7).*c_C16_TesA_AcACP - P.k7_1f(7).*c_TesA.*c_C16_AcACP;
d_C18_AcACP = P.k7_1r(8).*c_C18_TesA_AcACP - P.k7_1f(8).*c_TesA.*c_C18_AcACP;

% C2n (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP = P.k7_1r(9).*c_C20_TesA_AcACP - P.k7_1f(9).*c_TesA.*c_C20_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C12_AcACP_un = P.k7_1r(5).*c_C12_TesA_AcACP_un - P.k7_1f(5).*c_TesA.*c_C12_AcACP_un;
d_C14_AcACP_un = P.k7_1r(6).*c_C14_TesA_AcACP_un - P.k7_1f(6).*c_TesA.*c_C14_AcACP_un;
d_C16_AcACP_un = P.k7_1r(7).*c_C16_TesA_AcACP_un - P.k7_1f(7).*c_TesA.*c_C16_AcACP_un;
d_C18_AcACP_un = P.k7_1r(8).*c_C18_TesA_AcACP_un - P.k7_1f(8).*c_TesA.*c_C18_AcACP_un;

% C2n:1 (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP_un = P.k7_1r(9).*c_C20_TesA_AcACP_un - P.k7_1f(9).*c_TesA.*c_C20_AcACP_un;

% Fatty Acids
d_C4_FA   = P.kcat7(1).*c_C4_TesA_AcACP;
d_C6_FA   = P.kcat7(2).*c_C6_TesA_AcACP;
d_C8_FA   = P.kcat7(3).*c_C8_TesA_AcACP;
d_C10_FA = P.kcat7(4).*c_C10_TesA_AcACP;
d_C12_FA = P.kcat7(5).*c_C12_TesA_AcACP;
d_C14_FA = P.kcat7(6).*c_C14_TesA_AcACP;
d_C16_FA = P.kcat7(7).*c_C16_TesA_AcACP;
d_C18_FA = P.kcat7(8).*c_C18_TesA_AcACP;
d_C20_FA = P.kcat7(9).*c_C20_TesA_AcACP;

% Fatty Acids (unsaturated)
d_C12_FA_un = P.kcat7(5).*c_C12_TesA_AcACP_un;
d_C14_FA_un = P.kcat7(6).*c_C14_TesA_AcACP_un;
d_C16_FA_un = P.kcat7(7).*c_C16_TesA_AcACP_un;
d_C18_FA_un = P.kcat7(8).*c_C18_TesA_AcACP_un;
d_C20_FA_un = P.kcat7(9).*c_C20_TesA_AcACP_un;

% C2n (n=2:10) TesA-Acyl-ACPs
d_C4_TesA_AcACP   = P.k7_1f(1).*c_TesA.*c_C4_AcACP   - P.k7_1r(1).*c_C4_TesA_AcACP  - P.kcat7(1).*c_C4_TesA_AcACP;
d_C6_TesA_AcACP   = P.k7_1f(2).*c_TesA.*c_C6_AcACP   - P.k7_1r(2).*c_C6_TesA_AcACP  - P.kcat7(2).*c_C6_TesA_AcACP;
d_C8_TesA_AcACP   = P.k7_1f(3).*c_TesA.*c_C8_AcACP   - P.k7_1r(3).*c_C8_TesA_AcACP  - P.kcat7(3).*c_C8_TesA_AcACP;
d_C10_TesA_AcACP = P.k7_1f(4).*c_TesA.*c_C10_AcACP - P.k7_1r(4).*c_C10_TesA_AcACP - P.kcat7(4).*c_C10_TesA_AcACP;
d_C12_TesA_AcACP = P.k7_1f(5).*c_TesA.*c_C12_AcACP - P.k7_1r(5).*c_C12_TesA_AcACP - P.kcat7(5).*c_C12_TesA_AcACP;
d_C14_TesA_AcACP = P.k7_1f(6).*c_TesA.*c_C14_AcACP - P.k7_1r(6).*c_C14_TesA_AcACP - P.kcat7(6).*c_C14_TesA_AcACP;
d_C16_TesA_AcACP = P.k7_1f(7).*c_TesA.*c_C16_AcACP - P.k7_1r(7).*c_C16_TesA_AcACP - P.kcat7(7).*c_C16_TesA_AcACP;
d_C18_TesA_AcACP = P.k7_1f(8).*c_TesA.*c_C18_AcACP - P.k7_1r(8).*c_C18_TesA_AcACP - P.kcat7(8).*c_C18_TesA_AcACP;
d_C20_TesA_AcACP = P.k7_1f(9).*c_TesA.*c_C20_AcACP - P.k7_1r(9).*c_C20_TesA_AcACP - P.kcat7(9).*c_C20_TesA_AcACP;

% C2n:1 (n=6:10) TesA-Acyl-ACPs
d_C12_TesA_AcACP_un = P.k7_1f(5).*c_TesA.*c_C12_AcACP_un - P.k7_1r(5).*c_C12_TesA_AcACP_un - P.kcat7(5).*c_C12_TesA_AcACP_un;
d_C14_TesA_AcACP_un = P.k7_1f(6).*c_TesA.*c_C14_AcACP_un - P.k7_1r(6).*c_C14_TesA_AcACP_un - P.kcat7(6).*c_C14_TesA_AcACP_un;
d_C16_TesA_AcACP_un = P.k7_1f(7).*c_TesA.*c_C16_AcACP_un - P.k7_1r(7).*c_C16_TesA_AcACP_un - P.kcat7(7).*c_C16_TesA_AcACP_un;
d_C18_TesA_AcACP_un = P.k7_1f(8).*c_TesA.*c_C18_AcACP_un - P.k7_1r(8).*c_C18_TesA_AcACP_un - P.kcat7(8).*c_C18_TesA_AcACP_un;
d_C20_TesA_AcACP_un = P.k7_1f(9).*c_TesA.*c_C20_AcACP_un - P.k7_1r(9).*c_C20_TesA_AcACP_un - P.kcat7(9).*c_C20_TesA_AcACP_un;

% TesA-ACP
d_TesA_ACP = P.k7_inh_f.*c_TesA.*c_ACP - P.k7_inh_r.*c_TesA_ACP;


dcdt = [d_ACP;d_C4_AcACP;d_C6_AcACP;d_C8_AcACP;d_C10_AcACP;d_C12_AcACP;...
    d_C14_AcACP;d_C16_AcACP;d_C18_AcACP;d_C20_AcACP;d_C12_AcACP_un;...
    d_C14_AcACP_un;d_C16_AcACP_un;d_C18_AcACP_un;d_C20_AcACP_un;d_C4_FA;...
    d_C6_FA;d_C8_FA;d_C10_FA;d_C12_FA;d_C14_FA;d_C16_FA;d_C18_FA;d_C20_FA;...
    d_C12_FA_un;d_C14_FA_un;d_C16_FA_un;d_C18_FA_un;d_C20_FA_un;d_C4_TesA_AcACP;...
    d_C6_TesA_AcACP;d_C8_TesA_AcACP;d_C10_TesA_AcACP;d_C12_TesA_AcACP;...
    d_C14_TesA_AcACP;d_C16_TesA_AcACP;d_C18_TesA_AcACP;d_C20_TesA_AcACP;...
    d_C12_TesA_AcACP_un;d_C14_TesA_AcACP_un;d_C16_TesA_AcACP_un;...
    d_C18_TesA_AcACP_un;d_C20_TesA_AcACP_un;d_TesA_ACP];

end
