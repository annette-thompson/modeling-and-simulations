% Give access to all necessary folders

my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects';
cd(my_dir)
addpath(genpath(my_dir))

%% Variables

% Run variable code
S = set_vars();

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% FabI Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ACP', 'c_NADH', 'c_C4_EnAcACP', 'c_C6_EnAcACP', 'c_C8_EnAcACP', 'c_C10_EnAcACP',...
    'c_C12_EnAcACP', 'c_C14_EnAcACP', 'c_C16_EnAcACP', 'c_C18_EnAcACP', 'c_C20_EnAcACP',...
    'c_C12_EnAcACP_un', 'c_C14_EnAcACP_un', 'c_C16_EnAcACP_un', 'c_C18_EnAcACP_un',...
    'c_C20_EnAcACP_un', 'c_C4_AcACP', 'c_C6_AcACP', 'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP',...
    'c_C14_AcACP', 'c_C16_AcACP', 'c_C18_AcACP', 'c_C20_AcACP', 'c_C12_AcACP_un',...
    'c_C14_AcACP_un', 'c_C16_AcACP_un', 'c_C18_AcACP_un', 'c_C20_AcACP_un', 'c_FabI_NADH',...
    'c_C4_FabI_NADH_EnAcACP', 'c_C6_FabI_NADH_EnAcACP', 'c_C8_FabI_NADH_EnAcACP',...
    'c_C10_FabI_NADH_EnAcACP', 'c_C12_FabI_NADH_EnAcACP', 'c_C14_FabI_NADH_EnAcACP',...
    'c_C16_FabI_NADH_EnAcACP', 'c_C18_FabI_NADH_EnAcACP', 'c_C20_FabI_NADH_EnAcACP',...
    'c_C12_FabI_NADH_EnAcACP_un', 'c_C14_FabI_NADH_EnAcACP_un',...
    'c_C16_FabI_NADH_EnAcACP_un', 'c_C18_FabI_NADH_EnAcACP_un',...
    'c_C20_FabI_NADH_EnAcACP_un', 'c_FabI_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(2) = 1300; % NADH
S.init_cond(3:16) = 0.7143; % Enol-acyl-ACPs - evenly distribute 10 ACP to all of them

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 0 1 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabI(t,c,P);
tic
[T_FabI,C_FabI] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabI] = Calc_Function(T_FabI,C_FabI,S);

[balance_conc_FabI, balances_FabI, total_conc_FabI, carbon_FabI] = mass_balance(C_FabI,P);


function dcdt = ODE_Function_FabI(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabI
c_FabI = P.FabItot - c_FabI_NADH...
    - c_C4_FabI_NADH_EnAcACP - c_C6_FabI_NADH_EnAcACP - c_C8_FabI_NADH_EnAcACP - c_C10_FabI_NADH_EnAcACP - c_C12_FabI_NADH_EnAcACP - c_C14_FabI_NADH_EnAcACP - c_C16_FabI_NADH_EnAcACP - c_C18_FabI_NADH_EnAcACP - c_C20_FabI_NADH_EnAcACP...
    - c_C12_FabI_NADH_EnAcACP_un - c_C14_FabI_NADH_EnAcACP_un - c_C16_FabI_NADH_EnAcACP_un - c_C18_FabI_NADH_EnAcACP_un - c_C20_FabI_NADH_EnAcACP_un - c_FabI_ACP;

% Set of differential equations
% ACP
d_ACP = P.k6_inh_r.*c_FabI_ACP    - P.k6_inh_f.*c_FabI.*c_ACP; 

% NADH
d_NADH = P.k6_1r(1).*c_FabI_NADH - P.k6_1f(1).*c_FabI.*c_NADH; 

% C2n (n=2:10) Enoyl-Acyl-ACPs (FabZ + FabA - FabI) 
d_C4_EnAcACP   = P.k6_2r(1).*c_C4_FabI_NADH_EnAcACP   - P.k6_2f(1).*c_FabI_NADH.*c_C4_EnAcACP;
d_C6_EnAcACP   = P.k6_2r(2).*c_C6_FabI_NADH_EnAcACP   - P.k6_2f(2).*c_FabI_NADH.*c_C6_EnAcACP;
d_C8_EnAcACP   = P.k6_2r(3).*c_C8_FabI_NADH_EnAcACP   - P.k6_2f(3).*c_FabI_NADH.*c_C8_EnAcACP;
d_C10_EnAcACP = P.k6_2r(4).*c_C10_FabI_NADH_EnAcACP - P.k6_2f(4).*c_FabI_NADH.*c_C10_EnAcACP;
d_C12_EnAcACP = P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP - P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP;
d_C14_EnAcACP = P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP - P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP;
d_C16_EnAcACP = P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP - P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP;
d_C18_EnAcACP = P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP - P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP;
d_C20_EnAcACP = P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP - P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP;

% C2n:1 (n=6:10) Enoyl-Acyl-ACPs  (FabZ + FabA - FabI)
d_C12_EnAcACP_un = P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP_un - P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP_un;
d_C14_EnAcACP_un = P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP_un - P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP_un;
d_C16_EnAcACP_un = P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP_un - P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP_un;
d_C18_EnAcACP_un = P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP_un - P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP_un;
d_C20_EnAcACP_un = P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP_un - P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP_un;

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C4_AcACP   = P.kcat6(1).*c_C4_FabI_NADH_EnAcACP;
d_C6_AcACP   = P.kcat6(2).*c_C6_FabI_NADH_EnAcACP;
d_C8_AcACP   = P.kcat6(3).*c_C8_FabI_NADH_EnAcACP;
d_C10_AcACP = P.kcat6(4).*c_C10_FabI_NADH_EnAcACP;
d_C12_AcACP = P.kcat6(5).*c_C12_FabI_NADH_EnAcACP;
d_C14_AcACP = P.kcat6(6).*c_C14_FabI_NADH_EnAcACP;
d_C16_AcACP = P.kcat6(7).*c_C16_FabI_NADH_EnAcACP;
d_C18_AcACP = P.kcat6(8).*c_C18_FabI_NADH_EnAcACP;

% C2n (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP = P.kcat6(9).*c_C20_FabI_NADH_EnAcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C12_AcACP_un = P.kcat6(5).*c_C12_FabI_NADH_EnAcACP_un;
d_C14_AcACP_un = P.kcat6(6).*c_C14_FabI_NADH_EnAcACP_un;
d_C16_AcACP_un = P.kcat6(7).*c_C16_FabI_NADH_EnAcACP_un;
d_C18_AcACP_un = P.kcat6(8).*c_C18_FabI_NADH_EnAcACP_un;

% C2n:1 (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP_un = P.kcat6(9).*c_C20_FabI_NADH_EnAcACP_un;

% FabI-NADH
d_FabI_NADH = P.k6_1f(1).*c_FabI.*c_NADH - P.k6_1r(1).*c_FabI_NADH...
 + P.k6_2r(1).*c_C4_FabI_NADH_EnAcACP - P.k6_2f(1).*c_FabI_NADH.*c_C4_EnAcACP...
 + P.k6_2r(2).*c_C6_FabI_NADH_EnAcACP - P.k6_2f(2).*c_FabI_NADH.*c_C6_EnAcACP...
 + P.k6_2r(3).*c_C8_FabI_NADH_EnAcACP - P.k6_2f(3).*c_FabI_NADH.*c_C8_EnAcACP...
 + P.k6_2r(4).*c_C10_FabI_NADH_EnAcACP - P.k6_2f(4).*c_FabI_NADH.*c_C10_EnAcACP...
 + P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP - P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP...
 + P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP - P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP...
 + P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP - P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP...
 + P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP - P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP...
 + P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP - P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP...
 + P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP_un - P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP_un...
 + P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP_un - P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP_un...
 + P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP_un - P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP_un...
 + P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP_un - P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP_un...
 + P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP_un - P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP_un;

% C2n (n=2:10) FabI-NADH-Enoyl-Acyl-ACPs
d_C4_FabI_NADH_EnAcACP   = P.k6_2f(1).*c_FabI_NADH.*c_C4_EnAcACP   - P.k6_2r(1).*c_C4_FabI_NADH_EnAcACP  - P.kcat6(1).*c_C4_FabI_NADH_EnAcACP;
d_C6_FabI_NADH_EnAcACP   = P.k6_2f(2).*c_FabI_NADH.*c_C6_EnAcACP   - P.k6_2r(2).*c_C6_FabI_NADH_EnAcACP  - P.kcat6(2).*c_C6_FabI_NADH_EnAcACP;
d_C8_FabI_NADH_EnAcACP   = P.k6_2f(3).*c_FabI_NADH.*c_C8_EnAcACP   - P.k6_2r(3).*c_C8_FabI_NADH_EnAcACP  - P.kcat6(3).*c_C8_FabI_NADH_EnAcACP;
d_C10_FabI_NADH_EnAcACP = P.k6_2f(4).*c_FabI_NADH.*c_C10_EnAcACP - P.k6_2r(4).*c_C10_FabI_NADH_EnAcACP - P.kcat6(4).*c_C10_FabI_NADH_EnAcACP;
d_C12_FabI_NADH_EnAcACP = P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP - P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP - P.kcat6(5).*c_C12_FabI_NADH_EnAcACP;
d_C14_FabI_NADH_EnAcACP = P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP - P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP - P.kcat6(6).*c_C14_FabI_NADH_EnAcACP;
d_C16_FabI_NADH_EnAcACP = P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP - P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP - P.kcat6(7).*c_C16_FabI_NADH_EnAcACP;
d_C18_FabI_NADH_EnAcACP = P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP - P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP - P.kcat6(8).*c_C18_FabI_NADH_EnAcACP;
d_C20_FabI_NADH_EnAcACP = P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP - P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP - P.kcat6(9).*c_C20_FabI_NADH_EnAcACP;

% C2n:1 (n=6:10) FabI-NADH-Enoyl-Acyl-ACPs
d_C12_FabI_NADH_EnAcACP_un = P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP_un - P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP_un - P.kcat6(5).*c_C12_FabI_NADH_EnAcACP_un;
d_C14_FabI_NADH_EnAcACP_un = P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP_un - P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP_un - P.kcat6(6).*c_C14_FabI_NADH_EnAcACP_un;
d_C16_FabI_NADH_EnAcACP_un = P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP_un - P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP_un - P.kcat6(7).*c_C16_FabI_NADH_EnAcACP_un;
d_C18_FabI_NADH_EnAcACP_un = P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP_un - P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP_un - P.kcat6(8).*c_C18_FabI_NADH_EnAcACP_un;
d_C20_FabI_NADH_EnAcACP_un = P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP_un - P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP_un - P.kcat6(9).*c_C20_FabI_NADH_EnAcACP_un;

% FabI-ACP
d_FabI_ACP = P.k6_inh_f.*c_FabI.*c_ACP - P.k6_inh_r.*c_FabI_ACP;


dcdt = [d_ACP;d_NADH;d_C4_EnAcACP;d_C6_EnAcACP;d_C8_EnAcACP;d_C10_EnAcACP;...
    d_C12_EnAcACP;d_C14_EnAcACP;d_C16_EnAcACP;d_C18_EnAcACP;d_C20_EnAcACP;...
    d_C12_EnAcACP_un;d_C14_EnAcACP_un;d_C16_EnAcACP_un;d_C18_EnAcACP_un;...
    d_C20_EnAcACP_un;d_C4_AcACP;d_C6_AcACP;d_C8_AcACP;d_C10_AcACP;d_C12_AcACP;...
    d_C14_AcACP;d_C16_AcACP;d_C18_AcACP;d_C20_AcACP;d_C12_AcACP_un;...
    d_C14_AcACP_un;d_C16_AcACP_un;d_C18_AcACP_un;d_C20_AcACP_un;d_FabI_NADH;...
    d_C4_FabI_NADH_EnAcACP;d_C6_FabI_NADH_EnAcACP;d_C8_FabI_NADH_EnAcACP;...
    d_C10_FabI_NADH_EnAcACP;d_C12_FabI_NADH_EnAcACP;d_C14_FabI_NADH_EnAcACP;...
    d_C16_FabI_NADH_EnAcACP;d_C18_FabI_NADH_EnAcACP;d_C20_FabI_NADH_EnAcACP;...
    d_C12_FabI_NADH_EnAcACP_un;d_C14_FabI_NADH_EnAcACP_un;...
    d_C16_FabI_NADH_EnAcACP_un;d_C18_FabI_NADH_EnAcACP_un;...
    d_C20_FabI_NADH_EnAcACP_un;d_FabI_ACP];

end
