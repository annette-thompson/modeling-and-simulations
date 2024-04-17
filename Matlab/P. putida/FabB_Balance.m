%% FabB Balance
% Need to have ran the beginning of the real code

S.labels = {'c_C2_AcCoA', 'c_ACP', 'c_CoA', 'c_C3_MalACP', 'c_C1_CO2', 'c_C4_BKeACP', 'c_C6_BKeACP',...
    'c_C8_BKeACP', 'c_C10_BKeACP', 'c_C12_BKeACP', 'c_C14_BKeACP', 'c_C16_BKeACP',...
    'c_C18_BKeACP', 'c_C20_BKeACP', 'c_C12_BKeACP_un', 'c_C14_BKeACP_un', 'c_C16_BKeACP_un',...
    'c_C18_BKeACP_un', 'c_C20_BKeACP_un', 'c_C10_cis3EnAcACP', 'c_C4_AcACP', 'c_C6_AcACP',...
    'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP', 'c_C14_AcACP', 'c_C16_AcACP', 'c_C18_AcACP',...
    'c_C12_AcACP_un', 'c_C14_AcACP_un', 'c_C16_AcACP_un', 'c_C18_AcACP_un', 'c_C4_FabB_AcACP',...
    'c_C6_FabB_AcACP', 'c_C8_FabB_AcACP', 'c_C10_FabB_AcACP', 'c_C12_FabB_AcACP',...
    'c_C14_FabB_AcACP', 'c_C16_FabB_AcACP', 'c_C18_FabB_AcACP', 'c_C12_FabB_AcACP_un',...
    'c_C14_FabB_AcACP_un', 'c_C16_FabB_AcACP_un', 'c_C18_FabB_AcACP_un', 'c_C4_FabB_Act',...
    'c_C6_FabB_Act', 'c_C8_FabB_Act', 'c_C10_FabB_Act', 'c_C12_FabB_Act', 'c_C14_FabB_Act',...
    'c_C16_FabB_Act', 'c_C18_FabB_Act', 'c_C12_FabB_Act_un', 'c_C14_FabB_Act_un',...
    'c_C16_FabB_Act_un', 'c_C18_FabB_Act_un', 'c_C7_FabB_Act_MalACP', 'c_C9_FabB_Act_MalACP',...
    'c_C11_FabB_Act_MalACP', 'c_C13_FabB_Act_MalACP', 'c_C15_FabB_Act_MalACP',...
    'c_C17_FabB_Act_MalACP', 'c_C19_FabB_Act_MalACP', 'c_C21_FabB_Act_MalACP',...
    'c_C15_FabB_Act_MalACP_un', 'c_C17_FabB_Act_MalACP_un', 'c_C19_FabB_Act_MalACP_un',...
    'c_C21_FabB_Act_MalACP_un', 'c_C10_FabB_cis3EnAcACP', 'c_C10_FabB_Act_cis3',...
    'c_C10_FabB_Act_cis3MalACP', 'c_FabB_ACP', 'c_C2_FabB_AcCoA', 'c_C2_FabB_Act',...
    'c_C5_FabB_Act_MalACP', 'c_C3_FabB_MalACP', 'c_C2_AcACP', 'c_C2_FabB_AcACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 1; % Acetyl-CoA
S.init_cond(4) = 5; % Malonyl-ACP
S.init_cond(20:32) = 0.7143; % Acyl-ACPs - evenly distribute 10 ACP to all of them (half to Malonyl-ACP)
S.init_cond(71) = 0.7143; % Acetyl-ACP

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 0 0 0 0 0 1]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabB(t,c,P);
tic
[T_FabB,C_FabB] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabB] = Calc_Function(T_FabB,C_FabB,S);

[balance_conc_FabB, balances_FabB, total_conc_FabB, carbon_FabB] = mass_balance(C_FabB,P);


function dcdt = ODE_Function_FabB(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabB
c_FabB = P.FabBtot - c_FabB_ACP - c_C2_FabB_AcCoA - c_C2_FabB_Act - c_C5_FabB_Act_MalACP - c_C3_FabB_MalACP - c_C2_FabB_AcACP...
    - c_C4_FabB_AcACP - c_C6_FabB_AcACP - c_C8_FabB_AcACP - c_C10_FabB_AcACP - c_C12_FabB_AcACP - c_C14_FabB_AcACP - c_C16_FabB_AcACP - c_C18_FabB_AcACP...
    - c_C12_FabB_AcACP_un - c_C14_FabB_AcACP_un - c_C16_FabB_AcACP_un - c_C18_FabB_AcACP_un...
    - c_C4_FabB_Act - c_C6_FabB_Act - c_C8_FabB_Act - c_C10_FabB_Act - c_C12_FabB_Act - c_C14_FabB_Act - c_C16_FabB_Act - c_C18_FabB_Act...
    - c_C12_FabB_Act_un - c_C14_FabB_Act_un - c_C16_FabB_Act_un - c_C18_FabB_Act_un...
    - c_C7_FabB_Act_MalACP - c_C9_FabB_Act_MalACP - c_C11_FabB_Act_MalACP - c_C13_FabB_Act_MalACP - c_C15_FabB_Act_MalACP - c_C17_FabB_Act_MalACP - c_C19_FabB_Act_MalACP - c_C21_FabB_Act_MalACP...
    - c_C15_FabB_Act_MalACP_un - c_C17_FabB_Act_MalACP_un - c_C19_FabB_Act_MalACP_un - c_C21_FabB_Act_MalACP_un...
    - c_C10_FabB_cis3EnAcACP - c_C10_FabB_Act_cis3 - c_C10_FabB_Act_cis3MalACP;

% Set of differential equations
% C2n (n=1:9)-CoA % changed 
d_C2_AcCoA = 0*P.k10_4r.*c_C2_FabB_AcCoA - 0*P.k10_4f.*c_FabB.*c_C2_AcCoA;

% ACP
d_ACP = 0*P.k10_2f(1).*c_C4_FabB_AcACP   - 0*P.k10_2r(1).*c_C4_FabB_Act.*c_ACP...
    + 0*P.k10_2f(2).*c_C6_FabB_AcACP   - 0*P.k10_2r(2).*c_C6_FabB_Act.*c_ACP...
    + 0*P.k10_2f(3).*c_C8_FabB_AcACP   - 0*P.k10_2r(3).*c_C8_FabB_Act.*c_ACP...
    + 0*P.k10_2f(4).*c_C10_FabB_AcACP - 0*P.k10_2r(4).*c_C10_FabB_Act.*c_ACP...
    + 0*P.k10_2f(5).*c_C12_FabB_AcACP - 0*P.k10_2r(5).*c_C12_FabB_Act.*c_ACP...
    + 0*P.k10_2f(6).*c_C14_FabB_AcACP - 0*P.k10_2r(6).*c_C14_FabB_Act.*c_ACP...
    + 0*P.k10_2f(7).*c_C16_FabB_AcACP - 0*P.k10_2r(7).*c_C16_FabB_Act.*c_ACP...
    + 0*P.k10_2f(8).*c_C18_FabB_AcACP - 0*P.k10_2r(8).*c_C18_FabB_Act.*c_ACP...
    + 0*P.k10_2f(5).*c_C12_FabB_AcACP_un - 0*P.k10_2r(5).*c_C12_FabB_Act_un.*c_ACP...
    + 0*P.k10_2f(6).*c_C14_FabB_AcACP_un - 0*P.k10_2r(6).*c_C14_FabB_Act_un.*c_ACP...
    + 0*P.k10_2f(7).*c_C16_FabB_AcACP_un - 0*P.k10_2r(7).*c_C16_FabB_Act_un.*c_ACP...
    + 0*P.k10_2f(8).*c_C18_FabB_AcACP_un - 0*P.k10_2r(8).*c_C18_FabB_Act_un.*c_ACP...
    + 0*P.k10_2f(4).*c_C10_FabB_cis3EnAcACP - 0*P.k10_2r(4).*c_C10_FabB_Act_cis3.*c_ACP...
    + 0*P.k10_9f.*c_C2_FabB_AcACP - 0*P.k10_9r.*c_C2_FabB_Act.*c_ACP; 

% CoA % changed
d_CoA = 0*P.k10_5f.*c_C2_FabB_AcCoA - 0*P.k10_5r.*c_C2_FabB_Act.*c_CoA;

% Malonyl-ACP % changed
d_C3_MalACP = P.k10_3r(1).*c_C7_FabB_Act_MalACP - P.k10_3f(1).*c_C4_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(2).*c_C9_FabB_Act_MalACP   - P.k10_3f(2).*c_C6_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(3).*c_C11_FabB_Act_MalACP   - P.k10_3f(3).*c_C8_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(4).*c_C13_FabB_Act_MalACP - P.k10_3f(4).*c_C10_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(5).*c_C15_FabB_Act_MalACP - P.k10_3f(5).*c_C12_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(6).*c_C17_FabB_Act_MalACP - P.k10_3f(6).*c_C14_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(7).*c_C19_FabB_Act_MalACP - P.k10_3f(7).*c_C16_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(8).*c_C21_FabB_Act_MalACP - P.k10_3f(8).*c_C18_FabB_Act.*c_C3_MalACP...
    + P.k10_3r(5).*c_C15_FabB_Act_MalACP_un - P.k10_3f(5).*c_C12_FabB_Act_un.*c_C3_MalACP...
    + P.k10_3r(6).*c_C17_FabB_Act_MalACP_un - P.k10_3f(6).*c_C14_FabB_Act_un.*c_C3_MalACP...
    + P.k10_3r(7).*c_C19_FabB_Act_MalACP_un - P.k10_3f(7).*c_C16_FabB_Act_un.*c_C3_MalACP...
    + P.k10_3r(8).*c_C21_FabB_Act_MalACP_un - P.k10_3f(8).*c_C18_FabB_Act_un.*c_C3_MalACP...
    + P.k10_3r(4).*c_C10_FabB_Act_cis3MalACP - P.k10_3f(4).*c_C10_FabB_Act_cis3.*c_C3_MalACP...
    + 0*P.k10_7r.*c_C3_FabB_MalACP - 0*P.k10_7f.*c_FabB.*c_C3_MalACP;

% CO2 % changed
d_C1_CO2 = P.kcat10(1).*c_C7_FabB_Act_MalACP + P.kcat10(2).*c_C9_FabB_Act_MalACP + P.kcat10(3).*c_C11_FabB_Act_MalACP + P.kcat10(4).*c_C13_FabB_Act_MalACP...
    + P.kcat10(5).*c_C15_FabB_Act_MalACP + P.kcat10(6).*c_C17_FabB_Act_MalACP + P.kcat10(7).*c_C19_FabB_Act_MalACP + P.kcat10(8).*c_C21_FabB_Act_MalACP...
    + 0*P.kcat10_un(4).*c_C10_FabB_Act_cis3MalACP + 0*P.kcat10_un(5).*c_C15_FabB_Act_MalACP_un + 0*P.kcat10_un(6).*c_C17_FabB_Act_MalACP_un + 0*P.kcat10_un(7).*c_C19_FabB_Act_MalACP_un + 0*P.kcat10_un(8).*c_C21_FabB_Act_MalACP_un...
    + 0*P.kcat10_H.*c_C5_FabB_Act_MalACP + 0*P.kcat10_CO2.*c_C3_FabB_MalACP;

% C2n (n=2:10) B-ketoacyl-ACPs (FabH + FabF + FabB - FabG) % changed
d_C4_BKeACP    = 0*P.kcat10_H.*c_C5_FabB_Act_MalACP;
d_C6_BKeACP    = P.kcat10(1).*c_C7_FabB_Act_MalACP;
d_C8_BKeACP    = P.kcat10(2).*c_C9_FabB_Act_MalACP;
d_C10_BKeACP  = P.kcat10(3).*c_C11_FabB_Act_MalACP;
d_C12_BKeACP  = P.kcat10(4).*c_C13_FabB_Act_MalACP;
d_C14_BKeACP  = P.kcat10(5).*c_C15_FabB_Act_MalACP;
d_C16_BKeACP  = P.kcat10(6).*c_C17_FabB_Act_MalACP;
d_C18_BKeACP  = P.kcat10(7).*c_C19_FabB_Act_MalACP;
d_C20_BKeACP  = P.kcat10(8).*c_C21_FabB_Act_MalACP;

% C2n:1 (n=6:10) B-ketoacyl-ACPs (FabF + FabB - FabG)
d_C12_BKeACP_un = 0*P.kcat10_un(4).*c_C10_FabB_Act_cis3MalACP;
d_C14_BKeACP_un = 0*P.kcat10_un(5).*c_C15_FabB_Act_MalACP_un;
d_C16_BKeACP_un = 0*P.kcat10_un(6).*c_C17_FabB_Act_MalACP_un;
d_C18_BKeACP_un = 0*P.kcat10_un(7).*c_C19_FabB_Act_MalACP_un;
d_C20_BKeACP_un = 0*P.kcat10_un(8).*c_C21_FabB_Act_MalACP_un;

% C10 cis-3-Enoyl-Acyl-ACP (FabA - FabB)
d_C10_cis3EnAcACP = 0*P.k10_1r(4).*c_C10_FabB_cis3EnAcACP - 0*P.k10_1f(4).*c_FabB.*c_C10_cis3EnAcACP;

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C4_AcACP   = 0*P.k10_1r(1).*c_C4_FabB_AcACP   - 0*P.k10_1f(1).*c_FabB.*c_C4_AcACP;
d_C6_AcACP   = 0*P.k10_1r(2).*c_C6_FabB_AcACP   - 0*P.k10_1f(2).*c_FabB.*c_C6_AcACP;
d_C8_AcACP   = 0*P.k10_1r(3).*c_C8_FabB_AcACP   - 0*P.k10_1f(3).*c_FabB.*c_C8_AcACP;
d_C10_AcACP = 0*P.k10_1r(4).*c_C10_FabB_AcACP - 0*P.k10_1f(4).*c_FabB.*c_C10_AcACP;
d_C12_AcACP = 0*P.k10_1r(5).*c_C12_FabB_AcACP - 0*P.k10_1f(5).*c_FabB.*c_C12_AcACP;
d_C14_AcACP = 0*P.k10_1r(6).*c_C14_FabB_AcACP - 0*P.k10_1f(6).*c_FabB.*c_C14_AcACP;
d_C16_AcACP = 0*P.k10_1r(7).*c_C16_FabB_AcACP - 0*P.k10_1f(7).*c_FabB.*c_C16_AcACP;
d_C18_AcACP = 0*P.k10_1r(8).*c_C18_FabB_AcACP - 0*P.k10_1f(8).*c_FabB.*c_C18_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C12_AcACP_un = 0*P.k10_1r(5).*c_C12_FabB_AcACP_un - 0*P.k10_1f(5).*c_FabB.*c_C12_AcACP_un;
d_C14_AcACP_un = 0*P.k10_1r(6).*c_C14_FabB_AcACP_un - 0*P.k10_1f(6).*c_FabB.*c_C14_AcACP_un;
d_C16_AcACP_un = 0*P.k10_1r(7).*c_C16_FabB_AcACP_un - 0*P.k10_1f(7).*c_FabB.*c_C16_AcACP_un;
d_C18_AcACP_un = 0*P.k10_1r(8).*c_C18_FabB_AcACP_un - 0*P.k10_1f(8).*c_FabB.*c_C18_AcACP_un;

% C2n (n=2:9) FabB-Acyl-ACPs
d_C4_FabB_AcACP   = 0*P.k10_1f(1).*c_FabB.*c_C4_AcACP   - 0*P.k10_1r(1).*c_C4_FabB_AcACP  + 0*P.k10_2r(1).*c_C4_FabB_Act.*c_ACP   - 0*P.k10_2f(1).*c_C4_FabB_AcACP;
d_C6_FabB_AcACP   = 0*P.k10_1f(2).*c_FabB.*c_C6_AcACP   - 0*P.k10_1r(2).*c_C6_FabB_AcACP  + 0*P.k10_2r(2).*c_C6_FabB_Act.*c_ACP   - 0*P.k10_2f(2).*c_C6_FabB_AcACP;
d_C8_FabB_AcACP   = 0*P.k10_1f(3).*c_FabB.*c_C8_AcACP   - 0*P.k10_1r(3).*c_C8_FabB_AcACP  + 0*P.k10_2r(3).*c_C8_FabB_Act.*c_ACP   - 0*P.k10_2f(3).*c_C8_FabB_AcACP;
d_C10_FabB_AcACP = 0*P.k10_1f(4).*c_FabB.*c_C10_AcACP - 0*P.k10_1r(4).*c_C10_FabB_AcACP + 0*P.k10_2r(4).*c_C10_FabB_Act.*c_ACP - 0*P.k10_2f(4).*c_C10_FabB_AcACP;
d_C12_FabB_AcACP = 0*P.k10_1f(5).*c_FabB.*c_C12_AcACP - 0*P.k10_1r(5).*c_C12_FabB_AcACP + 0*P.k10_2r(5).*c_C12_FabB_Act.*c_ACP - 0*P.k10_2f(5).*c_C12_FabB_AcACP;
d_C14_FabB_AcACP = 0*P.k10_1f(6).*c_FabB.*c_C14_AcACP - 0*P.k10_1r(6).*c_C14_FabB_AcACP + 0*P.k10_2r(6).*c_C14_FabB_Act.*c_ACP - 0*P.k10_2f(6).*c_C14_FabB_AcACP;
d_C16_FabB_AcACP = 0*P.k10_1f(7).*c_FabB.*c_C16_AcACP - 0*P.k10_1r(7).*c_C16_FabB_AcACP + 0*P.k10_2r(7).*c_C16_FabB_Act.*c_ACP - 0*P.k10_2f(7).*c_C16_FabB_AcACP;
d_C18_FabB_AcACP = 0*P.k10_1f(8).*c_FabB.*c_C18_AcACP - 0*P.k10_1r(8).*c_C18_FabB_AcACP + 0*P.k10_2r(8).*c_C18_FabB_Act.*c_ACP - 0*P.k10_2f(8).*c_C18_FabB_AcACP;

% C2n:1 (n=6:9) FabB-Acyl-ACPs
d_C12_FabB_AcACP_un = 0*P.k10_1f(5).*c_FabB.*c_C12_AcACP_un - 0*P.k10_1r(5).*c_C12_FabB_AcACP_un + 0*P.k10_2r(5).*c_C12_FabB_Act_un.*c_ACP - 0*P.k10_2f(5).*c_C12_FabB_AcACP_un;
d_C14_FabB_AcACP_un = 0*P.k10_1f(6).*c_FabB.*c_C14_AcACP_un - 0*P.k10_1r(6).*c_C14_FabB_AcACP_un + 0*P.k10_2r(6).*c_C14_FabB_Act_un.*c_ACP - 0*P.k10_2f(6).*c_C14_FabB_AcACP_un;
d_C16_FabB_AcACP_un = 0*P.k10_1f(7).*c_FabB.*c_C16_AcACP_un - 0*P.k10_1r(7).*c_C16_FabB_AcACP_un + 0*P.k10_2r(7).*c_C16_FabB_Act_un.*c_ACP - 0*P.k10_2f(7).*c_C16_FabB_AcACP_un;
d_C18_FabB_AcACP_un = 0*P.k10_1f(8).*c_FabB.*c_C18_AcACP_un - 0*P.k10_1r(8).*c_C18_FabB_AcACP_un + 0*P.k10_2r(8).*c_C18_FabB_Act_un.*c_ACP - 0*P.k10_2f(8).*c_C18_FabB_AcACP_un;

% C2n (n=2:9) FabB*
d_C4_FabB_Act   = 0*P.k10_2f(1).*c_C4_FabB_AcACP   - 0*P.k10_2r(1).*c_C4_FabB_Act.*c_ACP  + P.k10_3r(1).*c_C7_FabB_Act_MalACP   - P.k10_3f(1).*c_C4_FabB_Act.*c_C3_MalACP;
d_C6_FabB_Act   = 0*P.k10_2f(2).*c_C6_FabB_AcACP   - 0*P.k10_2r(2).*c_C6_FabB_Act.*c_ACP  + P.k10_3r(2).*c_C9_FabB_Act_MalACP   - P.k10_3f(2).*c_C6_FabB_Act.*c_C3_MalACP;
d_C8_FabB_Act   = 0*P.k10_2f(3).*c_C8_FabB_AcACP   - 0*P.k10_2r(3).*c_C8_FabB_Act.*c_ACP  + P.k10_3r(3).*c_C11_FabB_Act_MalACP - P.k10_3f(3).*c_C8_FabB_Act.*c_C3_MalACP;
d_C10_FabB_Act = 0*P.k10_2f(4).*c_C10_FabB_AcACP - 0*P.k10_2r(4).*c_C10_FabB_Act.*c_ACP + P.k10_3r(4).*c_C13_FabB_Act_MalACP - P.k10_3f(4).*c_C10_FabB_Act.*c_C3_MalACP;
d_C12_FabB_Act = 0*P.k10_2f(5).*c_C12_FabB_AcACP - 0*P.k10_2r(5).*c_C12_FabB_Act.*c_ACP + P.k10_3r(5).*c_C15_FabB_Act_MalACP - P.k10_3f(5).*c_C12_FabB_Act.*c_C3_MalACP;
d_C14_FabB_Act = 0*P.k10_2f(6).*c_C14_FabB_AcACP - 0*P.k10_2r(6).*c_C14_FabB_Act.*c_ACP + P.k10_3r(6).*c_C17_FabB_Act_MalACP - P.k10_3f(6).*c_C14_FabB_Act.*c_C3_MalACP;
d_C16_FabB_Act = 0*P.k10_2f(7).*c_C16_FabB_AcACP - 0*P.k10_2r(7).*c_C16_FabB_Act.*c_ACP + P.k10_3r(7).*c_C19_FabB_Act_MalACP - P.k10_3f(7).*c_C16_FabB_Act.*c_C3_MalACP;
d_C18_FabB_Act = 0*P.k10_2f(8).*c_C18_FabB_AcACP - 0*P.k10_2r(8).*c_C18_FabB_Act.*c_ACP + P.k10_3r(8).*c_C21_FabB_Act_MalACP - P.k10_3f(8).*c_C18_FabB_Act.*c_C3_MalACP;

% C2n:1 (n=6:9) FabB*
d_C12_FabB_Act_un = 0*P.k10_2f(5).*c_C12_FabB_AcACP_un - 0*P.k10_2r(5).*c_C12_FabB_Act_un.*c_ACP + P.k10_3r(5).*c_C15_FabB_Act_MalACP_un - P.k10_3f(5).*c_C12_FabB_Act_un.*c_C3_MalACP;
d_C14_FabB_Act_un = 0*P.k10_2f(6).*c_C14_FabB_AcACP_un - 0*P.k10_2r(6).*c_C14_FabB_Act_un.*c_ACP + P.k10_3r(6).*c_C17_FabB_Act_MalACP_un - P.k10_3f(6).*c_C14_FabB_Act_un.*c_C3_MalACP;
d_C16_FabB_Act_un = 0*P.k10_2f(7).*c_C16_FabB_AcACP_un - 0*P.k10_2r(7).*c_C16_FabB_Act_un.*c_ACP + P.k10_3r(7).*c_C19_FabB_Act_MalACP_un - P.k10_3f(7).*c_C16_FabB_Act_un.*c_C3_MalACP;
d_C18_FabB_Act_un = 0*P.k10_2f(8).*c_C18_FabB_AcACP_un - 0*P.k10_2r(8).*c_C18_FabB_Act_un.*c_ACP + P.k10_3r(8).*c_C21_FabB_Act_MalACP_un - P.k10_3f(8).*c_C18_FabB_Act_un.*c_C3_MalACP;

% C2n (n=2:9) FabB*-Malonyl-ACPs
d_C7_FabB_Act_MalACP   = P.k10_3f(1).*c_C4_FabB_Act.*c_C3_MalACP   - P.k10_3r(1).*c_C7_FabB_Act_MalACP  - P.kcat10(1).*c_C7_FabB_Act_MalACP;
d_C9_FabB_Act_MalACP   = P.k10_3f(2).*c_C6_FabB_Act.*c_C3_MalACP   - P.k10_3r(2).*c_C9_FabB_Act_MalACP  - P.kcat10(2).*c_C9_FabB_Act_MalACP;
d_C11_FabB_Act_MalACP = P.k10_3f(3).*c_C8_FabB_Act.*c_C3_MalACP   - P.k10_3r(3).*c_C11_FabB_Act_MalACP  - P.kcat10(3).*c_C11_FabB_Act_MalACP;
d_C13_FabB_Act_MalACP = P.k10_3f(4).*c_C10_FabB_Act.*c_C3_MalACP - P.k10_3r(4).*c_C13_FabB_Act_MalACP - P.kcat10(4).*c_C13_FabB_Act_MalACP;
d_C15_FabB_Act_MalACP = P.k10_3f(5).*c_C12_FabB_Act.*c_C3_MalACP - P.k10_3r(5).*c_C15_FabB_Act_MalACP - P.kcat10(5).*c_C15_FabB_Act_MalACP;
d_C17_FabB_Act_MalACP = P.k10_3f(6).*c_C14_FabB_Act.*c_C3_MalACP - P.k10_3r(6).*c_C17_FabB_Act_MalACP - P.kcat10(6).*c_C17_FabB_Act_MalACP;
d_C19_FabB_Act_MalACP = P.k10_3f(7).*c_C16_FabB_Act.*c_C3_MalACP - P.k10_3r(7).*c_C19_FabB_Act_MalACP - P.kcat10(7).*c_C19_FabB_Act_MalACP;
d_C21_FabB_Act_MalACP = P.k10_3f(8).*c_C18_FabB_Act.*c_C3_MalACP - P.k10_3r(8).*c_C21_FabB_Act_MalACP - P.kcat10(8).*c_C21_FabB_Act_MalACP;

% C2n:1 (n=6:9) FabB*-Malonyl-ACPs
d_C15_FabB_Act_MalACP_un = P.k10_3f(5).*c_C12_FabB_Act_un.*c_C3_MalACP - P.k10_3r(5).*c_C15_FabB_Act_MalACP_un - 0*P.kcat10_un(5).*c_C15_FabB_Act_MalACP_un;
d_C17_FabB_Act_MalACP_un = P.k10_3f(6).*c_C14_FabB_Act_un.*c_C3_MalACP - P.k10_3r(6).*c_C17_FabB_Act_MalACP_un - 0*P.kcat10_un(6).*c_C17_FabB_Act_MalACP_un;
d_C19_FabB_Act_MalACP_un = P.k10_3f(7).*c_C16_FabB_Act_un.*c_C3_MalACP - P.k10_3r(7).*c_C19_FabB_Act_MalACP_un - 0*P.kcat10_un(7).*c_C19_FabB_Act_MalACP_un;
d_C21_FabB_Act_MalACP_un = P.k10_3f(8).*c_C18_FabB_Act_un.*c_C3_MalACP - P.k10_3r(8).*c_C21_FabB_Act_MalACP_un - 0*P.kcat10_un(8).*c_C21_FabB_Act_MalACP_un;

% FabB-(C10 cis-3-Enoyl-Acyl-ACP)
d_C10_FabB_cis3EnAcACP = 0*P.k10_1f(4).*c_FabB.*c_C10_cis3EnAcACP - 0*P.k10_1r(4).*c_C10_FabB_cis3EnAcACP + 0*P.k10_2r(4).*c_C10_FabB_Act_cis3.*c_ACP - 0*P.k10_2f(4).*c_C10_FabB_cis3EnAcACP;

% C10 cis-3-FabB*
d_C10_FabB_Act_cis3 = 0*P.k10_2f(4).*c_C10_FabB_cis3EnAcACP - 0*P.k10_2r(4).*c_C10_FabB_Act_cis3.*c_ACP + P.k10_3r(4).*c_C10_FabB_Act_cis3MalACP - P.k10_3f(4).*c_C10_FabB_Act_cis3.*c_C3_MalACP;

% C10 cis-3-FabB*-Malonyl-ACP
d_C10_FabB_Act_cis3MalACP = P.k10_3f(4).*c_C10_FabB_Act_cis3.*c_C3_MalACP - P.k10_3r(4).*c_C10_FabB_Act_cis3MalACP - 0*P.kcat10_un(4).*c_C10_FabB_Act_cis3MalACP;

% FabB-ACP
d_FabB_ACP = P.k10_inh_f.*c_FabB.*c_ACP - P.k10_inh_r.*c_FabB_ACP;

% Giving FabB FabH-like activity
% FabB-Acetyl-CoA
d_C2_FabB_AcCoA = 0*P.k10_4f.*c_FabB.*c_C2_AcCoA - 0*P.k10_4r.*c_C2_FabB_AcCoA + 0*P.k10_5r.*c_C2_FabB_Act.*c_CoA - 0*P.k10_5f.*c_C2_FabB_AcCoA;

% FabB*
d_C2_FabB_Act = 0*P.k10_5f.*c_C2_FabB_AcCoA - 0*P.k10_5r.*c_C2_FabB_Act.*c_CoA + 0*P.k10_6r.*c_C5_FabB_Act_MalACP - 0*P.k10_6f.*c_C2_FabB_Act.*c_C3_MalACP + 0*P.k10_9f.*c_C2_FabB_AcACP - 0*P.k10_9r.*c_C2_FabB_Act.*c_ACP;

% FabB*-Malonyl-ACP
d_C5_FabB_Act_MalACP = 0*P.k10_6f.*c_C2_FabB_Act.*c_C3_MalACP - 0*P.k10_6r.*c_C5_FabB_Act_MalACP - 0*P.kcat10_H.*c_C5_FabB_Act_MalACP; 

% FabB and FabF decarboxylating mACP to form aACP and reacting with it to form activated enzyme (initiation)
% FabB-Malonyl-ACP
d_C3_FabB_MalACP = 0*P.k10_7f.*c_FabB.*c_C3_MalACP - 0*P.k10_7r.*c_C3_FabB_MalACP - 0*P.kcat10_CO2.*c_C3_FabB_MalACP;

% Acetyl-ACP
d_C2_AcACP = 0*P.k10_8r.*c_C2_FabB_AcACP - 0*P.k10_8f.*c_FabB.*c_C2_AcACP; 

% FabB-Acetyl-ACP
d_C2_FabB_AcACP = 0*P.kcat10_CO2.*c_C3_FabB_MalACP + 0*P.k10_8f.*c_FabB.*c_C2_AcACP - 0*P.k10_8r.*c_C2_FabB_AcACP + 0*P.k10_9r.*c_C2_FabB_Act.*c_ACP - 0*P.k10_9f.*c_C2_FabB_AcACP; 


dcdt = [d_C2_AcCoA;d_ACP;d_CoA;d_C3_MalACP;d_C1_CO2;d_C4_BKeACP;d_C6_BKeACP;...
    d_C8_BKeACP;d_C10_BKeACP;d_C12_BKeACP;d_C14_BKeACP;d_C16_BKeACP;...
    d_C18_BKeACP;d_C20_BKeACP;d_C12_BKeACP_un;d_C14_BKeACP_un;d_C16_BKeACP_un;...
    d_C18_BKeACP_un;d_C20_BKeACP_un;d_C10_cis3EnAcACP;d_C4_AcACP;d_C6_AcACP;...
    d_C8_AcACP;d_C10_AcACP;d_C12_AcACP;d_C14_AcACP;d_C16_AcACP;d_C18_AcACP;...
    d_C12_AcACP_un;d_C14_AcACP_un;d_C16_AcACP_un;d_C18_AcACP_un;d_C4_FabB_AcACP;...
    d_C6_FabB_AcACP;d_C8_FabB_AcACP;d_C10_FabB_AcACP;d_C12_FabB_AcACP;...
    d_C14_FabB_AcACP;d_C16_FabB_AcACP;d_C18_FabB_AcACP;d_C12_FabB_AcACP_un;...
    d_C14_FabB_AcACP_un;d_C16_FabB_AcACP_un;d_C18_FabB_AcACP_un;d_C4_FabB_Act;...
    d_C6_FabB_Act;d_C8_FabB_Act;d_C10_FabB_Act;d_C12_FabB_Act;d_C14_FabB_Act;...
    d_C16_FabB_Act;d_C18_FabB_Act;d_C12_FabB_Act_un;d_C14_FabB_Act_un;...
    d_C16_FabB_Act_un;d_C18_FabB_Act_un;d_C7_FabB_Act_MalACP;d_C9_FabB_Act_MalACP;...
    d_C11_FabB_Act_MalACP;d_C13_FabB_Act_MalACP;d_C15_FabB_Act_MalACP;...
    d_C17_FabB_Act_MalACP;d_C19_FabB_Act_MalACP;d_C21_FabB_Act_MalACP;...
    d_C15_FabB_Act_MalACP_un;d_C17_FabB_Act_MalACP_un;d_C19_FabB_Act_MalACP_un;...
    d_C21_FabB_Act_MalACP_un;d_C10_FabB_cis3EnAcACP;d_C10_FabB_Act_cis3;...
    d_C10_FabB_Act_cis3MalACP;d_FabB_ACP;d_C2_FabB_AcCoA;d_C2_FabB_Act;...
    d_C5_FabB_Act_MalACP;d_C3_FabB_MalACP;d_C2_AcACP;d_C2_FabB_AcACP];

end
