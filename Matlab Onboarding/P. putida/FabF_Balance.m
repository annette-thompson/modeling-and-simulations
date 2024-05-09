%% FabF Balance
% Need to have ran the beginning of the real code

S.labels = {'c_C2_AcCoA', 'c_ACP', 'c_CoA', 'c_C3_MalACP', 'c_C1_CO2', 'c_C4_BKeACP', 'c_C6_BKeACP',...
    'c_C8_BKeACP', 'c_C10_BKeACP', 'c_C12_BKeACP', 'c_C14_BKeACP', 'c_C16_BKeACP',...
    'c_C18_BKeACP', 'c_C20_BKeACP', 'c_C14_BKeACP_un', 'c_C16_BKeACP_un', 'c_C18_BKeACP_un',...
    'c_C20_BKeACP_un', 'c_C4_AcACP', 'c_C6_AcACP', 'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP',...
    'c_C14_AcACP', 'c_C16_AcACP', 'c_C18_AcACP', 'c_C12_AcACP_un', 'c_C14_AcACP_un',...
    'c_C16_AcACP_un', 'c_C18_AcACP_un', 'c_C4_FabF_AcACP', 'c_C6_FabF_AcACP',...
    'c_C8_FabF_AcACP', 'c_C10_FabF_AcACP', 'c_C12_FabF_AcACP', 'c_C14_FabF_AcACP',...
    'c_C16_FabF_AcACP', 'c_C18_FabF_AcACP', 'c_C12_FabF_AcACP_un', 'c_C14_FabF_AcACP_un',...
    'c_C16_FabF_AcACP_un', 'c_C18_FabF_AcACP_un', 'c_C4_FabF_Act', 'c_C6_FabF_Act',...
    'c_C8_FabF_Act', 'c_C10_FabF_Act', 'c_C12_FabF_Act', 'c_C14_FabF_Act', 'c_C16_FabF_Act',...
    'c_C18_FabF_Act', 'c_C12_FabF_Act_un', 'c_C14_FabF_Act_un', 'c_C16_FabF_Act_un',...
    'c_C18_FabF_Act_un', 'c_C7_FabF_Act_MalACP', 'c_C9_FabF_Act_MalACP',...
    'c_C11_FabF_Act_MalACP', 'c_C13_FabF_Act_MalACP', 'c_C15_FabF_Act_MalACP',...
    'c_C17_FabF_Act_MalACP', 'c_C19_FabF_Act_MalACP', 'c_C21_FabF_Act_MalACP',...
    'c_C15_FabF_Act_MalACP_un', 'c_C17_FabF_Act_MalACP_un', 'c_C19_FabF_Act_MalACP_un',...
    'c_C21_FabF_Act_MalACP_un', 'c_FabF_ACP', 'c_C2_FabF_AcCoA', 'c_C2_FabF_Act',...
    'c_C5_FabF_Act_MalACP', 'c_C2_AcACP', 'c_C3_FabF_MalACP', 'c_C2_FabF_AcACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(4) = 5; % Malonyl-ACP
S.init_cond(19:30) = 0.3846; % Acyl-ACPs - evenly distribute 10 ACP to all of them (half to Malonyl-ACP)
S.init_cond(71) = 0.3846; % Acetyl-ACP

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 0 0 0 1 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabF(t,c,P);
tic
[T_FabF,C_FabF] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabF] = Calc_Function(T_FabF,C_FabF,S);

[balance_conc_FabF, balances_FabF, total_conc_FabF, carbon_FabF] = mass_balance(C_FabF,P);


function dcdt = ODE_Function_FabF(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabF
c_FabF = P.FabFtot - c_FabF_ACP - c_C2_FabF_AcCoA - c_C2_FabF_Act - c_C3_FabF_MalACP - c_C2_FabF_AcACP...
    - c_C4_FabF_AcACP - c_C6_FabF_AcACP - c_C8_FabF_AcACP - c_C10_FabF_AcACP - c_C12_FabF_AcACP - c_C14_FabF_AcACP - c_C16_FabF_AcACP - c_C18_FabF_AcACP...
    - c_C12_FabF_AcACP_un - c_C14_FabF_AcACP_un - c_C16_FabF_AcACP_un - c_C18_FabF_AcACP_un...
    - c_C4_FabF_Act - c_C6_FabF_Act - c_C8_FabF_Act - c_C10_FabF_Act - c_C12_FabF_Act - c_C14_FabF_Act - c_C16_FabF_Act - c_C18_FabF_Act...
    - c_C12_FabF_Act_un - c_C14_FabF_Act_un - c_C16_FabF_Act_un - c_C18_FabF_Act_un...
    - c_C5_FabF_Act_MalACP - c_C7_FabF_Act_MalACP - c_C9_FabF_Act_MalACP - c_C11_FabF_Act_MalACP - c_C13_FabF_Act_MalACP - c_C15_FabF_Act_MalACP - c_C17_FabF_Act_MalACP - c_C19_FabF_Act_MalACP - c_C21_FabF_Act_MalACP...
    - c_C15_FabF_Act_MalACP_un - c_C17_FabF_Act_MalACP_un - c_C19_FabF_Act_MalACP_un - c_C21_FabF_Act_MalACP_un;

% Set of differential equations
% C2n (n=1:9)-CoA % changed 
d_C2_AcCoA = P.k8_4r.*c_C2_FabF_AcCoA - P.k8_4f.*c_FabF.*c_C2_AcCoA;

d_ACP = P.k8_2f(1).*c_C4_FabF_AcACP    - P.k8_2r(1).*c_C4_FabF_Act.*c_ACP...
    + P.k8_2f(2).*c_C6_FabF_AcACP    - P.k8_2r(2).*c_C6_FabF_Act.*c_ACP...
    + P.k8_2f(3).*c_C8_FabF_AcACP    - P.k8_2r(3).*c_C8_FabF_Act.*c_ACP...
    + P.k8_2f(4).*c_C10_FabF_AcACP  - P.k8_2r(4).*c_C10_FabF_Act.*c_ACP...
    + P.k8_2f(5).*c_C12_FabF_AcACP  - P.k8_2r(5).*c_C12_FabF_Act.*c_ACP...
    + P.k8_2f(6).*c_C14_FabF_AcACP  - P.k8_2r(6).*c_C14_FabF_Act.*c_ACP...
    + P.k8_2f(7).*c_C16_FabF_AcACP  - P.k8_2r(7).*c_C16_FabF_Act.*c_ACP...
    + P.k8_2f(8).*c_C18_FabF_AcACP  - P.k8_2r(8).*c_C18_FabF_Act.*c_ACP...
    + P.k8_2f(5).*c_C12_FabF_AcACP_un  - P.k8_2r(5).*c_C12_FabF_Act_un.*c_ACP...
    + P.k8_2f(6).*c_C14_FabF_AcACP_un  - P.k8_2r(6).*c_C14_FabF_Act_un.*c_ACP...
    + P.k8_2f(7).*c_C16_FabF_AcACP_un  - P.k8_2r(7).*c_C16_FabF_Act_un.*c_ACP...
    + P.k8_2f(8).*c_C18_FabF_AcACP_un  - P.k8_2r(8).*c_C18_FabF_Act_un.*c_ACP...
    + P.k8_inh_r.*c_FabF_ACP   - P.k8_inh_f.*c_FabF.*c_ACP...
    + P.k8_9f.*c_C2_FabF_AcACP   - P.k8_9r.*c_C2_FabF_Act.*c_ACP; 

% CoA % changed
d_CoA = P.k8_5f.*c_C2_FabF_AcCoA   - P.k8_5r*c_C2_FabF_Act.*c_CoA;

% Malonyl-ACP % changed
d_C3_MalACP = P.k8_3r(1).*c_C7_FabF_Act_MalACP   - P.k8_3f(1).*c_C4_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(2).*c_C9_FabF_Act_MalACP   - P.k8_3f(2).*c_C6_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(3).*c_C11_FabF_Act_MalACP   - P.k8_3f(3).*c_C8_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(4).*c_C13_FabF_Act_MalACP - P.k8_3f(4).*c_C10_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(5).*c_C15_FabF_Act_MalACP - P.k8_3f(5).*c_C12_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(6).*c_C17_FabF_Act_MalACP - P.k8_3f(6).*c_C14_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(7).*c_C19_FabF_Act_MalACP - P.k8_3f(7).*c_C16_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(8).*c_C21_FabF_Act_MalACP - P.k8_3f(8).*c_C18_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(5).*c_C15_FabF_Act_MalACP_un - P.k8_3f(5).*c_C12_FabF_Act_un.*c_C3_MalACP...
    + P.k8_3r(6).*c_C17_FabF_Act_MalACP_un - P.k8_3f(6).*c_C14_FabF_Act_un.*c_C3_MalACP...
    + P.k8_3r(7).*c_C19_FabF_Act_MalACP_un - P.k8_3f(7).*c_C16_FabF_Act_un.*c_C3_MalACP...
    + P.k8_3r(8).*c_C21_FabF_Act_MalACP_un - P.k8_3f(8).*c_C18_FabF_Act_un.*c_C3_MalACP...
    + P.k8_6r.*c_C5_FabF_Act_MalACP   - P.k8_6f.*c_C2_FabF_Act.*c_C3_MalACP...
    + P.k8_7r.*c_C3_FabF_MalACP   - P.k8_7f.*c_FabF.*c_C3_MalACP;

% CO2 % changed
d_C1_CO2 = P.kcat8(1).*c_C7_FabF_Act_MalACP + P.kcat8(2).*c_C9_FabF_Act_MalACP + P.kcat8(3).*c_C11_FabF_Act_MalACP + P.kcat8(4).*c_C13_FabF_Act_MalACP...
    + P.kcat8(5).*c_C15_FabF_Act_MalACP + P.kcat8(6).*c_C17_FabF_Act_MalACP + P.kcat8(7).*c_C19_FabF_Act_MalACP + P.kcat8(8).*c_C21_FabF_Act_MalACP...
    + P.kcat8_un(5).*c_C15_FabF_Act_MalACP_un + P.kcat8_un(6).*c_C17_FabF_Act_MalACP_un + P.kcat8_un(7).*c_C19_FabF_Act_MalACP_un + P.kcat8_un(8).*c_C21_FabF_Act_MalACP_un...
    + P.kcat8_H.*c_C5_FabF_Act_MalACP + P.kcat8_CO2.*c_C3_FabF_MalACP;

% C2n (n=2:10) B-ketoacyl-ACPs (FabH + FabF + FabB - FabG) % changed
d_C4_BKeACP    = P.kcat8_H.*c_C5_FabF_Act_MalACP;
d_C6_BKeACP    = P.kcat8(1).*c_C7_FabF_Act_MalACP;
d_C8_BKeACP    = P.kcat8(2).*c_C9_FabF_Act_MalACP;
d_C10_BKeACP  = P.kcat8(3).*c_C11_FabF_Act_MalACP;
d_C12_BKeACP  = P.kcat8(4).*c_C13_FabF_Act_MalACP;
d_C14_BKeACP  = P.kcat8(5).*c_C15_FabF_Act_MalACP;
d_C16_BKeACP  = P.kcat8(6).*c_C17_FabF_Act_MalACP;
d_C18_BKeACP  = P.kcat8(7).*c_C19_FabF_Act_MalACP;
d_C20_BKeACP  = P.kcat8(8).*c_C21_FabF_Act_MalACP;

% C2n:1 (n=6:10) B-ketoacyl-ACPs (FabF + FabB - FabG)
d_C14_BKeACP_un = P.kcat8_un(5).*c_C15_FabF_Act_MalACP_un;
d_C16_BKeACP_un = P.kcat8_un(6).*c_C17_FabF_Act_MalACP_un;
d_C18_BKeACP_un = P.kcat8_un(7).*c_C19_FabF_Act_MalACP_un;
d_C20_BKeACP_un = P.kcat8_un(8).*c_C21_FabF_Act_MalACP_un;

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C4_AcACP   = P.k8_1r(1).*c_C4_FabF_AcACP   - P.k8_1f(1).*c_FabF.*c_C4_AcACP;
d_C6_AcACP   = P.k8_1r(2).*c_C6_FabF_AcACP   - P.k8_1f(2).*c_FabF.*c_C6_AcACP;
d_C8_AcACP   = P.k8_1r(3).*c_C8_FabF_AcACP   - P.k8_1f(3).*c_FabF.*c_C8_AcACP;
d_C10_AcACP = P.k8_1r(4).*c_C10_FabF_AcACP - P.k8_1f(4).*c_FabF.*c_C10_AcACP;
d_C12_AcACP = P.k8_1r(5).*c_C12_FabF_AcACP - P.k8_1f(5).*c_FabF.*c_C12_AcACP;
d_C14_AcACP = P.k8_1r(6).*c_C14_FabF_AcACP - P.k8_1f(6).*c_FabF.*c_C14_AcACP;
d_C16_AcACP = P.k8_1r(7).*c_C16_FabF_AcACP - P.k8_1f(7).*c_FabF.*c_C16_AcACP;
d_C18_AcACP = P.k8_1r(8).*c_C18_FabF_AcACP - P.k8_1f(8).*c_FabF.*c_C18_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C12_AcACP_un = P.k8_1r(5).*c_C12_FabF_AcACP_un - P.k8_1f(5).*c_FabF.*c_C12_AcACP_un;
d_C14_AcACP_un = P.k8_1r(6).*c_C14_FabF_AcACP_un - P.k8_1f(6).*c_FabF.*c_C14_AcACP_un;
d_C16_AcACP_un = P.k8_1r(7).*c_C16_FabF_AcACP_un - P.k8_1f(7).*c_FabF.*c_C16_AcACP_un;
d_C18_AcACP_un = P.k8_1r(8).*c_C18_FabF_AcACP_un - P.k8_1f(8).*c_FabF.*c_C18_AcACP_un;

% C2n (n=2:9) FabF-Acyl-ACPs
d_C4_FabF_AcACP   = P.k8_1f(1).*c_FabF.*c_C4_AcACP   - P.k8_1r(1).*c_C4_FabF_AcACP   + P.k8_2r(1).*c_C4_FabF_Act.*c_ACP  - P.k8_2f(1).*c_C4_FabF_AcACP;
d_C6_FabF_AcACP   = P.k8_1f(2).*c_FabF.*c_C6_AcACP   - P.k8_1r(2).*c_C6_FabF_AcACP   + P.k8_2r(2).*c_C6_FabF_Act.*c_ACP  - P.k8_2f(2).*c_C6_FabF_AcACP;
d_C8_FabF_AcACP   = P.k8_1f(3).*c_FabF.*c_C8_AcACP   - P.k8_1r(3).*c_C8_FabF_AcACP   + P.k8_2r(3).*c_C8_FabF_Act.*c_ACP  - P.k8_2f(3).*c_C8_FabF_AcACP;
d_C10_FabF_AcACP = P.k8_1f(4).*c_FabF.*c_C10_AcACP - P.k8_1r(4).*c_C10_FabF_AcACP + P.k8_2r(4).*c_C10_FabF_Act.*c_ACP - P.k8_2f(4).*c_C10_FabF_AcACP;
d_C12_FabF_AcACP = P.k8_1f(5).*c_FabF.*c_C12_AcACP - P.k8_1r(5).*c_C12_FabF_AcACP + P.k8_2r(5).*c_C12_FabF_Act.*c_ACP - P.k8_2f(5).*c_C12_FabF_AcACP;
d_C14_FabF_AcACP = P.k8_1f(6).*c_FabF.*c_C14_AcACP - P.k8_1r(6).*c_C14_FabF_AcACP + P.k8_2r(6).*c_C14_FabF_Act.*c_ACP - P.k8_2f(6).*c_C14_FabF_AcACP;
d_C16_FabF_AcACP = P.k8_1f(7).*c_FabF.*c_C16_AcACP - P.k8_1r(7).*c_C16_FabF_AcACP + P.k8_2r(7).*c_C16_FabF_Act.*c_ACP - P.k8_2f(7).*c_C16_FabF_AcACP;
d_C18_FabF_AcACP = P.k8_1f(8).*c_FabF.*c_C18_AcACP - P.k8_1r(8).*c_C18_FabF_AcACP + P.k8_2r(8).*c_C18_FabF_Act.*c_ACP - P.k8_2f(8).*c_C18_FabF_AcACP;

% C2n:1 (n=6:9) FabF-Acyl-ACPs
d_C12_FabF_AcACP_un = P.k8_1f(5).*c_FabF.*c_C12_AcACP_un - P.k8_1r(5).*c_C12_FabF_AcACP_un + P.k8_2r(5).*c_C12_FabF_Act_un.*c_ACP - P.k8_2f(5).*c_C12_FabF_AcACP_un;
d_C14_FabF_AcACP_un = P.k8_1f(6).*c_FabF.*c_C14_AcACP_un - P.k8_1r(6).*c_C14_FabF_AcACP_un + P.k8_2r(6).*c_C14_FabF_Act_un.*c_ACP - P.k8_2f(6).*c_C14_FabF_AcACP_un;
d_C16_FabF_AcACP_un = P.k8_1f(7).*c_FabF.*c_C16_AcACP_un - P.k8_1r(7).*c_C16_FabF_AcACP_un + P.k8_2r(7).*c_C16_FabF_Act_un.*c_ACP - P.k8_2f(7).*c_C16_FabF_AcACP_un;
d_C18_FabF_AcACP_un = P.k8_1f(8).*c_FabF.*c_C18_AcACP_un - P.k8_1r(8).*c_C18_FabF_AcACP_un + P.k8_2r(8).*c_C18_FabF_Act_un.*c_ACP - P.k8_2f(8).*c_C18_FabF_AcACP_un;

% C2n (n=2:9) FabF*
d_C4_FabF_Act   = P.k8_2f(1).*c_C4_FabF_AcACP   - P.k8_2r(1).*c_C4_FabF_Act.*c_ACP   + P.k8_3r(1).*c_C7_FabF_Act_MalACP  - P.k8_3f(1).*c_C4_FabF_Act.*c_C3_MalACP;
d_C6_FabF_Act   = P.k8_2f(2).*c_C6_FabF_AcACP   - P.k8_2r(2).*c_C6_FabF_Act.*c_ACP   + P.k8_3r(2).*c_C9_FabF_Act_MalACP  - P.k8_3f(2).*c_C6_FabF_Act.*c_C3_MalACP;
d_C8_FabF_Act   = P.k8_2f(3).*c_C8_FabF_AcACP   - P.k8_2r(3).*c_C8_FabF_Act.*c_ACP   + P.k8_3r(3).*c_C11_FabF_Act_MalACP - P.k8_3f(3).*c_C8_FabF_Act.*c_C3_MalACP;
d_C10_FabF_Act = P.k8_2f(4).*c_C10_FabF_AcACP - P.k8_2r(4).*c_C10_FabF_Act.*c_ACP + P.k8_3r(4).*c_C13_FabF_Act_MalACP - P.k8_3f(4).*c_C10_FabF_Act.*c_C3_MalACP;
d_C12_FabF_Act = P.k8_2f(5).*c_C12_FabF_AcACP - P.k8_2r(5).*c_C12_FabF_Act.*c_ACP + P.k8_3r(5).*c_C15_FabF_Act_MalACP - P.k8_3f(5).*c_C12_FabF_Act.*c_C3_MalACP;
d_C14_FabF_Act = P.k8_2f(6).*c_C14_FabF_AcACP - P.k8_2r(6).*c_C14_FabF_Act.*c_ACP + P.k8_3r(6).*c_C17_FabF_Act_MalACP - P.k8_3f(6).*c_C14_FabF_Act.*c_C3_MalACP;
d_C16_FabF_Act = P.k8_2f(7).*c_C16_FabF_AcACP - P.k8_2r(7).*c_C16_FabF_Act.*c_ACP + P.k8_3r(7).*c_C19_FabF_Act_MalACP - P.k8_3f(7).*c_C16_FabF_Act.*c_C3_MalACP;
d_C18_FabF_Act = P.k8_2f(8).*c_C18_FabF_AcACP - P.k8_2r(8).*c_C18_FabF_Act.*c_ACP + P.k8_3r(8).*c_C21_FabF_Act_MalACP - P.k8_3f(8).*c_C18_FabF_Act.*c_C3_MalACP;

% C2n:1 (n=6:9) FabF*
d_C12_FabF_Act_un = P.k8_2f(5).*c_C12_FabF_AcACP_un - P.k8_2r(5).*c_C12_FabF_Act_un.*c_ACP + P.k8_3r(5).*c_C15_FabF_Act_MalACP_un - P.k8_3f(5).*c_C12_FabF_Act_un.*c_C3_MalACP;
d_C14_FabF_Act_un = P.k8_2f(6).*c_C14_FabF_AcACP_un - P.k8_2r(6).*c_C14_FabF_Act_un.*c_ACP + P.k8_3r(6).*c_C17_FabF_Act_MalACP_un - P.k8_3f(6).*c_C14_FabF_Act_un.*c_C3_MalACP;
d_C16_FabF_Act_un = P.k8_2f(7).*c_C16_FabF_AcACP_un - P.k8_2r(7).*c_C16_FabF_Act_un.*c_ACP + P.k8_3r(7).*c_C19_FabF_Act_MalACP_un - P.k8_3f(7).*c_C16_FabF_Act_un.*c_C3_MalACP;
d_C18_FabF_Act_un = P.k8_2f(8).*c_C18_FabF_AcACP_un - P.k8_2r(8).*c_C18_FabF_Act_un.*c_ACP + P.k8_3r(8).*c_C21_FabF_Act_MalACP_un - P.k8_3f(8).*c_C18_FabF_Act_un.*c_C3_MalACP;

% C2n (n=2:9) FabF*-Malonyl-ACPs
d_C7_FabF_Act_MalACP   = P.k8_3f(1).*c_C4_FabF_Act.*c_C3_MalACP   - P.k8_3r(1).*c_C7_FabF_Act_MalACP  - P.kcat8(1).*c_C7_FabF_Act_MalACP;
d_C9_FabF_Act_MalACP   = P.k8_3f(2).*c_C6_FabF_Act.*c_C3_MalACP   - P.k8_3r(2).*c_C9_FabF_Act_MalACP  - P.kcat8(2).*c_C9_FabF_Act_MalACP;
d_C11_FabF_Act_MalACP = P.k8_3f(3).*c_C8_FabF_Act.*c_C3_MalACP   - P.k8_3r(3).*c_C11_FabF_Act_MalACP - P.kcat8(3).*c_C11_FabF_Act_MalACP;
d_C13_FabF_Act_MalACP = P.k8_3f(4).*c_C10_FabF_Act.*c_C3_MalACP - P.k8_3r(4).*c_C13_FabF_Act_MalACP - P.kcat8(4).*c_C13_FabF_Act_MalACP;
d_C15_FabF_Act_MalACP = P.k8_3f(5).*c_C12_FabF_Act.*c_C3_MalACP - P.k8_3r(5).*c_C15_FabF_Act_MalACP - P.kcat8(5).*c_C15_FabF_Act_MalACP;
d_C17_FabF_Act_MalACP = P.k8_3f(6).*c_C14_FabF_Act.*c_C3_MalACP - P.k8_3r(6).*c_C17_FabF_Act_MalACP - P.kcat8(6).*c_C17_FabF_Act_MalACP;
d_C19_FabF_Act_MalACP = P.k8_3f(7).*c_C16_FabF_Act.*c_C3_MalACP - P.k8_3r(7).*c_C19_FabF_Act_MalACP - P.kcat8(7).*c_C19_FabF_Act_MalACP;
d_C21_FabF_Act_MalACP = P.k8_3f(8).*c_C18_FabF_Act.*c_C3_MalACP - P.k8_3r(8).*c_C21_FabF_Act_MalACP - P.kcat8(8).*c_C21_FabF_Act_MalACP;

% C2n:1 (n=6:9) FabF*-Malonyl-ACPs
d_C15_FabF_Act_MalACP_un = P.k8_3f(5).*c_C12_FabF_Act_un.*c_C3_MalACP - P.k8_3r(5).*c_C15_FabF_Act_MalACP_un - P.kcat8_un(5).*c_C15_FabF_Act_MalACP_un;
d_C17_FabF_Act_MalACP_un = P.k8_3f(6).*c_C14_FabF_Act_un.*c_C3_MalACP - P.k8_3r(6).*c_C17_FabF_Act_MalACP_un - P.kcat8_un(6).*c_C17_FabF_Act_MalACP_un;
d_C19_FabF_Act_MalACP_un = P.k8_3f(7).*c_C16_FabF_Act_un.*c_C3_MalACP - P.k8_3r(7).*c_C19_FabF_Act_MalACP_un - P.kcat8_un(7).*c_C19_FabF_Act_MalACP_un;
d_C21_FabF_Act_MalACP_un = P.k8_3f(8).*c_C18_FabF_Act_un.*c_C3_MalACP - P.k8_3r(8).*c_C21_FabF_Act_MalACP_un - P.kcat8_un(8).*c_C21_FabF_Act_MalACP_un;

% FabF-ACP
d_FabF_ACP = P.k8_inh_f.*c_FabF.*c_ACP - P.k8_inh_r.*c_FabF_ACP;

% Giving FabF FabH-like activity
% FabF-Acetyl-CoA
d_C2_FabF_AcCoA = P.k8_4f.*c_FabF.*c_C2_AcCoA - P.k8_4r.*c_C2_FabF_AcCoA + P.k8_5r.*c_C2_FabF_Act.*c_CoA - P.k8_5f.*c_C2_FabF_AcCoA; 

% FabF*
d_C2_FabF_Act = P.k8_5f.*c_C2_FabF_AcCoA - P.k8_5r.*c_C2_FabF_Act.*c_CoA + P.k8_6r.*c_C5_FabF_Act_MalACP - P.k8_6f.*c_C2_FabF_Act.*c_C3_MalACP + P.k8_9f.*c_C2_FabF_AcACP - P.k8_9r.*c_C2_FabF_Act.*c_ACP; 

% FabF*-Malonyl-ACP
d_C5_FabF_Act_MalACP = P.k8_6f.*c_C2_FabF_Act.*c_C3_MalACP - P.k8_6r.*c_C5_FabF_Act_MalACP - P.kcat8_H.*c_C5_FabF_Act_MalACP; 

% Acetyl-ACP
d_C2_AcACP = P.k8_8r.*c_C2_FabF_AcACP - P.k8_8f.*c_FabF.*c_C2_AcACP; 

% FabF-Malonyl-ACP
d_C3_FabF_MalACP = P.k8_7f.*c_FabF.*c_C3_MalACP - P.k8_7r.*c_C3_FabF_MalACP - P.kcat8_CO2.*c_C3_FabF_MalACP; 

% FabF-Acetyl-ACP
d_C2_FabF_AcACP = P.kcat8_CO2.*c_C3_FabF_MalACP + P.k8_8f.*c_FabF.*c_C2_AcACP - P.k8_8r.*c_C2_FabF_AcACP + P.k8_9r.*c_C2_FabF_Act.*c_ACP - P.k8_9f.*c_C2_FabF_AcACP; 


dcdt = [d_C2_AcCoA;d_ACP;d_CoA;d_C3_MalACP;d_C1_CO2;d_C4_BKeACP;d_C6_BKeACP;...
    d_C8_BKeACP;d_C10_BKeACP;d_C12_BKeACP;d_C14_BKeACP;d_C16_BKeACP;...
    d_C18_BKeACP;d_C20_BKeACP;d_C14_BKeACP_un;d_C16_BKeACP_un;d_C18_BKeACP_un;...
    d_C20_BKeACP_un;d_C4_AcACP;d_C6_AcACP;d_C8_AcACP;d_C10_AcACP;d_C12_AcACP;...
    d_C14_AcACP;d_C16_AcACP;d_C18_AcACP;d_C12_AcACP_un;d_C14_AcACP_un;...
    d_C16_AcACP_un;d_C18_AcACP_un;d_C4_FabF_AcACP;d_C6_FabF_AcACP;...
    d_C8_FabF_AcACP;d_C10_FabF_AcACP;d_C12_FabF_AcACP;d_C14_FabF_AcACP;...
    d_C16_FabF_AcACP;d_C18_FabF_AcACP;d_C12_FabF_AcACP_un;d_C14_FabF_AcACP_un;...
    d_C16_FabF_AcACP_un;d_C18_FabF_AcACP_un;d_C4_FabF_Act;d_C6_FabF_Act;...
    d_C8_FabF_Act;d_C10_FabF_Act;d_C12_FabF_Act;d_C14_FabF_Act;d_C16_FabF_Act;...
    d_C18_FabF_Act;d_C12_FabF_Act_un;d_C14_FabF_Act_un;d_C16_FabF_Act_un;...
    d_C18_FabF_Act_un;d_C7_FabF_Act_MalACP;d_C9_FabF_Act_MalACP;...
    d_C11_FabF_Act_MalACP;d_C13_FabF_Act_MalACP;d_C15_FabF_Act_MalACP;...
    d_C17_FabF_Act_MalACP;d_C19_FabF_Act_MalACP;d_C21_FabF_Act_MalACP;...
    d_C15_FabF_Act_MalACP_un;d_C17_FabF_Act_MalACP_un;d_C19_FabF_Act_MalACP_un;...
    d_C21_FabF_Act_MalACP_un;d_FabF_ACP;d_C2_FabF_AcCoA;d_C2_FabF_Act;...
    d_C5_FabF_Act_MalACP;d_C2_AcACP;d_C3_FabF_MalACP;d_C2_FabF_AcACP];

end
