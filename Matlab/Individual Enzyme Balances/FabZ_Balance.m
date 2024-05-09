%% FabZ Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ACP', 'c_C4_BHyAcACP', 'c_C6_BHyAcACP', 'c_C8_BHyAcACP', 'c_C10_BHyAcACP',...
    'c_C12_BHyAcACP', 'c_C14_BHyAcACP', 'c_C16_BHyAcACP', 'c_C18_BHyAcACP',...
    'c_C20_BHyAcACP', 'c_C12_BHyAcACP_un', 'c_C14_BHyAcACP_un', 'c_C16_BHyAcACP_un',...
    'c_C18_BHyAcACP_un', 'c_C20_BHyAcACP_un', 'c_C4_EnAcACP', 'c_C6_EnAcACP', 'c_C8_EnAcACP',...
    'c_C10_EnAcACP', 'c_C12_EnAcACP', 'c_C14_EnAcACP', 'c_C16_EnAcACP', 'c_C18_EnAcACP',...
    'c_C20_EnAcACP', 'c_C12_EnAcACP_un', 'c_C14_EnAcACP_un', 'c_C16_EnAcACP_un',...
    'c_C18_EnAcACP_un', 'c_C20_EnAcACP_un', 'c_C4_FabZ_BHyAcACP', 'c_C6_FabZ_BHyAcACP',...
    'c_C8_FabZ_BHyAcACP', 'c_C10_FabZ_BHyAcACP', 'c_C12_FabZ_BHyAcACP',...
    'c_C14_FabZ_BHyAcACP', 'c_C16_FabZ_BHyAcACP', 'c_C18_FabZ_BHyAcACP',...
    'c_C20_FabZ_BHyAcACP', 'c_C12_FabZ_BHyAcACP_un', 'c_C14_FabZ_BHyAcACP_un',...
    'c_C16_FabZ_BHyAcACP_un', 'c_C18_FabZ_BHyAcACP_un', 'c_C20_FabZ_BHyAcACP_un',...
    'c_C4_FabZ_EnAcACP', 'c_C6_FabZ_EnAcACP', 'c_C8_FabZ_EnAcACP', 'c_C10_FabZ_EnAcACP',...
    'c_C12_FabZ_EnAcACP', 'c_C14_FabZ_EnAcACP', 'c_C16_FabZ_EnAcACP', 'c_C18_FabZ_EnAcACP',...
    'c_C20_FabZ_EnAcACP', 'c_C12_FabZ_EnAcACP_un', 'c_C14_FabZ_EnAcACP_un',...
    'c_C16_FabZ_EnAcACP_un', 'c_C18_FabZ_EnAcACP_un', 'c_C20_FabZ_EnAcACP_un', 'c_FabZ_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(2:15) = 0.7143; % Beta-hydroxy-acyl-ACPs - evenly distribute 10 ACP to all of them

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 0 1 0 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabZ(t,c,P);
tic
[T_FabZ,C_FabZ] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabZ] = Calc_Function(T_FabZ,C_FabZ,S);

[balance_conc_FabZ, balances_FabZ, total_conc_FabZ, carbon_FabZ] = mass_balance(C_FabZ,P);


function dcdt = ODE_Function_FabZ(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabZ
c_FabZ = P.FabZtot - c_FabZ_ACP...
    - c_C4_FabZ_BHyAcACP - c_C6_FabZ_BHyAcACP - c_C8_FabZ_BHyAcACP - c_C10_FabZ_BHyAcACP - c_C12_FabZ_BHyAcACP - c_C14_FabZ_BHyAcACP - c_C16_FabZ_BHyAcACP - c_C18_FabZ_BHyAcACP - c_C20_FabZ_BHyAcACP...
    - c_C12_FabZ_BHyAcACP_un - c_C14_FabZ_BHyAcACP_un - c_C16_FabZ_BHyAcACP_un - c_C18_FabZ_BHyAcACP_un - c_C20_FabZ_BHyAcACP_un...
    - c_C4_FabZ_EnAcACP - c_C6_FabZ_EnAcACP - c_C8_FabZ_EnAcACP - c_C10_FabZ_EnAcACP - c_C12_FabZ_EnAcACP - c_C14_FabZ_EnAcACP - c_C16_FabZ_EnAcACP - c_C18_FabZ_EnAcACP - c_C20_FabZ_EnAcACP...
    - c_C12_FabZ_EnAcACP_un - c_C14_FabZ_EnAcACP_un - c_C16_FabZ_EnAcACP_un - c_C18_FabZ_EnAcACP_un - c_C20_FabZ_EnAcACP_un;

% Set of differential equations
% ACP
d_ACP = P.k5_inh_r.*c_FabZ_ACP   - P.k5_inh_f.*c_FabZ.*c_ACP; 

% C2n (n=2:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C4_BHyAcACP   = P.k5_1r(1).*c_C4_FabZ_BHyAcACP  - P.k5_1f(1).*c_FabZ.*c_C4_BHyAcACP;
d_C6_BHyAcACP   = P.k5_1r(2).*c_C6_FabZ_BHyAcACP  - P.k5_1f(2).*c_FabZ.*c_C6_BHyAcACP;
d_C8_BHyAcACP   = P.k5_1r(3).*c_C8_FabZ_BHyAcACP  - P.k5_1f(3).*c_FabZ.*c_C8_BHyAcACP;
d_C10_BHyAcACP = P.k5_1r(4).*c_C10_FabZ_BHyAcACP - P.k5_1f(4).*c_FabZ.*c_C10_BHyAcACP;
d_C12_BHyAcACP = P.k5_1r(5).*c_C12_FabZ_BHyAcACP - P.k5_1f(5).*c_FabZ.*c_C12_BHyAcACP;
d_C14_BHyAcACP = P.k5_1r(6).*c_C14_FabZ_BHyAcACP - P.k5_1f(6).*c_FabZ.*c_C14_BHyAcACP;
d_C16_BHyAcACP = P.k5_1r(7).*c_C16_FabZ_BHyAcACP - P.k5_1f(7).*c_FabZ.*c_C16_BHyAcACP;
d_C18_BHyAcACP = P.k5_1r(8).*c_C18_FabZ_BHyAcACP - P.k5_1f(8).*c_FabZ.*c_C18_BHyAcACP;
d_C20_BHyAcACP = P.k5_1r(9).*c_C20_FabZ_BHyAcACP - P.k5_1f(9).*c_FabZ.*c_C20_BHyAcACP;

% C2n:1 (n=6:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C12_BHyAcACP_un = P.k5_1r(5).*c_C12_FabZ_BHyAcACP_un - P.k5_1f(5).*c_FabZ.*c_C12_BHyAcACP_un;
d_C14_BHyAcACP_un = P.k5_1r(6).*c_C14_FabZ_BHyAcACP_un - P.k5_1f(6).*c_FabZ.*c_C14_BHyAcACP_un;
d_C16_BHyAcACP_un = P.k5_1r(7).*c_C16_FabZ_BHyAcACP_un - P.k5_1f(7).*c_FabZ.*c_C16_BHyAcACP_un;
d_C18_BHyAcACP_un = P.k5_1r(8).*c_C18_FabZ_BHyAcACP_un - P.k5_1f(8).*c_FabZ.*c_C18_BHyAcACP_un;
d_C20_BHyAcACP_un = P.k5_1r(9).*c_C20_FabZ_BHyAcACP_un - P.k5_1f(9).*c_FabZ.*c_C20_BHyAcACP_un;

% C2n (n=2:10) Enoyl-Acyl-ACPs (FabZ + FabA - FabI) 
d_C4_EnAcACP   = P.k5_3f(1).*c_C4_FabZ_EnAcACP   - P.k5_3r(1).*c_FabZ.*c_C4_EnAcACP;
d_C6_EnAcACP   = P.k5_3f(2).*c_C6_FabZ_EnAcACP   - P.k5_3r(2).*c_FabZ.*c_C6_EnAcACP;
d_C8_EnAcACP   = P.k5_3f(3).*c_C8_FabZ_EnAcACP   - P.k5_3r(3).*c_FabZ.*c_C8_EnAcACP;
d_C10_EnAcACP = P.k5_3f(4).*c_C10_FabZ_EnAcACP - P.k5_3r(4).*c_FabZ.*c_C10_EnAcACP;
d_C12_EnAcACP = P.k5_3f(5).*c_C12_FabZ_EnAcACP - P.k5_3r(5).*c_FabZ.*c_C12_EnAcACP;
d_C14_EnAcACP = P.k5_3f(6).*c_C14_FabZ_EnAcACP - P.k5_3r(6).*c_FabZ.*c_C14_EnAcACP;
d_C16_EnAcACP = P.k5_3f(7).*c_C16_FabZ_EnAcACP - P.k5_3r(7).*c_FabZ.*c_C16_EnAcACP;
d_C18_EnAcACP = P.k5_3f(8).*c_C18_FabZ_EnAcACP - P.k5_3r(8).*c_FabZ.*c_C18_EnAcACP;
d_C20_EnAcACP = P.k5_3f(9).*c_C20_FabZ_EnAcACP - P.k5_3r(9).*c_FabZ.*c_C20_EnAcACP;

% C2n:1 (n=6:10) Enoyl-Acyl-ACPs  (FabZ + FabA - FabI)
d_C12_EnAcACP_un = P.k5_3f(5).*c_C12_FabZ_EnAcACP_un - P.k5_3r(5).*c_FabZ.*c_C12_EnAcACP_un;
d_C14_EnAcACP_un = P.k5_3f(6).*c_C14_FabZ_EnAcACP_un - P.k5_3r(6).*c_FabZ.*c_C14_EnAcACP_un;
d_C16_EnAcACP_un = P.k5_3f(7).*c_C16_FabZ_EnAcACP_un - P.k5_3r(7).*c_FabZ.*c_C16_EnAcACP_un;
d_C18_EnAcACP_un = P.k5_3f(8).*c_C18_FabZ_EnAcACP_un - P.k5_3r(8).*c_FabZ.*c_C18_EnAcACP_un;
d_C20_EnAcACP_un = P.k5_3f(9).*c_C20_FabZ_EnAcACP_un - P.k5_3r(9).*c_FabZ.*c_C20_EnAcACP_un;

% C2n (n=2:10) FabZ-B-hydroxy-acyl-ACPs
d_C4_FabZ_BHyAcACP   = P.k5_1f(1).*c_FabZ.*c_C4_BHyAcACP   - P.k5_1r(1).*c_C4_FabZ_BHyAcACP   + P.k5_2r(1).*c_C4_FabZ_EnAcACP  - P.kcat5(1).*c_C4_FabZ_BHyAcACP;
d_C6_FabZ_BHyAcACP   = P.k5_1f(2).*c_FabZ.*c_C6_BHyAcACP   - P.k5_1r(2).*c_C6_FabZ_BHyAcACP   + P.k5_2r(2).*c_C6_FabZ_EnAcACP  - P.kcat5(2).*c_C6_FabZ_BHyAcACP;
d_C8_FabZ_BHyAcACP   = P.k5_1f(3).*c_FabZ.*c_C8_BHyAcACP   - P.k5_1r(3).*c_C8_FabZ_BHyAcACP   + P.k5_2r(3).*c_C8_FabZ_EnAcACP  - P.kcat5(3).*c_C8_FabZ_BHyAcACP;
d_C10_FabZ_BHyAcACP = P.k5_1f(4).*c_FabZ.*c_C10_BHyAcACP - P.k5_1r(4).*c_C10_FabZ_BHyAcACP + P.k5_2r(4).*c_C10_FabZ_EnAcACP - P.kcat5(4).*c_C10_FabZ_BHyAcACP;
d_C12_FabZ_BHyAcACP = P.k5_1f(5).*c_FabZ.*c_C12_BHyAcACP - P.k5_1r(5).*c_C12_FabZ_BHyAcACP + P.k5_2r(5).*c_C12_FabZ_EnAcACP - P.kcat5(5).*c_C12_FabZ_BHyAcACP;
d_C14_FabZ_BHyAcACP = P.k5_1f(6).*c_FabZ.*c_C14_BHyAcACP - P.k5_1r(6).*c_C14_FabZ_BHyAcACP + P.k5_2r(6).*c_C14_FabZ_EnAcACP - P.kcat5(6).*c_C14_FabZ_BHyAcACP;
d_C16_FabZ_BHyAcACP = P.k5_1f(7).*c_FabZ.*c_C16_BHyAcACP - P.k5_1r(7).*c_C16_FabZ_BHyAcACP + P.k5_2r(7).*c_C16_FabZ_EnAcACP - P.kcat5(7).*c_C16_FabZ_BHyAcACP;
d_C18_FabZ_BHyAcACP = P.k5_1f(8).*c_FabZ.*c_C18_BHyAcACP - P.k5_1r(8).*c_C18_FabZ_BHyAcACP + P.k5_2r(8).*c_C18_FabZ_EnAcACP - P.kcat5(8).*c_C18_FabZ_BHyAcACP;
d_C20_FabZ_BHyAcACP = P.k5_1f(9).*c_FabZ.*c_C20_BHyAcACP - P.k5_1r(9).*c_C20_FabZ_BHyAcACP + P.k5_2r(9).*c_C20_FabZ_EnAcACP - P.kcat5(9).*c_C20_FabZ_BHyAcACP;

% C2n:1 (n=6:10) FabZ-B-hydroxy-acyl-ACPs
d_C12_FabZ_BHyAcACP_un = P.k5_1f(5).*c_FabZ.*c_C12_BHyAcACP_un - P.k5_1r(5).*c_C12_FabZ_BHyAcACP_un + P.k5_2r(5).*c_C12_FabZ_EnAcACP_un - P.kcat5(5).*c_C12_FabZ_BHyAcACP_un;
d_C14_FabZ_BHyAcACP_un = P.k5_1f(6).*c_FabZ.*c_C14_BHyAcACP_un - P.k5_1r(6).*c_C14_FabZ_BHyAcACP_un + P.k5_2r(6).*c_C14_FabZ_EnAcACP_un - P.kcat5(6).*c_C14_FabZ_BHyAcACP_un;
d_C16_FabZ_BHyAcACP_un = P.k5_1f(7).*c_FabZ.*c_C16_BHyAcACP_un - P.k5_1r(7).*c_C16_FabZ_BHyAcACP_un + P.k5_2r(7).*c_C16_FabZ_EnAcACP_un - P.kcat5(7).*c_C16_FabZ_BHyAcACP_un;
d_C18_FabZ_BHyAcACP_un = P.k5_1f(8).*c_FabZ.*c_C18_BHyAcACP_un - P.k5_1r(8).*c_C18_FabZ_BHyAcACP_un + P.k5_2r(8).*c_C18_FabZ_EnAcACP_un - P.kcat5(8).*c_C18_FabZ_BHyAcACP_un;
d_C20_FabZ_BHyAcACP_un = P.k5_1f(9).*c_FabZ.*c_C20_BHyAcACP_un - P.k5_1r(9).*c_C20_FabZ_BHyAcACP_un + P.k5_2r(9).*c_C20_FabZ_EnAcACP_un - P.kcat5(9).*c_C20_FabZ_BHyAcACP_un;

% C2n (n=2:10) FabZ-Enoyl-Acyl-ACPs
d_C4_FabZ_EnAcACP   = P.kcat5(1).*c_C4_FabZ_BHyAcACP   - P.k5_2r(1).*c_C4_FabZ_EnAcACP   + P.k5_3r(1).*c_FabZ.*c_C4_EnAcACP  - P.k5_3f(1).*c_C4_FabZ_EnAcACP;
d_C6_FabZ_EnAcACP   = P.kcat5(2).*c_C6_FabZ_BHyAcACP   - P.k5_2r(2).*c_C6_FabZ_EnAcACP   + P.k5_3r(2).*c_FabZ.*c_C6_EnAcACP  - P.k5_3f(2).*c_C6_FabZ_EnAcACP;
d_C8_FabZ_EnAcACP   = P.kcat5(3).*c_C8_FabZ_BHyAcACP   - P.k5_2r(3).*c_C8_FabZ_EnAcACP   + P.k5_3r(3).*c_FabZ.*c_C8_EnAcACP  - P.k5_3f(3).*c_C8_FabZ_EnAcACP;
d_C10_FabZ_EnAcACP = P.kcat5(4).*c_C10_FabZ_BHyAcACP - P.k5_2r(4).*c_C10_FabZ_EnAcACP + P.k5_3r(4).*c_FabZ.*c_C10_EnAcACP - P.k5_3f(4).*c_C10_FabZ_EnAcACP;
d_C12_FabZ_EnAcACP = P.kcat5(5).*c_C12_FabZ_BHyAcACP - P.k5_2r(5).*c_C12_FabZ_EnAcACP + P.k5_3r(5).*c_FabZ.*c_C12_EnAcACP - P.k5_3f(5).*c_C12_FabZ_EnAcACP;
d_C14_FabZ_EnAcACP = P.kcat5(6).*c_C14_FabZ_BHyAcACP - P.k5_2r(6).*c_C14_FabZ_EnAcACP + P.k5_3r(6).*c_FabZ.*c_C14_EnAcACP - P.k5_3f(6).*c_C14_FabZ_EnAcACP;
d_C16_FabZ_EnAcACP = P.kcat5(7).*c_C16_FabZ_BHyAcACP - P.k5_2r(7).*c_C16_FabZ_EnAcACP + P.k5_3r(7).*c_FabZ.*c_C16_EnAcACP - P.k5_3f(7).*c_C16_FabZ_EnAcACP;
d_C18_FabZ_EnAcACP = P.kcat5(8).*c_C18_FabZ_BHyAcACP - P.k5_2r(8).*c_C18_FabZ_EnAcACP + P.k5_3r(8).*c_FabZ.*c_C18_EnAcACP - P.k5_3f(8).*c_C18_FabZ_EnAcACP;
d_C20_FabZ_EnAcACP = P.kcat5(9).*c_C20_FabZ_BHyAcACP - P.k5_2r(9).*c_C20_FabZ_EnAcACP + P.k5_3r(9).*c_FabZ.*c_C20_EnAcACP - P.k5_3f(9).*c_C20_FabZ_EnAcACP;

% C2n:1 (n=6:10) FabZ-Enoyl-Acyl-ACPs
d_C12_FabZ_EnAcACP_un = P.kcat5(5).*c_C12_FabZ_BHyAcACP_un - P.k5_2r(5).*c_C12_FabZ_EnAcACP_un + P.k5_3r(5).*c_FabZ.*c_C12_EnAcACP_un - P.k5_3f(5).*c_C12_FabZ_EnAcACP_un;
d_C14_FabZ_EnAcACP_un = P.kcat5(6).*c_C14_FabZ_BHyAcACP_un - P.k5_2r(6).*c_C14_FabZ_EnAcACP_un + P.k5_3r(6).*c_FabZ.*c_C14_EnAcACP_un - P.k5_3f(6).*c_C14_FabZ_EnAcACP_un;
d_C16_FabZ_EnAcACP_un = P.kcat5(7).*c_C16_FabZ_BHyAcACP_un - P.k5_2r(7).*c_C16_FabZ_EnAcACP_un + P.k5_3r(7).*c_FabZ.*c_C16_EnAcACP_un - P.k5_3f(7).*c_C16_FabZ_EnAcACP_un;
d_C18_FabZ_EnAcACP_un = P.kcat5(8).*c_C18_FabZ_BHyAcACP_un - P.k5_2r(8).*c_C18_FabZ_EnAcACP_un + P.k5_3r(8).*c_FabZ.*c_C18_EnAcACP_un - P.k5_3f(8).*c_C18_FabZ_EnAcACP_un;
d_C20_FabZ_EnAcACP_un = P.kcat5(9).*c_C20_FabZ_BHyAcACP_un - P.k5_2r(9).*c_C20_FabZ_EnAcACP_un + P.k5_3r(9).*c_FabZ.*c_C20_EnAcACP_un - P.k5_3f(9).*c_C20_FabZ_EnAcACP_un;

% FabZ-ACP
d_FabZ_ACP = P.k5_inh_f.*c_FabZ.*c_ACP - P.k5_inh_r.*c_FabZ_ACP;


dcdt = [d_ACP;d_C4_BHyAcACP;d_C6_BHyAcACP;d_C8_BHyAcACP;d_C10_BHyAcACP;...
    d_C12_BHyAcACP;d_C14_BHyAcACP;d_C16_BHyAcACP;d_C18_BHyAcACP;...
    d_C20_BHyAcACP;d_C12_BHyAcACP_un;d_C14_BHyAcACP_un;d_C16_BHyAcACP_un;...
    d_C18_BHyAcACP_un;d_C20_BHyAcACP_un;d_C4_EnAcACP;d_C6_EnAcACP;d_C8_EnAcACP;...
    d_C10_EnAcACP;d_C12_EnAcACP;d_C14_EnAcACP;d_C16_EnAcACP;d_C18_EnAcACP;...
    d_C20_EnAcACP;d_C12_EnAcACP_un;d_C14_EnAcACP_un;d_C16_EnAcACP_un;...
    d_C18_EnAcACP_un;d_C20_EnAcACP_un;d_C4_FabZ_BHyAcACP;d_C6_FabZ_BHyAcACP;...
    d_C8_FabZ_BHyAcACP;d_C10_FabZ_BHyAcACP;d_C12_FabZ_BHyAcACP;...
    d_C14_FabZ_BHyAcACP;d_C16_FabZ_BHyAcACP;d_C18_FabZ_BHyAcACP;...
    d_C20_FabZ_BHyAcACP;d_C12_FabZ_BHyAcACP_un;d_C14_FabZ_BHyAcACP_un;...
    d_C16_FabZ_BHyAcACP_un;d_C18_FabZ_BHyAcACP_un;d_C20_FabZ_BHyAcACP_un;...
    d_C4_FabZ_EnAcACP;d_C6_FabZ_EnAcACP;d_C8_FabZ_EnAcACP;d_C10_FabZ_EnAcACP;...
    d_C12_FabZ_EnAcACP;d_C14_FabZ_EnAcACP;d_C16_FabZ_EnAcACP;d_C18_FabZ_EnAcACP;...
    d_C20_FabZ_EnAcACP;d_C12_FabZ_EnAcACP_un;d_C14_FabZ_EnAcACP_un;...
    d_C16_FabZ_EnAcACP_un;d_C18_FabZ_EnAcACP_un;d_C20_FabZ_EnAcACP_un;d_FabZ_ACP];

end
