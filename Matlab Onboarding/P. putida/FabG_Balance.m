%% FabG Balance
% Need to have ran the beginning of the real code

S.labels = {'c_ACP', 'c_NADPH', 'c_C4_BKeACP', 'c_C6_BKeACP', 'c_C8_BKeACP', 'c_C10_BKeACP',...
    'c_C12_BKeACP', 'c_C14_BKeACP', 'c_C16_BKeACP', 'c_C18_BKeACP', 'c_C20_BKeACP',...
    'c_C12_BKeACP_un', 'c_C14_BKeACP_un', 'c_C16_BKeACP_un', 'c_C18_BKeACP_un',...
    'c_C20_BKeACP_un', 'c_C4_BHyAcACP', 'c_C6_BHyAcACP', 'c_C8_BHyAcACP', 'c_C10_BHyAcACP',...
    'c_C12_BHyAcACP', 'c_C14_BHyAcACP', 'c_C16_BHyAcACP', 'c_C18_BHyAcACP',...
    'c_C20_BHyAcACP', 'c_C12_BHyAcACP_un', 'c_C14_BHyAcACP_un', 'c_C16_BHyAcACP_un',...
    'c_C18_BHyAcACP_un', 'c_C20_BHyAcACP_un', 'c_FabG_NADPH', 'c_C4_FabG_NADPH_BKeACP',...
    'c_C6_FabG_NADPH_BKeACP', 'c_C8_FabG_NADPH_BKeACP', 'c_C10_FabG_NADPH_BKeACP',...
    'c_C12_FabG_NADPH_BKeACP', 'c_C14_FabG_NADPH_BKeACP', 'c_C16_FabG_NADPH_BKeACP',...
    'c_C18_FabG_NADPH_BKeACP', 'c_C20_FabG_NADPH_BKeACP',...
    'c_C12_FabG_NADPH_BKeACP_un', 'c_C14_FabG_NADPH_BKeACP_un',...
    'c_C16_FabG_NADPH_BKeACP_un', 'c_C18_FabG_NADPH_BKeACP_un',...
    'c_C20_FabG_NADPH_BKeACP_un', 'c_FabG_ACP'};

S.range = [0 150]; %2.5 mins (initial rate)

S.num = length(S.labels); %how many diff eqs there are

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(2) = 1300; % NADPH
S.init_cond(3:16) = 0.7143; % Beta-ketoacyl-ACPs - evenly distribute 10 ACP to all of them

S.kcat_scaling_fabG = EC_kcat4_scaling;

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 0 0 1 0 0 0 0 0 0]; 

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function_FabG(t,c,P);
tic
[T_FabG,C_FabG] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_FabG] = Calc_Function(T_FabG,C_FabG,S);

[balance_conc_FabG, balances_FabG, total_conc_FabG, carbon_FabG] = mass_balance(C_FabG,P);


function dcdt = ODE_Function_FabG(t,c,P)

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% FabG
c_FabG = P.FabGtot - c_FabG_NADPH...
    - c_C4_FabG_NADPH_BKeACP - c_C6_FabG_NADPH_BKeACP - c_C8_FabG_NADPH_BKeACP - c_C10_FabG_NADPH_BKeACP - c_C12_FabG_NADPH_BKeACP - c_C14_FabG_NADPH_BKeACP - c_C16_FabG_NADPH_BKeACP - c_C18_FabG_NADPH_BKeACP - c_C20_FabG_NADPH_BKeACP...
    - c_C12_FabG_NADPH_BKeACP_un - c_C14_FabG_NADPH_BKeACP_un - c_C16_FabG_NADPH_BKeACP_un - c_C18_FabG_NADPH_BKeACP_un - c_C20_FabG_NADPH_BKeACP_un - c_FabG_ACP;

% Set of differential equations
% ACP
d_ACP = P.k4_inh_r.*c_FabG_ACP - P.k4_inh_f.*c_FabG.*c_ACP; 

% NADPH
d_NADPH = P.k4_1r(1).*c_FabG_NADPH - P.k4_1f(1).*c_FabG.*c_NADPH; 

% C2n (n=2:10) B-ketoacyl-ACPs (FabH + FabF + FabB - FabG) % changed
d_C4_BKeACP    = P.k4_2r(1).*c_C4_FabG_NADPH_BKeACP   - P.k4_2f(1).*c_FabG_NADPH.*c_C4_BKeACP;
d_C6_BKeACP    = P.k4_2r(2).*c_C6_FabG_NADPH_BKeACP   - P.k4_2f(2).*c_FabG_NADPH.*c_C6_BKeACP;
d_C8_BKeACP    = P.k4_2r(3).*c_C8_FabG_NADPH_BKeACP   - P.k4_2f(3).*c_FabG_NADPH.*c_C8_BKeACP;
d_C10_BKeACP  = P.k4_2r(4).*c_C10_FabG_NADPH_BKeACP - P.k4_2f(4).*c_FabG_NADPH.*c_C10_BKeACP;
d_C12_BKeACP  = P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP - P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP;
d_C14_BKeACP  = P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP - P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP;
d_C16_BKeACP  = P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP - P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP;
d_C18_BKeACP  = P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP - P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP;
d_C20_BKeACP  = P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP - P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP;

% C2n:1 (n=6:10) B-ketoacyl-ACPs (FabF + FabB - FabG)
d_C12_BKeACP_un = P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP_un - P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP_un;
d_C14_BKeACP_un = P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP_un - P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP_un;
d_C16_BKeACP_un = P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP_un - P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP_un;
d_C18_BKeACP_un = P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP_un - P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP_un;
d_C20_BKeACP_un = P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP_un - P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP_un;

% C2n (n=2:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C4_BHyAcACP   = P.kcat4(1).*c_C4_FabG_NADPH_BKeACP;
d_C6_BHyAcACP   = P.kcat4(2).*c_C6_FabG_NADPH_BKeACP;
d_C8_BHyAcACP   = P.kcat4(3).*c_C8_FabG_NADPH_BKeACP;
d_C10_BHyAcACP = P.kcat4(4).*c_C10_FabG_NADPH_BKeACP;
d_C12_BHyAcACP = P.kcat4(5).*c_C12_FabG_NADPH_BKeACP;
d_C14_BHyAcACP = P.kcat4(6).*c_C14_FabG_NADPH_BKeACP;
d_C16_BHyAcACP = P.kcat4(7).*c_C16_FabG_NADPH_BKeACP;
d_C18_BHyAcACP = P.kcat4(8).*c_C18_FabG_NADPH_BKeACP;
d_C20_BHyAcACP = P.kcat4(9).*c_C20_FabG_NADPH_BKeACP;

% C2n:1 (n=6:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C12_BHyAcACP_un = P.kcat4(5).*c_C12_FabG_NADPH_BKeACP_un;
d_C14_BHyAcACP_un = P.kcat4(6).*c_C14_FabG_NADPH_BKeACP_un;
d_C16_BHyAcACP_un = P.kcat4(7).*c_C16_FabG_NADPH_BKeACP_un;
d_C18_BHyAcACP_un = P.kcat4(8).*c_C18_FabG_NADPH_BKeACP_un;
d_C20_BHyAcACP_un = P.kcat4(9).*c_C20_FabG_NADPH_BKeACP_un;

% FabG-NADPH
d_FabG_NADPH = P.k4_1f(1).*c_FabG.*c_NADPH - P.k4_1r(1).*c_FabG_NADPH...
 + P.k4_2r(1).*c_C4_FabG_NADPH_BKeACP   - P.k4_2f(1).*c_FabG_NADPH.*c_C4_BKeACP...
 + P.k4_2r(2).*c_C6_FabG_NADPH_BKeACP   - P.k4_2f(2).*c_FabG_NADPH.*c_C6_BKeACP...
 + P.k4_2r(3).*c_C8_FabG_NADPH_BKeACP   - P.k4_2f(3).*c_FabG_NADPH.*c_C8_BKeACP...
 + P.k4_2r(4).*c_C10_FabG_NADPH_BKeACP - P.k4_2f(4).*c_FabG_NADPH.*c_C10_BKeACP...
 + P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP - P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP...
 + P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP - P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP...
 + P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP - P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP...
 + P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP - P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP...
 + P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP - P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP...
 + P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP_un - P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP_un...
 + P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP_un - P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP_un...
 + P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP_un - P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP_un...
 + P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP_un - P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP_un...
 + P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP_un - P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP_un; 

% C2n (n=2:10) FabG-NADPH-B-ketoacyl-ACPs
d_C4_FabG_NADPH_BKeACP   = P.k4_2f(1).*c_FabG_NADPH.*c_C4_BKeACP   - P.k4_2r(1).*c_C4_FabG_NADPH_BKeACP   - P.kcat4(1).*c_C4_FabG_NADPH_BKeACP;
d_C6_FabG_NADPH_BKeACP   = P.k4_2f(2).*c_FabG_NADPH.*c_C6_BKeACP   - P.k4_2r(2).*c_C6_FabG_NADPH_BKeACP   - P.kcat4(2).*c_C6_FabG_NADPH_BKeACP;
d_C8_FabG_NADPH_BKeACP   = P.k4_2f(3).*c_FabG_NADPH.*c_C8_BKeACP   - P.k4_2r(3).*c_C8_FabG_NADPH_BKeACP   - P.kcat4(3).*c_C8_FabG_NADPH_BKeACP;
d_C10_FabG_NADPH_BKeACP = P.k4_2f(4).*c_FabG_NADPH.*c_C10_BKeACP - P.k4_2r(4).*c_C10_FabG_NADPH_BKeACP - P.kcat4(4).*c_C10_FabG_NADPH_BKeACP;
d_C12_FabG_NADPH_BKeACP = P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP - P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP - P.kcat4(5).*c_C12_FabG_NADPH_BKeACP;
d_C14_FabG_NADPH_BKeACP = P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP - P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP - P.kcat4(6).*c_C14_FabG_NADPH_BKeACP;
d_C16_FabG_NADPH_BKeACP = P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP - P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP - P.kcat4(7).*c_C16_FabG_NADPH_BKeACP;
d_C18_FabG_NADPH_BKeACP = P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP - P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP - P.kcat4(8).*c_C18_FabG_NADPH_BKeACP;
d_C20_FabG_NADPH_BKeACP = P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP - P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP - P.kcat4(9).*c_C20_FabG_NADPH_BKeACP;

% C2n:1 (n=6:10) FabG-NADPH-B-ketoacyl-ACPs
d_C12_FabG_NADPH_BKeACP_un = P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP_un - P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP_un - P.kcat4(5).*c_C12_FabG_NADPH_BKeACP_un;
d_C14_FabG_NADPH_BKeACP_un = P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP_un - P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP_un - P.kcat4(6).*c_C14_FabG_NADPH_BKeACP_un;
d_C16_FabG_NADPH_BKeACP_un = P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP_un - P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP_un - P.kcat4(7).*c_C16_FabG_NADPH_BKeACP_un;
d_C18_FabG_NADPH_BKeACP_un = P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP_un - P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP_un - P.kcat4(8).*c_C18_FabG_NADPH_BKeACP_un;
d_C20_FabG_NADPH_BKeACP_un = P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP_un - P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP_un - P.kcat4(9).*c_C20_FabG_NADPH_BKeACP_un;

% FabG-ACP
d_FabG_ACP = P.k4_inh_f.*c_FabG.*c_ACP - P.k4_inh_r.*c_FabG_ACP;


dcdt = [d_ACP;d_NADPH;d_C4_BKeACP;d_C6_BKeACP;d_C8_BKeACP;d_C10_BKeACP;...
    d_C12_BKeACP;d_C14_BKeACP;d_C16_BKeACP;d_C18_BKeACP;d_C20_BKeACP;...
    d_C12_BKeACP_un;d_C14_BKeACP_un;d_C16_BKeACP_un;d_C18_BKeACP_un;...
    d_C20_BKeACP_un;d_C4_BHyAcACP;d_C6_BHyAcACP;d_C8_BHyAcACP;d_C10_BHyAcACP;...
    d_C12_BHyAcACP;d_C14_BHyAcACP;d_C16_BHyAcACP;d_C18_BHyAcACP;...
    d_C20_BHyAcACP;d_C12_BHyAcACP_un;d_C14_BHyAcACP_un;d_C16_BHyAcACP_un;...
    d_C18_BHyAcACP_un;d_C20_BHyAcACP_un;d_FabG_NADPH;d_C4_FabG_NADPH_BKeACP;...
    d_C6_FabG_NADPH_BKeACP;d_C8_FabG_NADPH_BKeACP;d_C10_FabG_NADPH_BKeACP;...
    d_C12_FabG_NADPH_BKeACP;d_C14_FabG_NADPH_BKeACP;d_C16_FabG_NADPH_BKeACP;...
    d_C18_FabG_NADPH_BKeACP;d_C20_FabG_NADPH_BKeACP;...
    d_C12_FabG_NADPH_BKeACP_un;d_C14_FabG_NADPH_BKeACP_un;...
    d_C16_FabG_NADPH_BKeACP_un;d_C18_FabG_NADPH_BKeACP_un;...
    d_C20_FabG_NADPH_BKeACP_un;d_FabG_ACP];

end
