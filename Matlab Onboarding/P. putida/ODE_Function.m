function dcdt = ODE_Function(t,c,P)
% Contains all the differential equations and enzyme balances that define
% the FAS model
% Input:
% t: time values (required as input for MATLAB ODE solver, sec)
% c: concentration values (all components and intermediates, uM)
% P: structure containing all kinetic parameters
% num: number of variables (y)
% Output:
% dcdt: values of differential equations for given 
% concentrations and kinetic parameters

%% ODEs

% Assign values to each variable
for var = 1:numel(P.labels)
    % Get the variable name
    var_name = P.labels{var};
    
    % Assign the value to the variable in the workspace
    eval([var_name ' = c(var, :);']);
end

% Enzyme concentration balance equations
% ACC
c_ACC = P.ACCtot - c_ACC_s1 - c_C1_ACC_s2 - c_C1_ACC_s3 - c_C3_ACC_s4;

% FabD
c_FabD = P.FabDtot - c_C3_FabD_MalCoA - c_C3_FabD_Act - c_C3_FabD_Act_ACP;

% FabH % changed
c_FabH = P.FabHtot - c_FabH_ACP...
    - c_C2_FabH_CoA - c_C4_FabH_CoA - c_C6_FabH_CoA - c_C8_FabH_CoA - c_C10_FabH_CoA - c_C12_FabH_CoA - c_C14_FabH_CoA - c_C16_FabH_CoA - c_C18_FabH_CoA...
    - c_C2_FabH_Act - c_C4_FabH_Act - c_C6_FabH_Act - c_C8_FabH_Act - c_C10_FabH_Act - c_C12_FabH_Act - c_C14_FabH_Act - c_C16_FabH_Act - c_C18_FabH_Act...
    - c_C5_FabH_Act_MalACP - c_C7_FabH_Act_MalACP - c_C9_FabH_Act_MalACP - c_C11_FabH_Act_MalACP - c_C13_FabH_Act_MalACP - c_C15_FabH_Act_MalACP - c_C17_FabH_Act_MalACP - c_C19_FabH_Act_MalACP - c_C21_FabH_Act_MalACP...
    - c_C4_FabH_AcACP - c_C6_FabH_AcACP - c_C8_FabH_AcACP - c_C10_FabH_AcACP - c_C12_FabH_AcACP - c_C14_FabH_AcACP - c_C16_FabH_AcACP - c_C18_FabH_AcACP - c_C20_FabH_AcACP...
    - c_C12_FabH_AcACP_un - c_C14_FabH_AcACP_un - c_C16_FabH_AcACP_un - c_C18_FabH_AcACP_un - c_C20_FabH_AcACP_un...
    - c_C6_FabH_Act_AcACP - c_C8_FabH_Act_AcACP - c_C10_FabH_Act_AcACP - c_C12_FabH_Act_AcACP - c_C14_FabH_Act_AcACP - c_C16_FabH_Act_AcACP - c_C18_FabH_Act_AcACP - c_C20_FabH_Act_AcACP - c_C22_FabH_Act_AcACP...
    - c_C14_FabH_Act_AcACP_un - c_C16_FabH_Act_AcACP_un - c_C18_FabH_Act_AcACP_un - c_C20_FabH_Act_AcACP_un - c_C22_FabH_Act_AcACP_un;

% FabG
c_FabG = P.FabGtot - c_FabG_NADPH...
    - c_C4_FabG_NADPH_BKeACP - c_C6_FabG_NADPH_BKeACP - c_C8_FabG_NADPH_BKeACP - c_C10_FabG_NADPH_BKeACP - c_C12_FabG_NADPH_BKeACP - c_C14_FabG_NADPH_BKeACP - c_C16_FabG_NADPH_BKeACP - c_C18_FabG_NADPH_BKeACP - c_C20_FabG_NADPH_BKeACP...
    - c_C12_FabG_NADPH_BKeACP_un - c_C14_FabG_NADPH_BKeACP_un - c_C16_FabG_NADPH_BKeACP_un - c_C18_FabG_NADPH_BKeACP_un - c_C20_FabG_NADPH_BKeACP_un - c_FabG_ACP;

% FabZ
c_FabZ = P.FabZtot - c_FabZ_ACP...
    - c_C4_FabZ_BHyAcACP - c_C6_FabZ_BHyAcACP - c_C8_FabZ_BHyAcACP - c_C10_FabZ_BHyAcACP - c_C12_FabZ_BHyAcACP - c_C14_FabZ_BHyAcACP - c_C16_FabZ_BHyAcACP - c_C18_FabZ_BHyAcACP - c_C20_FabZ_BHyAcACP...
    - c_C12_FabZ_BHyAcACP_un - c_C14_FabZ_BHyAcACP_un - c_C16_FabZ_BHyAcACP_un - c_C18_FabZ_BHyAcACP_un - c_C20_FabZ_BHyAcACP_un...
    - c_C4_FabZ_EnAcACP - c_C6_FabZ_EnAcACP - c_C8_FabZ_EnAcACP - c_C10_FabZ_EnAcACP - c_C12_FabZ_EnAcACP - c_C14_FabZ_EnAcACP - c_C16_FabZ_EnAcACP - c_C18_FabZ_EnAcACP - c_C20_FabZ_EnAcACP...
    - c_C12_FabZ_EnAcACP_un - c_C14_FabZ_EnAcACP_un - c_C16_FabZ_EnAcACP_un - c_C18_FabZ_EnAcACP_un - c_C20_FabZ_EnAcACP_un;

% FabI
c_FabI = P.FabItot - c_FabI_NADH...
    - c_C4_FabI_NADH_EnAcACP - c_C6_FabI_NADH_EnAcACP - c_C8_FabI_NADH_EnAcACP - c_C10_FabI_NADH_EnAcACP - c_C12_FabI_NADH_EnAcACP - c_C14_FabI_NADH_EnAcACP - c_C16_FabI_NADH_EnAcACP - c_C18_FabI_NADH_EnAcACP - c_C20_FabI_NADH_EnAcACP...
    - c_C12_FabI_NADH_EnAcACP_un - c_C14_FabI_NADH_EnAcACP_un - c_C16_FabI_NADH_EnAcACP_un - c_C18_FabI_NADH_EnAcACP_un - c_C20_FabI_NADH_EnAcACP_un - c_FabI_ACP;

% TesA
c_TesA = P.TesAtot- c_TesA_ACP...
    - c_C4_TesA_AcACP - c_C6_TesA_AcACP - c_C8_TesA_AcACP - c_C10_TesA_AcACP - c_C12_TesA_AcACP - c_C14_TesA_AcACP - c_C16_TesA_AcACP - c_C18_TesA_AcACP - c_C20_TesA_AcACP...
    - c_C12_TesA_AcACP_un - c_C14_TesA_AcACP_un - c_C16_TesA_AcACP_un - c_C18_TesA_AcACP_un - c_C20_TesA_AcACP_un;

% FabF
c_FabF = P.FabFtot - c_FabF_ACP - c_C2_FabF_AcCoA - c_C2_FabF_Act - c_C3_FabF_MalACP - c_C2_FabF_AcACP...
    - c_C4_FabF_AcACP - c_C6_FabF_AcACP - c_C8_FabF_AcACP - c_C10_FabF_AcACP - c_C12_FabF_AcACP - c_C14_FabF_AcACP - c_C16_FabF_AcACP - c_C18_FabF_AcACP...
    - c_C12_FabF_AcACP_un - c_C14_FabF_AcACP_un - c_C16_FabF_AcACP_un - c_C18_FabF_AcACP_un...
    - c_C4_FabF_Act - c_C6_FabF_Act - c_C8_FabF_Act - c_C10_FabF_Act - c_C12_FabF_Act - c_C14_FabF_Act - c_C16_FabF_Act - c_C18_FabF_Act...
    - c_C12_FabF_Act_un - c_C14_FabF_Act_un - c_C16_FabF_Act_un - c_C18_FabF_Act_un...
    - c_C5_FabF_Act_MalACP - c_C7_FabF_Act_MalACP - c_C9_FabF_Act_MalACP - c_C11_FabF_Act_MalACP - c_C13_FabF_Act_MalACP - c_C15_FabF_Act_MalACP - c_C17_FabF_Act_MalACP - c_C19_FabF_Act_MalACP - c_C21_FabF_Act_MalACP...
    - c_C15_FabF_Act_MalACP_un - c_C17_FabF_Act_MalACP_un - c_C19_FabF_Act_MalACP_un - c_C21_FabF_Act_MalACP_un;

% FabA
c_FabA = P.FabAtot - c_FabA_ACP - c_C10_FabA_cis3EnAcACP...
    - c_C4_FabA_BHyAcACP - c_C6_FabA_BHyAcACP - c_C8_FabA_BHyAcACP - c_C10_FabA_BHyAcACP - c_C12_FabA_BHyAcACP - c_C14_FabA_BHyAcACP - c_C16_FabA_BHyAcACP - c_C18_FabA_BHyAcACP - c_C20_FabA_BHyAcACP...
    - c_C12_FabA_BHyAcACP_un - c_C14_FabA_BHyAcACP_un - c_C16_FabA_BHyAcACP_un - c_C18_FabA_BHyAcACP_un - c_C20_FabA_BHyAcACP_un...
    - c_C4_FabA_EnAcACP - c_C6_FabA_EnAcACP - c_C8_FabA_EnAcACP - c_C10_FabA_EnAcACP - c_C12_FabA_EnAcACP - c_C14_FabA_EnAcACP - c_C16_FabA_EnAcACP - c_C18_FabA_EnAcACP - c_C20_FabA_EnAcACP...
    - c_C12_FabA_EnAcACP_un - c_C14_FabA_EnAcACP_un - c_C16_FabA_EnAcACP_un - c_C18_FabA_EnAcACP_un - c_C20_FabA_EnAcACP_un;

% FabB
c_FabB = P.FabBtot - c_FabB_ACP - c_C2_FabB_AcCoA - c_C2_FabB_Act - c_C5_FabB_Act_MalACP - c_C3_FabB_MalACP - c_C2_FabB_AcACP...
    - c_C4_FabB_AcACP - c_C6_FabB_AcACP - c_C8_FabB_AcACP - c_C10_FabB_AcACP - c_C12_FabB_AcACP - c_C14_FabB_AcACP - c_C16_FabB_AcACP - c_C18_FabB_AcACP...
    - c_C12_FabB_AcACP_un - c_C14_FabB_AcACP_un - c_C16_FabB_AcACP_un - c_C18_FabB_AcACP_un...
    - c_C4_FabB_Act - c_C6_FabB_Act - c_C8_FabB_Act - c_C10_FabB_Act - c_C12_FabB_Act - c_C14_FabB_Act - c_C16_FabB_Act - c_C18_FabB_Act...
    - c_C12_FabB_Act_un - c_C14_FabB_Act_un - c_C16_FabB_Act_un - c_C18_FabB_Act_un...
    - c_C7_FabB_Act_MalACP - c_C9_FabB_Act_MalACP - c_C11_FabB_Act_MalACP - c_C13_FabB_Act_MalACP - c_C15_FabB_Act_MalACP - c_C17_FabB_Act_MalACP - c_C19_FabB_Act_MalACP - c_C21_FabB_Act_MalACP...
    - c_C15_FabB_Act_MalACP_un - c_C17_FabB_Act_MalACP_un - c_C19_FabB_Act_MalACP_un - c_C21_FabB_Act_MalACP_un...
    - c_C10_FabB_cis3EnAcACP - c_C10_FabB_Act_cis3 - c_C13_FabB_Act_cis3MalACP;

% Set of differential equations
% ATP
d_ATP = P.k1_1r.*c_ACC_s1 - P.k1_1f.*c_ACC.*c_ATP;
%d_ATP = 0*c_ATP;

% Bicarbonate
d_C1_Bicarbonate = P.k1_2r.*c_C1_ACC_s2 - P.k1_2f.*c_ACC_s1.*c_C1_Bicarbonate;
%d_C1_Bicarbonate = 0*c_C1_Bicarbonate;

% C2n (n=1:9)-CoA % changed 
d_C2_AcCoA        = P.k3_1r(1).*c_C2_FabH_CoA   - P.k3_1f(1).*c_FabH.*c_C2_AcCoA + P.k1_3r.*c_C3_ACC_s4 - P.k1_3f.*c_C1_ACC_s3.*c_C2_AcCoA + P.k8_4r.*c_C2_FabF_AcCoA - P.k8_4f.*c_FabF.*c_C2_AcCoA + P.k10_4r.*c_C2_FabB_AcCoA - P.k10_4f.*c_FabB.*c_C2_AcCoA;
%d_C2_AcCoA        = P.k3_1r(1).*c_C2_FabH_CoA   - P.k3_1f(1).*c_FabH.*c_C2_AcCoA;
d_C4_SucCoA      = P.k3_1r(2).*c_C4_FabH_CoA   - P.k3_1f(2).*c_FabH.*c_C4_SucCoA;
d_C6_HexCoA      = P.k3_1r(3).*c_C6_FabH_CoA   - P.k3_1f(3).*c_FabH.*c_C6_HexCoA;
d_C8_OcCoA        = P.k3_1r(4).*c_C8_FabH_CoA   - P.k3_1f(4).*c_FabH.*c_C8_OcCoA;
d_C10_DecCoA    = P.k3_1r(5).*c_C10_FabH_CoA  - P.k3_1f(5).*c_FabH.*c_C10_DecCoA;
d_C12_LauCoA     = P.k3_1r(6).*c_C12_FabH_CoA - P.k3_1f(6).*c_FabH.*c_C12_LauCoA;
d_C14_EthCoA     = P.k3_1r(7).*c_C14_FabH_CoA  - P.k3_1f(7).*c_FabH.*c_C14_EthCoA;
d_C16_PalCoA      = P.k3_1r(8).*c_C16_FabH_CoA - P.k3_1f(8).*c_FabH.*c_C16_PalCoA;
d_C18_OcDecCoA = P.k3_1r(9).*c_C18_FabH_CoA - P.k3_1f(9).*c_FabH.*c_C18_OcDecCoA;

% ACP
d_ACP = P.k2_3r.*c_C3_FabD_Act_ACP - P.k2_3f.*c_C3_FabD_Act.*c_ACP...
    + P.kcat7(1).*c_C4_TesA_AcACP + P.kcat7(2).*c_C6_TesA_AcACP + P.kcat7(3).*c_C8_TesA_AcACP + P.kcat7(4).*c_C10_TesA_AcACP...
    + P.kcat7(5).*c_C12_TesA_AcACP + P.kcat7(6).*c_C14_TesA_AcACP + P.kcat7(7).*c_C16_TesA_AcACP + P.kcat7(8).*c_C18_TesA_AcACP + P.kcat7(9).*c_C20_TesA_AcACP...
    + P.kcat7(5).*c_C12_TesA_AcACP_un + P.kcat7(6).*c_C14_TesA_AcACP_un + P.kcat7(7).*c_C16_TesA_AcACP_un + P.kcat7(8).*c_C18_TesA_AcACP_un + P.kcat7(9).*c_C20_TesA_AcACP_un...
    + P.k8_2f(1).*c_C4_FabF_AcACP    - P.k8_2r(1).*c_C4_FabF_Act.*c_ACP...
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
    + P.k10_2f(1).*c_C4_FabB_AcACP   - P.k10_2r(1).*c_C4_FabB_Act.*c_ACP...
    + P.k10_2f(2).*c_C6_FabB_AcACP   - P.k10_2r(2).*c_C6_FabB_Act.*c_ACP...
    + P.k10_2f(3).*c_C8_FabB_AcACP   - P.k10_2r(3).*c_C8_FabB_Act.*c_ACP...
    + P.k10_2f(4).*c_C10_FabB_AcACP - P.k10_2r(4).*c_C10_FabB_Act.*c_ACP...
    + P.k10_2f(5).*c_C12_FabB_AcACP - P.k10_2r(5).*c_C12_FabB_Act.*c_ACP...
    + P.k10_2f(6).*c_C14_FabB_AcACP - P.k10_2r(6).*c_C14_FabB_Act.*c_ACP...
    + P.k10_2f(7).*c_C16_FabB_AcACP - P.k10_2r(7).*c_C16_FabB_Act.*c_ACP...
    + P.k10_2f(8).*c_C18_FabB_AcACP - P.k10_2r(8).*c_C18_FabB_Act.*c_ACP...
    + P.k10_2f(5).*c_C12_FabB_AcACP_un - P.k10_2r(5).*c_C12_FabB_Act_un.*c_ACP...
    + P.k10_2f(6).*c_C14_FabB_AcACP_un - P.k10_2r(6).*c_C14_FabB_Act_un.*c_ACP...
    + P.k10_2f(7).*c_C16_FabB_AcACP_un - P.k10_2r(7).*c_C16_FabB_Act_un.*c_ACP...
    + P.k10_2f(8).*c_C18_FabB_AcACP_un - P.k10_2r(8).*c_C18_FabB_Act_un.*c_ACP...
    + P.k10_2f(4).*c_C10_FabB_cis3EnAcACP - P.k10_2r(4).*c_C10_FabB_Act_cis3.*c_ACP...
    + P.k3_inh_r.*c_FabH_ACP   - P.k3_inh_f.*c_FabH.*c_ACP...
    + P.k4_inh_r.*c_FabG_ACP   - P.k4_inh_f.*c_FabG.*c_ACP...
    + P.k5_inh_r.*c_FabZ_ACP   - P.k5_inh_f.*c_FabZ.*c_ACP...
    + P.k6_inh_r.*c_FabI_ACP    - P.k6_inh_f.*c_FabI.*c_ACP...
    + P.k7_inh_r.*c_TesA_ACP   - P.k7_inh_f.*c_TesA.*c_ACP...
    + P.k8_inh_r.*c_FabF_ACP   - P.k8_inh_f.*c_FabF.*c_ACP...
    + P.k9_inh_r.*c_FabA_ACP   - P.k9_inh_f.*c_FabA.*c_ACP...
    + P.k10_inh_r.*c_FabB_ACP - P.k10_inh_f.*c_FabB.*c_ACP...
    + P.k8_9f.*c_C2_FabF_AcACP   - P.k8_9r.*c_C2_FabF_Act.*c_ACP...
    + P.k10_9f.*c_C2_FabB_AcACP - P.k10_9r.*c_C2_FabB_Act.*c_ACP; 

% NADPH
d_NADPH = P.k4_1r(1).*c_FabG_NADPH - P.k4_1f(1).*c_FabG.*c_NADPH; 

% % NADP+
% d_NADP = P.kcat4(1).*c_C4_FabG_NADPH_BKeACP + P.kcat4(2).*c_C6_FabG_NADPH_BKeACP + P.kcat4(3).*c_C8_FabG_NADPH_BKeACP + P.kcat4(4).*c_C10_FabG_NADPH_BKeACP...
%     + P.kcat4(5).*c_C12_FabG_NADPH_BKeACP + P.kcat4(6).*c_C14_FabG_NADPH_BKeACP + P.kcat4(7).*c_C16_FabG_NADPH_BKeACP + P.kcat4(8).*c_C18_FabG_NADPH_BKeACP + P.kcat4(9).*c_C20_FabG_NADPH_BKeACP; 

% NADH
d_NADH = P.k6_1r(1).*c_FabI_NADH - P.k6_1f(1).*c_FabI.*c_NADH; 

% NAD+
% d_NAD = P.kcat6(1).*c_C4_FabI_NADH_EnAcACP + P.kcat6(2).*c_C6_FabI_NADH_EnAcACP + P.kcat6(3).*c_C8_FabI_NADH_EnAcACP + P.kcat6(4).*c_C10_FabI_NADH_EnAcACP...
%     + P.kcat6(5).*c_C12_FabI_NADH_EnAcACP + P.kcat6(6).*c_C14_FabI_NADH_EnAcACP + P.kcat6(7).*c_C16_FabI_NADH_EnAcACP + P.kcat6(8).*c_C18_FabI_NADH_EnAcACP + P.kcat6(9).*c_C20_FabI_NADH_EnAcACP; 

% ADP
d_ADP = P.kcat1_2.*c_C3_ACC_s4;
%d_ADP = 0*c_ADP;

% Malonyl-CoA
d_C3_MalCoA = P.kcat1_2.*c_C3_ACC_s4 + P.k2_1r.*c_C3_FabD_MalCoA - P.k2_1f.*c_FabD.*c_C3_MalCoA; 

% CoA % changed
d_CoA = P.k2_2f.*c_C3_FabD_MalCoA - P.k2_2r.*c_C3_FabD_Act.*c_CoA...
    + P.k3_2f(1).*c_C2_FabH_CoA   - P.k3_2r(1).*c_C2_FabH_Act.*c_CoA...
    + P.k3_2f(2).*c_C4_FabH_CoA   - P.k3_2r(2).*c_C4_FabH_Act.*c_CoA...
    + P.k3_2f(3).*c_C6_FabH_CoA   - P.k3_2r(3).*c_C6_FabH_Act.*c_CoA...
    + P.k3_2f(4).*c_C8_FabH_CoA   - P.k3_2r(4).*c_C8_FabH_Act.*c_CoA...
    + P.k3_2f(5).*c_C10_FabH_CoA - P.k3_2r(5).*c_C10_FabH_Act.*c_CoA...
    + P.k3_2f(6).*c_C12_FabH_CoA - P.k3_2r(6).*c_C12_FabH_Act.*c_CoA...
    + P.k3_2f(7).*c_C14_FabH_CoA - P.k3_2r(7).*c_C14_FabH_Act.*c_CoA...
    + P.k3_2f(8).*c_C16_FabH_CoA - P.k3_2r(8).*c_C16_FabH_Act.*c_CoA...
    + P.k3_2f(9).*c_C18_FabH_CoA - P.k3_2r(9).*c_C18_FabH_Act.*c_CoA...
    + P.k8_5f.*c_C2_FabF_AcCoA   - P.k8_5r*c_C2_FabF_Act.*c_CoA...
    + P.k10_5f.*c_C2_FabB_AcCoA - P.k10_5r.*c_C2_FabB_Act.*c_CoA;

% Malonyl-ACP % changed
d_C3_MalACP = P.k2_4f.*c_C3_FabD_Act_ACP - P.k2_4r.*c_FabD.*c_C3_MalACP...
    + P.k3_3r(1).*c_C5_FabH_Act_MalACP   - P.k3_3f(1).*c_C2_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(2).*c_C7_FabH_Act_MalACP   - P.k3_3f(2).*c_C4_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(3).*c_C9_FabH_Act_MalACP   - P.k3_3f(3).*c_C6_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(4).*c_C11_FabH_Act_MalACP - P.k3_3f(4).*c_C8_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(5).*c_C13_FabH_Act_MalACP - P.k3_3f(5).*c_C10_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(6).*c_C15_FabH_Act_MalACP - P.k3_3f(6).*c_C12_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(7).*c_C17_FabH_Act_MalACP - P.k3_3f(7).*c_C14_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(8).*c_C19_FabH_Act_MalACP - P.k3_3f(8).*c_C16_FabH_Act.*c_C3_MalACP...
    + P.k3_3r(9).*c_C21_FabH_Act_MalACP - P.k3_3f(9).*c_C18_FabH_Act.*c_C3_MalACP...
    + P.k8_3r(1).*c_C7_FabF_Act_MalACP   - P.k8_3f(1).*c_C4_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(2).*c_C9_FabF_Act_MalACP   - P.k8_3f(2).*c_C6_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(3).*c_C11_FabF_Act_MalACP - P.k8_3f(3).*c_C8_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(4).*c_C13_FabF_Act_MalACP - P.k8_3f(4).*c_C10_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(5).*c_C15_FabF_Act_MalACP - P.k8_3f(5).*c_C12_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(6).*c_C17_FabF_Act_MalACP - P.k8_3f(6).*c_C14_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(7).*c_C19_FabF_Act_MalACP - P.k8_3f(7).*c_C16_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(8).*c_C21_FabF_Act_MalACP - P.k8_3f(8).*c_C18_FabF_Act.*c_C3_MalACP...
    + P.k8_3r(5).*c_C15_FabF_Act_MalACP_un - P.k8_3f(5).*c_C12_FabF_Act_un.*c_C3_MalACP...
    + P.k8_3r(6).*c_C17_FabF_Act_MalACP_un - P.k8_3f(6).*c_C14_FabF_Act_un.*c_C3_MalACP...
    + P.k8_3r(7).*c_C19_FabF_Act_MalACP_un - P.k8_3f(7).*c_C16_FabF_Act_un.*c_C3_MalACP...
    + P.k8_3r(8).*c_C21_FabF_Act_MalACP_un - P.k8_3f(8).*c_C18_FabF_Act_un.*c_C3_MalACP...
    + P.k10_3r(1).*c_C7_FabB_Act_MalACP   - P.k10_3f(1).*c_C4_FabB_Act.*c_C3_MalACP...
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
    + P.k10_3r(4).*c_C13_FabB_Act_cis3MalACP - P.k10_3f(4).*c_C10_FabB_Act_cis3.*c_C3_MalACP...
    + P.k8_6r.*c_C5_FabF_Act_MalACP   - P.k8_6f.*c_C2_FabF_Act.*c_C3_MalACP...
    + P.k10_6r.*c_C5_FabB_Act_MalACP - P.k10_6f.*c_C2_FabB_Act.*c_C3_MalACP...
    + P.k8_7r.*c_C3_FabF_MalACP   - P.k8_7f.*c_FabF.*c_C3_MalACP...
    + P.k10_7r.*c_C3_FabB_MalACP - P.k10_7f.*c_FabB.*c_C3_MalACP;

% CO2 % changed
d_C1_CO2 = P.kcat3(1).*c_C5_FabH_Act_MalACP + P.kcat3(2).*c_C7_FabH_Act_MalACP + P.kcat3(3).*c_C9_FabH_Act_MalACP + P.kcat3(4).*c_C11_FabH_Act_MalACP...
    + P.kcat3(5).*c_C13_FabH_Act_MalACP + P.kcat3(6).*c_C15_FabH_Act_MalACP + P.kcat3(7).*c_C17_FabH_Act_MalACP + P.kcat3(8).*c_C19_FabH_Act_MalACP + P.kcat3(9).*c_C21_FabH_Act_MalACP...
    + P.kcat8(1).*c_C7_FabF_Act_MalACP + P.kcat8(2).*c_C9_FabF_Act_MalACP + P.kcat8(3).*c_C11_FabF_Act_MalACP + P.kcat8(4).*c_C13_FabF_Act_MalACP...
    + P.kcat8(5).*c_C15_FabF_Act_MalACP + P.kcat8(6).*c_C17_FabF_Act_MalACP + P.kcat8(7).*c_C19_FabF_Act_MalACP + P.kcat8(8).*c_C21_FabF_Act_MalACP...
    + P.kcat8_un(5).*c_C15_FabF_Act_MalACP_un + P.kcat8_un(6).*c_C17_FabF_Act_MalACP_un + P.kcat8_un(7).*c_C19_FabF_Act_MalACP_un + P.kcat8_un(8).*c_C21_FabF_Act_MalACP_un...
    + P.kcat10(1).*c_C7_FabB_Act_MalACP + P.kcat10(2).*c_C9_FabB_Act_MalACP + P.kcat10(3).*c_C11_FabB_Act_MalACP + P.kcat10(4).*c_C13_FabB_Act_MalACP...
    + P.kcat10(5).*c_C15_FabB_Act_MalACP + P.kcat10(6).*c_C17_FabB_Act_MalACP + P.kcat10(7).*c_C19_FabB_Act_MalACP + P.kcat10(8).*c_C21_FabB_Act_MalACP...
    + P.kcat10_un(4).*c_C13_FabB_Act_cis3MalACP + P.kcat10_un(5).*c_C15_FabB_Act_MalACP_un + P.kcat10_un(6).*c_C17_FabB_Act_MalACP_un + P.kcat10_un(7).*c_C19_FabB_Act_MalACP_un + P.kcat10_un(8).*c_C21_FabB_Act_MalACP_un...
    + P.kcat8_H.*c_C5_FabF_Act_MalACP + P.kcat8_CO2.*c_C3_FabF_MalACP...
    + P.kcat10_H.*c_C5_FabB_Act_MalACP + P.kcat10_CO2.*c_C3_FabB_MalACP;

% C2n (n=2:10) B-ketoacyl-ACPs (FabH + FabF + FabB - FabG) % changed
d_C4_BKeACP    = P.kcat3(1).*c_C5_FabH_Act_MalACP   + P.kcat8_H.*c_C5_FabF_Act_MalACP   + P.kcat10_H.*c_C5_FabB_Act_MalACP         + P.k4_2r(1).*c_C4_FabG_NADPH_BKeACP   - P.k4_2f(1).*c_FabG_NADPH.*c_C4_BKeACP;
d_C6_BKeACP    = P.kcat3(2).*c_C7_FabH_Act_MalACP   + P.kcat8(1).*c_C7_FabF_Act_MalACP  + P.kcat10(1).*c_C7_FabB_Act_MalACP   + P.k4_2r(2).*c_C6_FabG_NADPH_BKeACP   - P.k4_2f(2).*c_FabG_NADPH.*c_C6_BKeACP;
d_C8_BKeACP    = P.kcat3(3).*c_C9_FabH_Act_MalACP   + P.kcat8(2).*c_C9_FabF_Act_MalACP  + P.kcat10(2).*c_C9_FabB_Act_MalACP   + P.k4_2r(3).*c_C8_FabG_NADPH_BKeACP   - P.k4_2f(3).*c_FabG_NADPH.*c_C8_BKeACP;
d_C10_BKeACP  = P.kcat3(4).*c_C11_FabH_Act_MalACP   + P.kcat8(3).*c_C11_FabF_Act_MalACP  + P.kcat10(3).*c_C11_FabB_Act_MalACP   + P.k4_2r(4).*c_C10_FabG_NADPH_BKeACP - P.k4_2f(4).*c_FabG_NADPH.*c_C10_BKeACP;
d_C12_BKeACP  = P.kcat3(5).*c_C13_FabH_Act_MalACP + P.kcat8(4).*c_C13_FabF_Act_MalACP + P.kcat10(4).*c_C13_FabB_Act_MalACP + P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP - P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP;
d_C14_BKeACP  = P.kcat3(6).*c_C15_FabH_Act_MalACP + P.kcat8(5).*c_C15_FabF_Act_MalACP + P.kcat10(5).*c_C15_FabB_Act_MalACP + P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP - P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP;
d_C16_BKeACP  = P.kcat3(7).*c_C17_FabH_Act_MalACP + P.kcat8(6).*c_C17_FabF_Act_MalACP + P.kcat10(6).*c_C17_FabB_Act_MalACP + P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP - P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP;
d_C18_BKeACP  = P.kcat3(8).*c_C19_FabH_Act_MalACP + P.kcat8(7).*c_C19_FabF_Act_MalACP + P.kcat10(7).*c_C19_FabB_Act_MalACP + P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP - P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP;
d_C20_BKeACP  = P.kcat3(9).*c_C21_FabH_Act_MalACP + P.kcat8(8).*c_C21_FabF_Act_MalACP + P.kcat10(8).*c_C21_FabB_Act_MalACP + P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP - P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP;

% C2n:1 (n=6:10) B-ketoacyl-ACPs (FabF + FabB - FabG)
d_C12_BKeACP_un =                                                                      P.kcat10_un(4).*c_C13_FabB_Act_cis3MalACP + P.k4_2r(5).*c_C12_FabG_NADPH_BKeACP_un - P.k4_2f(5).*c_FabG_NADPH.*c_C12_BKeACP_un;
d_C14_BKeACP_un = P.kcat8_un(5).*c_C15_FabF_Act_MalACP_un + P.kcat10_un(5).*c_C15_FabB_Act_MalACP_un + P.k4_2r(6).*c_C14_FabG_NADPH_BKeACP_un - P.k4_2f(6).*c_FabG_NADPH.*c_C14_BKeACP_un;
d_C16_BKeACP_un = P.kcat8_un(6).*c_C17_FabF_Act_MalACP_un + P.kcat10_un(6).*c_C17_FabB_Act_MalACP_un + P.k4_2r(7).*c_C16_FabG_NADPH_BKeACP_un - P.k4_2f(7).*c_FabG_NADPH.*c_C16_BKeACP_un;
d_C18_BKeACP_un = P.kcat8_un(7).*c_C19_FabF_Act_MalACP_un + P.kcat10_un(7).*c_C19_FabB_Act_MalACP_un + P.k4_2r(8).*c_C18_FabG_NADPH_BKeACP_un - P.k4_2f(8).*c_FabG_NADPH.*c_C18_BKeACP_un;
d_C20_BKeACP_un = P.kcat8_un(8).*c_C21_FabF_Act_MalACP_un + P.kcat10_un(8).*c_C21_FabB_Act_MalACP_un + P.k4_2r(9).*c_C20_FabG_NADPH_BKeACP_un - P.k4_2f(9).*c_FabG_NADPH.*c_C20_BKeACP_un;

% C2n (n=2:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C4_BHyAcACP   = P.kcat4(1).*c_C4_FabG_NADPH_BKeACP   + P.k5_1r(1).*c_C4_FabZ_BHyAcACP  - P.k5_1f(1).*c_FabZ.*c_C4_BHyAcACP   + P.k9_1r(1).*c_C4_FabA_BHyAcACP   - P.k9_1f(1).*c_FabA.*c_C4_BHyAcACP;
d_C6_BHyAcACP   = P.kcat4(2).*c_C6_FabG_NADPH_BKeACP   + P.k5_1r(2).*c_C6_FabZ_BHyAcACP  - P.k5_1f(2).*c_FabZ.*c_C6_BHyAcACP   + P.k9_1r(2).*c_C6_FabA_BHyAcACP   - P.k9_1f(2).*c_FabA.*c_C6_BHyAcACP;
d_C8_BHyAcACP   = P.kcat4(3).*c_C8_FabG_NADPH_BKeACP   + P.k5_1r(3).*c_C8_FabZ_BHyAcACP  - P.k5_1f(3).*c_FabZ.*c_C8_BHyAcACP   + P.k9_1r(3).*c_C8_FabA_BHyAcACP   - P.k9_1f(3).*c_FabA.*c_C8_BHyAcACP;
d_C10_BHyAcACP = P.kcat4(4).*c_C10_FabG_NADPH_BKeACP + P.k5_1r(4).*c_C10_FabZ_BHyAcACP - P.k5_1f(4).*c_FabZ.*c_C10_BHyAcACP + P.k9_1r(4).*c_C10_FabA_BHyAcACP - P.k9_1f(4).*c_FabA.*c_C10_BHyAcACP;
d_C12_BHyAcACP = P.kcat4(5).*c_C12_FabG_NADPH_BKeACP + P.k5_1r(5).*c_C12_FabZ_BHyAcACP - P.k5_1f(5).*c_FabZ.*c_C12_BHyAcACP + P.k9_1r(5).*c_C12_FabA_BHyAcACP - P.k9_1f(5).*c_FabA.*c_C12_BHyAcACP;
d_C14_BHyAcACP = P.kcat4(6).*c_C14_FabG_NADPH_BKeACP + P.k5_1r(6).*c_C14_FabZ_BHyAcACP - P.k5_1f(6).*c_FabZ.*c_C14_BHyAcACP + P.k9_1r(6).*c_C14_FabA_BHyAcACP - P.k9_1f(6).*c_FabA.*c_C14_BHyAcACP;
d_C16_BHyAcACP = P.kcat4(7).*c_C16_FabG_NADPH_BKeACP + P.k5_1r(7).*c_C16_FabZ_BHyAcACP - P.k5_1f(7).*c_FabZ.*c_C16_BHyAcACP + P.k9_1r(7).*c_C16_FabA_BHyAcACP - P.k9_1f(7).*c_FabA.*c_C16_BHyAcACP;
d_C18_BHyAcACP = P.kcat4(8).*c_C18_FabG_NADPH_BKeACP + P.k5_1r(8).*c_C18_FabZ_BHyAcACP - P.k5_1f(8).*c_FabZ.*c_C18_BHyAcACP + P.k9_1r(8).*c_C18_FabA_BHyAcACP - P.k9_1f(8).*c_FabA.*c_C18_BHyAcACP;
d_C20_BHyAcACP = P.kcat4(9).*c_C20_FabG_NADPH_BKeACP + P.k5_1r(9).*c_C20_FabZ_BHyAcACP - P.k5_1f(9).*c_FabZ.*c_C20_BHyAcACP + P.k9_1r(9).*c_C20_FabA_BHyAcACP - P.k9_1f(9).*c_FabA.*c_C20_BHyAcACP;

% C2n:1 (n=6:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C12_BHyAcACP_un = P.kcat4(5).*c_C12_FabG_NADPH_BKeACP_un + P.k5_1r(5).*c_C12_FabZ_BHyAcACP_un - P.k5_1f(5).*c_FabZ.*c_C12_BHyAcACP_un + P.k9_1r_un(5).*c_C12_FabA_BHyAcACP_un - P.k9_1f_un(5).*c_FabA.*c_C12_BHyAcACP_un;
d_C14_BHyAcACP_un = P.kcat4(6).*c_C14_FabG_NADPH_BKeACP_un + P.k5_1r(6).*c_C14_FabZ_BHyAcACP_un - P.k5_1f(6).*c_FabZ.*c_C14_BHyAcACP_un + P.k9_1r_un(6).*c_C14_FabA_BHyAcACP_un - P.k9_1f_un(6).*c_FabA.*c_C14_BHyAcACP_un;
d_C16_BHyAcACP_un = P.kcat4(7).*c_C16_FabG_NADPH_BKeACP_un + P.k5_1r(7).*c_C16_FabZ_BHyAcACP_un - P.k5_1f(7).*c_FabZ.*c_C16_BHyAcACP_un + P.k9_1r_un(7).*c_C16_FabA_BHyAcACP_un - P.k9_1f_un(7).*c_FabA.*c_C16_BHyAcACP_un;
d_C18_BHyAcACP_un = P.kcat4(8).*c_C18_FabG_NADPH_BKeACP_un + P.k5_1r(8).*c_C18_FabZ_BHyAcACP_un - P.k5_1f(8).*c_FabZ.*c_C18_BHyAcACP_un + P.k9_1r_un(8).*c_C18_FabA_BHyAcACP_un - P.k9_1f_un(8).*c_FabA.*c_C18_BHyAcACP_un;
d_C20_BHyAcACP_un = P.kcat4(9).*c_C20_FabG_NADPH_BKeACP_un + P.k5_1r(9).*c_C20_FabZ_BHyAcACP_un - P.k5_1f(9).*c_FabZ.*c_C20_BHyAcACP_un + P.k9_1r_un(9).*c_C20_FabA_BHyAcACP_un - P.k9_1f_un(9).*c_FabA.*c_C20_BHyAcACP_un;

% C2n (n=2:10) Enoyl-Acyl-ACPs (FabZ + FabA - FabI) 
d_C4_EnAcACP   = P.k5_3f(1).*c_C4_FabZ_EnAcACP   - P.k5_3r(1).*c_FabZ.*c_C4_EnAcACP  + P.k6_2r(1).*c_C4_FabI_NADH_EnAcACP   - P.k6_2f(1).*c_FabI_NADH.*c_C4_EnAcACP   + P.k9_3f(1).*c_C4_FabA_EnAcACP   - P.k9_3r(1).*c_FabA.*c_C4_EnAcACP;
d_C6_EnAcACP   = P.k5_3f(2).*c_C6_FabZ_EnAcACP   - P.k5_3r(2).*c_FabZ.*c_C6_EnAcACP  + P.k6_2r(2).*c_C6_FabI_NADH_EnAcACP   - P.k6_2f(2).*c_FabI_NADH.*c_C6_EnAcACP   + P.k9_3f(2).*c_C6_FabA_EnAcACP   - P.k9_3r(2).*c_FabA.*c_C6_EnAcACP;
d_C8_EnAcACP   = P.k5_3f(3).*c_C8_FabZ_EnAcACP   - P.k5_3r(3).*c_FabZ.*c_C8_EnAcACP  + P.k6_2r(3).*c_C8_FabI_NADH_EnAcACP   - P.k6_2f(3).*c_FabI_NADH.*c_C8_EnAcACP   + P.k9_3f(3).*c_C8_FabA_EnAcACP   - P.k9_3r(3).*c_FabA.*c_C8_EnAcACP;
d_C10_EnAcACP = P.k5_3f(4).*c_C10_FabZ_EnAcACP - P.k5_3r(4).*c_FabZ.*c_C10_EnAcACP + P.k6_2r(4).*c_C10_FabI_NADH_EnAcACP - P.k6_2f(4).*c_FabI_NADH.*c_C10_EnAcACP + P.k9_3f(4).*c_C10_FabA_EnAcACP - P.k9_3r(4).*c_FabA.*c_C10_EnAcACP;
d_C12_EnAcACP = P.k5_3f(5).*c_C12_FabZ_EnAcACP - P.k5_3r(5).*c_FabZ.*c_C12_EnAcACP + P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP - P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP + P.k9_3f(5).*c_C12_FabA_EnAcACP - P.k9_3r(5).*c_FabA.*c_C12_EnAcACP;
d_C14_EnAcACP = P.k5_3f(6).*c_C14_FabZ_EnAcACP - P.k5_3r(6).*c_FabZ.*c_C14_EnAcACP + P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP - P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP + P.k9_3f(6).*c_C14_FabA_EnAcACP - P.k9_3r(6).*c_FabA.*c_C14_EnAcACP;
d_C16_EnAcACP = P.k5_3f(7).*c_C16_FabZ_EnAcACP - P.k5_3r(7).*c_FabZ.*c_C16_EnAcACP + P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP - P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP + P.k9_3f(7).*c_C16_FabA_EnAcACP - P.k9_3r(7).*c_FabA.*c_C16_EnAcACP;
d_C18_EnAcACP = P.k5_3f(8).*c_C18_FabZ_EnAcACP - P.k5_3r(8).*c_FabZ.*c_C18_EnAcACP + P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP - P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP + P.k9_3f(8).*c_C18_FabA_EnAcACP - P.k9_3r(8).*c_FabA.*c_C18_EnAcACP;
d_C20_EnAcACP = P.k5_3f(9).*c_C20_FabZ_EnAcACP - P.k5_3r(9).*c_FabZ.*c_C20_EnAcACP + P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP - P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP + P.k9_3f(9).*c_C20_FabA_EnAcACP - P.k9_3r(9).*c_FabA.*c_C20_EnAcACP;

% C10 cis-3-Enoyl-Acyl-ACP (FabA - FabB)
d_C10_cis3EnAcACP = P.k9_3f_un(4).*c_C10_FabA_cis3EnAcACP - P.k9_3r_un(4).*c_FabA.*c_C10_cis3EnAcACP + P.k10_1r(4).*c_C10_FabB_cis3EnAcACP - P.k10_1f(4).*c_FabB.*c_C10_cis3EnAcACP;

% C2n:1 (n=6:10) Enoyl-Acyl-ACPs  (FabZ + FabA - FabI)
d_C12_EnAcACP_un = P.k5_3f(5).*c_C12_FabZ_EnAcACP_un - P.k5_3r(5).*c_FabZ.*c_C12_EnAcACP_un + P.k6_2r(5).*c_C12_FabI_NADH_EnAcACP_un - P.k6_2f(5).*c_FabI_NADH.*c_C12_EnAcACP_un + P.k9_3f(5).*c_C12_FabA_EnAcACP_un - P.k9_3r(5).*c_FabA.*c_C12_EnAcACP_un;
d_C14_EnAcACP_un = P.k5_3f(6).*c_C14_FabZ_EnAcACP_un - P.k5_3r(6).*c_FabZ.*c_C14_EnAcACP_un + P.k6_2r(6).*c_C14_FabI_NADH_EnAcACP_un - P.k6_2f(6).*c_FabI_NADH.*c_C14_EnAcACP_un + P.k9_3f(6).*c_C14_FabA_EnAcACP_un - P.k9_3r(6).*c_FabA.*c_C14_EnAcACP_un;
d_C16_EnAcACP_un = P.k5_3f(7).*c_C16_FabZ_EnAcACP_un - P.k5_3r(7).*c_FabZ.*c_C16_EnAcACP_un + P.k6_2r(7).*c_C16_FabI_NADH_EnAcACP_un - P.k6_2f(7).*c_FabI_NADH.*c_C16_EnAcACP_un + P.k9_3f(7).*c_C16_FabA_EnAcACP_un - P.k9_3r(7).*c_FabA.*c_C16_EnAcACP_un;
d_C18_EnAcACP_un = P.k5_3f(8).*c_C18_FabZ_EnAcACP_un - P.k5_3r(8).*c_FabZ.*c_C18_EnAcACP_un + P.k6_2r(8).*c_C18_FabI_NADH_EnAcACP_un - P.k6_2f(8).*c_FabI_NADH.*c_C18_EnAcACP_un + P.k9_3f(8).*c_C18_FabA_EnAcACP_un - P.k9_3r(8).*c_FabA.*c_C18_EnAcACP_un;
d_C20_EnAcACP_un = P.k5_3f(9).*c_C20_FabZ_EnAcACP_un - P.k5_3r(9).*c_FabZ.*c_C20_EnAcACP_un + P.k6_2r(9).*c_C20_FabI_NADH_EnAcACP_un - P.k6_2f(9).*c_FabI_NADH.*c_C20_EnAcACP_un + P.k9_3f(9).*c_C20_FabA_EnAcACP_un - P.k9_3r(9).*c_FabA.*c_C20_EnAcACP_un;

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C4_AcACP   = P.kcat6(1).*c_C4_FabI_NADH_EnAcACP   + P.k3_4r(1).*c_C4_FabH_AcACP  - P.k3_4f(1).*c_FabH.*c_C4_AcACP   + P.k3_5r(1).*c_C6_FabH_Act_AcACP   - P.k3_5f(1).*c_C2_FabH_Act.*c_C4_AcACP   + P.k7_1r(1).*c_C4_TesA_AcACP   + P.k8_1r(1).*c_C4_FabF_AcACP   + P.k10_1r(1).*c_C4_FabB_AcACP  - P.k7_1f(1).*c_TesA.*c_C4_AcACP   - P.k8_1f(1).*c_FabF.*c_C4_AcACP   - P.k10_1f(1).*c_FabB.*c_C4_AcACP;
d_C6_AcACP   = P.kcat6(2).*c_C6_FabI_NADH_EnAcACP   + P.k3_4r(2).*c_C6_FabH_AcACP  - P.k3_4f(2).*c_FabH.*c_C6_AcACP   + P.k3_5r(2).*c_C8_FabH_Act_AcACP   - P.k3_5f(2).*c_C2_FabH_Act.*c_C6_AcACP   + P.k7_1r(2).*c_C6_TesA_AcACP   + P.k8_1r(2).*c_C6_FabF_AcACP   + P.k10_1r(2).*c_C6_FabB_AcACP  - P.k7_1f(2).*c_TesA.*c_C6_AcACP   - P.k8_1f(2).*c_FabF.*c_C6_AcACP   - P.k10_1f(2).*c_FabB.*c_C6_AcACP;
d_C8_AcACP   = P.kcat6(3).*c_C8_FabI_NADH_EnAcACP   + P.k3_4r(3).*c_C8_FabH_AcACP  - P.k3_4f(3).*c_FabH.*c_C8_AcACP   + P.k3_5r(3).*c_C10_FabH_Act_AcACP   - P.k3_5f(3).*c_C2_FabH_Act.*c_C8_AcACP   + P.k7_1r(3).*c_C8_TesA_AcACP   + P.k8_1r(3).*c_C8_FabF_AcACP   + P.k10_1r(3).*c_C8_FabB_AcACP  - P.k7_1f(3).*c_TesA.*c_C8_AcACP   - P.k8_1f(3).*c_FabF.*c_C8_AcACP   - P.k10_1f(3).*c_FabB.*c_C8_AcACP;
d_C10_AcACP = P.kcat6(4).*c_C10_FabI_NADH_EnAcACP + P.k3_4r(4).*c_C10_FabH_AcACP - P.k3_4f(4).*c_FabH.*c_C10_AcACP + P.k3_5r(4).*c_C12_FabH_Act_AcACP - P.k3_5f(4).*c_C2_FabH_Act.*c_C10_AcACP + P.k7_1r(4).*c_C10_TesA_AcACP + P.k8_1r(4).*c_C10_FabF_AcACP + P.k10_1r(4).*c_C10_FabB_AcACP - P.k7_1f(4).*c_TesA.*c_C10_AcACP - P.k8_1f(4).*c_FabF.*c_C10_AcACP - P.k10_1f(4).*c_FabB.*c_C10_AcACP;
d_C12_AcACP = P.kcat6(5).*c_C12_FabI_NADH_EnAcACP + P.k3_4r(5).*c_C12_FabH_AcACP - P.k3_4f(5).*c_FabH.*c_C12_AcACP + P.k3_5r(5).*c_C14_FabH_Act_AcACP - P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP + P.k7_1r(5).*c_C12_TesA_AcACP + P.k8_1r(5).*c_C12_FabF_AcACP + P.k10_1r(5).*c_C12_FabB_AcACP - P.k7_1f(5).*c_TesA.*c_C12_AcACP - P.k8_1f(5).*c_FabF.*c_C12_AcACP - P.k10_1f(5).*c_FabB.*c_C12_AcACP;
d_C14_AcACP = P.kcat6(6).*c_C14_FabI_NADH_EnAcACP + P.k3_4r(6).*c_C14_FabH_AcACP - P.k3_4f(6).*c_FabH.*c_C14_AcACP + P.k3_5r(6).*c_C16_FabH_Act_AcACP - P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP + P.k7_1r(6).*c_C14_TesA_AcACP + P.k8_1r(6).*c_C14_FabF_AcACP + P.k10_1r(6).*c_C14_FabB_AcACP - P.k7_1f(6).*c_TesA.*c_C14_AcACP - P.k8_1f(6).*c_FabF.*c_C14_AcACP - P.k10_1f(6).*c_FabB.*c_C14_AcACP;
d_C16_AcACP = P.kcat6(7).*c_C16_FabI_NADH_EnAcACP + P.k3_4r(7).*c_C16_FabH_AcACP - P.k3_4f(7).*c_FabH.*c_C16_AcACP + P.k3_5r(7).*c_C18_FabH_Act_AcACP - P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP + P.k7_1r(7).*c_C16_TesA_AcACP + P.k8_1r(7).*c_C16_FabF_AcACP + P.k10_1r(7).*c_C16_FabB_AcACP - P.k7_1f(7).*c_TesA.*c_C16_AcACP - P.k8_1f(7).*c_FabF.*c_C16_AcACP - P.k10_1f(7).*c_FabB.*c_C16_AcACP;
d_C18_AcACP = P.kcat6(8).*c_C18_FabI_NADH_EnAcACP + P.k3_4r(8).*c_C18_FabH_AcACP - P.k3_4f(8).*c_FabH.*c_C18_AcACP + P.k3_5r(8).*c_C20_FabH_Act_AcACP - P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP + P.k7_1r(8).*c_C18_TesA_AcACP + P.k8_1r(8).*c_C18_FabF_AcACP + P.k10_1r(8).*c_C18_FabB_AcACP - P.k7_1f(8).*c_TesA.*c_C18_AcACP - P.k8_1f(8).*c_FabF.*c_C18_AcACP - P.k10_1f(8).*c_FabB.*c_C18_AcACP;

% C2n (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP = P.kcat6(9).*c_C20_FabI_NADH_EnAcACP + P.k3_4r(9).*c_C20_FabH_AcACP - P.k3_4f(9).*c_FabH.*c_C20_AcACP + P.k3_5r(9).*c_C22_FabH_Act_AcACP - P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP + P.k7_1r(9).*c_C20_TesA_AcACP - P.k7_1f(9).*c_TesA.*c_C20_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH) %no change
d_C12_AcACP_un = P.kcat6(5).*c_C12_FabI_NADH_EnAcACP_un + P.k3_4r(5).*c_C12_FabH_AcACP_un - P.k3_4f(5).*c_FabH.*c_C12_AcACP_un + P.k3_5r(5).*c_C14_FabH_Act_AcACP_un - P.k3_5f(5).*c_C2_FabH_Act.*c_C12_AcACP_un + P.k7_1r(5).*c_C12_TesA_AcACP_un - P.k7_1f(5).*c_TesA.*c_C12_AcACP_un + P.k8_1r(5).*c_C12_FabF_AcACP_un - P.k8_1f(5).*c_FabF.*c_C12_AcACP_un + P.k10_1r(5).*c_C12_FabB_AcACP_un - P.k10_1f(5).*c_FabB.*c_C12_AcACP_un;
d_C14_AcACP_un = P.kcat6(6).*c_C14_FabI_NADH_EnAcACP_un + P.k3_4r(6).*c_C14_FabH_AcACP_un - P.k3_4f(6).*c_FabH.*c_C14_AcACP_un + P.k3_5r(6).*c_C16_FabH_Act_AcACP_un - P.k3_5f(6).*c_C2_FabH_Act.*c_C14_AcACP_un + P.k7_1r(6).*c_C14_TesA_AcACP_un - P.k7_1f(6).*c_TesA.*c_C14_AcACP_un + P.k8_1r(6).*c_C14_FabF_AcACP_un - P.k8_1f(6).*c_FabF.*c_C14_AcACP_un + P.k10_1r(6).*c_C14_FabB_AcACP_un - P.k10_1f(6).*c_FabB.*c_C14_AcACP_un;
d_C16_AcACP_un = P.kcat6(7).*c_C16_FabI_NADH_EnAcACP_un + P.k3_4r(7).*c_C16_FabH_AcACP_un - P.k3_4f(7).*c_FabH.*c_C16_AcACP_un + P.k3_5r(7).*c_C18_FabH_Act_AcACP_un - P.k3_5f(7).*c_C2_FabH_Act.*c_C16_AcACP_un + P.k7_1r(7).*c_C16_TesA_AcACP_un - P.k7_1f(7).*c_TesA.*c_C16_AcACP_un + P.k8_1r(7).*c_C16_FabF_AcACP_un - P.k8_1f(7).*c_FabF.*c_C16_AcACP_un + P.k10_1r(7).*c_C16_FabB_AcACP_un - P.k10_1f(7).*c_FabB.*c_C16_AcACP_un;
d_C18_AcACP_un = P.kcat6(8).*c_C18_FabI_NADH_EnAcACP_un + P.k3_4r(8).*c_C18_FabH_AcACP_un - P.k3_4f(8).*c_FabH.*c_C18_AcACP_un + P.k3_5r(8).*c_C20_FabH_Act_AcACP_un - P.k3_5f(8).*c_C2_FabH_Act.*c_C18_AcACP_un + P.k7_1r(8).*c_C18_TesA_AcACP_un - P.k7_1f(8).*c_TesA.*c_C18_AcACP_un + P.k8_1r(8).*c_C18_FabF_AcACP_un - P.k8_1f(8).*c_FabF.*c_C18_AcACP_un + P.k10_1r(8).*c_C18_FabB_AcACP_un - P.k10_1f(8).*c_FabB.*c_C18_AcACP_un;

% C2n:1 (n=10) Acyl-ACPs (FabI - TesA - FabH) %no change
d_C20_AcACP_un = P.kcat6(9).*c_C20_FabI_NADH_EnAcACP_un + P.k7_1r(9).*c_C20_TesA_AcACP_un - P.k7_1f(9).*c_TesA.*c_C20_AcACP_un + P.k3_4r(9).*c_C20_FabH_AcACP_un - P.k3_4f(9).*c_FabH.*c_C20_AcACP_un + P.k3_5r(9).*c_C22_FabH_Act_AcACP_un - P.k3_5f(9).*c_C2_FabH_Act.*c_C20_AcACP_un;

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

% ACC Step 1
d_ACC_s1 = P.k1_1f.*c_ACC.*c_ATP - P.k1_1r.*c_ACC_s1 + P.k1_2r.*c_C1_ACC_s2 - P.k1_2f.*c_ACC_s1.*c_C1_Bicarbonate;
%d_ACC_s1 = 0*c_ACC_s1;

% ACC Step 2
d_C1_ACC_s2 = P.k1_2f.*c_ACC_s1.*c_C1_Bicarbonate - P.k1_2r.*c_C1_ACC_s2 - P.kcat1_1.*c_C1_ACC_s2;
%d_C1_ACC_s2 = 0*c_C1_ACC_s2;

% ACC Step 3
d_C1_ACC_s3 = P.kcat1_1.*c_C1_ACC_s2 + P.k1_3r.*c_C3_ACC_s4 - P.k1_3f.*c_C1_ACC_s3.*c_C2_AcCoA;
%d_C1_ACC_s3 = 0*c_C1_ACC_s3;

% ACC Step 4
d_C3_ACC_s4 = P.k1_3f.*c_C1_ACC_s3.*c_C2_AcCoA - P.k1_3r.*c_C3_ACC_s4 - P.kcat1_2.*c_C3_ACC_s4;
%d_C3_ACC_s4 = 0*c_C3_ACC_s4;

% FabD-Malonyl-CoA
d_C3_FabD_MalCoA = P.k2_1f.*c_FabD.*c_C3_MalCoA - P.k2_1r.*c_C3_FabD_MalCoA + P.k2_2r.*c_C3_FabD_Act.*c_CoA - P.k2_2f.*c_C3_FabD_MalCoA; 

% FabD*
d_C3_FabD_Act = P.k2_2f.*c_C3_FabD_MalCoA - P.k2_2r.*c_C3_FabD_Act.*c_CoA + P.k2_3r.*c_C3_FabD_Act_ACP - P.k2_3f.*c_C3_FabD_Act.*c_ACP;

% FabD*-ACP
d_C3_FabD_Act_ACP = P.k2_3f.*c_C3_FabD_Act.*c_ACP - P.k2_3r.*c_C3_FabD_Act_ACP + P.k2_4r.*c_FabD.*c_C3_MalACP - P.k2_4f.*c_C3_FabD_Act_ACP;

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
d_C8_FabF_Act   = P.k8_2f(3).*c_C8_FabF_AcACP   - P.k8_2r(3).*c_C8_FabF_Act.*c_ACP   + P.k8_3r(3).*c_C11_FabF_Act_MalACP  - P.k8_3f(3).*c_C8_FabF_Act.*c_C3_MalACP;
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
d_C11_FabF_Act_MalACP   = P.k8_3f(3).*c_C8_FabF_Act.*c_C3_MalACP   - P.k8_3r(3).*c_C11_FabF_Act_MalACP  - P.kcat8(3).*c_C11_FabF_Act_MalACP;
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

% C2n (n=2:9) FabB-Acyl-ACPs
d_C4_FabB_AcACP   = P.k10_1f(1).*c_FabB.*c_C4_AcACP   - P.k10_1r(1).*c_C4_FabB_AcACP  + P.k10_2r(1).*c_C4_FabB_Act.*c_ACP   - P.k10_2f(1).*c_C4_FabB_AcACP;
d_C6_FabB_AcACP   = P.k10_1f(2).*c_FabB.*c_C6_AcACP   - P.k10_1r(2).*c_C6_FabB_AcACP  + P.k10_2r(2).*c_C6_FabB_Act.*c_ACP   - P.k10_2f(2).*c_C6_FabB_AcACP;
d_C8_FabB_AcACP   = P.k10_1f(3).*c_FabB.*c_C8_AcACP   - P.k10_1r(3).*c_C8_FabB_AcACP  + P.k10_2r(3).*c_C8_FabB_Act.*c_ACP   - P.k10_2f(3).*c_C8_FabB_AcACP;
d_C10_FabB_AcACP = P.k10_1f(4).*c_FabB.*c_C10_AcACP - P.k10_1r(4).*c_C10_FabB_AcACP + P.k10_2r(4).*c_C10_FabB_Act.*c_ACP - P.k10_2f(4).*c_C10_FabB_AcACP;
d_C12_FabB_AcACP = P.k10_1f(5).*c_FabB.*c_C12_AcACP - P.k10_1r(5).*c_C12_FabB_AcACP + P.k10_2r(5).*c_C12_FabB_Act.*c_ACP - P.k10_2f(5).*c_C12_FabB_AcACP;
d_C14_FabB_AcACP = P.k10_1f(6).*c_FabB.*c_C14_AcACP - P.k10_1r(6).*c_C14_FabB_AcACP + P.k10_2r(6).*c_C14_FabB_Act.*c_ACP - P.k10_2f(6).*c_C14_FabB_AcACP;
d_C16_FabB_AcACP = P.k10_1f(7).*c_FabB.*c_C16_AcACP - P.k10_1r(7).*c_C16_FabB_AcACP + P.k10_2r(7).*c_C16_FabB_Act.*c_ACP - P.k10_2f(7).*c_C16_FabB_AcACP;
d_C18_FabB_AcACP = P.k10_1f(8).*c_FabB.*c_C18_AcACP - P.k10_1r(8).*c_C18_FabB_AcACP + P.k10_2r(8).*c_C18_FabB_Act.*c_ACP - P.k10_2f(8).*c_C18_FabB_AcACP;

% C2n:1 (n=6:9) FabB-Acyl-ACPs
d_C12_FabB_AcACP_un = P.k10_1f(5).*c_FabB.*c_C12_AcACP_un - P.k10_1r(5).*c_C12_FabB_AcACP_un + P.k10_2r(5).*c_C12_FabB_Act_un.*c_ACP - P.k10_2f(5).*c_C12_FabB_AcACP_un;
d_C14_FabB_AcACP_un = P.k10_1f(6).*c_FabB.*c_C14_AcACP_un - P.k10_1r(6).*c_C14_FabB_AcACP_un + P.k10_2r(6).*c_C14_FabB_Act_un.*c_ACP - P.k10_2f(6).*c_C14_FabB_AcACP_un;
d_C16_FabB_AcACP_un = P.k10_1f(7).*c_FabB.*c_C16_AcACP_un - P.k10_1r(7).*c_C16_FabB_AcACP_un + P.k10_2r(7).*c_C16_FabB_Act_un.*c_ACP - P.k10_2f(7).*c_C16_FabB_AcACP_un;
d_C18_FabB_AcACP_un = P.k10_1f(8).*c_FabB.*c_C18_AcACP_un - P.k10_1r(8).*c_C18_FabB_AcACP_un + P.k10_2r(8).*c_C18_FabB_Act_un.*c_ACP - P.k10_2f(8).*c_C18_FabB_AcACP_un;

% C2n (n=2:9) FabB*
d_C4_FabB_Act   = P.k10_2f(1).*c_C4_FabB_AcACP   - P.k10_2r(1).*c_C4_FabB_Act.*c_ACP  + P.k10_3r(1).*c_C7_FabB_Act_MalACP   - P.k10_3f(1).*c_C4_FabB_Act.*c_C3_MalACP;
d_C6_FabB_Act   = P.k10_2f(2).*c_C6_FabB_AcACP   - P.k10_2r(2).*c_C6_FabB_Act.*c_ACP  + P.k10_3r(2).*c_C9_FabB_Act_MalACP   - P.k10_3f(2).*c_C6_FabB_Act.*c_C3_MalACP;
d_C8_FabB_Act   = P.k10_2f(3).*c_C8_FabB_AcACP   - P.k10_2r(3).*c_C8_FabB_Act.*c_ACP  + P.k10_3r(3).*c_C11_FabB_Act_MalACP   - P.k10_3f(3).*c_C8_FabB_Act.*c_C3_MalACP;
d_C10_FabB_Act = P.k10_2f(4).*c_C10_FabB_AcACP - P.k10_2r(4).*c_C10_FabB_Act.*c_ACP + P.k10_3r(4).*c_C13_FabB_Act_MalACP - P.k10_3f(4).*c_C10_FabB_Act.*c_C3_MalACP;
d_C12_FabB_Act = P.k10_2f(5).*c_C12_FabB_AcACP - P.k10_2r(5).*c_C12_FabB_Act.*c_ACP + P.k10_3r(5).*c_C15_FabB_Act_MalACP - P.k10_3f(5).*c_C12_FabB_Act.*c_C3_MalACP;
d_C14_FabB_Act = P.k10_2f(6).*c_C14_FabB_AcACP - P.k10_2r(6).*c_C14_FabB_Act.*c_ACP + P.k10_3r(6).*c_C17_FabB_Act_MalACP - P.k10_3f(6).*c_C14_FabB_Act.*c_C3_MalACP;
d_C16_FabB_Act = P.k10_2f(7).*c_C16_FabB_AcACP - P.k10_2r(7).*c_C16_FabB_Act.*c_ACP + P.k10_3r(7).*c_C19_FabB_Act_MalACP - P.k10_3f(7).*c_C16_FabB_Act.*c_C3_MalACP;
d_C18_FabB_Act = P.k10_2f(8).*c_C18_FabB_AcACP - P.k10_2r(8).*c_C18_FabB_Act.*c_ACP + P.k10_3r(8).*c_C21_FabB_Act_MalACP - P.k10_3f(8).*c_C18_FabB_Act.*c_C3_MalACP;

% C2n:1 (n=6:9) FabB*
d_C12_FabB_Act_un = P.k10_2f(5).*c_C12_FabB_AcACP_un - P.k10_2r(5).*c_C12_FabB_Act_un.*c_ACP + P.k10_3r(5).*c_C15_FabB_Act_MalACP_un - P.k10_3f(5).*c_C12_FabB_Act_un.*c_C3_MalACP;
d_C14_FabB_Act_un = P.k10_2f(6).*c_C14_FabB_AcACP_un - P.k10_2r(6).*c_C14_FabB_Act_un.*c_ACP + P.k10_3r(6).*c_C17_FabB_Act_MalACP_un - P.k10_3f(6).*c_C14_FabB_Act_un.*c_C3_MalACP;
d_C16_FabB_Act_un = P.k10_2f(7).*c_C16_FabB_AcACP_un - P.k10_2r(7).*c_C16_FabB_Act_un.*c_ACP + P.k10_3r(7).*c_C19_FabB_Act_MalACP_un - P.k10_3f(7).*c_C16_FabB_Act_un.*c_C3_MalACP;
d_C18_FabB_Act_un = P.k10_2f(8).*c_C18_FabB_AcACP_un - P.k10_2r(8).*c_C18_FabB_Act_un.*c_ACP + P.k10_3r(8).*c_C21_FabB_Act_MalACP_un - P.k10_3f(8).*c_C18_FabB_Act_un.*c_C3_MalACP;

% C2n (n=2:9) FabB*-Malonyl-ACPs
d_C7_FabB_Act_MalACP   = P.k10_3f(1).*c_C4_FabB_Act.*c_C3_MalACP   - P.k10_3r(1).*c_C7_FabB_Act_MalACP  - P.kcat10(1).*c_C7_FabB_Act_MalACP;
d_C9_FabB_Act_MalACP   = P.k10_3f(2).*c_C6_FabB_Act.*c_C3_MalACP   - P.k10_3r(2).*c_C9_FabB_Act_MalACP  - P.kcat10(2).*c_C9_FabB_Act_MalACP;
d_C11_FabB_Act_MalACP   = P.k10_3f(3).*c_C8_FabB_Act.*c_C3_MalACP   - P.k10_3r(3).*c_C11_FabB_Act_MalACP  - P.kcat10(3).*c_C11_FabB_Act_MalACP;
d_C13_FabB_Act_MalACP = P.k10_3f(4).*c_C10_FabB_Act.*c_C3_MalACP - P.k10_3r(4).*c_C13_FabB_Act_MalACP - P.kcat10(4).*c_C13_FabB_Act_MalACP;
d_C15_FabB_Act_MalACP = P.k10_3f(5).*c_C12_FabB_Act.*c_C3_MalACP - P.k10_3r(5).*c_C15_FabB_Act_MalACP - P.kcat10(5).*c_C15_FabB_Act_MalACP;
d_C17_FabB_Act_MalACP = P.k10_3f(6).*c_C14_FabB_Act.*c_C3_MalACP - P.k10_3r(6).*c_C17_FabB_Act_MalACP - P.kcat10(6).*c_C17_FabB_Act_MalACP;
d_C19_FabB_Act_MalACP = P.k10_3f(7).*c_C16_FabB_Act.*c_C3_MalACP - P.k10_3r(7).*c_C19_FabB_Act_MalACP - P.kcat10(7).*c_C19_FabB_Act_MalACP;
d_C21_FabB_Act_MalACP = P.k10_3f(8).*c_C18_FabB_Act.*c_C3_MalACP - P.k10_3r(8).*c_C21_FabB_Act_MalACP - P.kcat10(8).*c_C21_FabB_Act_MalACP;

% C2n:1 (n=6:9) FabB*-Malonyl-ACPs
d_C15_FabB_Act_MalACP_un = P.k10_3f(5).*c_C12_FabB_Act_un.*c_C3_MalACP - P.k10_3r(5).*c_C15_FabB_Act_MalACP_un - P.kcat10_un(5).*c_C15_FabB_Act_MalACP_un;
d_C17_FabB_Act_MalACP_un = P.k10_3f(6).*c_C14_FabB_Act_un.*c_C3_MalACP - P.k10_3r(6).*c_C17_FabB_Act_MalACP_un - P.kcat10_un(6).*c_C17_FabB_Act_MalACP_un;
d_C19_FabB_Act_MalACP_un = P.k10_3f(7).*c_C16_FabB_Act_un.*c_C3_MalACP - P.k10_3r(7).*c_C19_FabB_Act_MalACP_un - P.kcat10_un(7).*c_C19_FabB_Act_MalACP_un;
d_C21_FabB_Act_MalACP_un = P.k10_3f(8).*c_C18_FabB_Act_un.*c_C3_MalACP - P.k10_3r(8).*c_C21_FabB_Act_MalACP_un - P.kcat10_un(8).*c_C21_FabB_Act_MalACP_un;

% FabB-(C10 cis-3-Enoyl-Acyl-ACP)
d_C10_FabB_cis3EnAcACP = P.k10_1f(4).*c_FabB.*c_C10_cis3EnAcACP - P.k10_1r(4).*c_C10_FabB_cis3EnAcACP + P.k10_2r(4).*c_C10_FabB_Act_cis3.*c_ACP - P.k10_2f(4).*c_C10_FabB_cis3EnAcACP;

% C10 cis-3-FabB*
d_C10_FabB_Act_cis3 = P.k10_2f(4).*c_C10_FabB_cis3EnAcACP - P.k10_2r(4).*c_C10_FabB_Act_cis3.*c_ACP + P.k10_3r(4).*c_C13_FabB_Act_cis3MalACP - P.k10_3f(4).*c_C10_FabB_Act_cis3.*c_C3_MalACP;

% C10 cis-3-FabB*-Malonyl-ACP
d_C13_FabB_Act_cis3MalACP = P.k10_3f(4).*c_C10_FabB_Act_cis3.*c_C3_MalACP - P.k10_3r(4).*c_C13_FabB_Act_cis3MalACP - P.kcat10_un(4).*c_C13_FabB_Act_cis3MalACP;

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

% TesA-ACP
d_TesA_ACP = P.k7_inh_f.*c_TesA.*c_ACP - P.k7_inh_r.*c_TesA_ACP;

% FabH-ACP
d_FabH_ACP = P.k3_inh_f.*c_FabH.*c_ACP - P.k3_inh_r.*c_FabH_ACP;

% FabG-ACP
d_FabG_ACP = P.k4_inh_f.*c_FabG.*c_ACP - P.k4_inh_r.*c_FabG_ACP;

% FabZ-ACP
d_FabZ_ACP = P.k5_inh_f.*c_FabZ.*c_ACP - P.k5_inh_r.*c_FabZ_ACP;

% FabI-ACP
d_FabI_ACP = P.k6_inh_f.*c_FabI.*c_ACP - P.k6_inh_r.*c_FabI_ACP;

% FabF-ACP
d_FabF_ACP = P.k8_inh_f.*c_FabF.*c_ACP - P.k8_inh_r.*c_FabF_ACP;

% FabA-ACP
d_FabA_ACP = P.k9_inh_f.*c_FabA.*c_ACP - P.k9_inh_r.*c_FabA_ACP;

% FabB-ACP
d_FabB_ACP = P.k10_inh_f.*c_FabB.*c_ACP - P.k10_inh_r.*c_FabB_ACP;


% Giving FabB FabH-like activity
% FabB-Acetyl-CoA
d_C2_FabB_AcCoA = P.k10_4f.*c_FabB.*c_C2_AcCoA - P.k10_4r.*c_C2_FabB_AcCoA + P.k10_5r.*c_C2_FabB_Act.*c_CoA - P.k10_5f.*c_C2_FabB_AcCoA;
%d_C2_FabB_AcCoA = 0*c_C2_FabB_AcCoA;

% FabB*
d_C2_FabB_Act = P.k10_5f.*c_C2_FabB_AcCoA - P.k10_5r.*c_C2_FabB_Act.*c_CoA + P.k10_6r.*c_C5_FabB_Act_MalACP - P.k10_6f.*c_C2_FabB_Act.*c_C3_MalACP + P.k10_9f.*c_C2_FabB_AcACP - P.k10_9r.*c_C2_FabB_Act.*c_ACP;
%d_C2_FabB_Act = 0*c_C2_FabB_Act;

% FabB*-Malonyl-ACP
d_C5_FabB_Act_MalACP = P.k10_6f.*c_C2_FabB_Act.*c_C3_MalACP - P.k10_6r.*c_C5_FabB_Act_MalACP - P.kcat10_H.*c_C5_FabB_Act_MalACP; 
%d_C5_FabB_Act_MalACP = 0*c_C5_FabB_Act_MalACP;

% Giving FabF FabH-like activity
% FabF-Acetyl-CoA
d_C2_FabF_AcCoA = P.k8_4f.*c_FabF.*c_C2_AcCoA - P.k8_4r.*c_C2_FabF_AcCoA + P.k8_5r.*c_C2_FabF_Act.*c_CoA - P.k8_5f.*c_C2_FabF_AcCoA; 
%d_C2_FabF_AcCoA = 0*c_C2_FabF_AcCoA;

% FabF*
d_C2_FabF_Act = P.k8_5f.*c_C2_FabF_AcCoA - P.k8_5r.*c_C2_FabF_Act.*c_CoA + P.k8_6r.*c_C5_FabF_Act_MalACP - P.k8_6f.*c_C2_FabF_Act.*c_C3_MalACP + P.k8_9f.*c_C2_FabF_AcACP - P.k8_9r.*c_C2_FabF_Act.*c_ACP; 
%d_C2_FabF_Act = 0*c_C2_FabF_Act;

% FabF*-Malonyl-ACP
d_C5_FabF_Act_MalACP = P.k8_6f.*c_C2_FabF_Act.*c_C3_MalACP - P.k8_6r.*c_C5_FabF_Act_MalACP - P.kcat8_H.*c_C5_FabF_Act_MalACP; 
%d_C5_FabF_Act_MalACP = 0*c_C5_FabF_Act_MalACP;

% FabB and FabF decarboxylating mACP to form aACP and reacting with it to form activated enzyme (initiation)
% FabB-Malonyl-ACP (e10p4)
d_C3_FabB_MalACP = P.k10_7f.*c_FabB.*c_C3_MalACP - P.k10_7r.*c_C3_FabB_MalACP - P.kcat10_CO2.*c_C3_FabB_MalACP;
%d_C3_FabB_MalACP = 0*c_C3_FabB_MalACP;

% Acetyl-ACP (T_2)
d_C2_AcACP = P.k8_8r.*c_C2_FabF_AcACP - P.k8_8f.*c_FabF.*c_C2_AcACP + P.k10_8r.*c_C2_FabB_AcACP - P.k10_8f.*c_FabB.*c_C2_AcACP; 
%d_C2_AcACP = 0*c_C2_AcACP;

% FabB-Acetyl-ACP (e10T_2)
d_C2_FabB_AcACP = P.kcat10_CO2.*c_C3_FabB_MalACP + P.k10_8f.*c_FabB.*c_C2_AcACP - P.k10_8r.*c_C2_FabB_AcACP + P.k10_9r.*c_C2_FabB_Act.*c_ACP - P.k10_9f.*c_C2_FabB_AcACP; 
%d_C2_FabB_AcACP = 0*c_C2_FabB_AcACP;

% FabF-Malonyl-ACP (e8p4)
d_C3_FabF_MalACP = P.k8_7f.*c_FabF.*c_C3_MalACP - P.k8_7r.*c_C3_FabF_MalACP - P.kcat8_CO2.*c_C3_FabF_MalACP; 
%d_C3_FabF_MalACP = 0*c_C3_FabF_MalACP;

% FabF-Acetyl-ACP (e8T_2)
d_C2_FabF_AcACP = P.kcat8_CO2.*c_C3_FabF_MalACP + P.k8_8f.*c_FabF.*c_C2_AcACP - P.k8_8r.*c_C2_FabF_AcACP + P.k8_9r.*c_C2_FabF_Act.*c_ACP - P.k8_9f.*c_C2_FabF_AcACP; 
%d_C2_FabF_AcACP = 0*c_C2_FabF_AcACP;


dcdt = [d_ATP; d_C1_Bicarbonate; d_C2_AcCoA; d_C4_SucCoA; d_C6_HexCoA; d_C8_OcCoA; d_C10_DecCoA; d_C12_LauCoA; d_C14_EthCoA; d_C16_PalCoA; d_C18_OcDecCoA; d_ACP; d_NADPH; d_NADH; d_ADP; ...
    d_C3_MalCoA; d_CoA; d_C3_MalACP; d_C1_CO2; d_C4_BKeACP; d_C6_BKeACP; d_C8_BKeACP; d_C10_BKeACP; d_C12_BKeACP; d_C14_BKeACP; d_C16_BKeACP; d_C18_BKeACP; ...
    d_C20_BKeACP; d_C12_BKeACP_un; d_C14_BKeACP_un; d_C16_BKeACP_un; d_C18_BKeACP_un; d_C20_BKeACP_un; d_C4_BHyAcACP; d_C6_BHyAcACP; d_C8_BHyAcACP; ...
    d_C10_BHyAcACP; d_C12_BHyAcACP; d_C14_BHyAcACP; d_C16_BHyAcACP; d_C18_BHyAcACP; d_C20_BHyAcACP; d_C12_BHyAcACP_un; d_C14_BHyAcACP_un; ...
    d_C16_BHyAcACP_un; d_C18_BHyAcACP_un; d_C20_BHyAcACP_un; d_C4_EnAcACP; d_C6_EnAcACP; d_C8_EnAcACP; d_C10_EnAcACP; d_C12_EnAcACP; d_C14_EnAcACP; ...
    d_C16_EnAcACP; d_C18_EnAcACP; d_C20_EnAcACP; d_C10_cis3EnAcACP; d_C12_EnAcACP_un; d_C14_EnAcACP_un; d_C16_EnAcACP_un; d_C18_EnAcACP_un; ...
    d_C20_EnAcACP_un; d_C4_AcACP; d_C6_AcACP; d_C8_AcACP; d_C10_AcACP; d_C12_AcACP; d_C14_AcACP; d_C16_AcACP; d_C18_AcACP; d_C20_AcACP; ...
    d_C12_AcACP_un; d_C14_AcACP_un; d_C16_AcACP_un; d_C18_AcACP_un; d_C20_AcACP_un; d_C4_FA; d_C6_FA; d_C8_FA; d_C10_FA; d_C12_FA; d_C14_FA; ...
    d_C16_FA; d_C18_FA; d_C20_FA; d_C12_FA_un; d_C14_FA_un; d_C16_FA_un; d_C18_FA_un; d_C20_FA_un; d_ACC_s1; d_C1_ACC_s2; d_C1_ACC_s3; d_C3_ACC_s4; d_C3_FabD_MalCoA; ...
    d_C3_FabD_Act; d_C3_FabD_Act_ACP; d_C2_FabH_CoA; d_C4_FabH_CoA; d_C6_FabH_CoA; d_C8_FabH_CoA; d_C10_FabH_CoA; d_C12_FabH_CoA; d_C14_FabH_CoA; ...
    d_C16_FabH_CoA; d_C18_FabH_CoA; d_C2_FabH_Act; d_C4_FabH_Act; d_C6_FabH_Act; d_C8_FabH_Act; d_C10_FabH_Act; d_C12_FabH_Act; d_C14_FabH_Act; ...
    d_C16_FabH_Act; d_C18_FabH_Act; d_C5_FabH_Act_MalACP; d_C7_FabH_Act_MalACP; d_C9_FabH_Act_MalACP; d_C11_FabH_Act_MalACP; d_C13_FabH_Act_MalACP; ...
    d_C15_FabH_Act_MalACP; d_C17_FabH_Act_MalACP; d_C19_FabH_Act_MalACP; d_C21_FabH_Act_MalACP; d_FabG_NADPH; d_C4_FabG_NADPH_BKeACP; ...
    d_C6_FabG_NADPH_BKeACP; d_C8_FabG_NADPH_BKeACP; d_C10_FabG_NADPH_BKeACP; d_C12_FabG_NADPH_BKeACP; d_C14_FabG_NADPH_BKeACP; ...
    d_C16_FabG_NADPH_BKeACP; d_C18_FabG_NADPH_BKeACP; d_C20_FabG_NADPH_BKeACP; d_C12_FabG_NADPH_BKeACP_un; d_C14_FabG_NADPH_BKeACP_un; ...
    d_C16_FabG_NADPH_BKeACP_un; d_C18_FabG_NADPH_BKeACP_un; d_C20_FabG_NADPH_BKeACP_un; d_C4_FabZ_BHyAcACP; d_C6_FabZ_BHyAcACP; ...
    d_C8_FabZ_BHyAcACP; d_C10_FabZ_BHyAcACP; d_C12_FabZ_BHyAcACP; d_C14_FabZ_BHyAcACP; d_C16_FabZ_BHyAcACP; d_C18_FabZ_BHyAcACP; ...
    d_C20_FabZ_BHyAcACP; d_C12_FabZ_BHyAcACP_un; d_C14_FabZ_BHyAcACP_un; d_C16_FabZ_BHyAcACP_un; d_C18_FabZ_BHyAcACP_un; d_C20_FabZ_BHyAcACP_un; ...
    d_C4_FabZ_EnAcACP; d_C6_FabZ_EnAcACP; d_C8_FabZ_EnAcACP; d_C10_FabZ_EnAcACP; d_C12_FabZ_EnAcACP; d_C14_FabZ_EnAcACP; d_C16_FabZ_EnAcACP; ...
    d_C18_FabZ_EnAcACP; d_C20_FabZ_EnAcACP; d_C12_FabZ_EnAcACP_un; d_C14_FabZ_EnAcACP_un; d_C16_FabZ_EnAcACP_un; d_C18_FabZ_EnAcACP_un; ...
    d_C20_FabZ_EnAcACP_un; d_C4_FabA_BHyAcACP; d_C6_FabA_BHyAcACP; d_C8_FabA_BHyAcACP; d_C10_FabA_BHyAcACP; d_C12_FabA_BHyAcACP; ...
    d_C14_FabA_BHyAcACP; d_C16_FabA_BHyAcACP; d_C18_FabA_BHyAcACP; d_C20_FabA_BHyAcACP; d_C12_FabA_BHyAcACP_un; d_C14_FabA_BHyAcACP_un; ...
    d_C16_FabA_BHyAcACP_un; d_C18_FabA_BHyAcACP_un; d_C20_FabA_BHyAcACP_un; d_C4_FabA_EnAcACP; d_C6_FabA_EnAcACP; d_C8_FabA_EnAcACP; d_C10_FabA_EnAcACP; ...
    d_C12_FabA_EnAcACP; d_C14_FabA_EnAcACP; d_C16_FabA_EnAcACP; d_C18_FabA_EnAcACP; d_C20_FabA_EnAcACP; d_C10_FabA_cis3EnAcACP; d_C12_FabA_EnAcACP_un; ...
    d_C14_FabA_EnAcACP_un; d_C16_FabA_EnAcACP_un; d_C18_FabA_EnAcACP_un; d_C20_FabA_EnAcACP_un; d_FabI_NADH; d_C4_FabI_NADH_EnAcACP; ...
    d_C6_FabI_NADH_EnAcACP; d_C8_FabI_NADH_EnAcACP; d_C10_FabI_NADH_EnAcACP; d_C12_FabI_NADH_EnAcACP; d_C14_FabI_NADH_EnAcACP; ...
    d_C16_FabI_NADH_EnAcACP; d_C18_FabI_NADH_EnAcACP; d_C20_FabI_NADH_EnAcACP; d_C12_FabI_NADH_EnAcACP_un; d_C14_FabI_NADH_EnAcACP_un; ...
    d_C16_FabI_NADH_EnAcACP_un; d_C18_FabI_NADH_EnAcACP_un; d_C20_FabI_NADH_EnAcACP_un; d_C4_TesA_AcACP; d_C6_TesA_AcACP; d_C8_TesA_AcACP; ...
    d_C10_TesA_AcACP; d_C12_TesA_AcACP; d_C14_TesA_AcACP; d_C16_TesA_AcACP; d_C18_TesA_AcACP; d_C20_TesA_AcACP; d_C12_TesA_AcACP_un; d_C14_TesA_AcACP_un; ...
    d_C16_TesA_AcACP_un; d_C18_TesA_AcACP_un; d_C20_TesA_AcACP_un; d_C4_FabF_AcACP; d_C6_FabF_AcACP; d_C8_FabF_AcACP; d_C10_FabF_AcACP; ...
    d_C12_FabF_AcACP; d_C14_FabF_AcACP; d_C16_FabF_AcACP; d_C18_FabF_AcACP; d_C12_FabF_AcACP_un; d_C14_FabF_AcACP_un; d_C16_FabF_AcACP_un; ...
    d_C18_FabF_AcACP_un; d_C4_FabF_Act; d_C6_FabF_Act; d_C8_FabF_Act; d_C10_FabF_Act; d_C12_FabF_Act; d_C14_FabF_Act; d_C16_FabF_Act; d_C18_FabF_Act; ...
    d_C12_FabF_Act_un; d_C14_FabF_Act_un; d_C16_FabF_Act_un; d_C18_FabF_Act_un; d_C7_FabF_Act_MalACP; d_C9_FabF_Act_MalACP; d_C11_FabF_Act_MalACP; ...
    d_C13_FabF_Act_MalACP; d_C15_FabF_Act_MalACP; d_C17_FabF_Act_MalACP; d_C19_FabF_Act_MalACP; d_C21_FabF_Act_MalACP; d_C15_FabF_Act_MalACP_un; ...
    d_C17_FabF_Act_MalACP_un; d_C19_FabF_Act_MalACP_un; d_C21_FabF_Act_MalACP_un; d_C4_FabB_AcACP; d_C6_FabB_AcACP; d_C8_FabB_AcACP; d_C10_FabB_AcACP; ...
    d_C12_FabB_AcACP; d_C14_FabB_AcACP; d_C16_FabB_AcACP; d_C18_FabB_AcACP; d_C12_FabB_AcACP_un; d_C14_FabB_AcACP_un; d_C16_FabB_AcACP_un; ...
    d_C18_FabB_AcACP_un; d_C4_FabB_Act; d_C6_FabB_Act; d_C8_FabB_Act; d_C10_FabB_Act; d_C12_FabB_Act; d_C14_FabB_Act; d_C16_FabB_Act; d_C18_FabB_Act; ...
    d_C12_FabB_Act_un; d_C14_FabB_Act_un; d_C16_FabB_Act_un; d_C18_FabB_Act_un; d_C7_FabB_Act_MalACP; d_C9_FabB_Act_MalACP; d_C11_FabB_Act_MalACP; ...
    d_C13_FabB_Act_MalACP; d_C15_FabB_Act_MalACP; d_C17_FabB_Act_MalACP; d_C19_FabB_Act_MalACP; d_C21_FabB_Act_MalACP; d_C15_FabB_Act_MalACP_un; ...
    d_C17_FabB_Act_MalACP_un; d_C19_FabB_Act_MalACP_un; d_C21_FabB_Act_MalACP_un; d_C10_FabB_cis3EnAcACP; d_C10_FabB_Act_cis3; d_C13_FabB_Act_cis3MalACP; ...
    d_C4_FabH_AcACP; d_C6_FabH_AcACP; d_C8_FabH_AcACP; d_C10_FabH_AcACP; d_C12_FabH_AcACP; d_C14_FabH_AcACP; d_C16_FabH_AcACP; d_C18_FabH_AcACP; ...
    d_C20_FabH_AcACP; d_C12_FabH_AcACP_un; d_C14_FabH_AcACP_un; d_C16_FabH_AcACP_un; d_C18_FabH_AcACP_un; d_C20_FabH_AcACP_un; d_C6_FabH_Act_AcACP; ...
    d_C8_FabH_Act_AcACP; d_C10_FabH_Act_AcACP; d_C12_FabH_Act_AcACP; d_C14_FabH_Act_AcACP; d_C16_FabH_Act_AcACP; d_C18_FabH_Act_AcACP; d_C20_FabH_Act_AcACP; ...
    d_C22_FabH_Act_AcACP; d_C14_FabH_Act_AcACP_un; d_C16_FabH_Act_AcACP_un; d_C18_FabH_Act_AcACP_un; d_C20_FabH_Act_AcACP_un; d_C22_FabH_Act_AcACP_un; ...
    d_TesA_ACP; d_FabH_ACP; d_FabG_ACP; d_FabZ_ACP; d_FabI_ACP; d_FabF_ACP; d_FabA_ACP; d_FabB_ACP; d_C2_FabB_AcCoA; d_C2_FabB_Act; d_C5_FabB_Act_MalACP; d_C2_FabF_AcCoA; ...
    d_C2_FabF_Act; d_C5_FabF_Act_MalACP; d_C3_FabB_MalACP; d_C2_AcACP; d_C2_FabB_AcACP; d_C3_FabF_MalACP; d_C2_FabF_AcACP];


