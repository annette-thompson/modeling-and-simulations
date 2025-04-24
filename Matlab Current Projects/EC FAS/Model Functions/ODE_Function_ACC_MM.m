%% ODE_Function_ACC_MM
function dcdt = ODE_Function_ACC_MM(t,conc,P)
% Contains all the differential equations and enzyme balances that define
% the FAS model with Michaelis Menton rate equations for ACC
% Input:
% t: time values (required as input for MATLAB ODE solver, sec)
% c: centration values (all components and intermediates, uM)
% P: structure containing all kinetic parameters
% Output:
% dcdt: values of differential equations for given 
% centrations and kinetic parameters

%% ODEs

% Preallocate structure
c = struct();

% Assign values dynamically
for i = 1:numel(P.labels)
    c.(P.labels{i}) = conc(i,:);
end

% ACC % changed for ACC
% BC (AccC)
c.ACC_C = P.ACC_Ctot - c.BC_ATP - c.C1_BC_ATP_HCO3 - c.C1_BC_Pi_HCO3 - c.C1_BC_Pi_HCO3_BCCP_Biotin;
% BCCP-Biotin (AccB-Biotin)
c.ACC_B = P.ACC_Btot - c.C1_BC_Pi_HCO3_BCCP_Biotin - c.C1_BCCP_Biotin_CO2 - c.C1_CT_BCCP_Biotin_CO2;
% CT (AccAD)
c.ACC_AD = P.ACC_ADtot - c.C1_CT_BCCP_Biotin_CO2 - c.C1_CT_Act - c.C3_CT_Act_AcCoA;

% FabD
c.FabD = P.FabDtot - c.C3_FabD_MalCoA - c.C3_FabD_Act - c.C3_FabD_Act_ACP;

% FabH
c.FabH = P.FabHtot - c.FabH_ACP...
 - c.C2_FabH_CoA - c.C4_FabH_CoA - c.C6_FabH_CoA - c.C8_FabH_CoA - c.C10_FabH_CoA - c.C12_FabH_CoA - c.C14_FabH_CoA - c.C16_FabH_CoA - c.C18_FabH_CoA...
 - c.C2_FabH_Act - c.C4_FabH_Act - c.C6_FabH_Act - c.C8_FabH_Act - c.C10_FabH_Act - c.C12_FabH_Act - c.C14_FabH_Act - c.C16_FabH_Act - c.C18_FabH_Act...
 - c.C5_FabH_Act_MalACP - c.C7_FabH_Act_MalACP - c.C9_FabH_Act_MalACP - c.C11_FabH_Act_MalACP - c.C13_FabH_Act_MalACP - c.C15_FabH_Act_MalACP - c.C17_FabH_Act_MalACP - c.C19_FabH_Act_MalACP - c.C21_FabH_Act_MalACP...
 - c.C4_FabH_AcACP - c.C6_FabH_AcACP - c.C8_FabH_AcACP - c.C10_FabH_AcACP - c.C12_FabH_AcACP - c.C14_FabH_AcACP - c.C16_FabH_AcACP - c.C18_FabH_AcACP - c.C20_FabH_AcACP...
 - c.C12_FabH_AcACP_un - c.C14_FabH_AcACP_un - c.C16_FabH_AcACP_un - c.C18_FabH_AcACP_un - c.C20_FabH_AcACP_un...
 - c.C6_FabH_Act_AcACP - c.C8_FabH_Act_AcACP - c.C10_FabH_Act_AcACP - c.C12_FabH_Act_AcACP - c.C14_FabH_Act_AcACP - c.C16_FabH_Act_AcACP - c.C18_FabH_Act_AcACP - c.C20_FabH_Act_AcACP - c.C22_FabH_Act_AcACP...
 - c.C14_FabH_Act_AcACP_un - c.C16_FabH_Act_AcACP_un - c.C18_FabH_Act_AcACP_un - c.C20_FabH_Act_AcACP_un - c.C22_FabH_Act_AcACP_un...
 - c.C3_FabH_MalACP - c.C2_FabH_AcACP;

% FabG
c.FabG = P.FabGtot - c.FabG_NADPH...
 - c.C4_FabG_NADPH_BKeACP - c.C6_FabG_NADPH_BKeACP - c.C8_FabG_NADPH_BKeACP - c.C10_FabG_NADPH_BKeACP - c.C12_FabG_NADPH_BKeACP - c.C14_FabG_NADPH_BKeACP - c.C16_FabG_NADPH_BKeACP - c.C18_FabG_NADPH_BKeACP - c.C20_FabG_NADPH_BKeACP...
 - c.C12_FabG_NADPH_BKeACP_un - c.C14_FabG_NADPH_BKeACP_un - c.C16_FabG_NADPH_BKeACP_un - c.C18_FabG_NADPH_BKeACP_un - c.C20_FabG_NADPH_BKeACP_un - c.FabG_ACP;

% FabZ
c.FabZ = P.FabZtot - c.FabZ_ACP...
 - c.C4_FabZ_BHyAcACP - c.C6_FabZ_BHyAcACP - c.C8_FabZ_BHyAcACP - c.C10_FabZ_BHyAcACP - c.C12_FabZ_BHyAcACP - c.C14_FabZ_BHyAcACP - c.C16_FabZ_BHyAcACP - c.C18_FabZ_BHyAcACP - c.C20_FabZ_BHyAcACP...
 - c.C12_FabZ_BHyAcACP_un - c.C14_FabZ_BHyAcACP_un - c.C16_FabZ_BHyAcACP_un - c.C18_FabZ_BHyAcACP_un - c.C20_FabZ_BHyAcACP_un...
 - c.C4_FabZ_EnAcACP - c.C6_FabZ_EnAcACP - c.C8_FabZ_EnAcACP - c.C10_FabZ_EnAcACP - c.C12_FabZ_EnAcACP - c.C14_FabZ_EnAcACP - c.C16_FabZ_EnAcACP - c.C18_FabZ_EnAcACP - c.C20_FabZ_EnAcACP...
 - c.C12_FabZ_EnAcACP_un - c.C14_FabZ_EnAcACP_un - c.C16_FabZ_EnAcACP_un - c.C18_FabZ_EnAcACP_un - c.C20_FabZ_EnAcACP_un;

% FabI
c.FabI = P.FabItot - c.FabI_NADH...
 - c.C4_FabI_NADH_EnAcACP - c.C6_FabI_NADH_EnAcACP - c.C8_FabI_NADH_EnAcACP - c.C10_FabI_NADH_EnAcACP - c.C12_FabI_NADH_EnAcACP - c.C14_FabI_NADH_EnAcACP - c.C16_FabI_NADH_EnAcACP - c.C18_FabI_NADH_EnAcACP - c.C20_FabI_NADH_EnAcACP...
 - c.C12_FabI_NADH_EnAcACP_un - c.C14_FabI_NADH_EnAcACP_un - c.C16_FabI_NADH_EnAcACP_un - c.C18_FabI_NADH_EnAcACP_un - c.C20_FabI_NADH_EnAcACP_un - c.FabI_ACP;

% TesA
c.TesA = P.TesAtot- c.TesA_ACP...
 - c.C4_TesA_AcACP - c.C6_TesA_AcACP - c.C8_TesA_AcACP - c.C10_TesA_AcACP - c.C12_TesA_AcACP - c.C14_TesA_AcACP - c.C16_TesA_AcACP - c.C18_TesA_AcACP - c.C20_TesA_AcACP...
 - c.C12_TesA_AcACP_un - c.C14_TesA_AcACP_un - c.C16_TesA_AcACP_un - c.C18_TesA_AcACP_un - c.C20_TesA_AcACP_un;

% FabF
c.FabF = P.FabFtot - c.FabF_ACP - c.C2_FabF_AcCoA - c.C2_FabF_Act - c.C3_FabF_MalACP - c.C2_FabF_AcACP...
 - c.C4_FabF_AcACP - c.C6_FabF_AcACP - c.C8_FabF_AcACP - c.C10_FabF_AcACP - c.C12_FabF_AcACP - c.C14_FabF_AcACP - c.C16_FabF_AcACP - c.C18_FabF_AcACP...
 - c.C12_FabF_AcACP_un - c.C14_FabF_AcACP_un - c.C16_FabF_AcACP_un - c.C18_FabF_AcACP_un...
 - c.C4_FabF_Act - c.C6_FabF_Act - c.C8_FabF_Act - c.C10_FabF_Act - c.C12_FabF_Act - c.C14_FabF_Act - c.C16_FabF_Act - c.C18_FabF_Act...
 - c.C12_FabF_Act_un - c.C14_FabF_Act_un - c.C16_FabF_Act_un - c.C18_FabF_Act_un...
 - c.C5_FabF_Act_MalACP - c.C7_FabF_Act_MalACP - c.C9_FabF_Act_MalACP - c.C11_FabF_Act_MalACP - c.C13_FabF_Act_MalACP - c.C15_FabF_Act_MalACP - c.C17_FabF_Act_MalACP - c.C19_FabF_Act_MalACP - c.C21_FabF_Act_MalACP...
 - c.C15_FabF_Act_MalACP_un - c.C17_FabF_Act_MalACP_un - c.C19_FabF_Act_MalACP_un - c.C21_FabF_Act_MalACP_un;

% FabA
c.FabA = P.FabAtot - c.FabA_ACP - c.C10_FabA_cis3EnAcACP...
 - c.C4_FabA_BHyAcACP - c.C6_FabA_BHyAcACP - c.C8_FabA_BHyAcACP - c.C10_FabA_BHyAcACP - c.C12_FabA_BHyAcACP - c.C14_FabA_BHyAcACP - c.C16_FabA_BHyAcACP - c.C18_FabA_BHyAcACP - c.C20_FabA_BHyAcACP...
 - c.C12_FabA_BHyAcACP_un - c.C14_FabA_BHyAcACP_un - c.C16_FabA_BHyAcACP_un - c.C18_FabA_BHyAcACP_un - c.C20_FabA_BHyAcACP_un...
 - c.C4_FabA_EnAcACP - c.C6_FabA_EnAcACP - c.C8_FabA_EnAcACP - c.C10_FabA_EnAcACP - c.C12_FabA_EnAcACP - c.C14_FabA_EnAcACP - c.C16_FabA_EnAcACP - c.C18_FabA_EnAcACP - c.C20_FabA_EnAcACP...
 - c.C12_FabA_EnAcACP_un - c.C14_FabA_EnAcACP_un - c.C16_FabA_EnAcACP_un - c.C18_FabA_EnAcACP_un - c.C20_FabA_EnAcACP_un;

% FabB
c.FabB = P.FabBtot - c.FabB_ACP - c.C2_FabB_AcCoA - c.C2_FabB_Act - c.C5_FabB_Act_MalACP - c.C3_FabB_MalACP - c.C2_FabB_AcACP...
 - c.C4_FabB_AcACP - c.C6_FabB_AcACP - c.C8_FabB_AcACP - c.C10_FabB_AcACP - c.C12_FabB_AcACP - c.C14_FabB_AcACP - c.C16_FabB_AcACP - c.C18_FabB_AcACP...
 - c.C12_FabB_AcACP_un - c.C14_FabB_AcACP_un - c.C16_FabB_AcACP_un - c.C18_FabB_AcACP_un...
 - c.C4_FabB_Act - c.C6_FabB_Act - c.C8_FabB_Act - c.C10_FabB_Act - c.C12_FabB_Act - c.C14_FabB_Act - c.C16_FabB_Act - c.C18_FabB_Act...
 - c.C12_FabB_Act_un - c.C14_FabB_Act_un - c.C16_FabB_Act_un - c.C18_FabB_Act_un...
 - c.C7_FabB_Act_MalACP - c.C9_FabB_Act_MalACP - c.C11_FabB_Act_MalACP - c.C13_FabB_Act_MalACP - c.C15_FabB_Act_MalACP - c.C17_FabB_Act_MalACP - c.C19_FabB_Act_MalACP - c.C21_FabB_Act_MalACP...
 - c.C15_FabB_Act_MalACP_un - c.C17_FabB_Act_MalACP_un - c.C19_FabB_Act_MalACP_un - c.C21_FabB_Act_MalACP_un...
 - c.C10_FabB_cis3EnAcACP - c.C10_FabB_Act_cis3 - c.C13_FabB_Act_cis3MalACP;

% Set of differential equations
% ATP % changed for ACC
d_ATP = -P.kcat1_1.*c.ACC_C.*c.ATP./(P.Km1_1 + c.ATP);

% Bicarbonate % changed for ACC
d_C1_Bicarbonate = -P.kcat1_2.*c.BC_ATP.*c.C1_Bicarbonate./(P.Km1_2 + c.C1_Bicarbonate);

% Acetyl-CoA (- ACC - FabH - FabF - FabB) % changed for ACC
d_C2_AcCoA = -P.kcat1_5.*c.C1_CT_Act.*c.C2_AcCoA./(P.Km1_5 + c.C2_AcCoA) + P.k3_1r(1).*c.C2_FabH_CoA - P.k3_1f(1).*c.FabH.*c.C2_AcCoA + P.k8_4r.*c.C2_FabF_AcCoA - P.k8_4f.*c.FabF.*c.C2_AcCoA + P.k10_4r.*c.C2_FabB_AcCoA - P.k10_4f.*c.FabB.*c.C2_AcCoA;

% C2n (n=2:9) Acyl-CoAs (- FabH)
d_C4_SucCoA = P.k3_1r(2).*c.C4_FabH_CoA - P.k3_1f(2).*c.FabH.*c.C4_SucCoA;
d_C6_HexCoA = P.k3_1r(3).*c.C6_FabH_CoA - P.k3_1f(3).*c.FabH.*c.C6_HexCoA;
d_C8_OcCoA = P.k3_1r(4).*c.C8_FabH_CoA - P.k3_1f(4).*c.FabH.*c.C8_OcCoA;
d_C10_DecCoA = P.k3_1r(5).*c.C10_FabH_CoA - P.k3_1f(5).*c.FabH.*c.C10_DecCoA;
d_C12_LauCoA = P.k3_1r(6).*c.C12_FabH_CoA - P.k3_1f(6).*c.FabH.*c.C12_LauCoA;
d_C14_EthCoA = P.k3_1r(7).*c.C14_FabH_CoA - P.k3_1f(7).*c.FabH.*c.C14_EthCoA;
d_C16_PalCoA = P.k3_1r(8).*c.C16_FabH_CoA - P.k3_1f(8).*c.FabH.*c.C16_PalCoA;
d_C18_OcDecCoA = P.k3_1r(9).*c.C18_FabH_CoA - P.k3_1f(9).*c.FabH.*c.C18_OcDecCoA;

% ACP (- FabD + TesA + FabF + FabB - Inhibition (H, G, Z, I, T, F, A, B))
d_ACP = P.k2_3r.*c.C3_FabD_Act_ACP - P.k2_3f.*c.C3_FabD_Act.*c.ACP...
 + P.kcat7(1).*c.C4_TesA_AcACP + P.kcat7(2).*c.C6_TesA_AcACP + P.kcat7(3).*c.C8_TesA_AcACP + P.kcat7(4).*c.C10_TesA_AcACP...
 + P.kcat7(5).*c.C12_TesA_AcACP + P.kcat7(6).*c.C14_TesA_AcACP + P.kcat7(7).*c.C16_TesA_AcACP + P.kcat7(8).*c.C18_TesA_AcACP + P.kcat7(9).*c.C20_TesA_AcACP...
 + P.kcat7(5).*c.C12_TesA_AcACP_un + P.kcat7(6).*c.C14_TesA_AcACP_un + P.kcat7(7).*c.C16_TesA_AcACP_un + P.kcat7(8).*c.C18_TesA_AcACP_un + P.kcat7(9).*c.C20_TesA_AcACP_un...
 + P.k8_2f(1).*c.C4_FabF_AcACP - P.k8_2r(1).*c.C4_FabF_Act.*c.ACP...
 + P.k8_2f(2).*c.C6_FabF_AcACP - P.k8_2r(2).*c.C6_FabF_Act.*c.ACP...
 + P.k8_2f(3).*c.C8_FabF_AcACP - P.k8_2r(3).*c.C8_FabF_Act.*c.ACP...
 + P.k8_2f(4).*c.C10_FabF_AcACP - P.k8_2r(4).*c.C10_FabF_Act.*c.ACP...
 + P.k8_2f(5).*c.C12_FabF_AcACP - P.k8_2r(5).*c.C12_FabF_Act.*c.ACP...
 + P.k8_2f(6).*c.C14_FabF_AcACP - P.k8_2r(6).*c.C14_FabF_Act.*c.ACP...
 + P.k8_2f(7).*c.C16_FabF_AcACP - P.k8_2r(7).*c.C16_FabF_Act.*c.ACP...
 + P.k8_2f(8).*c.C18_FabF_AcACP - P.k8_2r(8).*c.C18_FabF_Act.*c.ACP...
 + P.k8_2f(5).*c.C12_FabF_AcACP_un - P.k8_2r(5).*c.C12_FabF_Act_un.*c.ACP...
 + P.k8_2f(6).*c.C14_FabF_AcACP_un - P.k8_2r(6).*c.C14_FabF_Act_un.*c.ACP...
 + P.k8_2f(7).*c.C16_FabF_AcACP_un - P.k8_2r(7).*c.C16_FabF_Act_un.*c.ACP...
 + P.k8_2f(8).*c.C18_FabF_AcACP_un - P.k8_2r(8).*c.C18_FabF_Act_un.*c.ACP...
 + P.k10_2f(1).*c.C4_FabB_AcACP - P.k10_2r(1).*c.C4_FabB_Act.*c.ACP...
 + P.k10_2f(2).*c.C6_FabB_AcACP - P.k10_2r(2).*c.C6_FabB_Act.*c.ACP...
 + P.k10_2f(3).*c.C8_FabB_AcACP - P.k10_2r(3).*c.C8_FabB_Act.*c.ACP...
 + P.k10_2f(4).*c.C10_FabB_AcACP - P.k10_2r(4).*c.C10_FabB_Act.*c.ACP...
 + P.k10_2f(5).*c.C12_FabB_AcACP - P.k10_2r(5).*c.C12_FabB_Act.*c.ACP...
 + P.k10_2f(6).*c.C14_FabB_AcACP - P.k10_2r(6).*c.C14_FabB_Act.*c.ACP...
 + P.k10_2f(7).*c.C16_FabB_AcACP - P.k10_2r(7).*c.C16_FabB_Act.*c.ACP...
 + P.k10_2f(8).*c.C18_FabB_AcACP - P.k10_2r(8).*c.C18_FabB_Act.*c.ACP...
 + P.k10_2f(5).*c.C12_FabB_AcACP_un - P.k10_2r(5).*c.C12_FabB_Act_un.*c.ACP...
 + P.k10_2f(6).*c.C14_FabB_AcACP_un - P.k10_2r(6).*c.C14_FabB_Act_un.*c.ACP...
 + P.k10_2f(7).*c.C16_FabB_AcACP_un - P.k10_2r(7).*c.C16_FabB_Act_un.*c.ACP...
 + P.k10_2f(8).*c.C18_FabB_AcACP_un - P.k10_2r(8).*c.C18_FabB_Act_un.*c.ACP...
 + P.k10_2f(4).*c.C10_FabB_cis3EnAcACP - P.k10_2r(4).*c.C10_FabB_Act_cis3.*c.ACP...
 + P.k3_inh_r.*c.FabH_ACP - P.k3_inh_f.*c.FabH.*c.ACP...
 + P.k4_inh_r.*c.FabG_ACP - P.k4_inh_f.*c.FabG.*c.ACP...
 + P.k5_inh_r.*c.FabZ_ACP - P.k5_inh_f.*c.FabZ.*c.ACP...
 + P.k6_inh_r.*c.FabI_ACP - P.k6_inh_f.*c.FabI.*c.ACP...
 + P.k7_inh_r.*c.TesA_ACP - P.k7_inh_f.*c.TesA.*c.ACP...
 + P.k8_inh_r.*c.FabF_ACP - P.k8_inh_f.*c.FabF.*c.ACP...
 + P.k9_inh_r.*c.FabA_ACP - P.k9_inh_f.*c.FabA.*c.ACP...
 + P.k10_inh_r.*c.FabB_ACP - P.k10_inh_f.*c.FabB.*c.ACP...
 + P.k8_9f.*c.C2_FabF_AcACP - P.k8_9r.*c.C2_FabF_Act.*c.ACP...
 + P.k10_9f.*c.C2_FabB_AcACP - P.k10_9r.*c.C2_FabB_Act.*c.ACP...
 + P.k3_8f.*c.C2_FabH_AcACP - P.k3_8r.*c.C2_FabH_Act.*c.ACP;

% NADPH (- FabG)
d_NADPH = P.k4_1r(1).*c.FabG_NADPH - P.k4_1f(1).*c.FabG.*c.NADPH; 

% NADP+ (FabG)
d_NADP = P.kcat4(1).*c.C4_FabG_NADPH_BKeACP + P.kcat4(2).*c.C6_FabG_NADPH_BKeACP + P.kcat4(3).*c.C8_FabG_NADPH_BKeACP + P.kcat4(4).*c.C10_FabG_NADPH_BKeACP...
 + P.kcat4(5).*c.C12_FabG_NADPH_BKeACP + P.kcat4(6).*c.C14_FabG_NADPH_BKeACP + P.kcat4(7).*c.C16_FabG_NADPH_BKeACP + P.kcat4(8).*c.C18_FabG_NADPH_BKeACP + P.kcat4(9).*c.C20_FabG_NADPH_BKeACP...
 + P.kcat4(5).*c.C12_FabG_NADPH_BKeACP_un + P.kcat4(6).*c.C14_FabG_NADPH_BKeACP_un + P.kcat4(7).*c.C16_FabG_NADPH_BKeACP_un + P.kcat4(8).*c.C18_FabG_NADPH_BKeACP_un + P.kcat4(9).*c.C20_FabG_NADPH_BKeACP_un;

% NADH (- FabI)
d_NADH = P.k6_1r(1).*c.FabI_NADH - P.k6_1f(1).*c.FabI.*c.NADH; 

% NAD+ (FabI)
d_NAD = P.kcat6(1).*c.C4_FabI_NADH_EnAcACP + P.kcat6(2).*c.C6_FabI_NADH_EnAcACP + P.kcat6(3).*c.C8_FabI_NADH_EnAcACP + P.kcat6(4).*c.C10_FabI_NADH_EnAcACP...
 + P.kcat6(5).*c.C12_FabI_NADH_EnAcACP + P.kcat6(6).*c.C14_FabI_NADH_EnAcACP + P.kcat6(7).*c.C16_FabI_NADH_EnAcACP + P.kcat6(8).*c.C18_FabI_NADH_EnAcACP + P.kcat6(9).*c.C20_FabI_NADH_EnAcACP...
 + P.kcat6(5).*c.C12_FabI_NADH_EnAcACP_un + P.kcat6(6).*c.C14_FabI_NADH_EnAcACP_un + P.kcat6(7).*c.C16_FabI_NADH_EnAcACP_un + P.kcat6(8).*c.C18_FabI_NADH_EnAcACP_un + P.kcat6(9).*c.C20_FabI_NADH_EnAcACP_un; 

% ADP (ACC) % changed for ACC
d_ADP = P.kcat1_2.*c.BC_ATP.*c.C1_Bicarbonate./(P.Km1_2 + c.C1_Bicarbonate);

% Malonyl-CoA (ACC - FabD) % changed for ACC
d_C3_MalCoA = P.kcat1_5.*c.C1_CT_Act.*c.C2_AcCoA./(P.Km1_5 + c.C2_AcCoA) + P.k2_1r.*c.C3_FabD_MalCoA - P.k2_1f.*c.FabD.*c.C3_MalCoA; 

% CoA (FabD + FabH + FabF + FabB)
d_CoA = P.k2_2f.*c.C3_FabD_MalCoA - P.k2_2r.*c.C3_FabD_Act.*c.CoA...
 + P.k3_2f(1).*c.C2_FabH_CoA - P.k3_2r(1).*c.C2_FabH_Act.*c.CoA...
 + P.k3_2f(2).*c.C4_FabH_CoA - P.k3_2r(2).*c.C4_FabH_Act.*c.CoA...
 + P.k3_2f(3).*c.C6_FabH_CoA - P.k3_2r(3).*c.C6_FabH_Act.*c.CoA...
 + P.k3_2f(4).*c.C8_FabH_CoA - P.k3_2r(4).*c.C8_FabH_Act.*c.CoA...
 + P.k3_2f(5).*c.C10_FabH_CoA - P.k3_2r(5).*c.C10_FabH_Act.*c.CoA...
 + P.k3_2f(6).*c.C12_FabH_CoA - P.k3_2r(6).*c.C12_FabH_Act.*c.CoA...
 + P.k3_2f(7).*c.C14_FabH_CoA - P.k3_2r(7).*c.C14_FabH_Act.*c.CoA...
 + P.k3_2f(8).*c.C16_FabH_CoA - P.k3_2r(8).*c.C16_FabH_Act.*c.CoA...
 + P.k3_2f(9).*c.C18_FabH_CoA - P.k3_2r(9).*c.C18_FabH_Act.*c.CoA...
 + P.k8_5f.*c.C2_FabF_AcCoA - P.k8_5r*c.C2_FabF_Act.*c.CoA...
 + P.k10_5f.*c.C2_FabB_AcCoA - P.k10_5r.*c.C2_FabB_Act.*c.CoA;

% Malonyl-ACP (FabD - FabH - FabF - FabB)
d_C3_MalACP = P.k2_4f.*c.C3_FabD_Act_ACP - P.k2_4r.*c.FabD.*c.C3_MalACP...
 + P.k3_3r(1).*c.C5_FabH_Act_MalACP - P.k3_3f(1).*c.C2_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(2).*c.C7_FabH_Act_MalACP - P.k3_3f(2).*c.C4_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(3).*c.C9_FabH_Act_MalACP - P.k3_3f(3).*c.C6_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(4).*c.C11_FabH_Act_MalACP - P.k3_3f(4).*c.C8_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(5).*c.C13_FabH_Act_MalACP - P.k3_3f(5).*c.C10_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(6).*c.C15_FabH_Act_MalACP - P.k3_3f(6).*c.C12_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(7).*c.C17_FabH_Act_MalACP - P.k3_3f(7).*c.C14_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(8).*c.C19_FabH_Act_MalACP - P.k3_3f(8).*c.C16_FabH_Act.*c.C3_MalACP...
 + P.k3_3r(9).*c.C21_FabH_Act_MalACP - P.k3_3f(9).*c.C18_FabH_Act.*c.C3_MalACP...
 + P.k8_3r(1).*c.C7_FabF_Act_MalACP - P.k8_3f(1).*c.C4_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(2).*c.C9_FabF_Act_MalACP - P.k8_3f(2).*c.C6_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(3).*c.C11_FabF_Act_MalACP - P.k8_3f(3).*c.C8_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(4).*c.C13_FabF_Act_MalACP - P.k8_3f(4).*c.C10_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(5).*c.C15_FabF_Act_MalACP - P.k8_3f(5).*c.C12_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(6).*c.C17_FabF_Act_MalACP - P.k8_3f(6).*c.C14_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(7).*c.C19_FabF_Act_MalACP - P.k8_3f(7).*c.C16_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(8).*c.C21_FabF_Act_MalACP - P.k8_3f(8).*c.C18_FabF_Act.*c.C3_MalACP...
 + P.k8_3r(5).*c.C15_FabF_Act_MalACP_un - P.k8_3f(5).*c.C12_FabF_Act_un.*c.C3_MalACP...
 + P.k8_3r(6).*c.C17_FabF_Act_MalACP_un - P.k8_3f(6).*c.C14_FabF_Act_un.*c.C3_MalACP...
 + P.k8_3r(7).*c.C19_FabF_Act_MalACP_un - P.k8_3f(7).*c.C16_FabF_Act_un.*c.C3_MalACP...
 + P.k8_3r(8).*c.C21_FabF_Act_MalACP_un - P.k8_3f(8).*c.C18_FabF_Act_un.*c.C3_MalACP...
 + P.k10_3r(1).*c.C7_FabB_Act_MalACP - P.k10_3f(1).*c.C4_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(2).*c.C9_FabB_Act_MalACP - P.k10_3f(2).*c.C6_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(3).*c.C11_FabB_Act_MalACP - P.k10_3f(3).*c.C8_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(4).*c.C13_FabB_Act_MalACP - P.k10_3f(4).*c.C10_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(5).*c.C15_FabB_Act_MalACP - P.k10_3f(5).*c.C12_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(6).*c.C17_FabB_Act_MalACP - P.k10_3f(6).*c.C14_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(7).*c.C19_FabB_Act_MalACP - P.k10_3f(7).*c.C16_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(8).*c.C21_FabB_Act_MalACP - P.k10_3f(8).*c.C18_FabB_Act.*c.C3_MalACP...
 + P.k10_3r(5).*c.C15_FabB_Act_MalACP_un - P.k10_3f(5).*c.C12_FabB_Act_un.*c.C3_MalACP...
 + P.k10_3r(6).*c.C17_FabB_Act_MalACP_un - P.k10_3f(6).*c.C14_FabB_Act_un.*c.C3_MalACP...
 + P.k10_3r(7).*c.C19_FabB_Act_MalACP_un - P.k10_3f(7).*c.C16_FabB_Act_un.*c.C3_MalACP...
 + P.k10_3r(8).*c.C21_FabB_Act_MalACP_un - P.k10_3f(8).*c.C18_FabB_Act_un.*c.C3_MalACP...
 + P.k10_3r(4).*c.C13_FabB_Act_cis3MalACP - P.k10_3f(4).*c.C10_FabB_Act_cis3.*c.C3_MalACP...
 + P.k8_6r.*c.C5_FabF_Act_MalACP - P.k8_6f.*c.C2_FabF_Act.*c.C3_MalACP...
 + P.k10_6r.*c.C5_FabB_Act_MalACP - P.k10_6f.*c.C2_FabB_Act.*c.C3_MalACP...
 + P.k8_7r.*c.C3_FabF_MalACP - P.k8_7f.*c.FabF.*c.C3_MalACP...
 + P.k10_7r.*c.C3_FabB_MalACP - P.k10_7f.*c.FabB.*c.C3_MalACP...
 + P.k3_6r.*c.C3_FabH_MalACP - P.k3_6f.*c.FabH.*c.C3_MalACP;

% CO2 (FabH + FabF + FabB)
d_C1_CO2 = P.kcat3(1).*c.C5_FabH_Act_MalACP + P.kcat3(2).*c.C7_FabH_Act_MalACP + P.kcat3(3).*c.C9_FabH_Act_MalACP + P.kcat3(4).*c.C11_FabH_Act_MalACP...
 + P.kcat3(5).*c.C13_FabH_Act_MalACP + P.kcat3(6).*c.C15_FabH_Act_MalACP + P.kcat3(7).*c.C17_FabH_Act_MalACP + P.kcat3(8).*c.C19_FabH_Act_MalACP + P.kcat3(9).*c.C21_FabH_Act_MalACP...
 + P.kcat8(1).*c.C7_FabF_Act_MalACP + P.kcat8(2).*c.C9_FabF_Act_MalACP + P.kcat8(3).*c.C11_FabF_Act_MalACP + P.kcat8(4).*c.C13_FabF_Act_MalACP...
 + P.kcat8(5).*c.C15_FabF_Act_MalACP + P.kcat8(6).*c.C17_FabF_Act_MalACP + P.kcat8(7).*c.C19_FabF_Act_MalACP + P.kcat8(8).*c.C21_FabF_Act_MalACP...
 + P.kcat8_un(5).*c.C15_FabF_Act_MalACP_un + P.kcat8_un(6).*c.C17_FabF_Act_MalACP_un + P.kcat8_un(7).*c.C19_FabF_Act_MalACP_un + P.kcat8_un(8).*c.C21_FabF_Act_MalACP_un...
 + P.kcat10(1).*c.C7_FabB_Act_MalACP + P.kcat10(2).*c.C9_FabB_Act_MalACP + P.kcat10(3).*c.C11_FabB_Act_MalACP + P.kcat10(4).*c.C13_FabB_Act_MalACP...
 + P.kcat10(5).*c.C15_FabB_Act_MalACP + P.kcat10(6).*c.C17_FabB_Act_MalACP + P.kcat10(7).*c.C19_FabB_Act_MalACP + P.kcat10(8).*c.C21_FabB_Act_MalACP...
 + P.kcat10_un(4).*c.C13_FabB_Act_cis3MalACP + P.kcat10_un(5).*c.C15_FabB_Act_MalACP_un + P.kcat10_un(6).*c.C17_FabB_Act_MalACP_un + P.kcat10_un(7).*c.C19_FabB_Act_MalACP_un + P.kcat10_un(8).*c.C21_FabB_Act_MalACP_un...
 + P.kcat8_H.*c.C5_FabF_Act_MalACP + P.kcat8_CO2.*c.C3_FabF_MalACP...
 + P.kcat10_H.*c.C5_FabB_Act_MalACP + P.kcat10_CO2.*c.C3_FabB_MalACP...
 + P.kcat3_CO2.*c.C3_FabF_MalACP;

% C2n (n=2:10) B-ketoacyl-ACPs (FabH + FabF + FabB - FabG)
d_C4_BKeACP = P.kcat3(1).*c.C5_FabH_Act_MalACP + P.kcat8_H.*c.C5_FabF_Act_MalACP + P.kcat10_H.*c.C5_FabB_Act_MalACP + P.k4_2r(1).*c.C4_FabG_NADPH_BKeACP - P.k4_2f(1).*c.FabG_NADPH.*c.C4_BKeACP;
d_C6_BKeACP = P.kcat3(2).*c.C7_FabH_Act_MalACP + P.kcat8(1).*c.C7_FabF_Act_MalACP + P.kcat10(1).*c.C7_FabB_Act_MalACP + P.k4_2r(2).*c.C6_FabG_NADPH_BKeACP - P.k4_2f(2).*c.FabG_NADPH.*c.C6_BKeACP;
d_C8_BKeACP = P.kcat3(3).*c.C9_FabH_Act_MalACP + P.kcat8(2).*c.C9_FabF_Act_MalACP + P.kcat10(2).*c.C9_FabB_Act_MalACP + P.k4_2r(3).*c.C8_FabG_NADPH_BKeACP - P.k4_2f(3).*c.FabG_NADPH.*c.C8_BKeACP;
d_C10_BKeACP = P.kcat3(4).*c.C11_FabH_Act_MalACP + P.kcat8(3).*c.C11_FabF_Act_MalACP + P.kcat10(3).*c.C11_FabB_Act_MalACP + P.k4_2r(4).*c.C10_FabG_NADPH_BKeACP - P.k4_2f(4).*c.FabG_NADPH.*c.C10_BKeACP;
d_C12_BKeACP = P.kcat3(5).*c.C13_FabH_Act_MalACP + P.kcat8(4).*c.C13_FabF_Act_MalACP + P.kcat10(4).*c.C13_FabB_Act_MalACP + P.k4_2r(5).*c.C12_FabG_NADPH_BKeACP - P.k4_2f(5).*c.FabG_NADPH.*c.C12_BKeACP;
d_C14_BKeACP = P.kcat3(6).*c.C15_FabH_Act_MalACP + P.kcat8(5).*c.C15_FabF_Act_MalACP + P.kcat10(5).*c.C15_FabB_Act_MalACP + P.k4_2r(6).*c.C14_FabG_NADPH_BKeACP - P.k4_2f(6).*c.FabG_NADPH.*c.C14_BKeACP;
d_C16_BKeACP = P.kcat3(7).*c.C17_FabH_Act_MalACP + P.kcat8(6).*c.C17_FabF_Act_MalACP + P.kcat10(6).*c.C17_FabB_Act_MalACP + P.k4_2r(7).*c.C16_FabG_NADPH_BKeACP - P.k4_2f(7).*c.FabG_NADPH.*c.C16_BKeACP;
d_C18_BKeACP = P.kcat3(8).*c.C19_FabH_Act_MalACP + P.kcat8(7).*c.C19_FabF_Act_MalACP + P.kcat10(7).*c.C19_FabB_Act_MalACP + P.k4_2r(8).*c.C18_FabG_NADPH_BKeACP - P.k4_2f(8).*c.FabG_NADPH.*c.C18_BKeACP;
d_C20_BKeACP = P.kcat3(9).*c.C21_FabH_Act_MalACP + P.kcat8(8).*c.C21_FabF_Act_MalACP + P.kcat10(8).*c.C21_FabB_Act_MalACP + P.k4_2r(9).*c.C20_FabG_NADPH_BKeACP - P.k4_2f(9).*c.FabG_NADPH.*c.C20_BKeACP;

% C2n:1 (n=6:10) B-ketoacyl-ACPs (FabF + FabB - FabG)
d_C12_BKeACP_un =                                                   P.kcat10_un(4).*c.C13_FabB_Act_cis3MalACP + P.k4_2r(5).*c.C12_FabG_NADPH_BKeACP_un - P.k4_2f(5).*c.FabG_NADPH.*c.C12_BKeACP_un;
d_C14_BKeACP_un = P.kcat8_un(5).*c.C15_FabF_Act_MalACP_un + P.kcat10_un(5).*c.C15_FabB_Act_MalACP_un + P.k4_2r(6).*c.C14_FabG_NADPH_BKeACP_un - P.k4_2f(6).*c.FabG_NADPH.*c.C14_BKeACP_un;
d_C16_BKeACP_un = P.kcat8_un(6).*c.C17_FabF_Act_MalACP_un + P.kcat10_un(6).*c.C17_FabB_Act_MalACP_un + P.k4_2r(7).*c.C16_FabG_NADPH_BKeACP_un - P.k4_2f(7).*c.FabG_NADPH.*c.C16_BKeACP_un;
d_C18_BKeACP_un = P.kcat8_un(7).*c.C19_FabF_Act_MalACP_un + P.kcat10_un(7).*c.C19_FabB_Act_MalACP_un + P.k4_2r(8).*c.C18_FabG_NADPH_BKeACP_un - P.k4_2f(8).*c.FabG_NADPH.*c.C18_BKeACP_un;
d_C20_BKeACP_un = P.kcat8_un(8).*c.C21_FabF_Act_MalACP_un + P.kcat10_un(8).*c.C21_FabB_Act_MalACP_un + P.k4_2r(9).*c.C20_FabG_NADPH_BKeACP_un - P.k4_2f(9).*c.FabG_NADPH.*c.C20_BKeACP_un;

% C2n (n=2:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C4_BHyAcACP = P.kcat4(1).*c.C4_FabG_NADPH_BKeACP + P.k5_1r(1).*c.C4_FabZ_BHyAcACP - P.k5_1f(1).*c.FabZ.*c.C4_BHyAcACP + P.k9_1r(1).*c.C4_FabA_BHyAcACP - P.k9_1f(1).*c.FabA.*c.C4_BHyAcACP;
d_C6_BHyAcACP = P.kcat4(2).*c.C6_FabG_NADPH_BKeACP + P.k5_1r(2).*c.C6_FabZ_BHyAcACP - P.k5_1f(2).*c.FabZ.*c.C6_BHyAcACP + P.k9_1r(2).*c.C6_FabA_BHyAcACP - P.k9_1f(2).*c.FabA.*c.C6_BHyAcACP;
d_C8_BHyAcACP = P.kcat4(3).*c.C8_FabG_NADPH_BKeACP + P.k5_1r(3).*c.C8_FabZ_BHyAcACP - P.k5_1f(3).*c.FabZ.*c.C8_BHyAcACP + P.k9_1r(3).*c.C8_FabA_BHyAcACP - P.k9_1f(3).*c.FabA.*c.C8_BHyAcACP;
d_C10_BHyAcACP = P.kcat4(4).*c.C10_FabG_NADPH_BKeACP + P.k5_1r(4).*c.C10_FabZ_BHyAcACP - P.k5_1f(4).*c.FabZ.*c.C10_BHyAcACP + P.k9_1r(4).*c.C10_FabA_BHyAcACP - P.k9_1f(4).*c.FabA.*c.C10_BHyAcACP;
d_C12_BHyAcACP = P.kcat4(5).*c.C12_FabG_NADPH_BKeACP + P.k5_1r(5).*c.C12_FabZ_BHyAcACP - P.k5_1f(5).*c.FabZ.*c.C12_BHyAcACP + P.k9_1r(5).*c.C12_FabA_BHyAcACP - P.k9_1f(5).*c.FabA.*c.C12_BHyAcACP;
d_C14_BHyAcACP = P.kcat4(6).*c.C14_FabG_NADPH_BKeACP + P.k5_1r(6).*c.C14_FabZ_BHyAcACP - P.k5_1f(6).*c.FabZ.*c.C14_BHyAcACP + P.k9_1r(6).*c.C14_FabA_BHyAcACP - P.k9_1f(6).*c.FabA.*c.C14_BHyAcACP;
d_C16_BHyAcACP = P.kcat4(7).*c.C16_FabG_NADPH_BKeACP + P.k5_1r(7).*c.C16_FabZ_BHyAcACP - P.k5_1f(7).*c.FabZ.*c.C16_BHyAcACP + P.k9_1r(7).*c.C16_FabA_BHyAcACP - P.k9_1f(7).*c.FabA.*c.C16_BHyAcACP;
d_C18_BHyAcACP = P.kcat4(8).*c.C18_FabG_NADPH_BKeACP + P.k5_1r(8).*c.C18_FabZ_BHyAcACP - P.k5_1f(8).*c.FabZ.*c.C18_BHyAcACP + P.k9_1r(8).*c.C18_FabA_BHyAcACP - P.k9_1f(8).*c.FabA.*c.C18_BHyAcACP;
d_C20_BHyAcACP = P.kcat4(9).*c.C20_FabG_NADPH_BKeACP + P.k5_1r(9).*c.C20_FabZ_BHyAcACP - P.k5_1f(9).*c.FabZ.*c.C20_BHyAcACP + P.k9_1r(9).*c.C20_FabA_BHyAcACP - P.k9_1f(9).*c.FabA.*c.C20_BHyAcACP;

% C2n:1 (n=6:10) B-hydroxy-acyl-ACPs (FabG - FabZ - FabA)
d_C12_BHyAcACP_un = P.kcat4(5).*c.C12_FabG_NADPH_BKeACP_un + P.k5_1r(5).*c.C12_FabZ_BHyAcACP_un - P.k5_1f(5).*c.FabZ.*c.C12_BHyAcACP_un + P.k9_1r_un(5).*c.C12_FabA_BHyAcACP_un - P.k9_1f_un(5).*c.FabA.*c.C12_BHyAcACP_un;
d_C14_BHyAcACP_un = P.kcat4(6).*c.C14_FabG_NADPH_BKeACP_un + P.k5_1r(6).*c.C14_FabZ_BHyAcACP_un - P.k5_1f(6).*c.FabZ.*c.C14_BHyAcACP_un + P.k9_1r_un(6).*c.C14_FabA_BHyAcACP_un - P.k9_1f_un(6).*c.FabA.*c.C14_BHyAcACP_un;
d_C16_BHyAcACP_un = P.kcat4(7).*c.C16_FabG_NADPH_BKeACP_un + P.k5_1r(7).*c.C16_FabZ_BHyAcACP_un - P.k5_1f(7).*c.FabZ.*c.C16_BHyAcACP_un + P.k9_1r_un(7).*c.C16_FabA_BHyAcACP_un - P.k9_1f_un(7).*c.FabA.*c.C16_BHyAcACP_un;
d_C18_BHyAcACP_un = P.kcat4(8).*c.C18_FabG_NADPH_BKeACP_un + P.k5_1r(8).*c.C18_FabZ_BHyAcACP_un - P.k5_1f(8).*c.FabZ.*c.C18_BHyAcACP_un + P.k9_1r_un(8).*c.C18_FabA_BHyAcACP_un - P.k9_1f_un(8).*c.FabA.*c.C18_BHyAcACP_un;
d_C20_BHyAcACP_un = P.kcat4(9).*c.C20_FabG_NADPH_BKeACP_un + P.k5_1r(9).*c.C20_FabZ_BHyAcACP_un - P.k5_1f(9).*c.FabZ.*c.C20_BHyAcACP_un + P.k9_1r_un(9).*c.C20_FabA_BHyAcACP_un - P.k9_1f_un(9).*c.FabA.*c.C20_BHyAcACP_un;

% C2n (n=2:10) Enoyl-Acyl-ACPs (FabZ + FabA - FabI) 
d_C4_EnAcACP = P.k5_3f(1).*c.C4_FabZ_EnAcACP - P.k5_3r(1).*c.FabZ.*c.C4_EnAcACP + P.k6_2r(1).*c.C4_FabI_NADH_EnAcACP - P.k6_2f(1).*c.FabI_NADH.*c.C4_EnAcACP + P.k9_3f(1).*c.C4_FabA_EnAcACP - P.k9_3r(1).*c.FabA.*c.C4_EnAcACP;
d_C6_EnAcACP = P.k5_3f(2).*c.C6_FabZ_EnAcACP - P.k5_3r(2).*c.FabZ.*c.C6_EnAcACP + P.k6_2r(2).*c.C6_FabI_NADH_EnAcACP - P.k6_2f(2).*c.FabI_NADH.*c.C6_EnAcACP + P.k9_3f(2).*c.C6_FabA_EnAcACP - P.k9_3r(2).*c.FabA.*c.C6_EnAcACP;
d_C8_EnAcACP = P.k5_3f(3).*c.C8_FabZ_EnAcACP - P.k5_3r(3).*c.FabZ.*c.C8_EnAcACP + P.k6_2r(3).*c.C8_FabI_NADH_EnAcACP - P.k6_2f(3).*c.FabI_NADH.*c.C8_EnAcACP + P.k9_3f(3).*c.C8_FabA_EnAcACP - P.k9_3r(3).*c.FabA.*c.C8_EnAcACP;
d_C10_EnAcACP = P.k5_3f(4).*c.C10_FabZ_EnAcACP - P.k5_3r(4).*c.FabZ.*c.C10_EnAcACP + P.k6_2r(4).*c.C10_FabI_NADH_EnAcACP - P.k6_2f(4).*c.FabI_NADH.*c.C10_EnAcACP + P.k9_3f(4).*c.C10_FabA_EnAcACP - P.k9_3r(4).*c.FabA.*c.C10_EnAcACP;
d_C12_EnAcACP = P.k5_3f(5).*c.C12_FabZ_EnAcACP - P.k5_3r(5).*c.FabZ.*c.C12_EnAcACP + P.k6_2r(5).*c.C12_FabI_NADH_EnAcACP - P.k6_2f(5).*c.FabI_NADH.*c.C12_EnAcACP + P.k9_3f(5).*c.C12_FabA_EnAcACP - P.k9_3r(5).*c.FabA.*c.C12_EnAcACP;
d_C14_EnAcACP = P.k5_3f(6).*c.C14_FabZ_EnAcACP - P.k5_3r(6).*c.FabZ.*c.C14_EnAcACP + P.k6_2r(6).*c.C14_FabI_NADH_EnAcACP - P.k6_2f(6).*c.FabI_NADH.*c.C14_EnAcACP + P.k9_3f(6).*c.C14_FabA_EnAcACP - P.k9_3r(6).*c.FabA.*c.C14_EnAcACP;
d_C16_EnAcACP = P.k5_3f(7).*c.C16_FabZ_EnAcACP - P.k5_3r(7).*c.FabZ.*c.C16_EnAcACP + P.k6_2r(7).*c.C16_FabI_NADH_EnAcACP - P.k6_2f(7).*c.FabI_NADH.*c.C16_EnAcACP + P.k9_3f(7).*c.C16_FabA_EnAcACP - P.k9_3r(7).*c.FabA.*c.C16_EnAcACP;
d_C18_EnAcACP = P.k5_3f(8).*c.C18_FabZ_EnAcACP - P.k5_3r(8).*c.FabZ.*c.C18_EnAcACP + P.k6_2r(8).*c.C18_FabI_NADH_EnAcACP - P.k6_2f(8).*c.FabI_NADH.*c.C18_EnAcACP + P.k9_3f(8).*c.C18_FabA_EnAcACP - P.k9_3r(8).*c.FabA.*c.C18_EnAcACP;
d_C20_EnAcACP = P.k5_3f(9).*c.C20_FabZ_EnAcACP - P.k5_3r(9).*c.FabZ.*c.C20_EnAcACP + P.k6_2r(9).*c.C20_FabI_NADH_EnAcACP - P.k6_2f(9).*c.FabI_NADH.*c.C20_EnAcACP + P.k9_3f(9).*c.C20_FabA_EnAcACP - P.k9_3r(9).*c.FabA.*c.C20_EnAcACP;

% C10 cis-3-Enoyl-Acyl-ACP (FabA - FabB)
d_C10_cis3EnAcACP = P.k9_3f_un(4).*c.C10_FabA_cis3EnAcACP - P.k9_3r_un(4).*c.FabA.*c.C10_cis3EnAcACP + P.k10_1r(4).*c.C10_FabB_cis3EnAcACP - P.k10_1f(4).*c.FabB.*c.C10_cis3EnAcACP;

% C2n:1 (n=6:10) Enoyl-Acyl-ACPs (FabZ + FabA - FabI)
d_C12_EnAcACP_un = P.k5_3f(5).*c.C12_FabZ_EnAcACP_un - P.k5_3r(5).*c.FabZ.*c.C12_EnAcACP_un + P.k6_2r(5).*c.C12_FabI_NADH_EnAcACP_un - P.k6_2f(5).*c.FabI_NADH.*c.C12_EnAcACP_un + P.k9_3f(5).*c.C12_FabA_EnAcACP_un - P.k9_3r(5).*c.FabA.*c.C12_EnAcACP_un;
d_C14_EnAcACP_un = P.k5_3f(6).*c.C14_FabZ_EnAcACP_un - P.k5_3r(6).*c.FabZ.*c.C14_EnAcACP_un + P.k6_2r(6).*c.C14_FabI_NADH_EnAcACP_un - P.k6_2f(6).*c.FabI_NADH.*c.C14_EnAcACP_un + P.k9_3f(6).*c.C14_FabA_EnAcACP_un - P.k9_3r(6).*c.FabA.*c.C14_EnAcACP_un;
d_C16_EnAcACP_un = P.k5_3f(7).*c.C16_FabZ_EnAcACP_un - P.k5_3r(7).*c.FabZ.*c.C16_EnAcACP_un + P.k6_2r(7).*c.C16_FabI_NADH_EnAcACP_un - P.k6_2f(7).*c.FabI_NADH.*c.C16_EnAcACP_un + P.k9_3f(7).*c.C16_FabA_EnAcACP_un - P.k9_3r(7).*c.FabA.*c.C16_EnAcACP_un;
d_C18_EnAcACP_un = P.k5_3f(8).*c.C18_FabZ_EnAcACP_un - P.k5_3r(8).*c.FabZ.*c.C18_EnAcACP_un + P.k6_2r(8).*c.C18_FabI_NADH_EnAcACP_un - P.k6_2f(8).*c.FabI_NADH.*c.C18_EnAcACP_un + P.k9_3f(8).*c.C18_FabA_EnAcACP_un - P.k9_3r(8).*c.FabA.*c.C18_EnAcACP_un;
d_C20_EnAcACP_un = P.k5_3f(9).*c.C20_FabZ_EnAcACP_un - P.k5_3r(9).*c.FabZ.*c.C20_EnAcACP_un + P.k6_2r(9).*c.C20_FabI_NADH_EnAcACP_un - P.k6_2f(9).*c.FabI_NADH.*c.C20_EnAcACP_un + P.k9_3f(9).*c.C20_FabA_EnAcACP_un - P.k9_3r(9).*c.FabA.*c.C20_EnAcACP_un;

% C2n (n=2:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH)
d_C4_AcACP = P.kcat6(1).*c.C4_FabI_NADH_EnAcACP + P.k3_4r(1).*c.C4_FabH_AcACP - P.k3_4f(1).*c.FabH.*c.C4_AcACP + P.k3_5r(1).*c.C6_FabH_Act_AcACP - P.k3_5f(1).*c.C2_FabH_Act.*c.C4_AcACP + P.k7_1r(1).*c.C4_TesA_AcACP + P.k8_1r(1).*c.C4_FabF_AcACP + P.k10_1r(1).*c.C4_FabB_AcACP - P.k7_1f(1).*c.TesA.*c.C4_AcACP - P.k8_1f(1).*c.FabF.*c.C4_AcACP - P.k10_1f(1).*c.FabB.*c.C4_AcACP;
d_C6_AcACP = P.kcat6(2).*c.C6_FabI_NADH_EnAcACP + P.k3_4r(2).*c.C6_FabH_AcACP - P.k3_4f(2).*c.FabH.*c.C6_AcACP + P.k3_5r(2).*c.C8_FabH_Act_AcACP - P.k3_5f(2).*c.C2_FabH_Act.*c.C6_AcACP + P.k7_1r(2).*c.C6_TesA_AcACP + P.k8_1r(2).*c.C6_FabF_AcACP + P.k10_1r(2).*c.C6_FabB_AcACP - P.k7_1f(2).*c.TesA.*c.C6_AcACP - P.k8_1f(2).*c.FabF.*c.C6_AcACP - P.k10_1f(2).*c.FabB.*c.C6_AcACP;
d_C8_AcACP = P.kcat6(3).*c.C8_FabI_NADH_EnAcACP + P.k3_4r(3).*c.C8_FabH_AcACP - P.k3_4f(3).*c.FabH.*c.C8_AcACP + P.k3_5r(3).*c.C10_FabH_Act_AcACP - P.k3_5f(3).*c.C2_FabH_Act.*c.C8_AcACP + P.k7_1r(3).*c.C8_TesA_AcACP + P.k8_1r(3).*c.C8_FabF_AcACP + P.k10_1r(3).*c.C8_FabB_AcACP - P.k7_1f(3).*c.TesA.*c.C8_AcACP - P.k8_1f(3).*c.FabF.*c.C8_AcACP - P.k10_1f(3).*c.FabB.*c.C8_AcACP;
d_C10_AcACP = P.kcat6(4).*c.C10_FabI_NADH_EnAcACP + P.k3_4r(4).*c.C10_FabH_AcACP - P.k3_4f(4).*c.FabH.*c.C10_AcACP + P.k3_5r(4).*c.C12_FabH_Act_AcACP - P.k3_5f(4).*c.C2_FabH_Act.*c.C10_AcACP + P.k7_1r(4).*c.C10_TesA_AcACP + P.k8_1r(4).*c.C10_FabF_AcACP + P.k10_1r(4).*c.C10_FabB_AcACP - P.k7_1f(4).*c.TesA.*c.C10_AcACP - P.k8_1f(4).*c.FabF.*c.C10_AcACP - P.k10_1f(4).*c.FabB.*c.C10_AcACP;
d_C12_AcACP = P.kcat6(5).*c.C12_FabI_NADH_EnAcACP + P.k3_4r(5).*c.C12_FabH_AcACP - P.k3_4f(5).*c.FabH.*c.C12_AcACP + P.k3_5r(5).*c.C14_FabH_Act_AcACP - P.k3_5f(5).*c.C2_FabH_Act.*c.C12_AcACP + P.k7_1r(5).*c.C12_TesA_AcACP + P.k8_1r(5).*c.C12_FabF_AcACP + P.k10_1r(5).*c.C12_FabB_AcACP - P.k7_1f(5).*c.TesA.*c.C12_AcACP - P.k8_1f(5).*c.FabF.*c.C12_AcACP - P.k10_1f(5).*c.FabB.*c.C12_AcACP;
d_C14_AcACP = P.kcat6(6).*c.C14_FabI_NADH_EnAcACP + P.k3_4r(6).*c.C14_FabH_AcACP - P.k3_4f(6).*c.FabH.*c.C14_AcACP + P.k3_5r(6).*c.C16_FabH_Act_AcACP - P.k3_5f(6).*c.C2_FabH_Act.*c.C14_AcACP + P.k7_1r(6).*c.C14_TesA_AcACP + P.k8_1r(6).*c.C14_FabF_AcACP + P.k10_1r(6).*c.C14_FabB_AcACP - P.k7_1f(6).*c.TesA.*c.C14_AcACP - P.k8_1f(6).*c.FabF.*c.C14_AcACP - P.k10_1f(6).*c.FabB.*c.C14_AcACP;
d_C16_AcACP = P.kcat6(7).*c.C16_FabI_NADH_EnAcACP + P.k3_4r(7).*c.C16_FabH_AcACP - P.k3_4f(7).*c.FabH.*c.C16_AcACP + P.k3_5r(7).*c.C18_FabH_Act_AcACP - P.k3_5f(7).*c.C2_FabH_Act.*c.C16_AcACP + P.k7_1r(7).*c.C16_TesA_AcACP + P.k8_1r(7).*c.C16_FabF_AcACP + P.k10_1r(7).*c.C16_FabB_AcACP - P.k7_1f(7).*c.TesA.*c.C16_AcACP - P.k8_1f(7).*c.FabF.*c.C16_AcACP - P.k10_1f(7).*c.FabB.*c.C16_AcACP;
d_C18_AcACP = P.kcat6(8).*c.C18_FabI_NADH_EnAcACP + P.k3_4r(8).*c.C18_FabH_AcACP - P.k3_4f(8).*c.FabH.*c.C18_AcACP + P.k3_5r(8).*c.C20_FabH_Act_AcACP - P.k3_5f(8).*c.C2_FabH_Act.*c.C18_AcACP + P.k7_1r(8).*c.C18_TesA_AcACP + P.k8_1r(8).*c.C18_FabF_AcACP + P.k10_1r(8).*c.C18_FabB_AcACP - P.k7_1f(8).*c.TesA.*c.C18_AcACP - P.k8_1f(8).*c.FabF.*c.C18_AcACP - P.k10_1f(8).*c.FabB.*c.C18_AcACP;

% C2n (n=10) Acyl-ACPs (FabI - TesA - FabH)
d_C20_AcACP = P.kcat6(9).*c.C20_FabI_NADH_EnAcACP + P.k3_4r(9).*c.C20_FabH_AcACP - P.k3_4f(9).*c.FabH.*c.C20_AcACP + P.k3_5r(9).*c.C22_FabH_Act_AcACP - P.k3_5f(9).*c.C2_FabH_Act.*c.C20_AcACP + P.k7_1r(9).*c.C20_TesA_AcACP - P.k7_1f(9).*c.TesA.*c.C20_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (FabI - TesA - FabF - FabB - FabH)
d_C12_AcACP_un = P.kcat6(5).*c.C12_FabI_NADH_EnAcACP_un + P.k3_4r(5).*c.C12_FabH_AcACP_un - P.k3_4f(5).*c.FabH.*c.C12_AcACP_un + P.k3_5r(5).*c.C14_FabH_Act_AcACP_un - P.k3_5f(5).*c.C2_FabH_Act.*c.C12_AcACP_un + P.k7_1r(5).*c.C12_TesA_AcACP_un - P.k7_1f(5).*c.TesA.*c.C12_AcACP_un + P.k8_1r(5).*c.C12_FabF_AcACP_un - P.k8_1f(5).*c.FabF.*c.C12_AcACP_un + P.k10_1r(5).*c.C12_FabB_AcACP_un - P.k10_1f(5).*c.FabB.*c.C12_AcACP_un;
d_C14_AcACP_un = P.kcat6(6).*c.C14_FabI_NADH_EnAcACP_un + P.k3_4r(6).*c.C14_FabH_AcACP_un - P.k3_4f(6).*c.FabH.*c.C14_AcACP_un + P.k3_5r(6).*c.C16_FabH_Act_AcACP_un - P.k3_5f(6).*c.C2_FabH_Act.*c.C14_AcACP_un + P.k7_1r(6).*c.C14_TesA_AcACP_un - P.k7_1f(6).*c.TesA.*c.C14_AcACP_un + P.k8_1r(6).*c.C14_FabF_AcACP_un - P.k8_1f(6).*c.FabF.*c.C14_AcACP_un + P.k10_1r(6).*c.C14_FabB_AcACP_un - P.k10_1f(6).*c.FabB.*c.C14_AcACP_un;
d_C16_AcACP_un = P.kcat6(7).*c.C16_FabI_NADH_EnAcACP_un + P.k3_4r(7).*c.C16_FabH_AcACP_un - P.k3_4f(7).*c.FabH.*c.C16_AcACP_un + P.k3_5r(7).*c.C18_FabH_Act_AcACP_un - P.k3_5f(7).*c.C2_FabH_Act.*c.C16_AcACP_un + P.k7_1r(7).*c.C16_TesA_AcACP_un - P.k7_1f(7).*c.TesA.*c.C16_AcACP_un + P.k8_1r(7).*c.C16_FabF_AcACP_un - P.k8_1f(7).*c.FabF.*c.C16_AcACP_un + P.k10_1r(7).*c.C16_FabB_AcACP_un - P.k10_1f(7).*c.FabB.*c.C16_AcACP_un;
d_C18_AcACP_un = P.kcat6(8).*c.C18_FabI_NADH_EnAcACP_un + P.k3_4r(8).*c.C18_FabH_AcACP_un - P.k3_4f(8).*c.FabH.*c.C18_AcACP_un + P.k3_5r(8).*c.C20_FabH_Act_AcACP_un - P.k3_5f(8).*c.C2_FabH_Act.*c.C18_AcACP_un + P.k7_1r(8).*c.C18_TesA_AcACP_un - P.k7_1f(8).*c.TesA.*c.C18_AcACP_un + P.k8_1r(8).*c.C18_FabF_AcACP_un - P.k8_1f(8).*c.FabF.*c.C18_AcACP_un + P.k10_1r(8).*c.C18_FabB_AcACP_un - P.k10_1f(8).*c.FabB.*c.C18_AcACP_un;

% C2n:1 (n=10) Acyl-ACPs (FabI - TesA - FabH)
d_C20_AcACP_un = P.kcat6(9).*c.C20_FabI_NADH_EnAcACP_un + P.k7_1r(9).*c.C20_TesA_AcACP_un - P.k7_1f(9).*c.TesA.*c.C20_AcACP_un + P.k3_4r(9).*c.C20_FabH_AcACP_un - P.k3_4f(9).*c.FabH.*c.C20_AcACP_un + P.k3_5r(9).*c.C22_FabH_Act_AcACP_un - P.k3_5f(9).*c.C2_FabH_Act.*c.C20_AcACP_un;

% Fatty Acids (TesA)
d_C4_FA = P.kcat7(1).*c.C4_TesA_AcACP;
d_C6_FA = P.kcat7(2).*c.C6_TesA_AcACP;
d_C8_FA = P.kcat7(3).*c.C8_TesA_AcACP;
d_C10_FA = P.kcat7(4).*c.C10_TesA_AcACP;
d_C12_FA = P.kcat7(5).*c.C12_TesA_AcACP;
d_C14_FA = P.kcat7(6).*c.C14_TesA_AcACP;
d_C16_FA = P.kcat7(7).*c.C16_TesA_AcACP;
d_C18_FA = P.kcat7(8).*c.C18_TesA_AcACP;
d_C20_FA = P.kcat7(9).*c.C20_TesA_AcACP;

% Fatty Acids (unsaturated) (TesA)
d_C12_FA_un = P.kcat7(5).*c.C12_TesA_AcACP_un;
d_C14_FA_un = P.kcat7(6).*c.C14_TesA_AcACP_un;
d_C16_FA_un = P.kcat7(7).*c.C16_TesA_AcACP_un;
d_C18_FA_un = P.kcat7(8).*c.C18_TesA_AcACP_un;
d_C20_FA_un = P.kcat7(9).*c.C20_TesA_AcACP_un;

% BC-ATP % changed for ACC
d_BC_ATP = P.kcat1_1.*c.ACC_C.*c.ATP./(P.Km1_1 + c.ATP) - P.kcat1_2.*c.BC_ATP.*c.C1_Bicarbonate./(P.Km1_2 + c.C1_Bicarbonate);

% BC-ATP-HCO3 % changed for ACC
d_C1_BC_ATP_HCO3 = 0*c.ATP;

% BC-Pi-HCO3 % changed for ACC
d_C1_BC_Pi_HCO3 = P.kcat1_2.*c.BC_ATP.*c.C1_Bicarbonate./(P.Km1_2 + c.C1_Bicarbonate) - P.kcat1_3.*c.ACC_B.*c.C1_BC_Pi_HCO3./(P.Km1_3 + c.C1_BC_Pi_HCO3);

% BC-Pi-HCO3-BCCP-Biotin % changed for ACC
d_C1_BC_Pi_HCO3_BCCP_Biotin = 0*c.ATP;

% BCCP-Biotin-CO2 % changed for ACC
d_C1_BCCP_Biotin_CO2 = P.kcat1_3.*c.ACC_B.*c.C1_BC_Pi_HCO3./(P.Km1_3 + c.C1_BC_Pi_HCO3) - P.kcat1_4.*c.ACC_AD.*c.C1_BCCP_Biotin_CO2./(P.Km1_4 + c.C1_BCCP_Biotin_CO2);

% CT-BCCP-Biotin-CO2 % changed for ACC
d_C1_CT_BCCP_Biotin_CO2 = 0*c.ATP;

% CT* % changed for ACC
d_C1_CT_Act = P.kcat1_4.*c.ACC_AD.*c.C1_BCCP_Biotin_CO2./(P.Km1_4 + c.C1_BCCP_Biotin_CO2) - P.kcat1_5.*c.C1_CT_Act.*c.C2_AcCoA./(P.Km1_5 + c.C2_AcCoA);

% CT*-AcCoA % changed for ACC
d_C3_CT_Act_AcCoA = 0*c.ATP;

% FabD-Malonyl-CoA
d_C3_FabD_MalCoA = P.k2_1f.*c.FabD.*c.C3_MalCoA - P.k2_1r.*c.C3_FabD_MalCoA + P.k2_2r.*c.C3_FabD_Act.*c.CoA - P.k2_2f.*c.C3_FabD_MalCoA; 

% FabD*
d_C3_FabD_Act = P.k2_2f.*c.C3_FabD_MalCoA - P.k2_2r.*c.C3_FabD_Act.*c.CoA + P.k2_3r.*c.C3_FabD_Act_ACP - P.k2_3f.*c.C3_FabD_Act.*c.ACP;

% FabD*-ACP
d_C3_FabD_Act_ACP = P.k2_3f.*c.C3_FabD_Act.*c.ACP - P.k2_3r.*c.C3_FabD_Act_ACP + P.k2_4r.*c.FabD.*c.C3_MalACP - P.k2_4f.*c.C3_FabD_Act_ACP;

% C2n (n=1:9) FabH-CoA
d_C2_FabH_CoA = P.k3_1f(1).*c.FabH.*c.C2_AcCoA - P.k3_1r(1).*c.C2_FabH_CoA + P.k3_2r(1).*c.C2_FabH_Act.*c.CoA - P.k3_2f(1).*c.C2_FabH_CoA; 
d_C4_FabH_CoA = P.k3_1f(2).*c.FabH.*c.C4_SucCoA - P.k3_1r(2).*c.C4_FabH_CoA + P.k3_2r(2).*c.C4_FabH_Act.*c.CoA - P.k3_2f(2).*c.C4_FabH_CoA; 
d_C6_FabH_CoA = P.k3_1f(3).*c.FabH.*c.C6_HexCoA - P.k3_1r(3).*c.C6_FabH_CoA + P.k3_2r(3).*c.C6_FabH_Act.*c.CoA - P.k3_2f(3).*c.C6_FabH_CoA; 
d_C8_FabH_CoA = P.k3_1f(4).*c.FabH.*c.C8_OcCoA - P.k3_1r(4).*c.C8_FabH_CoA + P.k3_2r(4).*c.C8_FabH_Act.*c.CoA - P.k3_2f(4).*c.C8_FabH_CoA; 
d_C10_FabH_CoA = P.k3_1f(5).*c.FabH.*c.C10_DecCoA - P.k3_1r(5).*c.C10_FabH_CoA + P.k3_2r(5).*c.C10_FabH_Act.*c.CoA - P.k3_2f(5).*c.C10_FabH_CoA; 
d_C12_FabH_CoA = P.k3_1f(6).*c.FabH.*c.C12_LauCoA - P.k3_1r(6).*c.C12_FabH_CoA + P.k3_2r(6).*c.C12_FabH_Act.*c.CoA - P.k3_2f(6).*c.C12_FabH_CoA;
d_C14_FabH_CoA = P.k3_1f(7).*c.FabH.*c.C14_EthCoA - P.k3_1r(7).*c.C14_FabH_CoA + P.k3_2r(7).*c.C14_FabH_Act.*c.CoA - P.k3_2f(7).*c.C14_FabH_CoA; 
d_C16_FabH_CoA = P.k3_1f(8).*c.FabH.*c.C16_PalCoA - P.k3_1r(8).*c.C16_FabH_CoA + P.k3_2r(8).*c.C16_FabH_Act.*c.CoA - P.k3_2f(8).*c.C16_FabH_CoA; 
d_C18_FabH_CoA = P.k3_1f(9).*c.FabH.*c.C18_OcDecCoA - P.k3_1r(9).*c.C18_FabH_CoA + P.k3_2r(9).*c.C18_FabH_Act.*c.CoA - P.k3_2f(9).*c.C18_FabH_CoA; 

% C2n (n=1:9) FabH*
% making FabH* - using FabH* - inhibition from Acyl ACPs (only Acetyl-CoA derived FabH*)
d_C2_FabH_Act = P.k3_2f(1).*c.C2_FabH_CoA - P.k3_2r(1).*c.C2_FabH_Act.*c.CoA... 
 + P.k3_3r(1).*c.C5_FabH_Act_MalACP - P.k3_3f(1).*c.C2_FabH_Act.*c.C3_MalACP...
 + P.k3_8f.*c.C2_FabH_AcACP - P.k3_8r.*c.C2_FabH_Act.*c.ACP...
 + P.k3_5r(1).*c.C6_FabH_Act_AcACP - P.k3_5f(1).*c.C2_FabH_Act.*c.C4_AcACP...
 + P.k3_5r(2).*c.C8_FabH_Act_AcACP - P.k3_5f(2).*c.C2_FabH_Act.*c.C6_AcACP...
 + P.k3_5r(3).*c.C10_FabH_Act_AcACP - P.k3_5f(3).*c.C2_FabH_Act.*c.C8_AcACP...
 + P.k3_5r(4).*c.C12_FabH_Act_AcACP - P.k3_5f(4).*c.C2_FabH_Act.*c.C10_AcACP...
 + P.k3_5r(5).*c.C14_FabH_Act_AcACP - P.k3_5f(5).*c.C2_FabH_Act.*c.C12_AcACP...
 + P.k3_5r(6).*c.C16_FabH_Act_AcACP - P.k3_5f(6).*c.C2_FabH_Act.*c.C14_AcACP...
 + P.k3_5r(7).*c.C18_FabH_Act_AcACP - P.k3_5f(7).*c.C2_FabH_Act.*c.C16_AcACP...
 + P.k3_5r(8).*c.C20_FabH_Act_AcACP - P.k3_5f(8).*c.C2_FabH_Act.*c.C18_AcACP...
 + P.k3_5r(9).*c.C22_FabH_Act_AcACP - P.k3_5f(9).*c.C2_FabH_Act.*c.C20_AcACP...
 + P.k3_5r(5).*c.C14_FabH_Act_AcACP_un - P.k3_5f(5).*c.C2_FabH_Act.*c.C12_AcACP_un...
 + P.k3_5r(6).*c.C16_FabH_Act_AcACP_un - P.k3_5f(6).*c.C2_FabH_Act.*c.C14_AcACP_un...
 + P.k3_5r(7).*c.C18_FabH_Act_AcACP_un - P.k3_5f(7).*c.C2_FabH_Act.*c.C16_AcACP_un...
 + P.k3_5r(8).*c.C20_FabH_Act_AcACP_un - P.k3_5f(8).*c.C2_FabH_Act.*c.C18_AcACP_un...
 + P.k3_5r(9).*c.C22_FabH_Act_AcACP_un - P.k3_5f(9).*c.C2_FabH_Act.*c.C20_AcACP_un;
d_C4_FabH_Act = P.k3_2f(2).*c.C4_FabH_CoA - P.k3_2r(2).*c.C4_FabH_Act.*c.CoA + P.k3_3r(2).*c.C7_FabH_Act_MalACP - P.k3_3f(2).*c.C4_FabH_Act.*c.C3_MalACP;
d_C6_FabH_Act = P.k3_2f(3).*c.C6_FabH_CoA - P.k3_2r(3).*c.C6_FabH_Act.*c.CoA + P.k3_3r(3).*c.C9_FabH_Act_MalACP - P.k3_3f(3).*c.C6_FabH_Act.*c.C3_MalACP;
d_C8_FabH_Act = P.k3_2f(4).*c.C8_FabH_CoA - P.k3_2r(4).*c.C8_FabH_Act.*c.CoA + P.k3_3r(4).*c.C11_FabH_Act_MalACP - P.k3_3f(4).*c.C8_FabH_Act.*c.C3_MalACP;
d_C10_FabH_Act = P.k3_2f(5).*c.C10_FabH_CoA - P.k3_2r(5).*c.C10_FabH_Act.*c.CoA + P.k3_3r(5).*c.C13_FabH_Act_MalACP - P.k3_3f(5).*c.C10_FabH_Act.*c.C3_MalACP;
d_C12_FabH_Act = P.k3_2f(6).*c.C12_FabH_CoA - P.k3_2r(6).*c.C12_FabH_Act.*c.CoA + P.k3_3r(6).*c.C15_FabH_Act_MalACP - P.k3_3f(6).*c.C12_FabH_Act.*c.C3_MalACP;
d_C14_FabH_Act = P.k3_2f(7).*c.C14_FabH_CoA - P.k3_2r(7).*c.C14_FabH_Act.*c.CoA + P.k3_3r(7).*c.C17_FabH_Act_MalACP - P.k3_3f(7).*c.C14_FabH_Act.*c.C3_MalACP;
d_C16_FabH_Act = P.k3_2f(8).*c.C16_FabH_CoA - P.k3_2r(8).*c.C16_FabH_Act.*c.CoA + P.k3_3r(8).*c.C19_FabH_Act_MalACP - P.k3_3f(8).*c.C16_FabH_Act.*c.C3_MalACP;
d_C18_FabH_Act = P.k3_2f(9).*c.C18_FabH_CoA - P.k3_2r(9).*c.C18_FabH_Act.*c.CoA + P.k3_3r(9).*c.C21_FabH_Act_MalACP - P.k3_3f(9).*c.C18_FabH_Act.*c.C3_MalACP;

% C2n (n=1:9) FabH*-Malonyl-ACP
d_C5_FabH_Act_MalACP = P.k3_3f(1).*c.C2_FabH_Act.*c.C3_MalACP - P.k3_3r(1).*c.C5_FabH_Act_MalACP - P.kcat3(1).*c.C5_FabH_Act_MalACP; 
d_C7_FabH_Act_MalACP = P.k3_3f(2).*c.C4_FabH_Act.*c.C3_MalACP - P.k3_3r(2).*c.C7_FabH_Act_MalACP - P.kcat3(2).*c.C7_FabH_Act_MalACP; 
d_C9_FabH_Act_MalACP = P.k3_3f(3).*c.C6_FabH_Act.*c.C3_MalACP - P.k3_3r(3).*c.C9_FabH_Act_MalACP - P.kcat3(3).*c.C9_FabH_Act_MalACP; 
d_C11_FabH_Act_MalACP = P.k3_3f(4).*c.C8_FabH_Act.*c.C3_MalACP - P.k3_3r(4).*c.C11_FabH_Act_MalACP - P.kcat3(4).*c.C11_FabH_Act_MalACP; 
d_C13_FabH_Act_MalACP = P.k3_3f(5).*c.C10_FabH_Act.*c.C3_MalACP - P.k3_3r(5).*c.C13_FabH_Act_MalACP - P.kcat3(5).*c.C13_FabH_Act_MalACP; 
d_C15_FabH_Act_MalACP = P.k3_3f(6).*c.C12_FabH_Act.*c.C3_MalACP - P.k3_3r(6).*c.C15_FabH_Act_MalACP - P.kcat3(6).*c.C15_FabH_Act_MalACP; 
d_C17_FabH_Act_MalACP = P.k3_3f(7).*c.C14_FabH_Act.*c.C3_MalACP - P.k3_3r(7).*c.C17_FabH_Act_MalACP - P.kcat3(7).*c.C17_FabH_Act_MalACP; 
d_C19_FabH_Act_MalACP = P.k3_3f(8).*c.C16_FabH_Act.*c.C3_MalACP - P.k3_3r(8).*c.C19_FabH_Act_MalACP - P.kcat3(8).*c.C19_FabH_Act_MalACP; 
d_C21_FabH_Act_MalACP = P.k3_3f(9).*c.C18_FabH_Act.*c.C3_MalACP - P.k3_3r(9).*c.C21_FabH_Act_MalACP - P.kcat3(9).*c.C21_FabH_Act_MalACP; 

% FabG-NADPH
d_FabG_NADPH = P.k4_1f(1).*c.FabG.*c.NADPH - P.k4_1r(1).*c.FabG_NADPH...
 + P.k4_2r(1).*c.C4_FabG_NADPH_BKeACP - P.k4_2f(1).*c.FabG_NADPH.*c.C4_BKeACP...
 + P.k4_2r(2).*c.C6_FabG_NADPH_BKeACP - P.k4_2f(2).*c.FabG_NADPH.*c.C6_BKeACP...
 + P.k4_2r(3).*c.C8_FabG_NADPH_BKeACP - P.k4_2f(3).*c.FabG_NADPH.*c.C8_BKeACP...
 + P.k4_2r(4).*c.C10_FabG_NADPH_BKeACP - P.k4_2f(4).*c.FabG_NADPH.*c.C10_BKeACP...
 + P.k4_2r(5).*c.C12_FabG_NADPH_BKeACP - P.k4_2f(5).*c.FabG_NADPH.*c.C12_BKeACP...
 + P.k4_2r(6).*c.C14_FabG_NADPH_BKeACP - P.k4_2f(6).*c.FabG_NADPH.*c.C14_BKeACP...
 + P.k4_2r(7).*c.C16_FabG_NADPH_BKeACP - P.k4_2f(7).*c.FabG_NADPH.*c.C16_BKeACP...
 + P.k4_2r(8).*c.C18_FabG_NADPH_BKeACP - P.k4_2f(8).*c.FabG_NADPH.*c.C18_BKeACP...
 + P.k4_2r(9).*c.C20_FabG_NADPH_BKeACP - P.k4_2f(9).*c.FabG_NADPH.*c.C20_BKeACP...
 + P.k4_2r(5).*c.C12_FabG_NADPH_BKeACP_un - P.k4_2f(5).*c.FabG_NADPH.*c.C12_BKeACP_un...
 + P.k4_2r(6).*c.C14_FabG_NADPH_BKeACP_un - P.k4_2f(6).*c.FabG_NADPH.*c.C14_BKeACP_un...
 + P.k4_2r(7).*c.C16_FabG_NADPH_BKeACP_un - P.k4_2f(7).*c.FabG_NADPH.*c.C16_BKeACP_un...
 + P.k4_2r(8).*c.C18_FabG_NADPH_BKeACP_un - P.k4_2f(8).*c.FabG_NADPH.*c.C18_BKeACP_un...
 + P.k4_2r(9).*c.C20_FabG_NADPH_BKeACP_un - P.k4_2f(9).*c.FabG_NADPH.*c.C20_BKeACP_un; 

% C2n (n=2:10) FabG-NADPH-B-ketoacyl-ACPs
d_C4_FabG_NADPH_BKeACP = P.k4_2f(1).*c.FabG_NADPH.*c.C4_BKeACP - P.k4_2r(1).*c.C4_FabG_NADPH_BKeACP - P.kcat4(1).*c.C4_FabG_NADPH_BKeACP;
d_C6_FabG_NADPH_BKeACP = P.k4_2f(2).*c.FabG_NADPH.*c.C6_BKeACP - P.k4_2r(2).*c.C6_FabG_NADPH_BKeACP - P.kcat4(2).*c.C6_FabG_NADPH_BKeACP;
d_C8_FabG_NADPH_BKeACP = P.k4_2f(3).*c.FabG_NADPH.*c.C8_BKeACP - P.k4_2r(3).*c.C8_FabG_NADPH_BKeACP - P.kcat4(3).*c.C8_FabG_NADPH_BKeACP;
d_C10_FabG_NADPH_BKeACP = P.k4_2f(4).*c.FabG_NADPH.*c.C10_BKeACP - P.k4_2r(4).*c.C10_FabG_NADPH_BKeACP - P.kcat4(4).*c.C10_FabG_NADPH_BKeACP;
d_C12_FabG_NADPH_BKeACP = P.k4_2f(5).*c.FabG_NADPH.*c.C12_BKeACP - P.k4_2r(5).*c.C12_FabG_NADPH_BKeACP - P.kcat4(5).*c.C12_FabG_NADPH_BKeACP;
d_C14_FabG_NADPH_BKeACP = P.k4_2f(6).*c.FabG_NADPH.*c.C14_BKeACP - P.k4_2r(6).*c.C14_FabG_NADPH_BKeACP - P.kcat4(6).*c.C14_FabG_NADPH_BKeACP;
d_C16_FabG_NADPH_BKeACP = P.k4_2f(7).*c.FabG_NADPH.*c.C16_BKeACP - P.k4_2r(7).*c.C16_FabG_NADPH_BKeACP - P.kcat4(7).*c.C16_FabG_NADPH_BKeACP;
d_C18_FabG_NADPH_BKeACP = P.k4_2f(8).*c.FabG_NADPH.*c.C18_BKeACP - P.k4_2r(8).*c.C18_FabG_NADPH_BKeACP - P.kcat4(8).*c.C18_FabG_NADPH_BKeACP;
d_C20_FabG_NADPH_BKeACP = P.k4_2f(9).*c.FabG_NADPH.*c.C20_BKeACP - P.k4_2r(9).*c.C20_FabG_NADPH_BKeACP - P.kcat4(9).*c.C20_FabG_NADPH_BKeACP;

% C2n:1 (n=6:10) FabG-NADPH-B-ketoacyl-ACPs
d_C12_FabG_NADPH_BKeACP_un = P.k4_2f(5).*c.FabG_NADPH.*c.C12_BKeACP_un - P.k4_2r(5).*c.C12_FabG_NADPH_BKeACP_un - P.kcat4(5).*c.C12_FabG_NADPH_BKeACP_un;
d_C14_FabG_NADPH_BKeACP_un = P.k4_2f(6).*c.FabG_NADPH.*c.C14_BKeACP_un - P.k4_2r(6).*c.C14_FabG_NADPH_BKeACP_un - P.kcat4(6).*c.C14_FabG_NADPH_BKeACP_un;
d_C16_FabG_NADPH_BKeACP_un = P.k4_2f(7).*c.FabG_NADPH.*c.C16_BKeACP_un - P.k4_2r(7).*c.C16_FabG_NADPH_BKeACP_un - P.kcat4(7).*c.C16_FabG_NADPH_BKeACP_un;
d_C18_FabG_NADPH_BKeACP_un = P.k4_2f(8).*c.FabG_NADPH.*c.C18_BKeACP_un - P.k4_2r(8).*c.C18_FabG_NADPH_BKeACP_un - P.kcat4(8).*c.C18_FabG_NADPH_BKeACP_un;
d_C20_FabG_NADPH_BKeACP_un = P.k4_2f(9).*c.FabG_NADPH.*c.C20_BKeACP_un - P.k4_2r(9).*c.C20_FabG_NADPH_BKeACP_un - P.kcat4(9).*c.C20_FabG_NADPH_BKeACP_un;

% C2n (n=2:10) FabZ-B-hydroxy-acyl-ACPs
d_C4_FabZ_BHyAcACP = P.k5_1f(1).*c.FabZ.*c.C4_BHyAcACP - P.k5_1r(1).*c.C4_FabZ_BHyAcACP + P.k5_2r(1).*c.C4_FabZ_EnAcACP - P.kcat5(1).*c.C4_FabZ_BHyAcACP;
d_C6_FabZ_BHyAcACP = P.k5_1f(2).*c.FabZ.*c.C6_BHyAcACP - P.k5_1r(2).*c.C6_FabZ_BHyAcACP + P.k5_2r(2).*c.C6_FabZ_EnAcACP - P.kcat5(2).*c.C6_FabZ_BHyAcACP;
d_C8_FabZ_BHyAcACP = P.k5_1f(3).*c.FabZ.*c.C8_BHyAcACP - P.k5_1r(3).*c.C8_FabZ_BHyAcACP + P.k5_2r(3).*c.C8_FabZ_EnAcACP - P.kcat5(3).*c.C8_FabZ_BHyAcACP;
d_C10_FabZ_BHyAcACP = P.k5_1f(4).*c.FabZ.*c.C10_BHyAcACP - P.k5_1r(4).*c.C10_FabZ_BHyAcACP + P.k5_2r(4).*c.C10_FabZ_EnAcACP - P.kcat5(4).*c.C10_FabZ_BHyAcACP;
d_C12_FabZ_BHyAcACP = P.k5_1f(5).*c.FabZ.*c.C12_BHyAcACP - P.k5_1r(5).*c.C12_FabZ_BHyAcACP + P.k5_2r(5).*c.C12_FabZ_EnAcACP - P.kcat5(5).*c.C12_FabZ_BHyAcACP;
d_C14_FabZ_BHyAcACP = P.k5_1f(6).*c.FabZ.*c.C14_BHyAcACP - P.k5_1r(6).*c.C14_FabZ_BHyAcACP + P.k5_2r(6).*c.C14_FabZ_EnAcACP - P.kcat5(6).*c.C14_FabZ_BHyAcACP;
d_C16_FabZ_BHyAcACP = P.k5_1f(7).*c.FabZ.*c.C16_BHyAcACP - P.k5_1r(7).*c.C16_FabZ_BHyAcACP + P.k5_2r(7).*c.C16_FabZ_EnAcACP - P.kcat5(7).*c.C16_FabZ_BHyAcACP;
d_C18_FabZ_BHyAcACP = P.k5_1f(8).*c.FabZ.*c.C18_BHyAcACP - P.k5_1r(8).*c.C18_FabZ_BHyAcACP + P.k5_2r(8).*c.C18_FabZ_EnAcACP - P.kcat5(8).*c.C18_FabZ_BHyAcACP;
d_C20_FabZ_BHyAcACP = P.k5_1f(9).*c.FabZ.*c.C20_BHyAcACP - P.k5_1r(9).*c.C20_FabZ_BHyAcACP + P.k5_2r(9).*c.C20_FabZ_EnAcACP - P.kcat5(9).*c.C20_FabZ_BHyAcACP;

% C2n:1 (n=6:10) FabZ-B-hydroxy-acyl-ACPs
d_C12_FabZ_BHyAcACP_un = P.k5_1f(5).*c.FabZ.*c.C12_BHyAcACP_un - P.k5_1r(5).*c.C12_FabZ_BHyAcACP_un + P.k5_2r(5).*c.C12_FabZ_EnAcACP_un - P.kcat5(5).*c.C12_FabZ_BHyAcACP_un;
d_C14_FabZ_BHyAcACP_un = P.k5_1f(6).*c.FabZ.*c.C14_BHyAcACP_un - P.k5_1r(6).*c.C14_FabZ_BHyAcACP_un + P.k5_2r(6).*c.C14_FabZ_EnAcACP_un - P.kcat5(6).*c.C14_FabZ_BHyAcACP_un;
d_C16_FabZ_BHyAcACP_un = P.k5_1f(7).*c.FabZ.*c.C16_BHyAcACP_un - P.k5_1r(7).*c.C16_FabZ_BHyAcACP_un + P.k5_2r(7).*c.C16_FabZ_EnAcACP_un - P.kcat5(7).*c.C16_FabZ_BHyAcACP_un;
d_C18_FabZ_BHyAcACP_un = P.k5_1f(8).*c.FabZ.*c.C18_BHyAcACP_un - P.k5_1r(8).*c.C18_FabZ_BHyAcACP_un + P.k5_2r(8).*c.C18_FabZ_EnAcACP_un - P.kcat5(8).*c.C18_FabZ_BHyAcACP_un;
d_C20_FabZ_BHyAcACP_un = P.k5_1f(9).*c.FabZ.*c.C20_BHyAcACP_un - P.k5_1r(9).*c.C20_FabZ_BHyAcACP_un + P.k5_2r(9).*c.C20_FabZ_EnAcACP_un - P.kcat5(9).*c.C20_FabZ_BHyAcACP_un;

% C2n (n=2:10) FabZ-Enoyl-Acyl-ACPs
d_C4_FabZ_EnAcACP = P.kcat5(1).*c.C4_FabZ_BHyAcACP - P.k5_2r(1).*c.C4_FabZ_EnAcACP + P.k5_3r(1).*c.FabZ.*c.C4_EnAcACP - P.k5_3f(1).*c.C4_FabZ_EnAcACP;
d_C6_FabZ_EnAcACP = P.kcat5(2).*c.C6_FabZ_BHyAcACP - P.k5_2r(2).*c.C6_FabZ_EnAcACP + P.k5_3r(2).*c.FabZ.*c.C6_EnAcACP - P.k5_3f(2).*c.C6_FabZ_EnAcACP;
d_C8_FabZ_EnAcACP = P.kcat5(3).*c.C8_FabZ_BHyAcACP - P.k5_2r(3).*c.C8_FabZ_EnAcACP + P.k5_3r(3).*c.FabZ.*c.C8_EnAcACP - P.k5_3f(3).*c.C8_FabZ_EnAcACP;
d_C10_FabZ_EnAcACP = P.kcat5(4).*c.C10_FabZ_BHyAcACP - P.k5_2r(4).*c.C10_FabZ_EnAcACP + P.k5_3r(4).*c.FabZ.*c.C10_EnAcACP - P.k5_3f(4).*c.C10_FabZ_EnAcACP;
d_C12_FabZ_EnAcACP = P.kcat5(5).*c.C12_FabZ_BHyAcACP - P.k5_2r(5).*c.C12_FabZ_EnAcACP + P.k5_3r(5).*c.FabZ.*c.C12_EnAcACP - P.k5_3f(5).*c.C12_FabZ_EnAcACP;
d_C14_FabZ_EnAcACP = P.kcat5(6).*c.C14_FabZ_BHyAcACP - P.k5_2r(6).*c.C14_FabZ_EnAcACP + P.k5_3r(6).*c.FabZ.*c.C14_EnAcACP - P.k5_3f(6).*c.C14_FabZ_EnAcACP;
d_C16_FabZ_EnAcACP = P.kcat5(7).*c.C16_FabZ_BHyAcACP - P.k5_2r(7).*c.C16_FabZ_EnAcACP + P.k5_3r(7).*c.FabZ.*c.C16_EnAcACP - P.k5_3f(7).*c.C16_FabZ_EnAcACP;
d_C18_FabZ_EnAcACP = P.kcat5(8).*c.C18_FabZ_BHyAcACP - P.k5_2r(8).*c.C18_FabZ_EnAcACP + P.k5_3r(8).*c.FabZ.*c.C18_EnAcACP - P.k5_3f(8).*c.C18_FabZ_EnAcACP;
d_C20_FabZ_EnAcACP = P.kcat5(9).*c.C20_FabZ_BHyAcACP - P.k5_2r(9).*c.C20_FabZ_EnAcACP + P.k5_3r(9).*c.FabZ.*c.C20_EnAcACP - P.k5_3f(9).*c.C20_FabZ_EnAcACP;

% C2n:1 (n=6:10) FabZ-Enoyl-Acyl-ACPs
d_C12_FabZ_EnAcACP_un = P.kcat5(5).*c.C12_FabZ_BHyAcACP_un - P.k5_2r(5).*c.C12_FabZ_EnAcACP_un + P.k5_3r(5).*c.FabZ.*c.C12_EnAcACP_un - P.k5_3f(5).*c.C12_FabZ_EnAcACP_un;
d_C14_FabZ_EnAcACP_un = P.kcat5(6).*c.C14_FabZ_BHyAcACP_un - P.k5_2r(6).*c.C14_FabZ_EnAcACP_un + P.k5_3r(6).*c.FabZ.*c.C14_EnAcACP_un - P.k5_3f(6).*c.C14_FabZ_EnAcACP_un;
d_C16_FabZ_EnAcACP_un = P.kcat5(7).*c.C16_FabZ_BHyAcACP_un - P.k5_2r(7).*c.C16_FabZ_EnAcACP_un + P.k5_3r(7).*c.FabZ.*c.C16_EnAcACP_un - P.k5_3f(7).*c.C16_FabZ_EnAcACP_un;
d_C18_FabZ_EnAcACP_un = P.kcat5(8).*c.C18_FabZ_BHyAcACP_un - P.k5_2r(8).*c.C18_FabZ_EnAcACP_un + P.k5_3r(8).*c.FabZ.*c.C18_EnAcACP_un - P.k5_3f(8).*c.C18_FabZ_EnAcACP_un;
d_C20_FabZ_EnAcACP_un = P.kcat5(9).*c.C20_FabZ_BHyAcACP_un - P.k5_2r(9).*c.C20_FabZ_EnAcACP_un + P.k5_3r(9).*c.FabZ.*c.C20_EnAcACP_un - P.k5_3f(9).*c.C20_FabZ_EnAcACP_un;

% C2n (n=2:10) FabA-B-hydroxy-acyl-ACPs
d_C4_FabA_BHyAcACP = P.k9_1f(1).*c.FabA.*c.C4_BHyAcACP - P.k9_1r(1).*c.C4_FabA_BHyAcACP + P.k9_2r(1).*c.C4_FabA_EnAcACP - P.kcat9(1).*c.C4_FabA_BHyAcACP;
d_C6_FabA_BHyAcACP = P.k9_1f(2).*c.FabA.*c.C6_BHyAcACP - P.k9_1r(2).*c.C6_FabA_BHyAcACP + P.k9_2r(2).*c.C6_FabA_EnAcACP - P.kcat9(2).*c.C6_FabA_BHyAcACP;
d_C8_FabA_BHyAcACP = P.k9_1f(3).*c.FabA.*c.C8_BHyAcACP - P.k9_1r(3).*c.C8_FabA_BHyAcACP + P.k9_2r(3).*c.C8_FabA_EnAcACP - P.kcat9(3).*c.C8_FabA_BHyAcACP;
d_C10_FabA_BHyAcACP = P.k9_1f(4).*c.FabA.*c.C10_BHyAcACP - P.k9_1r(4).*c.C10_FabA_BHyAcACP + P.k9_2r(4).*c.C10_FabA_EnAcACP - P.kcat9(4).*c.C10_FabA_BHyAcACP;
d_C12_FabA_BHyAcACP = P.k9_1f(5).*c.FabA.*c.C12_BHyAcACP - P.k9_1r(5).*c.C12_FabA_BHyAcACP + P.k9_2r(5).*c.C12_FabA_EnAcACP - P.kcat9(5).*c.C12_FabA_BHyAcACP;
d_C14_FabA_BHyAcACP = P.k9_1f(6).*c.FabA.*c.C14_BHyAcACP - P.k9_1r(6).*c.C14_FabA_BHyAcACP + P.k9_2r(6).*c.C14_FabA_EnAcACP - P.kcat9(6).*c.C14_FabA_BHyAcACP;
d_C16_FabA_BHyAcACP = P.k9_1f(7).*c.FabA.*c.C16_BHyAcACP - P.k9_1r(7).*c.C16_FabA_BHyAcACP + P.k9_2r(7).*c.C16_FabA_EnAcACP - P.kcat9(7).*c.C16_FabA_BHyAcACP;
d_C18_FabA_BHyAcACP = P.k9_1f(8).*c.FabA.*c.C18_BHyAcACP - P.k9_1r(8).*c.C18_FabA_BHyAcACP + P.k9_2r(8).*c.C18_FabA_EnAcACP - P.kcat9(8).*c.C18_FabA_BHyAcACP;
d_C20_FabA_BHyAcACP = P.k9_1f(9).*c.FabA.*c.C20_BHyAcACP - P.k9_1r(9).*c.C20_FabA_BHyAcACP + P.k9_2r(9).*c.C20_FabA_EnAcACP - P.kcat9(9).*c.C20_FabA_BHyAcACP;

% C2n:1 (n=6:10) FabA-B-hydroxy-acyl-ACPs
d_C12_FabA_BHyAcACP_un = P.k9_1f_un(5).*c.FabA.*c.C12_BHyAcACP_un - P.k9_1r_un(5).*c.C12_FabA_BHyAcACP_un + P.k9_2r(5).*c.C12_FabA_EnAcACP_un - P.kcat9(5).*c.C12_FabA_BHyAcACP_un;
d_C14_FabA_BHyAcACP_un = P.k9_1f_un(6).*c.FabA.*c.C14_BHyAcACP_un - P.k9_1r_un(6).*c.C14_FabA_BHyAcACP_un + P.k9_2r(6).*c.C14_FabA_EnAcACP_un - P.kcat9(6).*c.C14_FabA_BHyAcACP_un;
d_C16_FabA_BHyAcACP_un = P.k9_1f_un(7).*c.FabA.*c.C16_BHyAcACP_un - P.k9_1r_un(7).*c.C16_FabA_BHyAcACP_un + P.k9_2r(7).*c.C16_FabA_EnAcACP_un - P.kcat9(7).*c.C16_FabA_BHyAcACP_un;
d_C18_FabA_BHyAcACP_un = P.k9_1f_un(8).*c.FabA.*c.C18_BHyAcACP_un - P.k9_1r_un(8).*c.C18_FabA_BHyAcACP_un + P.k9_2r(8).*c.C18_FabA_EnAcACP_un - P.kcat9(8).*c.C18_FabA_BHyAcACP_un;
d_C20_FabA_BHyAcACP_un = P.k9_1f_un(9).*c.FabA.*c.C20_BHyAcACP_un - P.k9_1r_un(9).*c.C20_FabA_BHyAcACP_un + P.k9_2r(9).*c.C20_FabA_EnAcACP_un - P.kcat9(9).*c.C20_FabA_BHyAcACP_un;

% C2n (n=2:10) FabA-Enoyl-Acyl-ACPs
d_C4_FabA_EnAcACP = P.kcat9(1).*c.C4_FabA_BHyAcACP - P.k9_2r(1).*c.C4_FabA_EnAcACP + P.k9_3r(1).*c.FabA.*c.C4_EnAcACP - P.k9_3f(1).*c.C4_FabA_EnAcACP;
d_C6_FabA_EnAcACP = P.kcat9(2).*c.C6_FabA_BHyAcACP - P.k9_2r(2).*c.C6_FabA_EnAcACP + P.k9_3r(2).*c.FabA.*c.C6_EnAcACP - P.k9_3f(2).*c.C6_FabA_EnAcACP;
d_C8_FabA_EnAcACP = P.kcat9(3).*c.C8_FabA_BHyAcACP - P.k9_2r(3).*c.C8_FabA_EnAcACP + P.k9_3r(3).*c.FabA.*c.C8_EnAcACP - P.k9_3f(3).*c.C8_FabA_EnAcACP;
d_C10_FabA_EnAcACP = P.kcat9(4).*c.C10_FabA_BHyAcACP - P.k9_2r(4).*c.C10_FabA_EnAcACP + P.k9_3r(4).*c.FabA.*c.C10_EnAcACP - P.k9_3f(4).*c.C10_FabA_EnAcACP + P.k9_2r_un(4).*c.C10_FabA_cis3EnAcACP - P.kcat9_un(4).*c.C10_FabA_EnAcACP;
d_C12_FabA_EnAcACP = P.kcat9(5).*c.C12_FabA_BHyAcACP - P.k9_2r(5).*c.C12_FabA_EnAcACP + P.k9_3r(5).*c.FabA.*c.C12_EnAcACP - P.k9_3f(5).*c.C12_FabA_EnAcACP;
d_C14_FabA_EnAcACP = P.kcat9(6).*c.C14_FabA_BHyAcACP - P.k9_2r(6).*c.C14_FabA_EnAcACP + P.k9_3r(6).*c.FabA.*c.C14_EnAcACP - P.k9_3f(6).*c.C14_FabA_EnAcACP;
d_C16_FabA_EnAcACP = P.kcat9(7).*c.C16_FabA_BHyAcACP - P.k9_2r(7).*c.C16_FabA_EnAcACP + P.k9_3r(7).*c.FabA.*c.C16_EnAcACP - P.k9_3f(7).*c.C16_FabA_EnAcACP;
d_C18_FabA_EnAcACP = P.kcat9(8).*c.C18_FabA_BHyAcACP - P.k9_2r(8).*c.C18_FabA_EnAcACP + P.k9_3r(8).*c.FabA.*c.C18_EnAcACP - P.k9_3f(8).*c.C18_FabA_EnAcACP;
d_C20_FabA_EnAcACP = P.kcat9(9).*c.C20_FabA_BHyAcACP - P.k9_2r(9).*c.C20_FabA_EnAcACP + P.k9_3r(9).*c.FabA.*c.C20_EnAcACP - P.k9_3f(9).*c.C20_FabA_EnAcACP;

% FabA-C10 cis-3-Enoyl-Acyl-ACP
d_C10_FabA_cis3EnAcACP = P.kcat9_un(4).*c.C10_FabA_EnAcACP - P.k9_2r_un(4).*c.C10_FabA_cis3EnAcACP + P.k9_3r_un(4).*c.FabA.*c.C10_cis3EnAcACP - P.k9_3f_un(4).*c.C10_FabA_cis3EnAcACP;

% C2n:1 (n=6:10) FabA-Enoyl-Acyl-ACPs
d_C12_FabA_EnAcACP_un = P.kcat9(5).*c.C12_FabA_BHyAcACP_un - P.k9_2r(5).*c.C12_FabA_EnAcACP_un + P.k9_3r(5).*c.FabA.*c.C12_EnAcACP_un - P.k9_3f(5).*c.C12_FabA_EnAcACP_un;
d_C14_FabA_EnAcACP_un = P.kcat9(6).*c.C14_FabA_BHyAcACP_un - P.k9_2r(6).*c.C14_FabA_EnAcACP_un + P.k9_3r(6).*c.FabA.*c.C14_EnAcACP_un - P.k9_3f(6).*c.C14_FabA_EnAcACP_un;
d_C16_FabA_EnAcACP_un = P.kcat9(7).*c.C16_FabA_BHyAcACP_un - P.k9_2r(7).*c.C16_FabA_EnAcACP_un + P.k9_3r(7).*c.FabA.*c.C16_EnAcACP_un - P.k9_3f(7).*c.C16_FabA_EnAcACP_un;
d_C18_FabA_EnAcACP_un = P.kcat9(8).*c.C18_FabA_BHyAcACP_un - P.k9_2r(8).*c.C18_FabA_EnAcACP_un + P.k9_3r(8).*c.FabA.*c.C18_EnAcACP_un - P.k9_3f(8).*c.C18_FabA_EnAcACP_un;
d_C20_FabA_EnAcACP_un = P.kcat9(9).*c.C20_FabA_BHyAcACP_un - P.k9_2r(9).*c.C20_FabA_EnAcACP_un + P.k9_3r(9).*c.FabA.*c.C20_EnAcACP_un - P.k9_3f(9).*c.C20_FabA_EnAcACP_un;

% FabI-NADH
d_FabI_NADH = P.k6_1f(1).*c.FabI.*c.NADH - P.k6_1r(1).*c.FabI_NADH...
 + P.k6_2r(1).*c.C4_FabI_NADH_EnAcACP - P.k6_2f(1).*c.FabI_NADH.*c.C4_EnAcACP...
 + P.k6_2r(2).*c.C6_FabI_NADH_EnAcACP - P.k6_2f(2).*c.FabI_NADH.*c.C6_EnAcACP...
 + P.k6_2r(3).*c.C8_FabI_NADH_EnAcACP - P.k6_2f(3).*c.FabI_NADH.*c.C8_EnAcACP...
 + P.k6_2r(4).*c.C10_FabI_NADH_EnAcACP - P.k6_2f(4).*c.FabI_NADH.*c.C10_EnAcACP...
 + P.k6_2r(5).*c.C12_FabI_NADH_EnAcACP - P.k6_2f(5).*c.FabI_NADH.*c.C12_EnAcACP...
 + P.k6_2r(6).*c.C14_FabI_NADH_EnAcACP - P.k6_2f(6).*c.FabI_NADH.*c.C14_EnAcACP...
 + P.k6_2r(7).*c.C16_FabI_NADH_EnAcACP - P.k6_2f(7).*c.FabI_NADH.*c.C16_EnAcACP...
 + P.k6_2r(8).*c.C18_FabI_NADH_EnAcACP - P.k6_2f(8).*c.FabI_NADH.*c.C18_EnAcACP...
 + P.k6_2r(9).*c.C20_FabI_NADH_EnAcACP - P.k6_2f(9).*c.FabI_NADH.*c.C20_EnAcACP...
 + P.k6_2r(5).*c.C12_FabI_NADH_EnAcACP_un - P.k6_2f(5).*c.FabI_NADH.*c.C12_EnAcACP_un...
 + P.k6_2r(6).*c.C14_FabI_NADH_EnAcACP_un - P.k6_2f(6).*c.FabI_NADH.*c.C14_EnAcACP_un...
 + P.k6_2r(7).*c.C16_FabI_NADH_EnAcACP_un - P.k6_2f(7).*c.FabI_NADH.*c.C16_EnAcACP_un...
 + P.k6_2r(8).*c.C18_FabI_NADH_EnAcACP_un - P.k6_2f(8).*c.FabI_NADH.*c.C18_EnAcACP_un...
 + P.k6_2r(9).*c.C20_FabI_NADH_EnAcACP_un - P.k6_2f(9).*c.FabI_NADH.*c.C20_EnAcACP_un;

% C2n (n=2:10) FabI-NADH-Enoyl-Acyl-ACPs
d_C4_FabI_NADH_EnAcACP = P.k6_2f(1).*c.FabI_NADH.*c.C4_EnAcACP - P.k6_2r(1).*c.C4_FabI_NADH_EnAcACP - P.kcat6(1).*c.C4_FabI_NADH_EnAcACP;
d_C6_FabI_NADH_EnAcACP = P.k6_2f(2).*c.FabI_NADH.*c.C6_EnAcACP - P.k6_2r(2).*c.C6_FabI_NADH_EnAcACP - P.kcat6(2).*c.C6_FabI_NADH_EnAcACP;
d_C8_FabI_NADH_EnAcACP = P.k6_2f(3).*c.FabI_NADH.*c.C8_EnAcACP - P.k6_2r(3).*c.C8_FabI_NADH_EnAcACP - P.kcat6(3).*c.C8_FabI_NADH_EnAcACP;
d_C10_FabI_NADH_EnAcACP = P.k6_2f(4).*c.FabI_NADH.*c.C10_EnAcACP - P.k6_2r(4).*c.C10_FabI_NADH_EnAcACP - P.kcat6(4).*c.C10_FabI_NADH_EnAcACP;
d_C12_FabI_NADH_EnAcACP = P.k6_2f(5).*c.FabI_NADH.*c.C12_EnAcACP - P.k6_2r(5).*c.C12_FabI_NADH_EnAcACP - P.kcat6(5).*c.C12_FabI_NADH_EnAcACP;
d_C14_FabI_NADH_EnAcACP = P.k6_2f(6).*c.FabI_NADH.*c.C14_EnAcACP - P.k6_2r(6).*c.C14_FabI_NADH_EnAcACP - P.kcat6(6).*c.C14_FabI_NADH_EnAcACP;
d_C16_FabI_NADH_EnAcACP = P.k6_2f(7).*c.FabI_NADH.*c.C16_EnAcACP - P.k6_2r(7).*c.C16_FabI_NADH_EnAcACP - P.kcat6(7).*c.C16_FabI_NADH_EnAcACP;
d_C18_FabI_NADH_EnAcACP = P.k6_2f(8).*c.FabI_NADH.*c.C18_EnAcACP - P.k6_2r(8).*c.C18_FabI_NADH_EnAcACP - P.kcat6(8).*c.C18_FabI_NADH_EnAcACP;
d_C20_FabI_NADH_EnAcACP = P.k6_2f(9).*c.FabI_NADH.*c.C20_EnAcACP - P.k6_2r(9).*c.C20_FabI_NADH_EnAcACP - P.kcat6(9).*c.C20_FabI_NADH_EnAcACP;

% C2n:1 (n=6:10) FabI-NADH-Enoyl-Acyl-ACPs
d_C12_FabI_NADH_EnAcACP_un = P.k6_2f(5).*c.FabI_NADH.*c.C12_EnAcACP_un - P.k6_2r(5).*c.C12_FabI_NADH_EnAcACP_un - P.kcat6(5).*c.C12_FabI_NADH_EnAcACP_un;
d_C14_FabI_NADH_EnAcACP_un = P.k6_2f(6).*c.FabI_NADH.*c.C14_EnAcACP_un - P.k6_2r(6).*c.C14_FabI_NADH_EnAcACP_un - P.kcat6(6).*c.C14_FabI_NADH_EnAcACP_un;
d_C16_FabI_NADH_EnAcACP_un = P.k6_2f(7).*c.FabI_NADH.*c.C16_EnAcACP_un - P.k6_2r(7).*c.C16_FabI_NADH_EnAcACP_un - P.kcat6(7).*c.C16_FabI_NADH_EnAcACP_un;
d_C18_FabI_NADH_EnAcACP_un = P.k6_2f(8).*c.FabI_NADH.*c.C18_EnAcACP_un - P.k6_2r(8).*c.C18_FabI_NADH_EnAcACP_un - P.kcat6(8).*c.C18_FabI_NADH_EnAcACP_un;
d_C20_FabI_NADH_EnAcACP_un = P.k6_2f(9).*c.FabI_NADH.*c.C20_EnAcACP_un - P.k6_2r(9).*c.C20_FabI_NADH_EnAcACP_un - P.kcat6(9).*c.C20_FabI_NADH_EnAcACP_un;

% C2n (n=2:10) TesA-Acyl-ACPs
d_C4_TesA_AcACP = P.k7_1f(1).*c.TesA.*c.C4_AcACP - P.k7_1r(1).*c.C4_TesA_AcACP - P.kcat7(1).*c.C4_TesA_AcACP;
d_C6_TesA_AcACP = P.k7_1f(2).*c.TesA.*c.C6_AcACP - P.k7_1r(2).*c.C6_TesA_AcACP - P.kcat7(2).*c.C6_TesA_AcACP;
d_C8_TesA_AcACP = P.k7_1f(3).*c.TesA.*c.C8_AcACP - P.k7_1r(3).*c.C8_TesA_AcACP - P.kcat7(3).*c.C8_TesA_AcACP;
d_C10_TesA_AcACP = P.k7_1f(4).*c.TesA.*c.C10_AcACP - P.k7_1r(4).*c.C10_TesA_AcACP - P.kcat7(4).*c.C10_TesA_AcACP;
d_C12_TesA_AcACP = P.k7_1f(5).*c.TesA.*c.C12_AcACP - P.k7_1r(5).*c.C12_TesA_AcACP - P.kcat7(5).*c.C12_TesA_AcACP;
d_C14_TesA_AcACP = P.k7_1f(6).*c.TesA.*c.C14_AcACP - P.k7_1r(6).*c.C14_TesA_AcACP - P.kcat7(6).*c.C14_TesA_AcACP;
d_C16_TesA_AcACP = P.k7_1f(7).*c.TesA.*c.C16_AcACP - P.k7_1r(7).*c.C16_TesA_AcACP - P.kcat7(7).*c.C16_TesA_AcACP;
d_C18_TesA_AcACP = P.k7_1f(8).*c.TesA.*c.C18_AcACP - P.k7_1r(8).*c.C18_TesA_AcACP - P.kcat7(8).*c.C18_TesA_AcACP;
d_C20_TesA_AcACP = P.k7_1f(9).*c.TesA.*c.C20_AcACP - P.k7_1r(9).*c.C20_TesA_AcACP - P.kcat7(9).*c.C20_TesA_AcACP;

% C2n:1 (n=6:10) TesA-Acyl-ACPs
d_C12_TesA_AcACP_un = P.k7_1f(5).*c.TesA.*c.C12_AcACP_un - P.k7_1r(5).*c.C12_TesA_AcACP_un - P.kcat7(5).*c.C12_TesA_AcACP_un;
d_C14_TesA_AcACP_un = P.k7_1f(6).*c.TesA.*c.C14_AcACP_un - P.k7_1r(6).*c.C14_TesA_AcACP_un - P.kcat7(6).*c.C14_TesA_AcACP_un;
d_C16_TesA_AcACP_un = P.k7_1f(7).*c.TesA.*c.C16_AcACP_un - P.k7_1r(7).*c.C16_TesA_AcACP_un - P.kcat7(7).*c.C16_TesA_AcACP_un;
d_C18_TesA_AcACP_un = P.k7_1f(8).*c.TesA.*c.C18_AcACP_un - P.k7_1r(8).*c.C18_TesA_AcACP_un - P.kcat7(8).*c.C18_TesA_AcACP_un;
d_C20_TesA_AcACP_un = P.k7_1f(9).*c.TesA.*c.C20_AcACP_un - P.k7_1r(9).*c.C20_TesA_AcACP_un - P.kcat7(9).*c.C20_TesA_AcACP_un;

% C2n (n=2:9) FabF-Acyl-ACPs
d_C4_FabF_AcACP = P.k8_1f(1).*c.FabF.*c.C4_AcACP - P.k8_1r(1).*c.C4_FabF_AcACP + P.k8_2r(1).*c.C4_FabF_Act.*c.ACP - P.k8_2f(1).*c.C4_FabF_AcACP;
d_C6_FabF_AcACP = P.k8_1f(2).*c.FabF.*c.C6_AcACP - P.k8_1r(2).*c.C6_FabF_AcACP + P.k8_2r(2).*c.C6_FabF_Act.*c.ACP - P.k8_2f(2).*c.C6_FabF_AcACP;
d_C8_FabF_AcACP = P.k8_1f(3).*c.FabF.*c.C8_AcACP - P.k8_1r(3).*c.C8_FabF_AcACP + P.k8_2r(3).*c.C8_FabF_Act.*c.ACP - P.k8_2f(3).*c.C8_FabF_AcACP;
d_C10_FabF_AcACP = P.k8_1f(4).*c.FabF.*c.C10_AcACP - P.k8_1r(4).*c.C10_FabF_AcACP + P.k8_2r(4).*c.C10_FabF_Act.*c.ACP - P.k8_2f(4).*c.C10_FabF_AcACP;
d_C12_FabF_AcACP = P.k8_1f(5).*c.FabF.*c.C12_AcACP - P.k8_1r(5).*c.C12_FabF_AcACP + P.k8_2r(5).*c.C12_FabF_Act.*c.ACP - P.k8_2f(5).*c.C12_FabF_AcACP;
d_C14_FabF_AcACP = P.k8_1f(6).*c.FabF.*c.C14_AcACP - P.k8_1r(6).*c.C14_FabF_AcACP + P.k8_2r(6).*c.C14_FabF_Act.*c.ACP - P.k8_2f(6).*c.C14_FabF_AcACP;
d_C16_FabF_AcACP = P.k8_1f(7).*c.FabF.*c.C16_AcACP - P.k8_1r(7).*c.C16_FabF_AcACP + P.k8_2r(7).*c.C16_FabF_Act.*c.ACP - P.k8_2f(7).*c.C16_FabF_AcACP;
d_C18_FabF_AcACP = P.k8_1f(8).*c.FabF.*c.C18_AcACP - P.k8_1r(8).*c.C18_FabF_AcACP + P.k8_2r(8).*c.C18_FabF_Act.*c.ACP - P.k8_2f(8).*c.C18_FabF_AcACP;

% C2n:1 (n=6:9) FabF-Acyl-ACPs
d_C12_FabF_AcACP_un = P.k8_1f(5).*c.FabF.*c.C12_AcACP_un - P.k8_1r(5).*c.C12_FabF_AcACP_un + P.k8_2r(5).*c.C12_FabF_Act_un.*c.ACP - P.k8_2f(5).*c.C12_FabF_AcACP_un;
d_C14_FabF_AcACP_un = P.k8_1f(6).*c.FabF.*c.C14_AcACP_un - P.k8_1r(6).*c.C14_FabF_AcACP_un + P.k8_2r(6).*c.C14_FabF_Act_un.*c.ACP - P.k8_2f(6).*c.C14_FabF_AcACP_un;
d_C16_FabF_AcACP_un = P.k8_1f(7).*c.FabF.*c.C16_AcACP_un - P.k8_1r(7).*c.C16_FabF_AcACP_un + P.k8_2r(7).*c.C16_FabF_Act_un.*c.ACP - P.k8_2f(7).*c.C16_FabF_AcACP_un;
d_C18_FabF_AcACP_un = P.k8_1f(8).*c.FabF.*c.C18_AcACP_un - P.k8_1r(8).*c.C18_FabF_AcACP_un + P.k8_2r(8).*c.C18_FabF_Act_un.*c.ACP - P.k8_2f(8).*c.C18_FabF_AcACP_un;

% C2n (n=2:9) FabF*
d_C4_FabF_Act = P.k8_2f(1).*c.C4_FabF_AcACP - P.k8_2r(1).*c.C4_FabF_Act.*c.ACP + P.k8_3r(1).*c.C7_FabF_Act_MalACP - P.k8_3f(1).*c.C4_FabF_Act.*c.C3_MalACP;
d_C6_FabF_Act = P.k8_2f(2).*c.C6_FabF_AcACP - P.k8_2r(2).*c.C6_FabF_Act.*c.ACP + P.k8_3r(2).*c.C9_FabF_Act_MalACP - P.k8_3f(2).*c.C6_FabF_Act.*c.C3_MalACP;
d_C8_FabF_Act = P.k8_2f(3).*c.C8_FabF_AcACP - P.k8_2r(3).*c.C8_FabF_Act.*c.ACP + P.k8_3r(3).*c.C11_FabF_Act_MalACP - P.k8_3f(3).*c.C8_FabF_Act.*c.C3_MalACP;
d_C10_FabF_Act = P.k8_2f(4).*c.C10_FabF_AcACP - P.k8_2r(4).*c.C10_FabF_Act.*c.ACP + P.k8_3r(4).*c.C13_FabF_Act_MalACP - P.k8_3f(4).*c.C10_FabF_Act.*c.C3_MalACP;
d_C12_FabF_Act = P.k8_2f(5).*c.C12_FabF_AcACP - P.k8_2r(5).*c.C12_FabF_Act.*c.ACP + P.k8_3r(5).*c.C15_FabF_Act_MalACP - P.k8_3f(5).*c.C12_FabF_Act.*c.C3_MalACP;
d_C14_FabF_Act = P.k8_2f(6).*c.C14_FabF_AcACP - P.k8_2r(6).*c.C14_FabF_Act.*c.ACP + P.k8_3r(6).*c.C17_FabF_Act_MalACP - P.k8_3f(6).*c.C14_FabF_Act.*c.C3_MalACP;
d_C16_FabF_Act = P.k8_2f(7).*c.C16_FabF_AcACP - P.k8_2r(7).*c.C16_FabF_Act.*c.ACP + P.k8_3r(7).*c.C19_FabF_Act_MalACP - P.k8_3f(7).*c.C16_FabF_Act.*c.C3_MalACP;
d_C18_FabF_Act = P.k8_2f(8).*c.C18_FabF_AcACP - P.k8_2r(8).*c.C18_FabF_Act.*c.ACP + P.k8_3r(8).*c.C21_FabF_Act_MalACP - P.k8_3f(8).*c.C18_FabF_Act.*c.C3_MalACP;

% C2n:1 (n=6:9) FabF*
d_C12_FabF_Act_un = P.k8_2f(5).*c.C12_FabF_AcACP_un - P.k8_2r(5).*c.C12_FabF_Act_un.*c.ACP + P.k8_3r(5).*c.C15_FabF_Act_MalACP_un - P.k8_3f(5).*c.C12_FabF_Act_un.*c.C3_MalACP;
d_C14_FabF_Act_un = P.k8_2f(6).*c.C14_FabF_AcACP_un - P.k8_2r(6).*c.C14_FabF_Act_un.*c.ACP + P.k8_3r(6).*c.C17_FabF_Act_MalACP_un - P.k8_3f(6).*c.C14_FabF_Act_un.*c.C3_MalACP;
d_C16_FabF_Act_un = P.k8_2f(7).*c.C16_FabF_AcACP_un - P.k8_2r(7).*c.C16_FabF_Act_un.*c.ACP + P.k8_3r(7).*c.C19_FabF_Act_MalACP_un - P.k8_3f(7).*c.C16_FabF_Act_un.*c.C3_MalACP;
d_C18_FabF_Act_un = P.k8_2f(8).*c.C18_FabF_AcACP_un - P.k8_2r(8).*c.C18_FabF_Act_un.*c.ACP + P.k8_3r(8).*c.C21_FabF_Act_MalACP_un - P.k8_3f(8).*c.C18_FabF_Act_un.*c.C3_MalACP;

% C2n (n=2:9) FabF*-Malonyl-ACPs
d_C7_FabF_Act_MalACP = P.k8_3f(1).*c.C4_FabF_Act.*c.C3_MalACP - P.k8_3r(1).*c.C7_FabF_Act_MalACP - P.kcat8(1).*c.C7_FabF_Act_MalACP;
d_C9_FabF_Act_MalACP = P.k8_3f(2).*c.C6_FabF_Act.*c.C3_MalACP - P.k8_3r(2).*c.C9_FabF_Act_MalACP - P.kcat8(2).*c.C9_FabF_Act_MalACP;
d_C11_FabF_Act_MalACP = P.k8_3f(3).*c.C8_FabF_Act.*c.C3_MalACP - P.k8_3r(3).*c.C11_FabF_Act_MalACP - P.kcat8(3).*c.C11_FabF_Act_MalACP;
d_C13_FabF_Act_MalACP = P.k8_3f(4).*c.C10_FabF_Act.*c.C3_MalACP - P.k8_3r(4).*c.C13_FabF_Act_MalACP - P.kcat8(4).*c.C13_FabF_Act_MalACP;
d_C15_FabF_Act_MalACP = P.k8_3f(5).*c.C12_FabF_Act.*c.C3_MalACP - P.k8_3r(5).*c.C15_FabF_Act_MalACP - P.kcat8(5).*c.C15_FabF_Act_MalACP;
d_C17_FabF_Act_MalACP = P.k8_3f(6).*c.C14_FabF_Act.*c.C3_MalACP - P.k8_3r(6).*c.C17_FabF_Act_MalACP - P.kcat8(6).*c.C17_FabF_Act_MalACP;
d_C19_FabF_Act_MalACP = P.k8_3f(7).*c.C16_FabF_Act.*c.C3_MalACP - P.k8_3r(7).*c.C19_FabF_Act_MalACP - P.kcat8(7).*c.C19_FabF_Act_MalACP;
d_C21_FabF_Act_MalACP = P.k8_3f(8).*c.C18_FabF_Act.*c.C3_MalACP - P.k8_3r(8).*c.C21_FabF_Act_MalACP - P.kcat8(8).*c.C21_FabF_Act_MalACP;

% C2n:1 (n=6:9) FabF*-Malonyl-ACPs
d_C15_FabF_Act_MalACP_un = P.k8_3f(5).*c.C12_FabF_Act_un.*c.C3_MalACP - P.k8_3r(5).*c.C15_FabF_Act_MalACP_un - P.kcat8_un(5).*c.C15_FabF_Act_MalACP_un;
d_C17_FabF_Act_MalACP_un = P.k8_3f(6).*c.C14_FabF_Act_un.*c.C3_MalACP - P.k8_3r(6).*c.C17_FabF_Act_MalACP_un - P.kcat8_un(6).*c.C17_FabF_Act_MalACP_un;
d_C19_FabF_Act_MalACP_un = P.k8_3f(7).*c.C16_FabF_Act_un.*c.C3_MalACP - P.k8_3r(7).*c.C19_FabF_Act_MalACP_un - P.kcat8_un(7).*c.C19_FabF_Act_MalACP_un;
d_C21_FabF_Act_MalACP_un = P.k8_3f(8).*c.C18_FabF_Act_un.*c.C3_MalACP - P.k8_3r(8).*c.C21_FabF_Act_MalACP_un - P.kcat8_un(8).*c.C21_FabF_Act_MalACP_un;

% C2n (n=2:9) FabB-Acyl-ACPs
d_C4_FabB_AcACP = P.k10_1f(1).*c.FabB.*c.C4_AcACP - P.k10_1r(1).*c.C4_FabB_AcACP + P.k10_2r(1).*c.C4_FabB_Act.*c.ACP - P.k10_2f(1).*c.C4_FabB_AcACP;
d_C6_FabB_AcACP = P.k10_1f(2).*c.FabB.*c.C6_AcACP - P.k10_1r(2).*c.C6_FabB_AcACP + P.k10_2r(2).*c.C6_FabB_Act.*c.ACP - P.k10_2f(2).*c.C6_FabB_AcACP;
d_C8_FabB_AcACP = P.k10_1f(3).*c.FabB.*c.C8_AcACP - P.k10_1r(3).*c.C8_FabB_AcACP + P.k10_2r(3).*c.C8_FabB_Act.*c.ACP - P.k10_2f(3).*c.C8_FabB_AcACP;
d_C10_FabB_AcACP = P.k10_1f(4).*c.FabB.*c.C10_AcACP - P.k10_1r(4).*c.C10_FabB_AcACP + P.k10_2r(4).*c.C10_FabB_Act.*c.ACP - P.k10_2f(4).*c.C10_FabB_AcACP;
d_C12_FabB_AcACP = P.k10_1f(5).*c.FabB.*c.C12_AcACP - P.k10_1r(5).*c.C12_FabB_AcACP + P.k10_2r(5).*c.C12_FabB_Act.*c.ACP - P.k10_2f(5).*c.C12_FabB_AcACP;
d_C14_FabB_AcACP = P.k10_1f(6).*c.FabB.*c.C14_AcACP - P.k10_1r(6).*c.C14_FabB_AcACP + P.k10_2r(6).*c.C14_FabB_Act.*c.ACP - P.k10_2f(6).*c.C14_FabB_AcACP;
d_C16_FabB_AcACP = P.k10_1f(7).*c.FabB.*c.C16_AcACP - P.k10_1r(7).*c.C16_FabB_AcACP + P.k10_2r(7).*c.C16_FabB_Act.*c.ACP - P.k10_2f(7).*c.C16_FabB_AcACP;
d_C18_FabB_AcACP = P.k10_1f(8).*c.FabB.*c.C18_AcACP - P.k10_1r(8).*c.C18_FabB_AcACP + P.k10_2r(8).*c.C18_FabB_Act.*c.ACP - P.k10_2f(8).*c.C18_FabB_AcACP;

% C2n:1 (n=6:9) FabB-Acyl-ACPs
d_C12_FabB_AcACP_un = P.k10_1f(5).*c.FabB.*c.C12_AcACP_un - P.k10_1r(5).*c.C12_FabB_AcACP_un + P.k10_2r(5).*c.C12_FabB_Act_un.*c.ACP - P.k10_2f(5).*c.C12_FabB_AcACP_un;
d_C14_FabB_AcACP_un = P.k10_1f(6).*c.FabB.*c.C14_AcACP_un - P.k10_1r(6).*c.C14_FabB_AcACP_un + P.k10_2r(6).*c.C14_FabB_Act_un.*c.ACP - P.k10_2f(6).*c.C14_FabB_AcACP_un;
d_C16_FabB_AcACP_un = P.k10_1f(7).*c.FabB.*c.C16_AcACP_un - P.k10_1r(7).*c.C16_FabB_AcACP_un + P.k10_2r(7).*c.C16_FabB_Act_un.*c.ACP - P.k10_2f(7).*c.C16_FabB_AcACP_un;
d_C18_FabB_AcACP_un = P.k10_1f(8).*c.FabB.*c.C18_AcACP_un - P.k10_1r(8).*c.C18_FabB_AcACP_un + P.k10_2r(8).*c.C18_FabB_Act_un.*c.ACP - P.k10_2f(8).*c.C18_FabB_AcACP_un;

% C2n (n=2:9) FabB*
d_C4_FabB_Act = P.k10_2f(1).*c.C4_FabB_AcACP - P.k10_2r(1).*c.C4_FabB_Act.*c.ACP + P.k10_3r(1).*c.C7_FabB_Act_MalACP - P.k10_3f(1).*c.C4_FabB_Act.*c.C3_MalACP;
d_C6_FabB_Act = P.k10_2f(2).*c.C6_FabB_AcACP - P.k10_2r(2).*c.C6_FabB_Act.*c.ACP + P.k10_3r(2).*c.C9_FabB_Act_MalACP - P.k10_3f(2).*c.C6_FabB_Act.*c.C3_MalACP;
d_C8_FabB_Act = P.k10_2f(3).*c.C8_FabB_AcACP - P.k10_2r(3).*c.C8_FabB_Act.*c.ACP + P.k10_3r(3).*c.C11_FabB_Act_MalACP - P.k10_3f(3).*c.C8_FabB_Act.*c.C3_MalACP;
d_C10_FabB_Act = P.k10_2f(4).*c.C10_FabB_AcACP - P.k10_2r(4).*c.C10_FabB_Act.*c.ACP + P.k10_3r(4).*c.C13_FabB_Act_MalACP - P.k10_3f(4).*c.C10_FabB_Act.*c.C3_MalACP;
d_C12_FabB_Act = P.k10_2f(5).*c.C12_FabB_AcACP - P.k10_2r(5).*c.C12_FabB_Act.*c.ACP + P.k10_3r(5).*c.C15_FabB_Act_MalACP - P.k10_3f(5).*c.C12_FabB_Act.*c.C3_MalACP;
d_C14_FabB_Act = P.k10_2f(6).*c.C14_FabB_AcACP - P.k10_2r(6).*c.C14_FabB_Act.*c.ACP + P.k10_3r(6).*c.C17_FabB_Act_MalACP - P.k10_3f(6).*c.C14_FabB_Act.*c.C3_MalACP;
d_C16_FabB_Act = P.k10_2f(7).*c.C16_FabB_AcACP - P.k10_2r(7).*c.C16_FabB_Act.*c.ACP + P.k10_3r(7).*c.C19_FabB_Act_MalACP - P.k10_3f(7).*c.C16_FabB_Act.*c.C3_MalACP;
d_C18_FabB_Act = P.k10_2f(8).*c.C18_FabB_AcACP - P.k10_2r(8).*c.C18_FabB_Act.*c.ACP + P.k10_3r(8).*c.C21_FabB_Act_MalACP - P.k10_3f(8).*c.C18_FabB_Act.*c.C3_MalACP;

% C2n:1 (n=6:9) FabB*
d_C12_FabB_Act_un = P.k10_2f(5).*c.C12_FabB_AcACP_un - P.k10_2r(5).*c.C12_FabB_Act_un.*c.ACP + P.k10_3r(5).*c.C15_FabB_Act_MalACP_un - P.k10_3f(5).*c.C12_FabB_Act_un.*c.C3_MalACP;
d_C14_FabB_Act_un = P.k10_2f(6).*c.C14_FabB_AcACP_un - P.k10_2r(6).*c.C14_FabB_Act_un.*c.ACP + P.k10_3r(6).*c.C17_FabB_Act_MalACP_un - P.k10_3f(6).*c.C14_FabB_Act_un.*c.C3_MalACP;
d_C16_FabB_Act_un = P.k10_2f(7).*c.C16_FabB_AcACP_un - P.k10_2r(7).*c.C16_FabB_Act_un.*c.ACP + P.k10_3r(7).*c.C19_FabB_Act_MalACP_un - P.k10_3f(7).*c.C16_FabB_Act_un.*c.C3_MalACP;
d_C18_FabB_Act_un = P.k10_2f(8).*c.C18_FabB_AcACP_un - P.k10_2r(8).*c.C18_FabB_Act_un.*c.ACP + P.k10_3r(8).*c.C21_FabB_Act_MalACP_un - P.k10_3f(8).*c.C18_FabB_Act_un.*c.C3_MalACP;

% C2n (n=2:9) FabB*-Malonyl-ACPs
d_C7_FabB_Act_MalACP = P.k10_3f(1).*c.C4_FabB_Act.*c.C3_MalACP - P.k10_3r(1).*c.C7_FabB_Act_MalACP - P.kcat10(1).*c.C7_FabB_Act_MalACP;
d_C9_FabB_Act_MalACP = P.k10_3f(2).*c.C6_FabB_Act.*c.C3_MalACP - P.k10_3r(2).*c.C9_FabB_Act_MalACP - P.kcat10(2).*c.C9_FabB_Act_MalACP;
d_C11_FabB_Act_MalACP = P.k10_3f(3).*c.C8_FabB_Act.*c.C3_MalACP - P.k10_3r(3).*c.C11_FabB_Act_MalACP - P.kcat10(3).*c.C11_FabB_Act_MalACP;
d_C13_FabB_Act_MalACP = P.k10_3f(4).*c.C10_FabB_Act.*c.C3_MalACP - P.k10_3r(4).*c.C13_FabB_Act_MalACP - P.kcat10(4).*c.C13_FabB_Act_MalACP;
d_C15_FabB_Act_MalACP = P.k10_3f(5).*c.C12_FabB_Act.*c.C3_MalACP - P.k10_3r(5).*c.C15_FabB_Act_MalACP - P.kcat10(5).*c.C15_FabB_Act_MalACP;
d_C17_FabB_Act_MalACP = P.k10_3f(6).*c.C14_FabB_Act.*c.C3_MalACP - P.k10_3r(6).*c.C17_FabB_Act_MalACP - P.kcat10(6).*c.C17_FabB_Act_MalACP;
d_C19_FabB_Act_MalACP = P.k10_3f(7).*c.C16_FabB_Act.*c.C3_MalACP - P.k10_3r(7).*c.C19_FabB_Act_MalACP - P.kcat10(7).*c.C19_FabB_Act_MalACP;
d_C21_FabB_Act_MalACP = P.k10_3f(8).*c.C18_FabB_Act.*c.C3_MalACP - P.k10_3r(8).*c.C21_FabB_Act_MalACP - P.kcat10(8).*c.C21_FabB_Act_MalACP;

% C2n:1 (n=6:9) FabB*-Malonyl-ACPs
d_C15_FabB_Act_MalACP_un = P.k10_3f(5).*c.C12_FabB_Act_un.*c.C3_MalACP - P.k10_3r(5).*c.C15_FabB_Act_MalACP_un - P.kcat10_un(5).*c.C15_FabB_Act_MalACP_un;
d_C17_FabB_Act_MalACP_un = P.k10_3f(6).*c.C14_FabB_Act_un.*c.C3_MalACP - P.k10_3r(6).*c.C17_FabB_Act_MalACP_un - P.kcat10_un(6).*c.C17_FabB_Act_MalACP_un;
d_C19_FabB_Act_MalACP_un = P.k10_3f(7).*c.C16_FabB_Act_un.*c.C3_MalACP - P.k10_3r(7).*c.C19_FabB_Act_MalACP_un - P.kcat10_un(7).*c.C19_FabB_Act_MalACP_un;
d_C21_FabB_Act_MalACP_un = P.k10_3f(8).*c.C18_FabB_Act_un.*c.C3_MalACP - P.k10_3r(8).*c.C21_FabB_Act_MalACP_un - P.kcat10_un(8).*c.C21_FabB_Act_MalACP_un;

% FabB-(C10 cis-3-Enoyl-Acyl-ACP)
d_C10_FabB_cis3EnAcACP = P.k10_1f(4).*c.FabB.*c.C10_cis3EnAcACP - P.k10_1r(4).*c.C10_FabB_cis3EnAcACP + P.k10_2r(4).*c.C10_FabB_Act_cis3.*c.ACP - P.k10_2f(4).*c.C10_FabB_cis3EnAcACP;

% C10 cis-3-FabB*
d_C10_FabB_Act_cis3 = P.k10_2f(4).*c.C10_FabB_cis3EnAcACP - P.k10_2r(4).*c.C10_FabB_Act_cis3.*c.ACP + P.k10_3r(4).*c.C13_FabB_Act_cis3MalACP - P.k10_3f(4).*c.C10_FabB_Act_cis3.*c.C3_MalACP;

% C10 cis-3-FabB*-Malonyl-ACP
d_C13_FabB_Act_cis3MalACP = P.k10_3f(4).*c.C10_FabB_Act_cis3.*c.C3_MalACP - P.k10_3r(4).*c.C13_FabB_Act_cis3MalACP - P.kcat10_un(4).*c.C13_FabB_Act_cis3MalACP;

% C2n (n=2:10) FabH-Acyl-ACPs
d_C4_FabH_AcACP = P.k3_4f(1).*c.FabH.*c.C4_AcACP - P.k3_4r(1).*c.C4_FabH_AcACP;
d_C6_FabH_AcACP = P.k3_4f(2).*c.FabH.*c.C6_AcACP - P.k3_4r(2).*c.C6_FabH_AcACP;
d_C8_FabH_AcACP = P.k3_4f(3).*c.FabH.*c.C8_AcACP - P.k3_4r(3).*c.C8_FabH_AcACP;
d_C10_FabH_AcACP = P.k3_4f(4).*c.FabH.*c.C10_AcACP - P.k3_4r(4).*c.C10_FabH_AcACP;
d_C12_FabH_AcACP = P.k3_4f(5).*c.FabH.*c.C12_AcACP - P.k3_4r(5).*c.C12_FabH_AcACP;
d_C14_FabH_AcACP = P.k3_4f(6).*c.FabH.*c.C14_AcACP - P.k3_4r(6).*c.C14_FabH_AcACP;
d_C16_FabH_AcACP = P.k3_4f(7).*c.FabH.*c.C16_AcACP - P.k3_4r(7).*c.C16_FabH_AcACP;
d_C18_FabH_AcACP = P.k3_4f(8).*c.FabH.*c.C18_AcACP - P.k3_4r(8).*c.C18_FabH_AcACP;
d_C20_FabH_AcACP = P.k3_4f(9).*c.FabH.*c.C20_AcACP - P.k3_4r(9).*c.C20_FabH_AcACP;

% C2n:1 (n=6:10) FabH-Acyl-ACPs
d_C12_FabH_AcACP_un = P.k3_4f(5).*c.FabH.*c.C12_AcACP_un - P.k3_4r(5).*c.C12_FabH_AcACP_un;
d_C14_FabH_AcACP_un = P.k3_4f(6).*c.FabH.*c.C14_AcACP_un - P.k3_4r(6).*c.C14_FabH_AcACP_un;
d_C16_FabH_AcACP_un = P.k3_4f(7).*c.FabH.*c.C16_AcACP_un - P.k3_4r(7).*c.C16_FabH_AcACP_un;
d_C18_FabH_AcACP_un = P.k3_4f(8).*c.FabH.*c.C18_AcACP_un - P.k3_4r(8).*c.C18_FabH_AcACP_un;
d_C20_FabH_AcACP_un = P.k3_4f(9).*c.FabH.*c.C20_AcACP_un - P.k3_4r(9).*c.C20_FabH_AcACP_un;

% C2n (n=2:10) FabH*-Acyl-ACPs
d_C6_FabH_Act_AcACP = P.k3_5f(1).*c.C2_FabH_Act.*c.C4_AcACP - P.k3_5r(1).*c.C6_FabH_Act_AcACP;
d_C8_FabH_Act_AcACP = P.k3_5f(2).*c.C2_FabH_Act.*c.C6_AcACP - P.k3_5r(2).*c.C8_FabH_Act_AcACP;
d_C10_FabH_Act_AcACP = P.k3_5f(3).*c.C2_FabH_Act.*c.C8_AcACP - P.k3_5r(3).*c.C10_FabH_Act_AcACP;
d_C12_FabH_Act_AcACP = P.k3_5f(4).*c.C2_FabH_Act.*c.C10_AcACP - P.k3_5r(4).*c.C12_FabH_Act_AcACP;
d_C14_FabH_Act_AcACP = P.k3_5f(5).*c.C2_FabH_Act.*c.C12_AcACP - P.k3_5r(5).*c.C14_FabH_Act_AcACP;
d_C16_FabH_Act_AcACP = P.k3_5f(6).*c.C2_FabH_Act.*c.C14_AcACP - P.k3_5r(6).*c.C16_FabH_Act_AcACP;
d_C18_FabH_Act_AcACP = P.k3_5f(7).*c.C2_FabH_Act.*c.C16_AcACP - P.k3_5r(7).*c.C18_FabH_Act_AcACP;
d_C20_FabH_Act_AcACP = P.k3_5f(8).*c.C2_FabH_Act.*c.C18_AcACP - P.k3_5r(8).*c.C20_FabH_Act_AcACP;
d_C22_FabH_Act_AcACP = P.k3_5f(9).*c.C2_FabH_Act.*c.C20_AcACP - P.k3_5r(9).*c.C22_FabH_Act_AcACP;

% C2n:1 (n=6:10) FabH*-Acyl-ACPs
d_C14_FabH_Act_AcACP_un = P.k3_5f(5).*c.C2_FabH_Act.*c.C12_AcACP_un - P.k3_5r(5).*c.C14_FabH_Act_AcACP_un;
d_C16_FabH_Act_AcACP_un = P.k3_5f(6).*c.C2_FabH_Act.*c.C14_AcACP_un - P.k3_5r(6).*c.C16_FabH_Act_AcACP_un;
d_C18_FabH_Act_AcACP_un = P.k3_5f(7).*c.C2_FabH_Act.*c.C16_AcACP_un - P.k3_5r(7).*c.C18_FabH_Act_AcACP_un;
d_C20_FabH_Act_AcACP_un = P.k3_5f(8).*c.C2_FabH_Act.*c.C18_AcACP_un - P.k3_5r(8).*c.C20_FabH_Act_AcACP_un;
d_C22_FabH_Act_AcACP_un = P.k3_5f(9).*c.C2_FabH_Act.*c.C20_AcACP_un - P.k3_5r(9).*c.C22_FabH_Act_AcACP_un;

% TesA-ACP
d_TesA_ACP = P.k7_inh_f.*c.TesA.*c.ACP - P.k7_inh_r.*c.TesA_ACP;

% FabH-ACP
d_FabH_ACP = P.k3_inh_f.*c.FabH.*c.ACP - P.k3_inh_r.*c.FabH_ACP;

% FabG-ACP
d_FabG_ACP = P.k4_inh_f.*c.FabG.*c.ACP - P.k4_inh_r.*c.FabG_ACP;

% FabZ-ACP
d_FabZ_ACP = P.k5_inh_f.*c.FabZ.*c.ACP - P.k5_inh_r.*c.FabZ_ACP;

% FabI-ACP
d_FabI_ACP = P.k6_inh_f.*c.FabI.*c.ACP - P.k6_inh_r.*c.FabI_ACP;

% FabF-ACP
d_FabF_ACP = P.k8_inh_f.*c.FabF.*c.ACP - P.k8_inh_r.*c.FabF_ACP;

% FabA-ACP
d_FabA_ACP = P.k9_inh_f.*c.FabA.*c.ACP - P.k9_inh_r.*c.FabA_ACP;

% FabB-ACP
d_FabB_ACP = P.k10_inh_f.*c.FabB.*c.ACP - P.k10_inh_r.*c.FabB_ACP;


% Giving FabB FabH-like activity
% FabB-Acetyl-CoA
d_C2_FabB_AcCoA = P.k10_4f.*c.FabB.*c.C2_AcCoA - P.k10_4r.*c.C2_FabB_AcCoA + P.k10_5r.*c.C2_FabB_Act.*c.CoA - P.k10_5f.*c.C2_FabB_AcCoA;

% FabB*
d_C2_FabB_Act = P.k10_5f.*c.C2_FabB_AcCoA - P.k10_5r.*c.C2_FabB_Act.*c.CoA + P.k10_6r.*c.C5_FabB_Act_MalACP - P.k10_6f.*c.C2_FabB_Act.*c.C3_MalACP + P.k10_9f.*c.C2_FabB_AcACP - P.k10_9r.*c.C2_FabB_Act.*c.ACP;

% FabB*-Malonyl-ACP
d_C5_FabB_Act_MalACP = P.k10_6f.*c.C2_FabB_Act.*c.C3_MalACP - P.k10_6r.*c.C5_FabB_Act_MalACP - P.kcat10_H.*c.C5_FabB_Act_MalACP; 

% Giving FabF FabH-like activity
% FabF-Acetyl-CoA
d_C2_FabF_AcCoA = P.k8_4f.*c.FabF.*c.C2_AcCoA - P.k8_4r.*c.C2_FabF_AcCoA + P.k8_5r.*c.C2_FabF_Act.*c.CoA - P.k8_5f.*c.C2_FabF_AcCoA; 

% FabF*
d_C2_FabF_Act = P.k8_5f.*c.C2_FabF_AcCoA - P.k8_5r.*c.C2_FabF_Act.*c.CoA + P.k8_6r.*c.C5_FabF_Act_MalACP - P.k8_6f.*c.C2_FabF_Act.*c.C3_MalACP + P.k8_9f.*c.C2_FabF_AcACP - P.k8_9r.*c.C2_FabF_Act.*c.ACP; 

% FabF*-Malonyl-ACP
d_C5_FabF_Act_MalACP = P.k8_6f.*c.C2_FabF_Act.*c.C3_MalACP - P.k8_6r.*c.C5_FabF_Act_MalACP - P.kcat8_H.*c.C5_FabF_Act_MalACP; 

% FabB and FabF decarboxylating mACP to form aACP and reacting with it to form activated enzyme (initiation)
% FabB-Malonyl-ACP
d_C3_FabB_MalACP = P.k10_7f.*c.FabB.*c.C3_MalACP - P.k10_7r.*c.C3_FabB_MalACP - P.kcat10_CO2.*c.C3_FabB_MalACP;

% Acetyl-ACP (- FabF - FabB - FabH)
d_C2_AcACP = P.k8_8r.*c.C2_FabF_AcACP - P.k8_8f.*c.FabF.*c.C2_AcACP + P.k10_8r.*c.C2_FabB_AcACP - P.k10_8f.*c.FabB.*c.C2_AcACP + P.k3_7r.*c.C2_FabH_AcACP - P.k3_7f.*c.FabH.*c.C2_AcACP; 

% FabB-Acetyl-ACP
d_C2_FabB_AcACP = P.kcat10_CO2.*c.C3_FabB_MalACP + P.k10_8f.*c.FabB.*c.C2_AcACP - P.k10_8r.*c.C2_FabB_AcACP + P.k10_9r.*c.C2_FabB_Act.*c.ACP - P.k10_9f.*c.C2_FabB_AcACP; 

% FabF-Malonyl-ACP
d_C3_FabF_MalACP = P.k8_7f.*c.FabF.*c.C3_MalACP - P.k8_7r.*c.C3_FabF_MalACP - P.kcat8_CO2.*c.C3_FabF_MalACP; 

% FabF-Acetyl-ACP
d_C2_FabF_AcACP = P.kcat8_CO2.*c.C3_FabF_MalACP + P.k8_8f.*c.FabF.*c.C2_AcACP - P.k8_8r.*c.C2_FabF_AcACP + P.k8_9r.*c.C2_FabF_Act.*c.ACP - P.k8_9f.*c.C2_FabF_AcACP; 

% Giving FabH decarboxylative activity
% FabH-Acetyl-ACP
d_C2_FabH_AcACP = P.kcat3_CO2.*c.C3_FabH_MalACP + P.k3_7f.*c.FabH.*c.C2_AcACP - P.k3_7r.*c.C2_FabH_AcACP + P.k3_8r.*c.C2_FabH_Act.*c.ACP - P.k3_8f.*c.C2_FabH_AcACP;

% FabH-Malonyl-ACP
d_C3_FabH_MalACP = P.k3_6f.*c.FabH.*c.C3_MalACP - P.k3_6r.*c.C3_FabH_MalACP - P.kcat3_CO2.*c.C3_FabH_MalACP;

dcdt = [d_ATP; d_C1_Bicarbonate; d_C2_AcCoA; d_C4_SucCoA; d_C6_HexCoA; d_C8_OcCoA; d_C10_DecCoA; d_C12_LauCoA; d_C14_EthCoA; d_C16_PalCoA; d_C18_OcDecCoA; d_ACP;...
 d_NADPH; d_NADP; d_NADH; d_NAD; d_ADP; d_C3_MalCoA; d_CoA; d_C3_MalACP; d_C1_CO2; d_C4_BKeACP; d_C6_BKeACP; d_C8_BKeACP; d_C10_BKeACP; d_C12_BKeACP;...
 d_C14_BKeACP; d_C16_BKeACP; d_C18_BKeACP; d_C20_BKeACP; d_C12_BKeACP_un; d_C14_BKeACP_un; d_C16_BKeACP_un; d_C18_BKeACP_un; d_C20_BKeACP_un; d_C4_BHyAcACP;...
 d_C6_BHyAcACP; d_C8_BHyAcACP; d_C10_BHyAcACP; d_C12_BHyAcACP; d_C14_BHyAcACP; d_C16_BHyAcACP; d_C18_BHyAcACP; d_C20_BHyAcACP; d_C12_BHyAcACP_un;...
 d_C14_BHyAcACP_un; d_C16_BHyAcACP_un; d_C18_BHyAcACP_un; d_C20_BHyAcACP_un; d_C4_EnAcACP; d_C6_EnAcACP; d_C8_EnAcACP; d_C10_EnAcACP; d_C12_EnAcACP;...
 d_C14_EnAcACP; d_C16_EnAcACP; d_C18_EnAcACP; d_C20_EnAcACP; d_C10_cis3EnAcACP; d_C12_EnAcACP_un; d_C14_EnAcACP_un; d_C16_EnAcACP_un; d_C18_EnAcACP_un;...
 d_C20_EnAcACP_un; d_C4_AcACP; d_C6_AcACP; d_C8_AcACP; d_C10_AcACP; d_C12_AcACP; d_C14_AcACP; d_C16_AcACP; d_C18_AcACP; d_C20_AcACP;...
 d_C12_AcACP_un; d_C14_AcACP_un; d_C16_AcACP_un; d_C18_AcACP_un; d_C20_AcACP_un; d_C4_FA; d_C6_FA; d_C8_FA; d_C10_FA; d_C12_FA; d_C14_FA;...
 d_C16_FA; d_C18_FA; d_C20_FA; d_C12_FA_un; d_C14_FA_un; d_C16_FA_un; d_C18_FA_un; d_C20_FA_un; d_BC_ATP; d_C1_BC_ATP_HCO3; d_C1_BC_Pi_HCO3;...
 d_C1_BC_Pi_HCO3_BCCP_Biotin; d_C1_BCCP_Biotin_CO2; d_C1_CT_BCCP_Biotin_CO2; d_C1_CT_Act; d_C3_CT_Act_AcCoA; d_C3_FabD_MalCoA;...
 d_C3_FabD_Act; d_C3_FabD_Act_ACP; d_C2_FabH_CoA; d_C4_FabH_CoA; d_C6_FabH_CoA; d_C8_FabH_CoA; d_C10_FabH_CoA; d_C12_FabH_CoA; d_C14_FabH_CoA;...
 d_C16_FabH_CoA; d_C18_FabH_CoA; d_C2_FabH_Act; d_C4_FabH_Act; d_C6_FabH_Act; d_C8_FabH_Act; d_C10_FabH_Act; d_C12_FabH_Act; d_C14_FabH_Act;...
 d_C16_FabH_Act; d_C18_FabH_Act; d_C5_FabH_Act_MalACP; d_C7_FabH_Act_MalACP; d_C9_FabH_Act_MalACP; d_C11_FabH_Act_MalACP; d_C13_FabH_Act_MalACP;...
 d_C15_FabH_Act_MalACP; d_C17_FabH_Act_MalACP; d_C19_FabH_Act_MalACP; d_C21_FabH_Act_MalACP; d_FabG_NADPH; d_C4_FabG_NADPH_BKeACP;...
 d_C6_FabG_NADPH_BKeACP; d_C8_FabG_NADPH_BKeACP; d_C10_FabG_NADPH_BKeACP; d_C12_FabG_NADPH_BKeACP; d_C14_FabG_NADPH_BKeACP;...
 d_C16_FabG_NADPH_BKeACP; d_C18_FabG_NADPH_BKeACP; d_C20_FabG_NADPH_BKeACP; d_C12_FabG_NADPH_BKeACP_un; d_C14_FabG_NADPH_BKeACP_un;...
 d_C16_FabG_NADPH_BKeACP_un; d_C18_FabG_NADPH_BKeACP_un; d_C20_FabG_NADPH_BKeACP_un; d_C4_FabZ_BHyAcACP; d_C6_FabZ_BHyAcACP;...
 d_C8_FabZ_BHyAcACP; d_C10_FabZ_BHyAcACP; d_C12_FabZ_BHyAcACP; d_C14_FabZ_BHyAcACP; d_C16_FabZ_BHyAcACP; d_C18_FabZ_BHyAcACP;...
 d_C20_FabZ_BHyAcACP; d_C12_FabZ_BHyAcACP_un; d_C14_FabZ_BHyAcACP_un; d_C16_FabZ_BHyAcACP_un; d_C18_FabZ_BHyAcACP_un; d_C20_FabZ_BHyAcACP_un;...
 d_C4_FabZ_EnAcACP; d_C6_FabZ_EnAcACP; d_C8_FabZ_EnAcACP; d_C10_FabZ_EnAcACP; d_C12_FabZ_EnAcACP; d_C14_FabZ_EnAcACP; d_C16_FabZ_EnAcACP;...
 d_C18_FabZ_EnAcACP; d_C20_FabZ_EnAcACP; d_C12_FabZ_EnAcACP_un; d_C14_FabZ_EnAcACP_un; d_C16_FabZ_EnAcACP_un; d_C18_FabZ_EnAcACP_un;...
 d_C20_FabZ_EnAcACP_un; d_C4_FabA_BHyAcACP; d_C6_FabA_BHyAcACP; d_C8_FabA_BHyAcACP; d_C10_FabA_BHyAcACP; d_C12_FabA_BHyAcACP;...
 d_C14_FabA_BHyAcACP; d_C16_FabA_BHyAcACP; d_C18_FabA_BHyAcACP; d_C20_FabA_BHyAcACP; d_C12_FabA_BHyAcACP_un; d_C14_FabA_BHyAcACP_un;...
 d_C16_FabA_BHyAcACP_un; d_C18_FabA_BHyAcACP_un; d_C20_FabA_BHyAcACP_un; d_C4_FabA_EnAcACP; d_C6_FabA_EnAcACP; d_C8_FabA_EnAcACP; d_C10_FabA_EnAcACP;...
 d_C12_FabA_EnAcACP; d_C14_FabA_EnAcACP; d_C16_FabA_EnAcACP; d_C18_FabA_EnAcACP; d_C20_FabA_EnAcACP; d_C10_FabA_cis3EnAcACP; d_C12_FabA_EnAcACP_un;...
 d_C14_FabA_EnAcACP_un; d_C16_FabA_EnAcACP_un; d_C18_FabA_EnAcACP_un; d_C20_FabA_EnAcACP_un; d_FabI_NADH; d_C4_FabI_NADH_EnAcACP;...
 d_C6_FabI_NADH_EnAcACP; d_C8_FabI_NADH_EnAcACP; d_C10_FabI_NADH_EnAcACP; d_C12_FabI_NADH_EnAcACP; d_C14_FabI_NADH_EnAcACP;...
 d_C16_FabI_NADH_EnAcACP; d_C18_FabI_NADH_EnAcACP; d_C20_FabI_NADH_EnAcACP; d_C12_FabI_NADH_EnAcACP_un; d_C14_FabI_NADH_EnAcACP_un;...
 d_C16_FabI_NADH_EnAcACP_un; d_C18_FabI_NADH_EnAcACP_un; d_C20_FabI_NADH_EnAcACP_un; d_C4_TesA_AcACP; d_C6_TesA_AcACP; d_C8_TesA_AcACP;...
 d_C10_TesA_AcACP; d_C12_TesA_AcACP; d_C14_TesA_AcACP; d_C16_TesA_AcACP; d_C18_TesA_AcACP; d_C20_TesA_AcACP; d_C12_TesA_AcACP_un; d_C14_TesA_AcACP_un;...
 d_C16_TesA_AcACP_un; d_C18_TesA_AcACP_un; d_C20_TesA_AcACP_un; d_C4_FabF_AcACP; d_C6_FabF_AcACP; d_C8_FabF_AcACP; d_C10_FabF_AcACP;...
 d_C12_FabF_AcACP; d_C14_FabF_AcACP; d_C16_FabF_AcACP; d_C18_FabF_AcACP; d_C12_FabF_AcACP_un; d_C14_FabF_AcACP_un; d_C16_FabF_AcACP_un;...
 d_C18_FabF_AcACP_un; d_C4_FabF_Act; d_C6_FabF_Act; d_C8_FabF_Act; d_C10_FabF_Act; d_C12_FabF_Act; d_C14_FabF_Act; d_C16_FabF_Act; d_C18_FabF_Act;...
 d_C12_FabF_Act_un; d_C14_FabF_Act_un; d_C16_FabF_Act_un; d_C18_FabF_Act_un; d_C7_FabF_Act_MalACP; d_C9_FabF_Act_MalACP; d_C11_FabF_Act_MalACP;...
 d_C13_FabF_Act_MalACP; d_C15_FabF_Act_MalACP; d_C17_FabF_Act_MalACP; d_C19_FabF_Act_MalACP; d_C21_FabF_Act_MalACP; d_C15_FabF_Act_MalACP_un;...
 d_C17_FabF_Act_MalACP_un; d_C19_FabF_Act_MalACP_un; d_C21_FabF_Act_MalACP_un; d_C4_FabB_AcACP; d_C6_FabB_AcACP; d_C8_FabB_AcACP; d_C10_FabB_AcACP;...
 d_C12_FabB_AcACP; d_C14_FabB_AcACP; d_C16_FabB_AcACP; d_C18_FabB_AcACP; d_C12_FabB_AcACP_un; d_C14_FabB_AcACP_un; d_C16_FabB_AcACP_un;...
 d_C18_FabB_AcACP_un; d_C4_FabB_Act; d_C6_FabB_Act; d_C8_FabB_Act; d_C10_FabB_Act; d_C12_FabB_Act; d_C14_FabB_Act; d_C16_FabB_Act; d_C18_FabB_Act;...
 d_C12_FabB_Act_un; d_C14_FabB_Act_un; d_C16_FabB_Act_un; d_C18_FabB_Act_un; d_C7_FabB_Act_MalACP; d_C9_FabB_Act_MalACP; d_C11_FabB_Act_MalACP;...
 d_C13_FabB_Act_MalACP; d_C15_FabB_Act_MalACP; d_C17_FabB_Act_MalACP; d_C19_FabB_Act_MalACP; d_C21_FabB_Act_MalACP; d_C15_FabB_Act_MalACP_un;...
 d_C17_FabB_Act_MalACP_un; d_C19_FabB_Act_MalACP_un; d_C21_FabB_Act_MalACP_un; d_C10_FabB_cis3EnAcACP; d_C10_FabB_Act_cis3; d_C13_FabB_Act_cis3MalACP;...
 d_C4_FabH_AcACP; d_C6_FabH_AcACP; d_C8_FabH_AcACP; d_C10_FabH_AcACP; d_C12_FabH_AcACP; d_C14_FabH_AcACP; d_C16_FabH_AcACP; d_C18_FabH_AcACP;...
 d_C20_FabH_AcACP; d_C12_FabH_AcACP_un; d_C14_FabH_AcACP_un; d_C16_FabH_AcACP_un; d_C18_FabH_AcACP_un; d_C20_FabH_AcACP_un; d_C6_FabH_Act_AcACP;...
 d_C8_FabH_Act_AcACP; d_C10_FabH_Act_AcACP; d_C12_FabH_Act_AcACP; d_C14_FabH_Act_AcACP; d_C16_FabH_Act_AcACP; d_C18_FabH_Act_AcACP; d_C20_FabH_Act_AcACP;...
 d_C22_FabH_Act_AcACP; d_C14_FabH_Act_AcACP_un; d_C16_FabH_Act_AcACP_un; d_C18_FabH_Act_AcACP_un; d_C20_FabH_Act_AcACP_un; d_C22_FabH_Act_AcACP_un;...
 d_TesA_ACP; d_FabH_ACP; d_FabG_ACP; d_FabZ_ACP; d_FabI_ACP; d_FabF_ACP; d_FabA_ACP; d_FabB_ACP; d_C2_FabB_AcCoA; d_C2_FabB_Act; d_C5_FabB_Act_MalACP;...
 d_C2_FabF_AcCoA; d_C2_FabF_Act; d_C5_FabF_Act_MalACP; d_C3_FabB_MalACP; d_C2_AcACP; d_C2_FabB_AcACP; d_C3_FabF_MalACP; d_C2_FabF_AcACP;d_C2_FabH_AcACP;d_C3_FabH_MalACP];


