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

% ACC % changed for ACC
% BC (AccC)
c_ACC_C = P.ACC_Ctot - c_BC_ATP - c_C1_BC_ATP_HCO3 - c_C1_BC_Pi_HCO3 - c_C1_BC_Pi_HCO3_BCCP_Biotin;
% BCCP-Biotin (AccB-Biotin)
c_ACC_B = P.ACC_Btot - c_C1_BC_Pi_HCO3_BCCP_Biotin - c_C1_BCCP_Biotin_CO2 - c_C1_CT_BCCP_Biotin_CO2;
% CT (AccAD)
c_ACC_AD = P.ACC_ADtot - c_C1_CT_BCCP_Biotin_CO2 - c_C1_CT_Act - c_C3_CT_Act_AcCoA;

% MCMT
c_MCMT = P.MCMTtot - c_C3_MCMT_MalCoA - c_C3_MCMT_Act - c_C3_MCMT_Act_ACP;

% KASIII
c_KASIII = P.KASIIItot - c_KASIII_ACP...
    - c_C2_KASIII_CoA - c_C4_KASIII_CoA - c_C6_KASIII_CoA - c_C8_KASIII_CoA - c_C10_KASIII_CoA - c_C12_KASIII_CoA - c_C14_KASIII_CoA - c_C16_KASIII_CoA - c_C18_KASIII_CoA...
    - c_C2_KASIII_Act - c_C4_KASIII_Act - c_C6_KASIII_Act - c_C8_KASIII_Act - c_C10_KASIII_Act - c_C12_KASIII_Act - c_C14_KASIII_Act - c_C16_KASIII_Act - c_C18_KASIII_Act...
    - c_C5_KASIII_Act_MalACP - c_C7_KASIII_Act_MalACP - c_C9_KASIII_Act_MalACP - c_C11_KASIII_Act_MalACP - c_C13_KASIII_Act_MalACP - c_C15_KASIII_Act_MalACP - c_C17_KASIII_Act_MalACP - c_C19_KASIII_Act_MalACP - c_C21_KASIII_Act_MalACP...
    - c_C4_KASIII_AcACP - c_C6_KASIII_AcACP - c_C8_KASIII_AcACP - c_C10_KASIII_AcACP - c_C12_KASIII_AcACP - c_C14_KASIII_AcACP - c_C16_KASIII_AcACP - c_C18_KASIII_AcACP - c_C20_KASIII_AcACP...
    - c_C12_KASIII_AcACP_un - c_C14_KASIII_AcACP_un - c_C16_KASIII_AcACP_un - c_C18_KASIII_AcACP_un - c_C20_KASIII_AcACP_un...
    - c_C6_KASIII_Act_AcACP - c_C8_KASIII_Act_AcACP - c_C10_KASIII_Act_AcACP - c_C12_KASIII_Act_AcACP - c_C14_KASIII_Act_AcACP - c_C16_KASIII_Act_AcACP - c_C18_KASIII_Act_AcACP - c_C20_KASIII_Act_AcACP - c_C22_KASIII_Act_AcACP...
    - c_C14_KASIII_Act_AcACP_un - c_C16_KASIII_Act_AcACP_un - c_C18_KASIII_Act_AcACP_un - c_C20_KASIII_Act_AcACP_un - c_C22_KASIII_Act_AcACP_un;

% KAR
c_KAR = P.KARtot - c_KAR_NADPH...
    - c_C4_KAR_NADPH_BKeACP - c_C6_KAR_NADPH_BKeACP - c_C8_KAR_NADPH_BKeACP - c_C10_KAR_NADPH_BKeACP - c_C12_KAR_NADPH_BKeACP - c_C14_KAR_NADPH_BKeACP - c_C16_KAR_NADPH_BKeACP - c_C18_KAR_NADPH_BKeACP - c_C20_KAR_NADPH_BKeACP...
    - c_C12_KAR_NADPH_BKeACP_un - c_C14_KAR_NADPH_BKeACP_un - c_C16_KAR_NADPH_BKeACP_un - c_C18_KAR_NADPH_BKeACP_un - c_C20_KAR_NADPH_BKeACP_un - c_KAR_ACP;

% HAD
c_HAD = P.HADtot - c_HAD_ACP...
    - c_C4_HAD_BHyAcACP - c_C6_HAD_BHyAcACP - c_C8_HAD_BHyAcACP - c_C10_HAD_BHyAcACP - c_C12_HAD_BHyAcACP - c_C14_HAD_BHyAcACP - c_C16_HAD_BHyAcACP - c_C18_HAD_BHyAcACP - c_C20_HAD_BHyAcACP...
    - c_C12_HAD_BHyAcACP_un - c_C14_HAD_BHyAcACP_un - c_C16_HAD_BHyAcACP_un - c_C18_HAD_BHyAcACP_un - c_C20_HAD_BHyAcACP_un...
    - c_C4_HAD_EnAcACP - c_C6_HAD_EnAcACP - c_C8_HAD_EnAcACP - c_C10_HAD_EnAcACP - c_C12_HAD_EnAcACP - c_C14_HAD_EnAcACP - c_C16_HAD_EnAcACP - c_C18_HAD_EnAcACP - c_C20_HAD_EnAcACP...
    - c_C12_HAD_EnAcACP_un - c_C14_HAD_EnAcACP_un - c_C16_HAD_EnAcACP_un - c_C18_HAD_EnAcACP_un - c_C20_HAD_EnAcACP_un;

% ER
c_ER = P.ERtot - c_ER_NADH...
    - c_C4_ER_NADH_EnAcACP - c_C6_ER_NADH_EnAcACP - c_C8_ER_NADH_EnAcACP - c_C10_ER_NADH_EnAcACP - c_C12_ER_NADH_EnAcACP - c_C14_ER_NADH_EnAcACP - c_C16_ER_NADH_EnAcACP - c_C18_ER_NADH_EnAcACP - c_C20_ER_NADH_EnAcACP...
    - c_C12_ER_NADH_EnAcACP_un - c_C14_ER_NADH_EnAcACP_un - c_C16_ER_NADH_EnAcACP_un - c_C18_ER_NADH_EnAcACP_un - c_C20_ER_NADH_EnAcACP_un - c_ER_ACP;

% FatA
c_FatA = P.FatAtot- c_FatA_ACP...
    - c_C4_FatA_AcACP - c_C6_FatA_AcACP - c_C8_FatA_AcACP - c_C10_FatA_AcACP - c_C12_FatA_AcACP - c_C14_FatA_AcACP - c_C16_FatA_AcACP - c_C18_FatA_AcACP - c_C20_FatA_AcACP...
    - c_C12_FatA_AcACP_un - c_C14_FatA_AcACP_un - c_C16_FatA_AcACP_un - c_C18_FatA_AcACP_un - c_C20_FatA_AcACP_un;

% KASI
c_KASI = P.KASItot - c_KASI_ACP - c_C2_KASI_AcCoA - c_C2_KASI_Act - c_C3_KASI_MalACP - c_C2_KASI_AcACP...
    - c_C4_KASI_AcACP - c_C6_KASI_AcACP - c_C8_KASI_AcACP - c_C10_KASI_AcACP - c_C12_KASI_AcACP - c_C14_KASI_AcACP - c_C16_KASI_AcACP - c_C18_KASI_AcACP...
    - c_C12_KASI_AcACP_un - c_C14_KASI_AcACP_un - c_C16_KASI_AcACP_un - c_C18_KASI_AcACP_un...
    - c_C4_KASI_Act - c_C6_KASI_Act - c_C8_KASI_Act - c_C10_KASI_Act - c_C12_KASI_Act - c_C14_KASI_Act - c_C16_KASI_Act - c_C18_KASI_Act...
    - c_C12_KASI_Act_un - c_C14_KASI_Act_un - c_C16_KASI_Act_un - c_C18_KASI_Act_un...
    - c_C5_KASI_Act_MalACP - c_C7_KASI_Act_MalACP - c_C9_KASI_Act_MalACP - c_C11_KASI_Act_MalACP - c_C13_KASI_Act_MalACP - c_C15_KASI_Act_MalACP - c_C17_KASI_Act_MalACP - c_C19_KASI_Act_MalACP - c_C21_KASI_Act_MalACP...
    - c_C15_KASI_Act_MalACP_un - c_C17_KASI_Act_MalACP_un - c_C19_KASI_Act_MalACP_un - c_C21_KASI_Act_MalACP_un;

% SAD
c_SAD = P.SADtot - c_SAD_ACP - c_C10_SAD_cis3EnAcACP...
    - c_C4_SAD_BHyAcACP - c_C6_SAD_BHyAcACP - c_C8_SAD_BHyAcACP - c_C10_SAD_BHyAcACP - c_C12_SAD_BHyAcACP - c_C14_SAD_BHyAcACP - c_C16_SAD_BHyAcACP - c_C18_SAD_BHyAcACP - c_C20_SAD_BHyAcACP...
    - c_C12_SAD_BHyAcACP_un - c_C14_SAD_BHyAcACP_un - c_C16_SAD_BHyAcACP_un - c_C18_SAD_BHyAcACP_un - c_C20_SAD_BHyAcACP_un...
    - c_C4_SAD_EnAcACP - c_C6_SAD_EnAcACP - c_C8_SAD_EnAcACP - c_C10_SAD_EnAcACP - c_C12_SAD_EnAcACP - c_C14_SAD_EnAcACP - c_C16_SAD_EnAcACP - c_C18_SAD_EnAcACP - c_C20_SAD_EnAcACP...
    - c_C12_SAD_EnAcACP_un - c_C14_SAD_EnAcACP_un - c_C16_SAD_EnAcACP_un - c_C18_SAD_EnAcACP_un - c_C20_SAD_EnAcACP_un;

% KASII
c_KASII = P.KASIItot - c_KASII_ACP - c_C2_KASII_AcCoA - c_C2_KASII_Act - c_C5_KASII_Act_MalACP - c_C3_KASII_MalACP - c_C2_KASII_AcACP...
    - c_C4_KASII_AcACP - c_C6_KASII_AcACP - c_C8_KASII_AcACP - c_C10_KASII_AcACP - c_C12_KASII_AcACP - c_C14_KASII_AcACP - c_C16_KASII_AcACP - c_C18_KASII_AcACP...
    - c_C12_KASII_AcACP_un - c_C14_KASII_AcACP_un - c_C16_KASII_AcACP_un - c_C18_KASII_AcACP_un...
    - c_C4_KASII_Act - c_C6_KASII_Act - c_C8_KASII_Act - c_C10_KASII_Act - c_C12_KASII_Act - c_C14_KASII_Act - c_C16_KASII_Act - c_C18_KASII_Act...
    - c_C12_KASII_Act_un - c_C14_KASII_Act_un - c_C16_KASII_Act_un - c_C18_KASII_Act_un...
    - c_C7_KASII_Act_MalACP - c_C9_KASII_Act_MalACP - c_C11_KASII_Act_MalACP - c_C13_KASII_Act_MalACP - c_C15_KASII_Act_MalACP - c_C17_KASII_Act_MalACP - c_C19_KASII_Act_MalACP - c_C21_KASII_Act_MalACP...
    - c_C15_KASII_Act_MalACP_un - c_C17_KASII_Act_MalACP_un - c_C19_KASII_Act_MalACP_un - c_C21_KASII_Act_MalACP_un...
    - c_C10_KASII_cis3EnAcACP - c_C10_KASII_Act_cis3 - c_C13_KASII_Act_cis3MalACP;

% Set of differential equations
% ATP % changed for ACC
d_ATP = P.k1_1r.*c_BC_ATP - P.k1_1f.*c_ACC_C.*c_ATP;

% Bicarbonate % changed for ACC
d_C1_Bicarbonate = P.k1_2r.*c_C1_BC_ATP_HCO3 - P.k1_2f.*c_BC_ATP.*c_C1_Bicarbonate;

% C2n (n=1:9)-CoA % changed for ACC
d_C2_AcCoA = P.k1_5r.*c_C3_CT_Act_AcCoA - P.k1_5f.*c_C1_CT_Act.*c_C2_AcCoA + P.k3_1r(1).*c_C2_KASIII_CoA - P.k3_1f(1).*c_KASIII.*c_C2_AcCoA + P.k8_4r.*c_C2_KASI_AcCoA - P.k8_4f.*c_KASI.*c_C2_AcCoA;
d_C4_SucCoA = P.k3_1r(2).*c_C4_KASIII_CoA - P.k3_1f(2).*c_KASIII.*c_C4_SucCoA;
d_C6_HexCoA = P.k3_1r(3).*c_C6_KASIII_CoA - P.k3_1f(3).*c_KASIII.*c_C6_HexCoA;
d_C8_OcCoA = P.k3_1r(4).*c_C8_KASIII_CoA - P.k3_1f(4).*c_KASIII.*c_C8_OcCoA;
d_C10_DecCoA = P.k3_1r(5).*c_C10_KASIII_CoA - P.k3_1f(5).*c_KASIII.*c_C10_DecCoA;
d_C12_LauCoA = P.k3_1r(6).*c_C12_KASIII_CoA - P.k3_1f(6).*c_KASIII.*c_C12_LauCoA;
d_C14_EthCoA = P.k3_1r(7).*c_C14_KASIII_CoA - P.k3_1f(7).*c_KASIII.*c_C14_EthCoA;
d_C16_PalCoA = P.k3_1r(8).*c_C16_KASIII_CoA - P.k3_1f(8).*c_KASIII.*c_C16_PalCoA;
d_C18_OcDecCoA = 0*c_C18_OcDecCoA;

% ACP
d_ACP = P.k2_3r.*c_C3_MCMT_Act_ACP - P.k2_3f.*c_C3_MCMT_Act.*c_ACP...
    + P.kcat7(1).*c_C4_FatA_AcACP + P.kcat7(2).*c_C6_FatA_AcACP + P.kcat7(3).*c_C8_FatA_AcACP + P.kcat7(4).*c_C10_FatA_AcACP...
    + P.kcat7(5).*c_C12_FatA_AcACP + P.kcat7(6).*c_C14_FatA_AcACP + P.kcat7(7).*c_C16_FatA_AcACP...
    + P.k8_2f(1).*c_C4_KASI_AcACP - P.k8_2r(1).*c_C4_KASI_Act.*c_ACP...
    + P.k8_2f(2).*c_C6_KASI_AcACP - P.k8_2r(2).*c_C6_KASI_Act.*c_ACP...
    + P.k8_2f(3).*c_C8_KASI_AcACP - P.k8_2r(3).*c_C8_KASI_Act.*c_ACP...
    + P.k8_2f(4).*c_C10_KASI_AcACP - P.k8_2r(4).*c_C10_KASI_Act.*c_ACP...
    + P.k8_2f(5).*c_C12_KASI_AcACP - P.k8_2r(5).*c_C12_KASI_Act.*c_ACP...
    + P.k8_2f(6).*c_C14_KASI_AcACP - P.k8_2r(6).*c_C14_KASI_Act.*c_ACP...
    + P.k8_2f(7).*c_C16_KASI_AcACP - P.k8_2r(7).*c_C16_KASI_Act.*c_ACP...
    + P.k3_inh_r.*c_KASIII_ACP - P.k3_inh_f.*c_KASIII.*c_ACP...
    + P.k4_inh_r.*c_KAR_ACP - P.k4_inh_f.*c_KAR.*c_ACP...
    + P.k5_inh_r.*c_HAD_ACP - P.k5_inh_f.*c_HAD.*c_ACP...
    + P.k6_inh_r.*c_ER_ACP - P.k6_inh_f.*c_ER.*c_ACP...
    + P.k7_inh_r.*c_FatA_ACP - P.k7_inh_f.*c_FatA.*c_ACP...
    + P.k8_inh_r.*c_KASI_ACP - P.k8_inh_f.*c_KASI.*c_ACP...
    + P.k8_9f.*c_C2_KASI_AcACP - P.k8_9r.*c_C2_KASI_Act.*c_ACP; 

% NADPH
d_NADPH = P.k4_1r(1).*c_KAR_NADPH - P.k4_1f(1).*c_KAR.*c_NADPH; 

% NADP+
d_NADP = P.kcat4(1).*c_C4_KAR_NADPH_BKeACP + P.kcat4(2).*c_C6_KAR_NADPH_BKeACP + P.kcat4(3).*c_C8_KAR_NADPH_BKeACP + P.kcat4(4).*c_C10_KAR_NADPH_BKeACP...
    + P.kcat4(5).*c_C12_KAR_NADPH_BKeACP + P.kcat4(6).*c_C14_KAR_NADPH_BKeACP + P.kcat4(7).*c_C16_KAR_NADPH_BKeACP;

% NADH
d_NADH = P.k6_1r(1).*c_ER_NADH - P.k6_1f(1).*c_ER.*c_NADH; 

% NAD+
d_NAD = P.kcat6(1).*c_C4_ER_NADH_EnAcACP + P.kcat6(2).*c_C6_ER_NADH_EnAcACP + P.kcat6(3).*c_C8_ER_NADH_EnAcACP + P.kcat6(4).*c_C10_ER_NADH_EnAcACP...
    + P.kcat6(5).*c_C12_ER_NADH_EnAcACP + P.kcat6(6).*c_C14_ER_NADH_EnAcACP + P.kcat6(7).*c_C16_ER_NADH_EnAcACP; 

% ADP % changed for ACC
d_ADP = P.kcat1_1.*c_C1_BC_ATP_HCO3;

% Malonyl-CoA % changed for ACC
d_C3_MalCoA = P.kcat1_4.*c_C3_CT_Act_AcCoA + P.k2_1r.*c_C3_MCMT_MalCoA - P.k2_1f.*c_MCMT.*c_C3_MalCoA; 


% CoA
d_CoA = P.k2_2f.*c_C3_MCMT_MalCoA - P.k2_2r.*c_C3_MCMT_Act.*c_CoA...
    + P.k3_2f(1).*c_C2_KASIII_CoA - P.k3_2r(1).*c_C2_KASIII_Act.*c_CoA...
    + P.k3_2f(2).*c_C4_KASIII_CoA - P.k3_2r(2).*c_C4_KASIII_Act.*c_CoA...
    + P.k3_2f(3).*c_C6_KASIII_CoA - P.k3_2r(3).*c_C6_KASIII_Act.*c_CoA...
    + P.k3_2f(4).*c_C8_KASIII_CoA - P.k3_2r(4).*c_C8_KASIII_Act.*c_CoA...
    + P.k3_2f(5).*c_C10_KASIII_CoA - P.k3_2r(5).*c_C10_KASIII_Act.*c_CoA...
    + P.k3_2f(6).*c_C12_KASIII_CoA - P.k3_2r(6).*c_C12_KASIII_Act.*c_CoA...
    + P.k3_2f(7).*c_C14_KASIII_CoA - P.k3_2r(7).*c_C14_KASIII_Act.*c_CoA...
    + P.k3_2f(8).*c_C16_KASIII_CoA - P.k3_2r(8).*c_C16_KASIII_Act.*c_CoA...
    + P.k8_5f.*c_C2_KASI_AcCoA - P.k8_5r*c_C2_KASI_Act.*c_CoA;

% Malonyl-ACP
d_C3_MalACP = P.k2_4f.*c_C3_MCMT_Act_ACP - P.k2_4r.*c_MCMT.*c_C3_MalACP...
    + P.k3_3r(1).*c_C5_KASIII_Act_MalACP - P.k3_3f(1).*c_C2_KASIII_Act.*c_C3_MalACP...
    + P.k3_3r(2).*c_C7_KASIII_Act_MalACP - P.k3_3f(2).*c_C4_KASIII_Act.*c_C3_MalACP...
    + P.k3_3r(3).*c_C9_KASIII_Act_MalACP - P.k3_3f(3).*c_C6_KASIII_Act.*c_C3_MalACP...
    + P.k3_3r(4).*c_C11_KASIII_Act_MalACP - P.k3_3f(4).*c_C8_KASIII_Act.*c_C3_MalACP...
    + P.k3_3r(5).*c_C13_KASIII_Act_MalACP - P.k3_3f(5).*c_C10_KASIII_Act.*c_C3_MalACP...
    + P.k3_3r(6).*c_C15_KASIII_Act_MalACP - P.k3_3f(6).*c_C12_KASIII_Act.*c_C3_MalACP...
    + P.k3_3r(7).*c_C17_KASIII_Act_MalACP - P.k3_3f(7).*c_C14_KASIII_Act.*c_C3_MalACP...
    + P.k8_3r(1).*c_C7_KASI_Act_MalACP - P.k8_3f(1).*c_C4_KASI_Act.*c_C3_MalACP...
    + P.k8_3r(2).*c_C9_KASI_Act_MalACP - P.k8_3f(2).*c_C6_KASI_Act.*c_C3_MalACP...
    + P.k8_3r(3).*c_C11_KASI_Act_MalACP - P.k8_3f(3).*c_C8_KASI_Act.*c_C3_MalACP...
    + P.k8_3r(4).*c_C13_KASI_Act_MalACP - P.k8_3f(4).*c_C10_KASI_Act.*c_C3_MalACP...
    + P.k8_3r(5).*c_C15_KASI_Act_MalACP - P.k8_3f(5).*c_C12_KASI_Act.*c_C3_MalACP...
    + P.k8_3r(6).*c_C17_KASI_Act_MalACP - P.k8_3f(6).*c_C14_KASI_Act.*c_C3_MalACP...
    + P.k8_6r.*c_C5_KASI_Act_MalACP - P.k8_6f.*c_C2_KASI_Act.*c_C3_MalACP...
    + P.k8_7r.*c_C3_KASI_MalACP - P.k8_7f.*c_KASI.*c_C3_MalACP;

% CO2
d_C1_CO2 = P.kcat3(1).*c_C5_KASIII_Act_MalACP + P.kcat3(2).*c_C7_KASIII_Act_MalACP + P.kcat3(3).*c_C9_KASIII_Act_MalACP + P.kcat3(4).*c_C11_KASIII_Act_MalACP...
    + P.kcat3(5).*c_C13_KASIII_Act_MalACP + P.kcat3(6).*c_C15_KASIII_Act_MalACP + P.kcat3(7).*c_C17_KASIII_Act_MalACP...
    + P.kcat8(1).*c_C7_KASI_Act_MalACP + P.kcat8(2).*c_C9_KASI_Act_MalACP + P.kcat8(3).*c_C11_KASI_Act_MalACP + P.kcat8(4).*c_C13_KASI_Act_MalACP...
    + P.kcat8(5).*c_C15_KASI_Act_MalACP + P.kcat8(6).*c_C17_KASI_Act_MalACP...
    + P.kcat8_H.*c_C5_KASI_Act_MalACP + P.kcat8_CO2.*c_C3_KASI_MalACP;

% C2n (n=2:10) B-ketoacyl-ACPs (KASIII + KASI + KASII - KAR)
d_C4_BKeACP = P.kcat3(1).*c_C5_KASIII_Act_MalACP + P.kcat8_H.*c_C5_KASI_Act_MalACP + P.k4_2r(1).*c_C4_KAR_NADPH_BKeACP - P.k4_2f(1).*c_KAR_NADPH.*c_C4_BKeACP;
d_C6_BKeACP = P.kcat3(2).*c_C7_KASIII_Act_MalACP + P.kcat8(1).*c_C7_KASI_Act_MalACP + P.k4_2r(2).*c_C6_KAR_NADPH_BKeACP - P.k4_2f(2).*c_KAR_NADPH.*c_C6_BKeACP;
d_C8_BKeACP = P.kcat3(3).*c_C9_KASIII_Act_MalACP + P.kcat8(2).*c_C9_KASI_Act_MalACP + P.k4_2r(3).*c_C8_KAR_NADPH_BKeACP - P.k4_2f(3).*c_KAR_NADPH.*c_C8_BKeACP;
d_C10_BKeACP = P.kcat3(4).*c_C11_KASIII_Act_MalACP + P.kcat8(3).*c_C11_KASI_Act_MalACP + P.k4_2r(4).*c_C10_KAR_NADPH_BKeACP - P.k4_2f(4).*c_KAR_NADPH.*c_C10_BKeACP;
d_C12_BKeACP = P.kcat3(5).*c_C13_KASIII_Act_MalACP + P.kcat8(4).*c_C13_KASI_Act_MalACP + P.k4_2r(5).*c_C12_KAR_NADPH_BKeACP - P.k4_2f(5).*c_KAR_NADPH.*c_C12_BKeACP;
d_C14_BKeACP = P.kcat3(6).*c_C15_KASIII_Act_MalACP + P.kcat8(5).*c_C15_KASI_Act_MalACP + P.k4_2r(6).*c_C14_KAR_NADPH_BKeACP - P.k4_2f(6).*c_KAR_NADPH.*c_C14_BKeACP;
d_C16_BKeACP = P.kcat3(7).*c_C17_KASIII_Act_MalACP + P.kcat8(6).*c_C17_KASI_Act_MalACP + P.k4_2r(7).*c_C16_KAR_NADPH_BKeACP - P.k4_2f(7).*c_KAR_NADPH.*c_C16_BKeACP;
d_C18_BKeACP = 0*c_C18_BKeACP;
d_C20_BKeACP = 0*c_C20_BKeACP;

% C2n:1 (n=6:10) B-ketoacyl-ACPs (KASI + KASII - KAR)
d_C12_BKeACP_un = 0*c_C12_BKeACP_un;
d_C14_BKeACP_un = 0*c_C14_BKeACP_un;
d_C16_BKeACP_un = 0*c_C16_BKeACP_un;
d_C18_BKeACP_un = 0*c_C18_BKeACP_un;
d_C20_BKeACP_un = 0*c_C20_BKeACP_un;

% C2n (n=2:10) B-hydroxy-acyl-ACPs (KAR - HAD - SAD)
d_C4_BHyAcACP = P.kcat4(1).*c_C4_KAR_NADPH_BKeACP + P.k5_1r(1).*c_C4_HAD_BHyAcACP - P.k5_1f(1).*c_HAD.*c_C4_BHyAcACP;
d_C6_BHyAcACP = P.kcat4(2).*c_C6_KAR_NADPH_BKeACP + P.k5_1r(2).*c_C6_HAD_BHyAcACP - P.k5_1f(2).*c_HAD.*c_C6_BHyAcACP;
d_C8_BHyAcACP = P.kcat4(3).*c_C8_KAR_NADPH_BKeACP + P.k5_1r(3).*c_C8_HAD_BHyAcACP - P.k5_1f(3).*c_HAD.*c_C8_BHyAcACP;
d_C10_BHyAcACP = P.kcat4(4).*c_C10_KAR_NADPH_BKeACP + P.k5_1r(4).*c_C10_HAD_BHyAcACP - P.k5_1f(4).*c_HAD.*c_C10_BHyAcACP;
d_C12_BHyAcACP = P.kcat4(5).*c_C12_KAR_NADPH_BKeACP + P.k5_1r(5).*c_C12_HAD_BHyAcACP - P.k5_1f(5).*c_HAD.*c_C12_BHyAcACP;
d_C14_BHyAcACP = P.kcat4(6).*c_C14_KAR_NADPH_BKeACP + P.k5_1r(6).*c_C14_HAD_BHyAcACP - P.k5_1f(6).*c_HAD.*c_C14_BHyAcACP;
d_C16_BHyAcACP = P.kcat4(7).*c_C16_KAR_NADPH_BKeACP + P.k5_1r(7).*c_C16_HAD_BHyAcACP - P.k5_1f(7).*c_HAD.*c_C16_BHyAcACP;
d_C18_BHyAcACP = 0*c_C18_BHyAcACP;
d_C20_BHyAcACP = 0*c_C20_BHyAcACP;

% C2n:1 (n=6:10) B-hydroxy-acyl-ACPs (KAR - HAD - SAD)
d_C12_BHyAcACP_un = 0*c_C12_BHyAcACP_un;
d_C14_BHyAcACP_un = 0*c_C14_BHyAcACP_un;
d_C16_BHyAcACP_un = 0*c_C16_BHyAcACP_un;
d_C18_BHyAcACP_un = 0*c_C18_BHyAcACP_un;
d_C20_BHyAcACP_un = 0*c_C20_BHyAcACP_un;

% C2n (n=2:10) Enoyl-Acyl-ACPs (HAD + SAD - ER) 
d_C4_EnAcACP = P.k5_3f(1).*c_C4_HAD_EnAcACP - P.k5_3r(1).*c_HAD.*c_C4_EnAcACP + P.k6_2r(1).*c_C4_ER_NADH_EnAcACP - P.k6_2f(1).*c_ER_NADH.*c_C4_EnAcACP;
d_C6_EnAcACP = P.k5_3f(2).*c_C6_HAD_EnAcACP - P.k5_3r(2).*c_HAD.*c_C6_EnAcACP + P.k6_2r(2).*c_C6_ER_NADH_EnAcACP - P.k6_2f(2).*c_ER_NADH.*c_C6_EnAcACP;
d_C8_EnAcACP = P.k5_3f(3).*c_C8_HAD_EnAcACP - P.k5_3r(3).*c_HAD.*c_C8_EnAcACP + P.k6_2r(3).*c_C8_ER_NADH_EnAcACP - P.k6_2f(3).*c_ER_NADH.*c_C8_EnAcACP;
d_C10_EnAcACP = P.k5_3f(4).*c_C10_HAD_EnAcACP - P.k5_3r(4).*c_HAD.*c_C10_EnAcACP + P.k6_2r(4).*c_C10_ER_NADH_EnAcACP - P.k6_2f(4).*c_ER_NADH.*c_C10_EnAcACP;
d_C12_EnAcACP = P.k5_3f(5).*c_C12_HAD_EnAcACP - P.k5_3r(5).*c_HAD.*c_C12_EnAcACP + P.k6_2r(5).*c_C12_ER_NADH_EnAcACP - P.k6_2f(5).*c_ER_NADH.*c_C12_EnAcACP;
d_C14_EnAcACP = P.k5_3f(6).*c_C14_HAD_EnAcACP - P.k5_3r(6).*c_HAD.*c_C14_EnAcACP + P.k6_2r(6).*c_C14_ER_NADH_EnAcACP - P.k6_2f(6).*c_ER_NADH.*c_C14_EnAcACP;
d_C16_EnAcACP = P.k5_3f(7).*c_C16_HAD_EnAcACP - P.k5_3r(7).*c_HAD.*c_C16_EnAcACP + P.k6_2r(7).*c_C16_ER_NADH_EnAcACP - P.k6_2f(7).*c_ER_NADH.*c_C16_EnAcACP;
d_C18_EnAcACP = 0*c_C18_EnAcACP;
d_C20_EnAcACP = 0*c_C20_EnAcACP;

% C10 cis-3-Enoyl-Acyl-ACP (SAD - KASII)
d_C10_cis3EnAcACP = 0*c_C10_cis3EnAcACP;

% C2n:1 (n=6:10) Enoyl-Acyl-ACPs (HAD + SAD - ER)
d_C12_EnAcACP_un = 0*c_C12_EnAcACP_un;
d_C14_EnAcACP_un = 0*c_C14_EnAcACP_un;
d_C16_EnAcACP_un = 0*c_C16_EnAcACP_un;
d_C18_EnAcACP_un = 0*c_C18_EnAcACP_un;
d_C20_EnAcACP_un = 0*c_C20_EnAcACP_un;

% C2n (n=2:9) Acyl-ACPs (ER - FatA - KASI - KASII - KASIII)
d_C4_AcACP = P.kcat6(1).*c_C4_ER_NADH_EnAcACP + P.k3_4r(1).*c_C4_KASIII_AcACP - P.k3_4f(1).*c_KASIII.*c_C4_AcACP + P.k3_5r(1).*c_C6_KASIII_Act_AcACP - P.k3_5f(1).*c_C2_KASIII_Act.*c_C4_AcACP + P.k7_1r(1).*c_C4_FatA_AcACP + P.k8_1r(1).*c_C4_KASI_AcACP - P.k7_1f(1).*c_FatA.*c_C4_AcACP - P.k8_1f(1).*c_KASI.*c_C4_AcACP;
d_C6_AcACP = P.kcat6(2).*c_C6_ER_NADH_EnAcACP + P.k3_4r(2).*c_C6_KASIII_AcACP - P.k3_4f(2).*c_KASIII.*c_C6_AcACP + P.k3_5r(2).*c_C8_KASIII_Act_AcACP - P.k3_5f(2).*c_C2_KASIII_Act.*c_C6_AcACP + P.k7_1r(2).*c_C6_FatA_AcACP + P.k8_1r(2).*c_C6_KASI_AcACP - P.k7_1f(2).*c_FatA.*c_C6_AcACP - P.k8_1f(2).*c_KASI.*c_C6_AcACP;
d_C8_AcACP = P.kcat6(3).*c_C8_ER_NADH_EnAcACP + P.k3_4r(3).*c_C8_KASIII_AcACP - P.k3_4f(3).*c_KASIII.*c_C8_AcACP + P.k3_5r(3).*c_C10_KASIII_Act_AcACP - P.k3_5f(3).*c_C2_KASIII_Act.*c_C8_AcACP + P.k7_1r(3).*c_C8_FatA_AcACP + P.k8_1r(3).*c_C8_KASI_AcACP - P.k7_1f(3).*c_FatA.*c_C8_AcACP - P.k8_1f(3).*c_KASI.*c_C8_AcACP;
d_C10_AcACP = P.kcat6(4).*c_C10_ER_NADH_EnAcACP + P.k3_4r(4).*c_C10_KASIII_AcACP - P.k3_4f(4).*c_KASIII.*c_C10_AcACP + P.k3_5r(4).*c_C12_KASIII_Act_AcACP - P.k3_5f(4).*c_C2_KASIII_Act.*c_C10_AcACP + P.k7_1r(4).*c_C10_FatA_AcACP + P.k8_1r(4).*c_C10_KASI_AcACP - P.k7_1f(4).*c_FatA.*c_C10_AcACP - P.k8_1f(4).*c_KASI.*c_C10_AcACP;
d_C12_AcACP = P.kcat6(5).*c_C12_ER_NADH_EnAcACP + P.k3_4r(5).*c_C12_KASIII_AcACP - P.k3_4f(5).*c_KASIII.*c_C12_AcACP + P.k3_5r(5).*c_C14_KASIII_Act_AcACP - P.k3_5f(5).*c_C2_KASIII_Act.*c_C12_AcACP + P.k7_1r(5).*c_C12_FatA_AcACP + P.k8_1r(5).*c_C12_KASI_AcACP - P.k7_1f(5).*c_FatA.*c_C12_AcACP - P.k8_1f(5).*c_KASI.*c_C12_AcACP;
d_C14_AcACP = P.kcat6(6).*c_C14_ER_NADH_EnAcACP + P.k3_4r(6).*c_C14_KASIII_AcACP - P.k3_4f(6).*c_KASIII.*c_C14_AcACP + P.k3_5r(6).*c_C16_KASIII_Act_AcACP - P.k3_5f(6).*c_C2_KASIII_Act.*c_C14_AcACP + P.k7_1r(6).*c_C14_FatA_AcACP + P.k8_1r(6).*c_C14_KASI_AcACP - P.k7_1f(6).*c_FatA.*c_C14_AcACP - P.k8_1f(6).*c_KASI.*c_C14_AcACP;
d_C16_AcACP = P.kcat6(7).*c_C16_ER_NADH_EnAcACP + P.k3_4r(7).*c_C16_KASIII_AcACP - P.k3_4f(7).*c_KASIII.*c_C16_AcACP + P.k7_1r(7).*c_C16_FatA_AcACP + P.k8_1r(7).*c_C16_KASI_AcACP - P.k7_1f(7).*c_FatA.*c_C16_AcACP - P.k8_1f(7).*c_KASI.*c_C16_AcACP;
d_C18_AcACP = 0*c_C18_AcACP;

% C2n (n=10) Acyl-ACPs (ER - FatA - KASIII)
d_C20_AcACP = 0*c_C20_AcACP;

% C2n:1 (n=6:9) Acyl-ACPs (ER - FatA - KASI - KASII - KASIII)
d_C12_AcACP_un = 0*c_C12_AcACP_un;
d_C14_AcACP_un = 0*c_C14_AcACP_un;
d_C16_AcACP_un = 0*c_C16_AcACP_un;
d_C18_AcACP_un = 0*c_C18_AcACP_un;

% C2n:1 (n=10) Acyl-ACPs (ER - FatA - KASIII)
d_C20_AcACP_un = 0*c_C20_AcACP_un;

% Fatty Acids
d_C4_FA = P.kcat7(1).*c_C4_FatA_AcACP;
d_C6_FA = P.kcat7(2).*c_C6_FatA_AcACP;
d_C8_FA = P.kcat7(3).*c_C8_FatA_AcACP;
d_C10_FA = P.kcat7(4).*c_C10_FatA_AcACP;
d_C12_FA = P.kcat7(5).*c_C12_FatA_AcACP;
d_C14_FA = P.kcat7(6).*c_C14_FatA_AcACP;
d_C16_FA = P.kcat7(7).*c_C16_FatA_AcACP;
d_C18_FA = 0*c_C18_FA;
d_C20_FA = 0*c_C20_FA;

% Fatty Acids (unsaturated)
d_C12_FA_un = 0*c_C12_FA_un;
d_C14_FA_un = 0*c_C14_FA_un;
d_C16_FA_un = 0*c_C16_FA_un;
d_C18_FA_un = 0*c_C18_FA_un;
d_C20_FA_un = 0*c_C20_FA_un;

% BC-ATP k1_1 % changed for ACC
d_BC_ATP = P.k1_1f.*c_ACC_C.*c_ATP - P.k1_1r.*c_BC_ATP + P.k1_2r.*c_C1_BC_ATP_HCO3 - P.k1_2f.*c_BC_ATP.*c_C1_Bicarbonate;

% BC-ATP-HCO3 k1_2 % changed for ACC
d_C1_BC_ATP_HCO3 = P.k1_2f.*c_BC_ATP.*c_C1_Bicarbonate - P.k1_2r.*c_C1_BC_ATP_HCO3 - P.kcat1_1.*c_C1_BC_ATP_HCO3;

% BC-Pi-HCO3 kcat1_1 % changed for ACC
d_C1_BC_Pi_HCO3 = P.kcat1_1.*c_C1_BC_ATP_HCO3 + P.k1_3r.*c_C1_BC_Pi_HCO3_BCCP_Biotin - P.k1_3f.*c_ACC_B.*c_C1_BC_Pi_HCO3;

% BC-Pi-HCO3-BCCP-Biotin k1_3 % changed for ACC
d_C1_BC_Pi_HCO3_BCCP_Biotin = P.k1_3f.*c_ACC_B.*c_C1_BC_Pi_HCO3 - P.k1_3r.*c_C1_BC_Pi_HCO3_BCCP_Biotin - P.kcat1_2.*c_C1_BC_Pi_HCO3_BCCP_Biotin;

% BCCP-Biotin-CO2 kcat1_2 % changed for ACC
d_C1_BCCP_Biotin_CO2 = P.kcat1_2.*c_C1_BC_Pi_HCO3_BCCP_Biotin + P.k1_4r.*c_C1_CT_BCCP_Biotin_CO2 - P.k1_4f.*c_ACC_AD.*c_C1_BCCP_Biotin_CO2;

% CT-BCCP-Biotin-CO2 k1_4 % changed for ACC
d_C1_CT_BCCP_Biotin_CO2 = P.k1_4f.*c_ACC_AD.*c_C1_BCCP_Biotin_CO2 - P.k1_4r.*c_C1_CT_BCCP_Biotin_CO2 - P.kcat1_3.*c_C1_CT_BCCP_Biotin_CO2;

% CT* kcat1_3 % changed for ACC
d_C1_CT_Act = P.kcat1_3.*c_C1_CT_BCCP_Biotin_CO2 + P.k1_5r.*c_C3_CT_Act_AcCoA - P.k1_5f.*c_C1_CT_Act.*c_C2_AcCoA;

% CT*-AcCoA k1_5 % changed for ACC
d_C3_CT_Act_AcCoA = P.k1_5f.*c_C1_CT_Act.*c_C2_AcCoA - P.k1_5r.*c_C3_CT_Act_AcCoA - P.kcat1_4.*c_C3_CT_Act_AcCoA;

% MCMT-Malonyl-CoA
d_C3_MCMT_MalCoA = P.k2_1f.*c_MCMT.*c_C3_MalCoA - P.k2_1r.*c_C3_MCMT_MalCoA + P.k2_2r.*c_C3_MCMT_Act.*c_CoA - P.k2_2f.*c_C3_MCMT_MalCoA; 

% MCMT*
d_C3_MCMT_Act = P.k2_2f.*c_C3_MCMT_MalCoA - P.k2_2r.*c_C3_MCMT_Act.*c_CoA + P.k2_3r.*c_C3_MCMT_Act_ACP - P.k2_3f.*c_C3_MCMT_Act.*c_ACP;

% MCMT*-ACP
d_C3_MCMT_Act_ACP = P.k2_3f.*c_C3_MCMT_Act.*c_ACP - P.k2_3r.*c_C3_MCMT_Act_ACP + P.k2_4r.*c_MCMT.*c_C3_MalACP - P.k2_4f.*c_C3_MCMT_Act_ACP;

% C2n (n=1:9) KASIII-CoA
d_C2_KASIII_CoA = P.k3_1f(1).*c_KASIII.*c_C2_AcCoA - P.k3_1r(1).*c_C2_KASIII_CoA + P.k3_2r(1).*c_C2_KASIII_Act.*c_CoA - P.k3_2f(1).*c_C2_KASIII_CoA; 
d_C4_KASIII_CoA = P.k3_1f(2).*c_KASIII.*c_C4_SucCoA - P.k3_1r(2).*c_C4_KASIII_CoA + P.k3_2r(2).*c_C4_KASIII_Act.*c_CoA - P.k3_2f(2).*c_C4_KASIII_CoA; 
d_C6_KASIII_CoA = P.k3_1f(3).*c_KASIII.*c_C6_HexCoA - P.k3_1r(3).*c_C6_KASIII_CoA + P.k3_2r(3).*c_C6_KASIII_Act.*c_CoA - P.k3_2f(3).*c_C6_KASIII_CoA; 
d_C8_KASIII_CoA = P.k3_1f(4).*c_KASIII.*c_C8_OcCoA - P.k3_1r(4).*c_C8_KASIII_CoA + P.k3_2r(4).*c_C8_KASIII_Act.*c_CoA - P.k3_2f(4).*c_C8_KASIII_CoA; 
d_C10_KASIII_CoA = P.k3_1f(5).*c_KASIII.*c_C10_DecCoA - P.k3_1r(5).*c_C10_KASIII_CoA + P.k3_2r(5).*c_C10_KASIII_Act.*c_CoA - P.k3_2f(5).*c_C10_KASIII_CoA; 
d_C12_KASIII_CoA = P.k3_1f(6).*c_KASIII.*c_C12_LauCoA - P.k3_1r(6).*c_C12_KASIII_CoA + P.k3_2r(6).*c_C12_KASIII_Act.*c_CoA - P.k3_2f(6).*c_C12_KASIII_CoA;
d_C14_KASIII_CoA = P.k3_1f(7).*c_KASIII.*c_C14_EthCoA - P.k3_1r(7).*c_C14_KASIII_CoA + P.k3_2r(7).*c_C14_KASIII_Act.*c_CoA - P.k3_2f(7).*c_C14_KASIII_CoA; 
d_C16_KASIII_CoA = P.k3_1f(8).*c_KASIII.*c_C16_PalCoA - P.k3_1r(8).*c_C16_KASIII_CoA + P.k3_2r(8).*c_C16_KASIII_Act.*c_CoA - P.k3_2f(8).*c_C16_KASIII_CoA; 
d_C18_KASIII_CoA = 0*c_C18_KASIII_CoA; 

% C2n (n=1:9) KASIII*
% making KASIII* - using KASIII* - inhibition from Acyl ACPs (only Acetyl-CoA derived KASIII*)
d_C2_KASIII_Act = P.k3_2f(1).*c_C2_KASIII_CoA - P.k3_2r(1).*c_C2_KASIII_Act.*c_CoA... 
 + P.k3_3r(1).*c_C5_KASIII_Act_MalACP - P.k3_3f(1).*c_C2_KASIII_Act.*c_C3_MalACP...
 + P.k3_5r(1).*c_C6_KASIII_Act_AcACP - P.k3_5f(1).*c_C2_KASIII_Act.*c_C4_AcACP...
 + P.k3_5r(2).*c_C8_KASIII_Act_AcACP - P.k3_5f(2).*c_C2_KASIII_Act.*c_C6_AcACP...
 + P.k3_5r(3).*c_C10_KASIII_Act_AcACP - P.k3_5f(3).*c_C2_KASIII_Act.*c_C8_AcACP...
 + P.k3_5r(4).*c_C12_KASIII_Act_AcACP - P.k3_5f(4).*c_C2_KASIII_Act.*c_C10_AcACP...
 + P.k3_5r(5).*c_C14_KASIII_Act_AcACP - P.k3_5f(5).*c_C2_KASIII_Act.*c_C12_AcACP...
 + P.k3_5r(6).*c_C16_KASIII_Act_AcACP - P.k3_5f(6).*c_C2_KASIII_Act.*c_C14_AcACP;
d_C4_KASIII_Act = P.k3_2f(2).*c_C4_KASIII_CoA - P.k3_2r(2).*c_C4_KASIII_Act.*c_CoA + P.k3_3r(2).*c_C7_KASIII_Act_MalACP - P.k3_3f(2).*c_C4_KASIII_Act.*c_C3_MalACP;
d_C6_KASIII_Act = P.k3_2f(3).*c_C6_KASIII_CoA - P.k3_2r(3).*c_C6_KASIII_Act.*c_CoA + P.k3_3r(3).*c_C9_KASIII_Act_MalACP - P.k3_3f(3).*c_C6_KASIII_Act.*c_C3_MalACP;
d_C8_KASIII_Act = P.k3_2f(4).*c_C8_KASIII_CoA - P.k3_2r(4).*c_C8_KASIII_Act.*c_CoA + P.k3_3r(4).*c_C11_KASIII_Act_MalACP - P.k3_3f(4).*c_C8_KASIII_Act.*c_C3_MalACP;
d_C10_KASIII_Act = P.k3_2f(5).*c_C10_KASIII_CoA - P.k3_2r(5).*c_C10_KASIII_Act.*c_CoA + P.k3_3r(5).*c_C13_KASIII_Act_MalACP - P.k3_3f(5).*c_C10_KASIII_Act.*c_C3_MalACP;
d_C12_KASIII_Act = P.k3_2f(6).*c_C12_KASIII_CoA - P.k3_2r(6).*c_C12_KASIII_Act.*c_CoA + P.k3_3r(6).*c_C15_KASIII_Act_MalACP - P.k3_3f(6).*c_C12_KASIII_Act.*c_C3_MalACP;
d_C14_KASIII_Act = P.k3_2f(7).*c_C14_KASIII_CoA - P.k3_2r(7).*c_C14_KASIII_Act.*c_CoA + P.k3_3r(7).*c_C17_KASIII_Act_MalACP - P.k3_3f(7).*c_C14_KASIII_Act.*c_C3_MalACP;
d_C16_KASIII_Act = P.k3_2f(8).*c_C16_KASIII_CoA - P.k3_2r(8).*c_C16_KASIII_Act.*c_CoA;
d_C18_KASIII_Act = 0*c_C18_KASIII_Act;

% C2n (n=1:9) KASIII*-Malonyl-ACP
d_C5_KASIII_Act_MalACP = P.k3_3f(1).*c_C2_KASIII_Act.*c_C3_MalACP - P.k3_3r(1).*c_C5_KASIII_Act_MalACP - P.kcat3(1).*c_C5_KASIII_Act_MalACP; 
d_C7_KASIII_Act_MalACP = P.k3_3f(2).*c_C4_KASIII_Act.*c_C3_MalACP - P.k3_3r(2).*c_C7_KASIII_Act_MalACP - P.kcat3(2).*c_C7_KASIII_Act_MalACP; 
d_C9_KASIII_Act_MalACP = P.k3_3f(3).*c_C6_KASIII_Act.*c_C3_MalACP - P.k3_3r(3).*c_C9_KASIII_Act_MalACP - P.kcat3(3).*c_C9_KASIII_Act_MalACP; 
d_C11_KASIII_Act_MalACP = P.k3_3f(4).*c_C8_KASIII_Act.*c_C3_MalACP - P.k3_3r(4).*c_C11_KASIII_Act_MalACP - P.kcat3(4).*c_C11_KASIII_Act_MalACP; 
d_C13_KASIII_Act_MalACP = P.k3_3f(5).*c_C10_KASIII_Act.*c_C3_MalACP - P.k3_3r(5).*c_C13_KASIII_Act_MalACP - P.kcat3(5).*c_C13_KASIII_Act_MalACP; 
d_C15_KASIII_Act_MalACP = P.k3_3f(6).*c_C12_KASIII_Act.*c_C3_MalACP - P.k3_3r(6).*c_C15_KASIII_Act_MalACP - P.kcat3(6).*c_C15_KASIII_Act_MalACP; 
d_C17_KASIII_Act_MalACP = P.k3_3f(7).*c_C14_KASIII_Act.*c_C3_MalACP - P.k3_3r(7).*c_C17_KASIII_Act_MalACP - P.kcat3(7).*c_C17_KASIII_Act_MalACP; 
d_C19_KASIII_Act_MalACP = 0*c_C19_KASIII_Act_MalACP; 
d_C21_KASIII_Act_MalACP = 0*c_C21_KASIII_Act_MalACP; 

% KAR-NADPH
d_KAR_NADPH = P.k4_1f(1).*c_KAR.*c_NADPH - P.k4_1r(1).*c_KAR_NADPH...
 + P.k4_2r(1).*c_C4_KAR_NADPH_BKeACP - P.k4_2f(1).*c_KAR_NADPH.*c_C4_BKeACP...
 + P.k4_2r(2).*c_C6_KAR_NADPH_BKeACP - P.k4_2f(2).*c_KAR_NADPH.*c_C6_BKeACP...
 + P.k4_2r(3).*c_C8_KAR_NADPH_BKeACP - P.k4_2f(3).*c_KAR_NADPH.*c_C8_BKeACP...
 + P.k4_2r(4).*c_C10_KAR_NADPH_BKeACP - P.k4_2f(4).*c_KAR_NADPH.*c_C10_BKeACP...
 + P.k4_2r(5).*c_C12_KAR_NADPH_BKeACP - P.k4_2f(5).*c_KAR_NADPH.*c_C12_BKeACP...
 + P.k4_2r(6).*c_C14_KAR_NADPH_BKeACP - P.k4_2f(6).*c_KAR_NADPH.*c_C14_BKeACP...
 + P.k4_2r(7).*c_C16_KAR_NADPH_BKeACP - P.k4_2f(7).*c_KAR_NADPH.*c_C16_BKeACP; 

% C2n (n=2:10) KAR-NADPH-B-ketoacyl-ACPs
d_C4_KAR_NADPH_BKeACP = P.k4_2f(1).*c_KAR_NADPH.*c_C4_BKeACP - P.k4_2r(1).*c_C4_KAR_NADPH_BKeACP - P.kcat4(1).*c_C4_KAR_NADPH_BKeACP;
d_C6_KAR_NADPH_BKeACP = P.k4_2f(2).*c_KAR_NADPH.*c_C6_BKeACP - P.k4_2r(2).*c_C6_KAR_NADPH_BKeACP - P.kcat4(2).*c_C6_KAR_NADPH_BKeACP;
d_C8_KAR_NADPH_BKeACP = P.k4_2f(3).*c_KAR_NADPH.*c_C8_BKeACP - P.k4_2r(3).*c_C8_KAR_NADPH_BKeACP - P.kcat4(3).*c_C8_KAR_NADPH_BKeACP;
d_C10_KAR_NADPH_BKeACP = P.k4_2f(4).*c_KAR_NADPH.*c_C10_BKeACP - P.k4_2r(4).*c_C10_KAR_NADPH_BKeACP - P.kcat4(4).*c_C10_KAR_NADPH_BKeACP;
d_C12_KAR_NADPH_BKeACP = P.k4_2f(5).*c_KAR_NADPH.*c_C12_BKeACP - P.k4_2r(5).*c_C12_KAR_NADPH_BKeACP - P.kcat4(5).*c_C12_KAR_NADPH_BKeACP;
d_C14_KAR_NADPH_BKeACP = P.k4_2f(6).*c_KAR_NADPH.*c_C14_BKeACP - P.k4_2r(6).*c_C14_KAR_NADPH_BKeACP - P.kcat4(6).*c_C14_KAR_NADPH_BKeACP;
d_C16_KAR_NADPH_BKeACP = P.k4_2f(7).*c_KAR_NADPH.*c_C16_BKeACP - P.k4_2r(7).*c_C16_KAR_NADPH_BKeACP - P.kcat4(7).*c_C16_KAR_NADPH_BKeACP;
d_C18_KAR_NADPH_BKeACP = 0*c_C18_KAR_NADPH_BKeACP;
d_C20_KAR_NADPH_BKeACP = 0*c_C20_KAR_NADPH_BKeACP;

% C2n:1 (n=6:10) KAR-NADPH-B-ketoacyl-ACPs
d_C12_KAR_NADPH_BKeACP_un = 0*c_C12_KAR_NADPH_BKeACP_un;
d_C14_KAR_NADPH_BKeACP_un = 0*c_C14_KAR_NADPH_BKeACP_un;
d_C16_KAR_NADPH_BKeACP_un = 0*c_C16_KAR_NADPH_BKeACP_un;
d_C18_KAR_NADPH_BKeACP_un = 0*c_C18_KAR_NADPH_BKeACP_un;
d_C20_KAR_NADPH_BKeACP_un = 0*c_C20_KAR_NADPH_BKeACP_un;

% C2n (n=2:10) HAD-B-hydroxy-acyl-ACPs
d_C4_HAD_BHyAcACP = P.k5_1f(1).*c_HAD.*c_C4_BHyAcACP - P.k5_1r(1).*c_C4_HAD_BHyAcACP + P.k5_2r(1).*c_C4_HAD_EnAcACP - P.kcat5(1).*c_C4_HAD_BHyAcACP;
d_C6_HAD_BHyAcACP = P.k5_1f(2).*c_HAD.*c_C6_BHyAcACP - P.k5_1r(2).*c_C6_HAD_BHyAcACP + P.k5_2r(2).*c_C6_HAD_EnAcACP - P.kcat5(2).*c_C6_HAD_BHyAcACP;
d_C8_HAD_BHyAcACP = P.k5_1f(3).*c_HAD.*c_C8_BHyAcACP - P.k5_1r(3).*c_C8_HAD_BHyAcACP + P.k5_2r(3).*c_C8_HAD_EnAcACP - P.kcat5(3).*c_C8_HAD_BHyAcACP;
d_C10_HAD_BHyAcACP = P.k5_1f(4).*c_HAD.*c_C10_BHyAcACP - P.k5_1r(4).*c_C10_HAD_BHyAcACP + P.k5_2r(4).*c_C10_HAD_EnAcACP - P.kcat5(4).*c_C10_HAD_BHyAcACP;
d_C12_HAD_BHyAcACP = P.k5_1f(5).*c_HAD.*c_C12_BHyAcACP - P.k5_1r(5).*c_C12_HAD_BHyAcACP + P.k5_2r(5).*c_C12_HAD_EnAcACP - P.kcat5(5).*c_C12_HAD_BHyAcACP;
d_C14_HAD_BHyAcACP = P.k5_1f(6).*c_HAD.*c_C14_BHyAcACP - P.k5_1r(6).*c_C14_HAD_BHyAcACP + P.k5_2r(6).*c_C14_HAD_EnAcACP - P.kcat5(6).*c_C14_HAD_BHyAcACP;
d_C16_HAD_BHyAcACP = P.k5_1f(7).*c_HAD.*c_C16_BHyAcACP - P.k5_1r(7).*c_C16_HAD_BHyAcACP + P.k5_2r(7).*c_C16_HAD_EnAcACP - P.kcat5(7).*c_C16_HAD_BHyAcACP;
d_C18_HAD_BHyAcACP = 0*c_C18_HAD_BHyAcACP;
d_C20_HAD_BHyAcACP = 0*c_C20_HAD_BHyAcACP;

% C2n:1 (n=6:10) HAD-B-hydroxy-acyl-ACPs
d_C12_HAD_BHyAcACP_un = 0*c_C12_HAD_BHyAcACP_un;
d_C14_HAD_BHyAcACP_un = 0*c_C14_HAD_BHyAcACP_un;
d_C16_HAD_BHyAcACP_un = 0*c_C16_HAD_BHyAcACP_un;
d_C18_HAD_BHyAcACP_un = 0*c_C18_HAD_BHyAcACP_un;
d_C20_HAD_BHyAcACP_un = 0*c_C20_HAD_BHyAcACP_un;

% C2n (n=2:10) HAD-Enoyl-Acyl-ACPs
d_C4_HAD_EnAcACP = P.kcat5(1).*c_C4_HAD_BHyAcACP - P.k5_2r(1).*c_C4_HAD_EnAcACP + P.k5_3r(1).*c_HAD.*c_C4_EnAcACP - P.k5_3f(1).*c_C4_HAD_EnAcACP;
d_C6_HAD_EnAcACP = P.kcat5(2).*c_C6_HAD_BHyAcACP - P.k5_2r(2).*c_C6_HAD_EnAcACP + P.k5_3r(2).*c_HAD.*c_C6_EnAcACP - P.k5_3f(2).*c_C6_HAD_EnAcACP;
d_C8_HAD_EnAcACP = P.kcat5(3).*c_C8_HAD_BHyAcACP - P.k5_2r(3).*c_C8_HAD_EnAcACP + P.k5_3r(3).*c_HAD.*c_C8_EnAcACP - P.k5_3f(3).*c_C8_HAD_EnAcACP;
d_C10_HAD_EnAcACP = P.kcat5(4).*c_C10_HAD_BHyAcACP - P.k5_2r(4).*c_C10_HAD_EnAcACP + P.k5_3r(4).*c_HAD.*c_C10_EnAcACP - P.k5_3f(4).*c_C10_HAD_EnAcACP;
d_C12_HAD_EnAcACP = P.kcat5(5).*c_C12_HAD_BHyAcACP - P.k5_2r(5).*c_C12_HAD_EnAcACP + P.k5_3r(5).*c_HAD.*c_C12_EnAcACP - P.k5_3f(5).*c_C12_HAD_EnAcACP;
d_C14_HAD_EnAcACP = P.kcat5(6).*c_C14_HAD_BHyAcACP - P.k5_2r(6).*c_C14_HAD_EnAcACP + P.k5_3r(6).*c_HAD.*c_C14_EnAcACP - P.k5_3f(6).*c_C14_HAD_EnAcACP;
d_C16_HAD_EnAcACP = P.kcat5(7).*c_C16_HAD_BHyAcACP - P.k5_2r(7).*c_C16_HAD_EnAcACP + P.k5_3r(7).*c_HAD.*c_C16_EnAcACP - P.k5_3f(7).*c_C16_HAD_EnAcACP;
d_C18_HAD_EnAcACP = 0*c_C18_HAD_EnAcACP;
d_C20_HAD_EnAcACP = 0*c_C20_HAD_EnAcACP;

% C2n:1 (n=6:10) HAD-Enoyl-Acyl-ACPs
d_C12_HAD_EnAcACP_un = 0*c_C12_HAD_EnAcACP_un;
d_C14_HAD_EnAcACP_un = 0*c_C14_HAD_EnAcACP_un;
d_C16_HAD_EnAcACP_un = 0*c_C16_HAD_EnAcACP_un;
d_C18_HAD_EnAcACP_un = 0*c_C18_HAD_EnAcACP_un;
d_C20_HAD_EnAcACP_un = 0*c_C20_HAD_EnAcACP_un;

% C2n (n=2:10) SAD-B-hydroxy-acyl-ACPs
d_C4_SAD_BHyAcACP = 0*c_C4_SAD_BHyAcACP;
d_C6_SAD_BHyAcACP = 0*c_C6_SAD_BHyAcACP;
d_C8_SAD_BHyAcACP = 0*c_C8_SAD_BHyAcACP;
d_C10_SAD_BHyAcACP = 0*c_C10_SAD_BHyAcACP;
d_C12_SAD_BHyAcACP = 0*c_C12_SAD_BHyAcACP;
d_C14_SAD_BHyAcACP = 0*c_C14_SAD_BHyAcACP;
d_C16_SAD_BHyAcACP = 0*c_C16_SAD_BHyAcACP;
d_C18_SAD_BHyAcACP = 0*c_C18_SAD_BHyAcACP;
d_C20_SAD_BHyAcACP = 0*c_C20_SAD_BHyAcACP;

% C2n:1 (n=6:10) SAD-B-hydroxy-acyl-ACPs
d_C12_SAD_BHyAcACP_un = 0*c_C12_SAD_BHyAcACP_un;
d_C14_SAD_BHyAcACP_un = 0*c_C14_SAD_BHyAcACP_un;
d_C16_SAD_BHyAcACP_un = 0*c_C16_SAD_BHyAcACP_un;
d_C18_SAD_BHyAcACP_un = 0*c_C18_SAD_BHyAcACP_un;
d_C20_SAD_BHyAcACP_un = 0*c_C20_SAD_BHyAcACP_un;

% C2n (n=2:10) SAD-Enoyl-Acyl-ACPs
d_C4_SAD_EnAcACP = 0*c_C4_SAD_EnAcACP;
d_C6_SAD_EnAcACP = 0*c_C6_SAD_EnAcACP;
d_C8_SAD_EnAcACP = 0*c_C8_SAD_EnAcACP;
d_C10_SAD_EnAcACP = 0*c_C10_SAD_EnAcACP;
d_C12_SAD_EnAcACP = 0*c_C12_SAD_EnAcACP;
d_C14_SAD_EnAcACP = 0*c_C14_SAD_EnAcACP;
d_C16_SAD_EnAcACP = 0*c_C16_SAD_EnAcACP;
d_C18_SAD_EnAcACP = 0*c_C18_SAD_EnAcACP;
d_C20_SAD_EnAcACP = 0*c_C20_SAD_EnAcACP;

% SAD-C10 cis-3-Enoyl-Acyl-ACP
d_C10_SAD_cis3EnAcACP = 0*c_C10_SAD_cis3EnAcACP;

% C2n:1 (n=6:10) SAD-Enoyl-Acyl-ACPs
d_C12_SAD_EnAcACP_un = 0*c_C12_SAD_EnAcACP_un;
d_C14_SAD_EnAcACP_un = 0*c_C14_SAD_EnAcACP_un;
d_C16_SAD_EnAcACP_un = 0*c_C16_SAD_EnAcACP_un;
d_C18_SAD_EnAcACP_un = 0*c_C18_SAD_EnAcACP_un;
d_C20_SAD_EnAcACP_un = 0*c_C20_SAD_EnAcACP_un;

% ER-NADH
d_ER_NADH = P.k6_1f(1).*c_ER.*c_NADH - P.k6_1r(1).*c_ER_NADH...
 + P.k6_2r(1).*c_C4_ER_NADH_EnAcACP - P.k6_2f(1).*c_ER_NADH.*c_C4_EnAcACP...
 + P.k6_2r(2).*c_C6_ER_NADH_EnAcACP - P.k6_2f(2).*c_ER_NADH.*c_C6_EnAcACP...
 + P.k6_2r(3).*c_C8_ER_NADH_EnAcACP - P.k6_2f(3).*c_ER_NADH.*c_C8_EnAcACP...
 + P.k6_2r(4).*c_C10_ER_NADH_EnAcACP - P.k6_2f(4).*c_ER_NADH.*c_C10_EnAcACP...
 + P.k6_2r(5).*c_C12_ER_NADH_EnAcACP - P.k6_2f(5).*c_ER_NADH.*c_C12_EnAcACP...
 + P.k6_2r(6).*c_C14_ER_NADH_EnAcACP - P.k6_2f(6).*c_ER_NADH.*c_C14_EnAcACP...
 + P.k6_2r(7).*c_C16_ER_NADH_EnAcACP - P.k6_2f(7).*c_ER_NADH.*c_C16_EnAcACP;

% C2n (n=2:10) ER-NADH-Enoyl-Acyl-ACPs
d_C4_ER_NADH_EnAcACP = P.k6_2f(1).*c_ER_NADH.*c_C4_EnAcACP - P.k6_2r(1).*c_C4_ER_NADH_EnAcACP - P.kcat6(1).*c_C4_ER_NADH_EnAcACP;
d_C6_ER_NADH_EnAcACP = P.k6_2f(2).*c_ER_NADH.*c_C6_EnAcACP - P.k6_2r(2).*c_C6_ER_NADH_EnAcACP - P.kcat6(2).*c_C6_ER_NADH_EnAcACP;
d_C8_ER_NADH_EnAcACP = P.k6_2f(3).*c_ER_NADH.*c_C8_EnAcACP - P.k6_2r(3).*c_C8_ER_NADH_EnAcACP - P.kcat6(3).*c_C8_ER_NADH_EnAcACP;
d_C10_ER_NADH_EnAcACP = P.k6_2f(4).*c_ER_NADH.*c_C10_EnAcACP - P.k6_2r(4).*c_C10_ER_NADH_EnAcACP - P.kcat6(4).*c_C10_ER_NADH_EnAcACP;
d_C12_ER_NADH_EnAcACP = P.k6_2f(5).*c_ER_NADH.*c_C12_EnAcACP - P.k6_2r(5).*c_C12_ER_NADH_EnAcACP - P.kcat6(5).*c_C12_ER_NADH_EnAcACP;
d_C14_ER_NADH_EnAcACP = P.k6_2f(6).*c_ER_NADH.*c_C14_EnAcACP - P.k6_2r(6).*c_C14_ER_NADH_EnAcACP - P.kcat6(6).*c_C14_ER_NADH_EnAcACP;
d_C16_ER_NADH_EnAcACP = P.k6_2f(7).*c_ER_NADH.*c_C16_EnAcACP - P.k6_2r(7).*c_C16_ER_NADH_EnAcACP - P.kcat6(7).*c_C16_ER_NADH_EnAcACP;
d_C18_ER_NADH_EnAcACP = 0*c_C18_ER_NADH_EnAcACP;
d_C20_ER_NADH_EnAcACP = 0*c_C20_ER_NADH_EnAcACP;

% C2n:1 (n=6:10) ER-NADH-Enoyl-Acyl-ACPs
d_C12_ER_NADH_EnAcACP_un = 0*c_C12_ER_NADH_EnAcACP_un;
d_C14_ER_NADH_EnAcACP_un = 0*c_C14_ER_NADH_EnAcACP_un;
d_C16_ER_NADH_EnAcACP_un = 0*c_C16_ER_NADH_EnAcACP_un;
d_C18_ER_NADH_EnAcACP_un = 0*c_C18_ER_NADH_EnAcACP_un;
d_C20_ER_NADH_EnAcACP_un = 0*c_C20_ER_NADH_EnAcACP_un;

% C2n (n=2:10) FatA-Acyl-ACPs
d_C4_FatA_AcACP = P.k7_1f(1).*c_FatA.*c_C4_AcACP - P.k7_1r(1).*c_C4_FatA_AcACP - P.kcat7(1).*c_C4_FatA_AcACP;
d_C6_FatA_AcACP = P.k7_1f(2).*c_FatA.*c_C6_AcACP - P.k7_1r(2).*c_C6_FatA_AcACP - P.kcat7(2).*c_C6_FatA_AcACP;
d_C8_FatA_AcACP = P.k7_1f(3).*c_FatA.*c_C8_AcACP - P.k7_1r(3).*c_C8_FatA_AcACP - P.kcat7(3).*c_C8_FatA_AcACP;
d_C10_FatA_AcACP = P.k7_1f(4).*c_FatA.*c_C10_AcACP - P.k7_1r(4).*c_C10_FatA_AcACP - P.kcat7(4).*c_C10_FatA_AcACP;
d_C12_FatA_AcACP = P.k7_1f(5).*c_FatA.*c_C12_AcACP - P.k7_1r(5).*c_C12_FatA_AcACP - P.kcat7(5).*c_C12_FatA_AcACP;
d_C14_FatA_AcACP = P.k7_1f(6).*c_FatA.*c_C14_AcACP - P.k7_1r(6).*c_C14_FatA_AcACP - P.kcat7(6).*c_C14_FatA_AcACP;
d_C16_FatA_AcACP = P.k7_1f(7).*c_FatA.*c_C16_AcACP - P.k7_1r(7).*c_C16_FatA_AcACP - P.kcat7(7).*c_C16_FatA_AcACP;
d_C18_FatA_AcACP = 0*c_C18_FatA_AcACP;
d_C20_FatA_AcACP = 0*c_C20_FatA_AcACP;

% C2n:1 (n=6:10) FatA-Acyl-ACPs
d_C12_FatA_AcACP_un = 0*c_C12_FatA_AcACP_un;
d_C14_FatA_AcACP_un = 0*c_C14_FatA_AcACP_un;
d_C16_FatA_AcACP_un = 0*c_C16_FatA_AcACP_un;
d_C18_FatA_AcACP_un = 0*c_C18_FatA_AcACP_un;
d_C20_FatA_AcACP_un = 0*c_C20_FatA_AcACP_un;

% C2n (n=2:9) KASI-Acyl-ACPs
d_C4_KASI_AcACP = P.k8_1f(1).*c_KASI.*c_C4_AcACP - P.k8_1r(1).*c_C4_KASI_AcACP + P.k8_2r(1).*c_C4_KASI_Act.*c_ACP - P.k8_2f(1).*c_C4_KASI_AcACP;
d_C6_KASI_AcACP = P.k8_1f(2).*c_KASI.*c_C6_AcACP - P.k8_1r(2).*c_C6_KASI_AcACP + P.k8_2r(2).*c_C6_KASI_Act.*c_ACP - P.k8_2f(2).*c_C6_KASI_AcACP;
d_C8_KASI_AcACP = P.k8_1f(3).*c_KASI.*c_C8_AcACP - P.k8_1r(3).*c_C8_KASI_AcACP + P.k8_2r(3).*c_C8_KASI_Act.*c_ACP - P.k8_2f(3).*c_C8_KASI_AcACP;
d_C10_KASI_AcACP = P.k8_1f(4).*c_KASI.*c_C10_AcACP - P.k8_1r(4).*c_C10_KASI_AcACP + P.k8_2r(4).*c_C10_KASI_Act.*c_ACP - P.k8_2f(4).*c_C10_KASI_AcACP;
d_C12_KASI_AcACP = P.k8_1f(5).*c_KASI.*c_C12_AcACP - P.k8_1r(5).*c_C12_KASI_AcACP + P.k8_2r(5).*c_C12_KASI_Act.*c_ACP - P.k8_2f(5).*c_C12_KASI_AcACP;
d_C14_KASI_AcACP = P.k8_1f(6).*c_KASI.*c_C14_AcACP - P.k8_1r(6).*c_C14_KASI_AcACP + P.k8_2r(6).*c_C14_KASI_Act.*c_ACP - P.k8_2f(6).*c_C14_KASI_AcACP;
d_C16_KASI_AcACP = P.k8_1f(7).*c_KASI.*c_C16_AcACP - P.k8_1r(7).*c_C16_KASI_AcACP + P.k8_2r(7).*c_C16_KASI_Act.*c_ACP - P.k8_2f(7).*c_C16_KASI_AcACP;
d_C18_KASI_AcACP = 0*c_C18_KASI_AcACP;

% C2n:1 (n=6:9) KASI-Acyl-ACPs
d_C12_KASI_AcACP_un = 0*c_C12_KASI_AcACP_un;
d_C14_KASI_AcACP_un = 0*c_C14_KASI_AcACP_un;
d_C16_KASI_AcACP_un = 0*c_C16_KASI_AcACP_un;
d_C18_KASI_AcACP_un = 0*c_C18_KASI_AcACP_un;

% C2n (n=2:9) KASI*
d_C4_KASI_Act = P.k8_2f(1).*c_C4_KASI_AcACP - P.k8_2r(1).*c_C4_KASI_Act.*c_ACP + P.k8_3r(1).*c_C7_KASI_Act_MalACP - P.k8_3f(1).*c_C4_KASI_Act.*c_C3_MalACP;
d_C6_KASI_Act = P.k8_2f(2).*c_C6_KASI_AcACP - P.k8_2r(2).*c_C6_KASI_Act.*c_ACP + P.k8_3r(2).*c_C9_KASI_Act_MalACP - P.k8_3f(2).*c_C6_KASI_Act.*c_C3_MalACP;
d_C8_KASI_Act = P.k8_2f(3).*c_C8_KASI_AcACP - P.k8_2r(3).*c_C8_KASI_Act.*c_ACP + P.k8_3r(3).*c_C11_KASI_Act_MalACP - P.k8_3f(3).*c_C8_KASI_Act.*c_C3_MalACP;
d_C10_KASI_Act = P.k8_2f(4).*c_C10_KASI_AcACP - P.k8_2r(4).*c_C10_KASI_Act.*c_ACP + P.k8_3r(4).*c_C13_KASI_Act_MalACP - P.k8_3f(4).*c_C10_KASI_Act.*c_C3_MalACP;
d_C12_KASI_Act = P.k8_2f(5).*c_C12_KASI_AcACP - P.k8_2r(5).*c_C12_KASI_Act.*c_ACP + P.k8_3r(5).*c_C15_KASI_Act_MalACP - P.k8_3f(5).*c_C12_KASI_Act.*c_C3_MalACP;
d_C14_KASI_Act = P.k8_2f(6).*c_C14_KASI_AcACP - P.k8_2r(6).*c_C14_KASI_Act.*c_ACP + P.k8_3r(6).*c_C17_KASI_Act_MalACP - P.k8_3f(6).*c_C14_KASI_Act.*c_C3_MalACP;
d_C16_KASI_Act = P.k8_2f(7).*c_C16_KASI_AcACP - P.k8_2r(7).*c_C16_KASI_Act.*c_ACP;
d_C18_KASI_Act = 0*c_C18_KASI_Act;

% C2n:1 (n=6:9) KASI*
d_C12_KASI_Act_un = 0*c_C12_KASI_Act_un;
d_C14_KASI_Act_un = 0*c_C14_KASI_Act_un;
d_C16_KASI_Act_un = 0*c_C16_KASI_Act_un;
d_C18_KASI_Act_un = 0*c_C18_KASI_Act_un;

% C2n (n=2:9) KASI*-Malonyl-ACPs
d_C7_KASI_Act_MalACP = P.k8_3f(1).*c_C4_KASI_Act.*c_C3_MalACP - P.k8_3r(1).*c_C7_KASI_Act_MalACP - P.kcat8(1).*c_C7_KASI_Act_MalACP;
d_C9_KASI_Act_MalACP = P.k8_3f(2).*c_C6_KASI_Act.*c_C3_MalACP - P.k8_3r(2).*c_C9_KASI_Act_MalACP - P.kcat8(2).*c_C9_KASI_Act_MalACP;
d_C11_KASI_Act_MalACP = P.k8_3f(3).*c_C8_KASI_Act.*c_C3_MalACP - P.k8_3r(3).*c_C11_KASI_Act_MalACP - P.kcat8(3).*c_C11_KASI_Act_MalACP;
d_C13_KASI_Act_MalACP = P.k8_3f(4).*c_C10_KASI_Act.*c_C3_MalACP - P.k8_3r(4).*c_C13_KASI_Act_MalACP - P.kcat8(4).*c_C13_KASI_Act_MalACP;
d_C15_KASI_Act_MalACP = P.k8_3f(5).*c_C12_KASI_Act.*c_C3_MalACP - P.k8_3r(5).*c_C15_KASI_Act_MalACP - P.kcat8(5).*c_C15_KASI_Act_MalACP;
d_C17_KASI_Act_MalACP = P.k8_3f(6).*c_C14_KASI_Act.*c_C3_MalACP - P.k8_3r(6).*c_C17_KASI_Act_MalACP - P.kcat8(6).*c_C17_KASI_Act_MalACP;
d_C19_KASI_Act_MalACP = 0*c_C19_KASI_Act_MalACP;
d_C21_KASI_Act_MalACP = 0*c_C21_KASI_Act_MalACP;

% C2n:1 (n=6:9) KASI*-Malonyl-ACPs
d_C15_KASI_Act_MalACP_un = 0*c_C15_KASI_Act_MalACP_un;
d_C17_KASI_Act_MalACP_un = 0*c_C17_KASI_Act_MalACP_un;
d_C19_KASI_Act_MalACP_un = 0*c_C19_KASI_Act_MalACP_un;
d_C21_KASI_Act_MalACP_un = 0*c_C21_KASI_Act_MalACP_un;

% C2n (n=2:9) KASII-Acyl-ACPs
d_C4_KASII_AcACP = 0*c_C4_KASII_AcACP;
d_C6_KASII_AcACP = 0*c_C6_KASII_AcACP;
d_C8_KASII_AcACP = 0*c_C8_KASII_AcACP;
d_C10_KASII_AcACP = 0*c_C10_KASII_AcACP;
d_C12_KASII_AcACP = 0*c_C12_KASII_AcACP;
d_C14_KASII_AcACP = 0*c_C14_KASII_AcACP;
d_C16_KASII_AcACP = 0*c_C16_KASII_AcACP;
d_C18_KASII_AcACP = 0*c_C18_KASII_AcACP;

% C2n:1 (n=6:9) KASII-Acyl-ACPs
d_C12_KASII_AcACP_un = 0*c_C12_KASII_AcACP_un;
d_C14_KASII_AcACP_un = 0*c_C14_KASII_AcACP_un;
d_C16_KASII_AcACP_un = 0*c_C16_KASII_AcACP_un;
d_C18_KASII_AcACP_un = 0*c_C18_KASII_AcACP_un;

% C2n (n=2:9) KASII*
d_C4_KASII_Act = 0*c_C4_KASII_Act;
d_C6_KASII_Act = 0*c_C6_KASII_Act;
d_C8_KASII_Act = 0*c_C8_KASII_Act;
d_C10_KASII_Act = 0*c_C10_KASII_Act;
d_C12_KASII_Act = 0*c_C12_KASII_Act;
d_C14_KASII_Act = 0*c_C14_KASII_Act;
d_C16_KASII_Act = 0*c_C16_KASII_Act;
d_C18_KASII_Act = 0*c_C18_KASII_Act;

% C2n:1 (n=6:9) KASII*
d_C12_KASII_Act_un = 0*c_C12_KASII_Act_un;
d_C14_KASII_Act_un = 0*c_C14_KASII_Act_un;
d_C16_KASII_Act_un = 0*c_C16_KASII_Act_un;
d_C18_KASII_Act_un = 0*c_C18_KASII_Act_un;

% C2n (n=2:9) KASII*-Malonyl-ACPs
d_C7_KASII_Act_MalACP = 0*c_C7_KASII_Act_MalACP;
d_C9_KASII_Act_MalACP = 0*c_C9_KASII_Act_MalACP;
d_C11_KASII_Act_MalACP = 0*c_C11_KASII_Act_MalACP;
d_C13_KASII_Act_MalACP = 0*c_C13_KASII_Act_MalACP;
d_C15_KASII_Act_MalACP = 0*c_C15_KASII_Act_MalACP;
d_C17_KASII_Act_MalACP = 0*c_C17_KASII_Act_MalACP;
d_C19_KASII_Act_MalACP = 0*c_C19_KASII_Act_MalACP;
d_C21_KASII_Act_MalACP = 0*c_C21_KASII_Act_MalACP;

% C2n:1 (n=6:9) KASII*-Malonyl-ACPs
d_C15_KASII_Act_MalACP_un = 0*c_C15_KASII_Act_MalACP_un;
d_C17_KASII_Act_MalACP_un = 0*c_C17_KASII_Act_MalACP_un;
d_C19_KASII_Act_MalACP_un = 0*c_C19_KASII_Act_MalACP_un;
d_C21_KASII_Act_MalACP_un = 0*c_C21_KASII_Act_MalACP_un;

% KASII-(C10 cis-3-Enoyl-Acyl-ACP)
d_C10_KASII_cis3EnAcACP = 0*c_C10_KASII_cis3EnAcACP;

% C10 cis-3-KASII*
d_C10_KASII_Act_cis3 = 0*c_C10_KASII_Act_cis3;

% C10 cis-3-KASII*-Malonyl-ACP
d_C13_KASII_Act_cis3MalACP = 0*c_C13_KASII_Act_cis3MalACP;

% C2n (n=2:10) KASIII-Acyl-ACPs
d_C4_KASIII_AcACP = P.k3_4f(1).*c_KASIII.*c_C4_AcACP - P.k3_4r(1).*c_C4_KASIII_AcACP;
d_C6_KASIII_AcACP = P.k3_4f(2).*c_KASIII.*c_C6_AcACP - P.k3_4r(2).*c_C6_KASIII_AcACP;
d_C8_KASIII_AcACP = P.k3_4f(3).*c_KASIII.*c_C8_AcACP - P.k3_4r(3).*c_C8_KASIII_AcACP;
d_C10_KASIII_AcACP = P.k3_4f(4).*c_KASIII.*c_C10_AcACP - P.k3_4r(4).*c_C10_KASIII_AcACP;
d_C12_KASIII_AcACP = P.k3_4f(5).*c_KASIII.*c_C12_AcACP - P.k3_4r(5).*c_C12_KASIII_AcACP;
d_C14_KASIII_AcACP = P.k3_4f(6).*c_KASIII.*c_C14_AcACP - P.k3_4r(6).*c_C14_KASIII_AcACP;
d_C16_KASIII_AcACP = P.k3_4f(7).*c_KASIII.*c_C16_AcACP - P.k3_4r(7).*c_C16_KASIII_AcACP;
d_C18_KASIII_AcACP = 0*c_C18_KASIII_AcACP;
d_C20_KASIII_AcACP = 0*c_C20_KASIII_AcACP;

% C2n:1 (n=6:10) KASIII-Acyl-ACPs
d_C12_KASIII_AcACP_un = 0*c_C12_KASIII_AcACP_un;
d_C14_KASIII_AcACP_un = 0*c_C14_KASIII_AcACP_un;
d_C16_KASIII_AcACP_un = 0*c_C16_KASIII_AcACP_un;
d_C18_KASIII_AcACP_un = 0*c_C18_KASIII_AcACP_un;
d_C20_KASIII_AcACP_un = 0*c_C20_KASIII_AcACP_un;

% C2n (n=2:10) KASIII*-Acyl-ACPs
d_C6_KASIII_Act_AcACP = P.k3_5f(1).*c_C2_KASIII_Act.*c_C4_AcACP - P.k3_5r(1).*c_C6_KASIII_Act_AcACP;
d_C8_KASIII_Act_AcACP = P.k3_5f(2).*c_C2_KASIII_Act.*c_C6_AcACP - P.k3_5r(2).*c_C8_KASIII_Act_AcACP;
d_C10_KASIII_Act_AcACP = P.k3_5f(3).*c_C2_KASIII_Act.*c_C8_AcACP - P.k3_5r(3).*c_C10_KASIII_Act_AcACP;
d_C12_KASIII_Act_AcACP = P.k3_5f(4).*c_C2_KASIII_Act.*c_C10_AcACP - P.k3_5r(4).*c_C12_KASIII_Act_AcACP;
d_C14_KASIII_Act_AcACP = P.k3_5f(5).*c_C2_KASIII_Act.*c_C12_AcACP - P.k3_5r(5).*c_C14_KASIII_Act_AcACP;
d_C16_KASIII_Act_AcACP = P.k3_5f(6).*c_C2_KASIII_Act.*c_C14_AcACP - P.k3_5r(6).*c_C16_KASIII_Act_AcACP;
d_C18_KASIII_Act_AcACP = 0*c_C18_KASIII_Act_AcACP;
d_C20_KASIII_Act_AcACP = 0*c_C20_KASIII_Act_AcACP;
d_C22_KASIII_Act_AcACP = 0*c_C22_KASIII_Act_AcACP;

% C2n:1 (n=6:10) KASIII*-Acyl-ACPs
d_C14_KASIII_Act_AcACP_un = 0*c_C14_KASIII_Act_AcACP_un;
d_C16_KASIII_Act_AcACP_un = 0*c_C16_KASIII_Act_AcACP_un;
d_C18_KASIII_Act_AcACP_un = 0*c_C18_KASIII_Act_AcACP_un;
d_C20_KASIII_Act_AcACP_un = 0*c_C20_KASIII_Act_AcACP_un;
d_C22_KASIII_Act_AcACP_un = 0*c_C22_KASIII_Act_AcACP_un;

% FatA-ACP
d_FatA_ACP = P.k7_inh_f.*c_FatA.*c_ACP - P.k7_inh_r.*c_FatA_ACP;

% KASIII-ACP
d_KASIII_ACP = P.k3_inh_f.*c_KASIII.*c_ACP - P.k3_inh_r.*c_KASIII_ACP;

% KAR-ACP
d_KAR_ACP = P.k4_inh_f.*c_KAR.*c_ACP - P.k4_inh_r.*c_KAR_ACP;

% HAD-ACP
d_HAD_ACP = P.k5_inh_f.*c_HAD.*c_ACP - P.k5_inh_r.*c_HAD_ACP;

% ER-ACP
d_ER_ACP = P.k6_inh_f.*c_ER.*c_ACP - P.k6_inh_r.*c_ER_ACP;

% KASI-ACP
d_KASI_ACP = P.k8_inh_f.*c_KASI.*c_ACP - P.k8_inh_r.*c_KASI_ACP;

% SAD-ACP
d_SAD_ACP = 0*c_SAD_ACP;

% KASII-ACP
d_KASII_ACP = 0*c_KASII_ACP;


% Giving KASII KASIII-like activity
% KASII-Acetyl-CoA
d_C2_KASII_AcCoA = 0*c_C2_KASII_AcCoA;

% KASII*
d_C2_KASII_Act = 0*c_C2_KASII_Act;

% KASII*-Malonyl-ACP
d_C5_KASII_Act_MalACP = 0*c_C5_KASII_Act_MalACP;

% Giving KASI KASIII-like activity
% KASI-Acetyl-CoA
d_C2_KASI_AcCoA = P.k8_4f.*c_KASI.*c_C2_AcCoA - P.k8_4r.*c_C2_KASI_AcCoA + P.k8_5r.*c_C2_KASI_Act.*c_CoA - P.k8_5f.*c_C2_KASI_AcCoA; 

% KASI*
d_C2_KASI_Act = P.k8_5f.*c_C2_KASI_AcCoA - P.k8_5r.*c_C2_KASI_Act.*c_CoA + P.k8_6r.*c_C5_KASI_Act_MalACP - P.k8_6f.*c_C2_KASI_Act.*c_C3_MalACP + P.k8_9f.*c_C2_KASI_AcACP - P.k8_9r.*c_C2_KASI_Act.*c_ACP; 

% KASI*-Malonyl-ACP
d_C5_KASI_Act_MalACP = P.k8_6f.*c_C2_KASI_Act.*c_C3_MalACP - P.k8_6r.*c_C5_KASI_Act_MalACP - P.kcat8_H.*c_C5_KASI_Act_MalACP; 

% KASII and KASI decarboxylating mACP to form aACP and reacting with it to form activated enzyme (initiation)
% KASII-Malonyl-ACP
d_C3_KASII_MalACP = 0*c_C3_KASII_MalACP;

% Acetyl-ACP
d_C2_AcACP = P.k8_8r.*c_C2_KASI_AcACP - P.k8_8f.*c_KASI.*c_C2_AcACP; 

% KASII-Acetyl-ACP
d_C2_KASII_AcACP = 0*c_C2_KASII_AcACP;

% KASI-Malonyl-ACP
d_C3_KASI_MalACP = P.k8_7f.*c_KASI.*c_C3_MalACP - P.k8_7r.*c_C3_KASI_MalACP - P.kcat8_CO2.*c_C3_KASI_MalACP; 

% KASI-Acetyl-ACP
d_C2_KASI_AcACP = P.kcat8_CO2.*c_C3_KASI_MalACP + P.k8_8f.*c_KASI.*c_C2_AcACP - P.k8_8r.*c_C2_KASI_AcACP + P.k8_9r.*c_C2_KASI_Act.*c_ACP - P.k8_9f.*c_C2_KASI_AcACP; 


dcdt = [d_ATP; d_C1_Bicarbonate; d_C2_AcCoA; d_C4_SucCoA; d_C6_HexCoA; d_C8_OcCoA; d_C10_DecCoA; d_C12_LauCoA; d_C14_EthCoA; d_C16_PalCoA; d_C18_OcDecCoA; d_ACP;...
    d_NADPH; d_NADP; d_NADH; d_NAD; d_ADP; d_C3_MalCoA; d_CoA; d_C3_MalACP; d_C1_CO2; d_C4_BKeACP; d_C6_BKeACP; d_C8_BKeACP; d_C10_BKeACP; d_C12_BKeACP;...
    d_C14_BKeACP; d_C16_BKeACP; d_C18_BKeACP; d_C20_BKeACP; d_C12_BKeACP_un; d_C14_BKeACP_un; d_C16_BKeACP_un; d_C18_BKeACP_un; d_C20_BKeACP_un; d_C4_BHyAcACP;...
    d_C6_BHyAcACP; d_C8_BHyAcACP; d_C10_BHyAcACP; d_C12_BHyAcACP; d_C14_BHyAcACP; d_C16_BHyAcACP; d_C18_BHyAcACP; d_C20_BHyAcACP; d_C12_BHyAcACP_un;...
    d_C14_BHyAcACP_un; d_C16_BHyAcACP_un; d_C18_BHyAcACP_un; d_C20_BHyAcACP_un; d_C4_EnAcACP; d_C6_EnAcACP; d_C8_EnAcACP; d_C10_EnAcACP; d_C12_EnAcACP;...
    d_C14_EnAcACP; d_C16_EnAcACP; d_C18_EnAcACP; d_C20_EnAcACP; d_C10_cis3EnAcACP; d_C12_EnAcACP_un; d_C14_EnAcACP_un; d_C16_EnAcACP_un; d_C18_EnAcACP_un;...
    d_C20_EnAcACP_un; d_C4_AcACP; d_C6_AcACP; d_C8_AcACP; d_C10_AcACP; d_C12_AcACP; d_C14_AcACP; d_C16_AcACP; d_C18_AcACP; d_C20_AcACP;...
    d_C12_AcACP_un; d_C14_AcACP_un; d_C16_AcACP_un; d_C18_AcACP_un; d_C20_AcACP_un; d_C4_FA; d_C6_FA; d_C8_FA; d_C10_FA; d_C12_FA; d_C14_FA;...
    d_C16_FA; d_C18_FA; d_C20_FA; d_C12_FA_un; d_C14_FA_un; d_C16_FA_un; d_C18_FA_un; d_C20_FA_un; d_BC_ATP; d_C1_BC_ATP_HCO3; d_C1_BC_Pi_HCO3;...
    d_C1_BC_Pi_HCO3_BCCP_Biotin; d_C1_BCCP_Biotin_CO2; d_C1_CT_BCCP_Biotin_CO2; d_C1_CT_Act; d_C3_CT_Act_AcCoA; d_C3_MCMT_MalCoA;...
    d_C3_MCMT_Act; d_C3_MCMT_Act_ACP; d_C2_KASIII_CoA; d_C4_KASIII_CoA; d_C6_KASIII_CoA; d_C8_KASIII_CoA; d_C10_KASIII_CoA; d_C12_KASIII_CoA; d_C14_KASIII_CoA;...
    d_C16_KASIII_CoA; d_C18_KASIII_CoA; d_C2_KASIII_Act; d_C4_KASIII_Act; d_C6_KASIII_Act; d_C8_KASIII_Act; d_C10_KASIII_Act; d_C12_KASIII_Act; d_C14_KASIII_Act;...
    d_C16_KASIII_Act; d_C18_KASIII_Act; d_C5_KASIII_Act_MalACP; d_C7_KASIII_Act_MalACP; d_C9_KASIII_Act_MalACP; d_C11_KASIII_Act_MalACP; d_C13_KASIII_Act_MalACP;...
    d_C15_KASIII_Act_MalACP; d_C17_KASIII_Act_MalACP; d_C19_KASIII_Act_MalACP; d_C21_KASIII_Act_MalACP; d_KAR_NADPH; d_C4_KAR_NADPH_BKeACP;...
    d_C6_KAR_NADPH_BKeACP; d_C8_KAR_NADPH_BKeACP; d_C10_KAR_NADPH_BKeACP; d_C12_KAR_NADPH_BKeACP; d_C14_KAR_NADPH_BKeACP;...
    d_C16_KAR_NADPH_BKeACP; d_C18_KAR_NADPH_BKeACP; d_C20_KAR_NADPH_BKeACP; d_C12_KAR_NADPH_BKeACP_un; d_C14_KAR_NADPH_BKeACP_un;...
    d_C16_KAR_NADPH_BKeACP_un; d_C18_KAR_NADPH_BKeACP_un; d_C20_KAR_NADPH_BKeACP_un; d_C4_HAD_BHyAcACP; d_C6_HAD_BHyAcACP;...
    d_C8_HAD_BHyAcACP; d_C10_HAD_BHyAcACP; d_C12_HAD_BHyAcACP; d_C14_HAD_BHyAcACP; d_C16_HAD_BHyAcACP; d_C18_HAD_BHyAcACP;...
    d_C20_HAD_BHyAcACP; d_C12_HAD_BHyAcACP_un; d_C14_HAD_BHyAcACP_un; d_C16_HAD_BHyAcACP_un; d_C18_HAD_BHyAcACP_un; d_C20_HAD_BHyAcACP_un;...
    d_C4_HAD_EnAcACP; d_C6_HAD_EnAcACP; d_C8_HAD_EnAcACP; d_C10_HAD_EnAcACP; d_C12_HAD_EnAcACP; d_C14_HAD_EnAcACP; d_C16_HAD_EnAcACP;...
    d_C18_HAD_EnAcACP; d_C20_HAD_EnAcACP; d_C12_HAD_EnAcACP_un; d_C14_HAD_EnAcACP_un; d_C16_HAD_EnAcACP_un; d_C18_HAD_EnAcACP_un;...
    d_C20_HAD_EnAcACP_un; d_C4_SAD_BHyAcACP; d_C6_SAD_BHyAcACP; d_C8_SAD_BHyAcACP; d_C10_SAD_BHyAcACP; d_C12_SAD_BHyAcACP;...
    d_C14_SAD_BHyAcACP; d_C16_SAD_BHyAcACP; d_C18_SAD_BHyAcACP; d_C20_SAD_BHyAcACP; d_C12_SAD_BHyAcACP_un; d_C14_SAD_BHyAcACP_un;...
    d_C16_SAD_BHyAcACP_un; d_C18_SAD_BHyAcACP_un; d_C20_SAD_BHyAcACP_un; d_C4_SAD_EnAcACP; d_C6_SAD_EnAcACP; d_C8_SAD_EnAcACP; d_C10_SAD_EnAcACP;...
    d_C12_SAD_EnAcACP; d_C14_SAD_EnAcACP; d_C16_SAD_EnAcACP; d_C18_SAD_EnAcACP; d_C20_SAD_EnAcACP; d_C10_SAD_cis3EnAcACP; d_C12_SAD_EnAcACP_un;...
    d_C14_SAD_EnAcACP_un; d_C16_SAD_EnAcACP_un; d_C18_SAD_EnAcACP_un; d_C20_SAD_EnAcACP_un; d_ER_NADH; d_C4_ER_NADH_EnAcACP;...
    d_C6_ER_NADH_EnAcACP; d_C8_ER_NADH_EnAcACP; d_C10_ER_NADH_EnAcACP; d_C12_ER_NADH_EnAcACP; d_C14_ER_NADH_EnAcACP;...
    d_C16_ER_NADH_EnAcACP; d_C18_ER_NADH_EnAcACP; d_C20_ER_NADH_EnAcACP; d_C12_ER_NADH_EnAcACP_un; d_C14_ER_NADH_EnAcACP_un;...
    d_C16_ER_NADH_EnAcACP_un; d_C18_ER_NADH_EnAcACP_un; d_C20_ER_NADH_EnAcACP_un; d_C4_FatA_AcACP; d_C6_FatA_AcACP; d_C8_FatA_AcACP;...
    d_C10_FatA_AcACP; d_C12_FatA_AcACP; d_C14_FatA_AcACP; d_C16_FatA_AcACP; d_C18_FatA_AcACP; d_C20_FatA_AcACP; d_C12_FatA_AcACP_un; d_C14_FatA_AcACP_un;...
    d_C16_FatA_AcACP_un; d_C18_FatA_AcACP_un; d_C20_FatA_AcACP_un; d_C4_KASI_AcACP; d_C6_KASI_AcACP; d_C8_KASI_AcACP; d_C10_KASI_AcACP;...
    d_C12_KASI_AcACP; d_C14_KASI_AcACP; d_C16_KASI_AcACP; d_C18_KASI_AcACP; d_C12_KASI_AcACP_un; d_C14_KASI_AcACP_un; d_C16_KASI_AcACP_un;...
    d_C18_KASI_AcACP_un; d_C4_KASI_Act; d_C6_KASI_Act; d_C8_KASI_Act; d_C10_KASI_Act; d_C12_KASI_Act; d_C14_KASI_Act; d_C16_KASI_Act; d_C18_KASI_Act;...
    d_C12_KASI_Act_un; d_C14_KASI_Act_un; d_C16_KASI_Act_un; d_C18_KASI_Act_un; d_C7_KASI_Act_MalACP; d_C9_KASI_Act_MalACP; d_C11_KASI_Act_MalACP;...
    d_C13_KASI_Act_MalACP; d_C15_KASI_Act_MalACP; d_C17_KASI_Act_MalACP; d_C19_KASI_Act_MalACP; d_C21_KASI_Act_MalACP; d_C15_KASI_Act_MalACP_un;...
    d_C17_KASI_Act_MalACP_un; d_C19_KASI_Act_MalACP_un; d_C21_KASI_Act_MalACP_un; d_C4_KASII_AcACP; d_C6_KASII_AcACP; d_C8_KASII_AcACP; d_C10_KASII_AcACP;...
    d_C12_KASII_AcACP; d_C14_KASII_AcACP; d_C16_KASII_AcACP; d_C18_KASII_AcACP; d_C12_KASII_AcACP_un; d_C14_KASII_AcACP_un; d_C16_KASII_AcACP_un;...
    d_C18_KASII_AcACP_un; d_C4_KASII_Act; d_C6_KASII_Act; d_C8_KASII_Act; d_C10_KASII_Act; d_C12_KASII_Act; d_C14_KASII_Act; d_C16_KASII_Act; d_C18_KASII_Act;...
    d_C12_KASII_Act_un; d_C14_KASII_Act_un; d_C16_KASII_Act_un; d_C18_KASII_Act_un; d_C7_KASII_Act_MalACP; d_C9_KASII_Act_MalACP; d_C11_KASII_Act_MalACP;...
    d_C13_KASII_Act_MalACP; d_C15_KASII_Act_MalACP; d_C17_KASII_Act_MalACP; d_C19_KASII_Act_MalACP; d_C21_KASII_Act_MalACP; d_C15_KASII_Act_MalACP_un;...
    d_C17_KASII_Act_MalACP_un; d_C19_KASII_Act_MalACP_un; d_C21_KASII_Act_MalACP_un; d_C10_KASII_cis3EnAcACP; d_C10_KASII_Act_cis3; d_C13_KASII_Act_cis3MalACP;...
    d_C4_KASIII_AcACP; d_C6_KASIII_AcACP; d_C8_KASIII_AcACP; d_C10_KASIII_AcACP; d_C12_KASIII_AcACP; d_C14_KASIII_AcACP; d_C16_KASIII_AcACP; d_C18_KASIII_AcACP;...
    d_C20_KASIII_AcACP; d_C12_KASIII_AcACP_un; d_C14_KASIII_AcACP_un; d_C16_KASIII_AcACP_un; d_C18_KASIII_AcACP_un; d_C20_KASIII_AcACP_un; d_C6_KASIII_Act_AcACP;...
    d_C8_KASIII_Act_AcACP; d_C10_KASIII_Act_AcACP; d_C12_KASIII_Act_AcACP; d_C14_KASIII_Act_AcACP; d_C16_KASIII_Act_AcACP; d_C18_KASIII_Act_AcACP; d_C20_KASIII_Act_AcACP;...
    d_C22_KASIII_Act_AcACP; d_C14_KASIII_Act_AcACP_un; d_C16_KASIII_Act_AcACP_un; d_C18_KASIII_Act_AcACP_un; d_C20_KASIII_Act_AcACP_un; d_C22_KASIII_Act_AcACP_un;...
    d_FatA_ACP; d_KASIII_ACP; d_KAR_ACP; d_HAD_ACP; d_ER_ACP; d_KASI_ACP; d_SAD_ACP; d_KASII_ACP; d_C2_KASII_AcCoA; d_C2_KASII_Act; d_C5_KASII_Act_MalACP;...
    d_C2_KASI_AcCoA; d_C2_KASI_Act; d_C5_KASI_Act_MalACP; d_C3_KASII_MalACP; d_C2_AcACP; d_C2_KASII_AcACP; d_C3_KASI_MalACP; d_C2_KASI_AcACP];

