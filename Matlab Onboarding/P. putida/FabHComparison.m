% clear all
% close all
% clc

%% Variables

% Set up structure to store variables
S = struct;

% p_vec = [a1 a2 a3 b1 c2 c3 c1 b2 b3 d1 d2 e f c4 x1 x2 x3 x4];
S.p_vec = [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
    4.645634978 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673 132.8499358...
    2180.050007 0.539756276	0.053673263	34.49718991	11.15058888];

% Potential FA lengths for calculating production profile
% Numbers repeated are saturated and unsaturated
% changed
S.FA_dist = [4,6,8,10,12,14,16,18,20,12,14,16,18,20];

% NO NAD+ NADP+ right now
% changed
S.labels = {'c_ATP', 'c_C1_Bicarbonate', 'c_C2_AcCoA', 'c_C4_SucCoA', 'c_C6_HexCoA', 'c_C8_OcCoA', 'c_C10_DecCoA', 'c_C12_LauCoA', 'c_C14_EthCoA', 'c_C16_PalCoA', 'c_C18_OcDecCoA', 'c_ACP', 'c_NADPH', 'c_NADH', 'c_ADP',...
    'c_C3_MalCoA', 'c_CoA', 'c_C3_MalACP', 'c_C1_CO2', 'c_C4_BKeACP', 'c_C6_BKeACP', 'c_C8_BKeACP', 'c_C10_BKeACP', 'c_C12_BKeACP', 'c_C14_BKeACP', 'c_C16_BKeACP', 'c_C18_BKeACP',...
    'c_C20_BKeACP', 'c_C12_BKeACP_un', 'c_C14_BKeACP_un', 'c_C16_BKeACP_un', 'c_C18_BKeACP_un', 'c_C20_BKeACP_un', 'c_C4_BHyAcACP', 'c_C6_BHyAcACP', 'c_C8_BHyAcACP',...
    'c_C10_BHyAcACP', 'c_C12_BHyAcACP', 'c_C14_BHyAcACP', 'c_C16_BHyAcACP', 'c_C18_BHyAcACP', 'c_C20_BHyAcACP', 'c_C12_BHyAcACP_un', 'c_C14_BHyAcACP_un',...
    'c_C16_BHyAcACP_un', 'c_C18_BHyAcACP_un', 'c_C20_BHyAcACP_un', 'c_C4_EnAcACP', 'c_C6_EnAcACP', 'c_C8_EnAcACP', 'c_C10_EnAcACP', 'c_C12_EnAcACP', 'c_C14_EnAcACP',...
    'c_C16_EnAcACP', 'c_C18_EnAcACP', 'c_C20_EnAcACP', 'c_C10_cis3EnAcACP', 'c_C12_EnAcACP_un', 'c_C14_EnAcACP_un', 'c_C16_EnAcACP_un', 'c_C18_EnAcACP_un',...
    'c_C20_EnAcACP_un', 'c_C4_AcACP', 'c_C6_AcACP', 'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP', 'c_C14_AcACP', 'c_C16_AcACP', 'c_C18_AcACP', 'c_C20_AcACP',...
    'c_C12_AcACP_un', 'c_C14_AcACP_un', 'c_C16_AcACP_un', 'c_C18_AcACP_un', 'c_C20_AcACP_un', 'c_C4_FA', 'c_C6_FA', 'c_C8_FA', 'c_C10_FA', 'c_C12_FA', 'c_C14_FA',...
    'c_C16_FA', 'c_C18_FA', 'c_C20_FA', 'c_C12_FA_un', 'c_C14_FA_un', 'c_C16_FA_un', 'c_C18_FA_un', 'c_C20_FA_un', 'c_ACC_s1', 'c_C1_ACC_s2', 'c_C1_ACC_s3', 'c_C3_ACC_s4', 'c_C3_FabD_MalCoA',...
    'c_C3_FabD_Act', 'c_C3_FabD_Act_ACP', 'c_C2_FabH_CoA', 'c_C4_FabH_CoA', 'c_C6_FabH_CoA', 'c_C8_FabH_CoA', 'c_C10_FabH_CoA', 'c_C12_FabH_CoA', 'c_C14_FabH_CoA',...
    'c_C16_FabH_CoA', 'c_C18_FabH_CoA', 'c_C2_FabH_Act', 'c_C4_FabH_Act', 'c_C6_FabH_Act', 'c_C8_FabH_Act', 'c_C10_FabH_Act', 'c_C12_FabH_Act', 'c_C14_FabH_Act',...
    'c_C16_FabH_Act', 'c_C18_FabH_Act', 'c_C5_FabH_Act_MalACP', 'c_C7_FabH_Act_MalACP', 'c_C9_FabH_Act_MalACP', 'c_C11_FabH_Act_MalACP', 'c_C13_FabH_Act_MalACP',...
    'c_C15_FabH_Act_MalACP', 'c_C17_FabH_Act_MalACP', 'c_C19_FabH_Act_MalACP', 'c_C21_FabH_Act_MalACP', 'c_FabG_NADPH', 'c_C4_FabG_NADPH_BKeACP',...
    'c_C6_FabG_NADPH_BKeACP', 'c_C8_FabG_NADPH_BKeACP', 'c_C10_FabG_NADPH_BKeACP', 'c_C12_FabG_NADPH_BKeACP', 'c_C14_FabG_NADPH_BKeACP',...
    'c_C16_FabG_NADPH_BKeACP', 'c_C18_FabG_NADPH_BKeACP', 'c_C20_FabG_NADPH_BKeACP', 'c_C12_FabG_NADPH_BKeACP_un', 'c_C14_FabG_NADPH_BKeACP_un',...
    'c_C16_FabG_NADPH_BKeACP_un', 'c_C18_FabG_NADPH_BKeACP_un', 'c_C20_FabG_NADPH_BKeACP_un', 'c_C4_FabZ_BHyAcACP', 'c_C6_FabZ_BHyAcACP',...
    'c_C8_FabZ_BHyAcACP', 'c_C10_FabZ_BHyAcACP', 'c_C12_FabZ_BHyAcACP', 'c_C14_FabZ_BHyAcACP', 'c_C16_FabZ_BHyAcACP', 'c_C18_FabZ_BHyAcACP',...
    'c_C20_FabZ_BHyAcACP', 'c_C12_FabZ_BHyAcACP_un', 'c_C14_FabZ_BHyAcACP_un', 'c_C16_FabZ_BHyAcACP_un', 'c_C18_FabZ_BHyAcACP_un', 'c_C20_FabZ_BHyAcACP_un',...
    'c_C4_FabZ_EnAcACP', 'c_C6_FabZ_EnAcACP', 'c_C8_FabZ_EnAcACP', 'c_C10_FabZ_EnAcACP', 'c_C12_FabZ_EnAcACP', 'c_C14_FabZ_EnAcACP', 'c_C16_FabZ_EnAcACP',...
    'c_C18_FabZ_EnAcACP', 'c_C20_FabZ_EnAcACP', 'c_C12_FabZ_EnAcACP_un', 'c_C14_FabZ_EnAcACP_un', 'c_C16_FabZ_EnAcACP_un', 'c_C18_FabZ_EnAcACP_un',...
    'c_C20_FabZ_EnAcACP_un', 'c_C4_FabA_BHyAcACP', 'c_C6_FabA_BHyAcACP', 'c_C8_FabA_BHyAcACP', 'c_C10_FabA_BHyAcACP', 'c_C12_FabA_BHyAcACP',...
    'c_C14_FabA_BHyAcACP', 'c_C16_FabA_BHyAcACP', 'c_C18_FabA_BHyAcACP', 'c_C20_FabA_BHyAcACP', 'c_C12_FabA_BHyAcACP_un', 'c_C14_FabA_BHyAcACP_un',...
    'c_C16_FabA_BHyAcACP_un', 'c_C18_FabA_BHyAcACP_un', 'c_C20_FabA_BHyAcACP_un', 'c_C4_FabA_EnAcACP', 'c_C6_FabA_EnAcACP', 'c_C8_FabA_EnAcACP', 'c_C10_FabA_EnAcACP',...
    'c_C12_FabA_EnAcACP', 'c_C14_FabA_EnAcACP', 'c_C16_FabA_EnAcACP', 'c_C18_FabA_EnAcACP', 'c_C20_FabA_EnAcACP', 'c_C10_FabA_cis3EnAcACP', 'c_C12_FabA_EnAcACP_un',...
    'c_C14_FabA_EnAcACP_un', 'c_C16_FabA_EnAcACP_un', 'c_C18_FabA_EnAcACP_un', 'c_C20_FabA_EnAcACP_un', 'c_FabI_NADH', 'c_C4_FabI_NADH_EnAcACP',...
    'c_C6_FabI_NADH_EnAcACP', 'c_C8_FabI_NADH_EnAcACP', 'c_C10_FabI_NADH_EnAcACP', 'c_C12_FabI_NADH_EnAcACP', 'c_C14_FabI_NADH_EnAcACP',...
    'c_C16_FabI_NADH_EnAcACP', 'c_C18_FabI_NADH_EnAcACP', 'c_C20_FabI_NADH_EnAcACP', 'c_C12_FabI_NADH_EnAcACP_un', 'c_C14_FabI_NADH_EnAcACP_un',...
    'c_C16_FabI_NADH_EnAcACP_un', 'c_C18_FabI_NADH_EnAcACP_un', 'c_C20_FabI_NADH_EnAcACP_un', 'c_C4_TesA_AcACP', 'c_C6_TesA_AcACP', 'c_C8_TesA_AcACP',...
    'c_C10_TesA_AcACP', 'c_C12_TesA_AcACP', 'c_C14_TesA_AcACP', 'c_C16_TesA_AcACP', 'c_C18_TesA_AcACP', 'c_C20_TesA_AcACP', 'c_C12_TesA_AcACP_un', 'c_C14_TesA_AcACP_un',...
    'c_C16_TesA_AcACP_un', 'c_C18_TesA_AcACP_un', 'c_C20_TesA_AcACP_un', 'c_C4_FabF_AcACP', 'c_C6_FabF_AcACP', 'c_C8_FabF_AcACP', 'c_C10_FabF_AcACP',...
    'c_C12_FabF_AcACP', 'c_C14_FabF_AcACP', 'c_C16_FabF_AcACP', 'c_C18_FabF_AcACP', 'c_C12_FabF_AcACP_un', 'c_C14_FabF_AcACP_un', 'c_C16_FabF_AcACP_un',...
    'c_C18_FabF_AcACP_un', 'c_C4_FabF_Act', 'c_C6_FabF_Act', 'c_C8_FabF_Act', 'c_C10_FabF_Act', 'c_C12_FabF_Act', 'c_C14_FabF_Act', 'c_C16_FabF_Act', 'c_C18_FabF_Act',...
    'c_C12_FabF_Act_un', 'c_C14_FabF_Act_un', 'c_C16_FabF_Act_un', 'c_C18_FabF_Act_un', 'c_C7_FabF_Act_MalACP', 'c_C9_FabF_Act_MalACP', 'c_C11_FabF_Act_MalACP',...
    'c_C13_FabF_Act_MalACP', 'c_C15_FabF_Act_MalACP', 'c_C17_FabF_Act_MalACP', 'c_C19_FabF_Act_MalACP', 'c_C21_FabF_Act_MalACP', 'c_C15_FabF_Act_MalACP_un',...
    'c_C17_FabF_Act_MalACP_un', 'c_C19_FabF_Act_MalACP_un', 'c_C21_FabF_Act_MalACP_un', 'c_C4_FabB_AcACP', 'c_C6_FabB_AcACP', 'c_C8_FabB_AcACP', 'c_C10_FabB_AcACP',...
    'c_C12_FabB_AcACP', 'c_C14_FabB_AcACP', 'c_C16_FabB_AcACP', 'c_C18_FabB_AcACP', 'c_C12_FabB_AcACP_un', 'c_C14_FabB_AcACP_un', 'c_C16_FabB_AcACP_un',...
    'c_C18_FabB_AcACP_un', 'c_C4_FabB_Act', 'c_C6_FabB_Act', 'c_C8_FabB_Act', 'c_C10_FabB_Act', 'c_C12_FabB_Act', 'c_C14_FabB_Act', 'c_C16_FabB_Act', 'c_C18_FabB_Act',...
    'c_C12_FabB_Act_un', 'c_C14_FabB_Act_un', 'c_C16_FabB_Act_un', 'c_C18_FabB_Act_un', 'c_C7_FabB_Act_MalACP', 'c_C9_FabB_Act_MalACP', 'c_C11_FabB_Act_MalACP',...
    'c_C13_FabB_Act_MalACP', 'c_C15_FabB_Act_MalACP', 'c_C17_FabB_Act_MalACP', 'c_C19_FabB_Act_MalACP', 'c_C21_FabB_Act_MalACP', 'c_C15_FabB_Act_MalACP_un',...
    'c_C17_FabB_Act_MalACP_un', 'c_C19_FabB_Act_MalACP_un', 'c_C21_FabB_Act_MalACP_un', 'c_C10_FabB_cis3EnAcACP', 'c_C10_FabB_Act_cis3', 'c_C13_FabB_Act_cis3MalACP',...
    'c_C4_FabH_AcACP', 'c_C6_FabH_AcACP', 'c_C8_FabH_AcACP', 'c_C10_FabH_AcACP', 'c_C12_FabH_AcACP', 'c_C14_FabH_AcACP', 'c_C16_FabH_AcACP', 'c_C18_FabH_AcACP',...
    'c_C20_FabH_AcACP', 'c_C12_FabH_AcACP_un', 'c_C14_FabH_AcACP_un', 'c_C16_FabH_AcACP_un', 'c_C18_FabH_AcACP_un', 'c_C20_FabH_AcACP_un', 'c_C6_FabH_Act_AcACP',...
    'c_C8_FabH_Act_AcACP', 'c_C10_FabH_Act_AcACP', 'c_C12_FabH_Act_AcACP', 'c_C14_FabH_Act_AcACP', 'c_C16_FabH_Act_AcACP', 'c_C18_FabH_Act_AcACP', 'c_C20_FabH_Act_AcACP',...
    'c_C22_FabH_Act_AcACP', 'c_C14_FabH_Act_AcACP_un', 'c_C16_FabH_Act_AcACP_un', 'c_C18_FabH_Act_AcACP_un', 'c_C20_FabH_Act_AcACP_un', 'c_C22_FabH_Act_AcACP_un',...
    'c_TesA_ACP', 'c_FabH_ACP', 'c_FabG_ACP', 'c_FabZ_ACP', 'c_FabI_ACP', 'c_FabF_ACP', 'c_FabA_ACP', 'c_FabB_ACP', 'c_C2_FabB_AcCoA', 'c_C2_FabB_Act', 'c_C5_FabB_Act_MalACP', 'c_C2_FabF_AcCoA',...
    'c_C2_FabF_Act', 'c_C5_FabF_Act_MalACP', 'c_C3_FabB_MalACP', 'c_C2_AcACP', 'c_C2_FabB_AcACP', 'c_C3_FabF_MalACP', 'c_C2_FabF_AcACP'};

S.num = length(S.labels); %how many diff eqs there are

% Set thioesterase fitting source
% 'Pf'; Pfleger group measurements
% 'R3M1';
% 'R3M4';
% 'Fox'; Fox group measurements
% 'Non-native'; For alternative thioesterases
S.TesA_fitting_source = 'Pf';

% Load parameter estimates
S.km_table = readtable('km_est.csv','ReadRowNames',true);
S.param_table = readtable('est_param.csv','ReadRowNames',true);

S.acp_bind = S.p_vec(12);%parameter "e"
S.ACP_inh = [S.acp_bind*2.41E-04,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,S.acp_bind*2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02];%ACP_inh lists kon,koff (in this order) for inhibitory ACP binding[FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB]
S.inhibition_kds = (1/S.acp_bind).*[4335.05,23.67;824.7,30.1;967.5,7.55;251.52,8.484;128.78,2.509];%(Acyl-ACP binding FabH Kd values in pairs (binding to FabH and FabH*), first two values are Kd for 4-12, subsequent values are 14-20)
S.inhibition_on_rates = [0.3088,1.552];%On rates of acyl-ACP binding FabH (binding to FabH and FabH*)
S.kd_fits = [S.p_vec(8),S.p_vec(9),S.p_vec(4)]; %[b2,b3,b1] Keq or Kd values that are fit
S.lin_int = S.p_vec(11);%parameter "d2"
S.lin_slope = S.p_vec(10);%parameter "d1" TesA linear free energy slope and intercept (as a function of chain length)
S.num_elong_steps = 9;%number of elongation steps

S.scaling_factor_elon = S.p_vec(2);%parameter "a2"
S.scaling_factor_init = S.p_vec(1);%parameter "a1"
S.scaling_factor_kcat = S.p_vec(5);%parameter "c2"
S.scaling_factor_kcat_init = S.p_vec(7);%parameter "c1"
S.scaling_factor_kcat_term = S.p_vec(6);%parameter "c3"
S.scaling_factor_term = S.p_vec(3);%parameter "a3"

S.kcat_scaling_fabI = [1,1,1,1,1,1,1,1,1];
S.kcat_scaling_fabH = [1,1,1,1,1,1,1,1,1];
S.kcat_scaling_fabG = [1,1,1,1,1,1,1,1,1];
S.kcat_scaling_fabA = [1,1,1,1,1,1,1,1,1];
S.kon_scaling_fabA = [0.0847,0.322,0.717,1,0.751,0.0847,0.0373,0.0373,0.0373];
S.scaling_factor_fabA_unsat = S.p_vec(13);%parameter "f"
S.kcat_scaling_fabA_unsat = [1,1,1,S.scaling_factor_fabA_unsat,1,1,1,1,1];
S.scaling_vector_fabA_unsat = [1,1,1,1,0.0358,0.0358,0.0358,0.0358,0.0358];%specificity of reaction with unsaturated acyl chains
S.scaling_factor_fabAZ_kcat = S.p_vec(14);%parameter "c4"
S.kcat_scaling_fabB = [0.855,0.855,0.975,0.967,1,0.125,0.0208,0.0208,.0208];
S.kcat_scaling_fabB_unsat = [1,1,1,1,1,0.125,0.0229,0.0229,0.0229];%specificity of reaction with unsaturated acyl chains
S.kcat_scaling_fabF = [0.914,0.914,0.901,1,0.9,0.289,0.0222,0.0222,0.0222];
S.kcat_scaling_fabF_unsat = [1,1,1,1,0.9,0.289,0.34,0.34,0.34];%specificity of reaction with unsaturated acyl chains
S.scaling_factor_fabF = S.p_vec(2);%parameter "a2" (option here to modify FabF scaling seperately)
S.kcat_scaling_fabZ = [1,1,1,1,1,1,1,1,1];
S.kon_scaling_fabZ = [0.469,1,0.296,0.372,0.2,0.0551,0.105,0.105,0.105];
S.scaling_factor_fabB_init = 0.1; % changed
S.scaling_factor_fabF_init = 1;

S.scaling_factor_kcat8_CO2 = S.p_vec(15);
S.scaling_factor_kcat10_CO2 = S.p_vec(16);
S.scaling_factor_aCoA_8 = S.p_vec(17);
S.scaling_factor_aCoA_10 = S.p_vec(18);

S.Pf_kcat_scaling = [0.0568,0.0509,0.1035,0.0158,0.25256,0.45819,1,1.221,1.5368];
S.Pf_kd_est_scaling = [473 293.9 52.986 14.79];
S.Pf_scaling = 0.519*14.79;
S.Fox_kcat_scaling = [0.47895,0.47895,0.361111111,0.110185185,1.212962963,1.453703704,1.444444444,3.064814815,3.064814815];
S.Fox_kd_est_scaling = [4815.2 4815.2 434 48];
S.Fox_scaling = 1.0416*48;
S.Non_kcat_scaling = [0,0,0,0,0,0,1,1,1];
S.Non_kd_est = [1,1,1,1,1,1,1,1,1];
S.R3M1_kcat_scaling = [0.0568 0.0509 0.25256 0.0158 1.5368 0.45819 0.25256 1.221 1.5368];
S.R3M1_kd_est = [56.91208755 35.36250007 0.92358933 1.779555549 0.093940844 0.521582239 0.92358933 0.166345313 0.093940844];
S.R3M4_kcat_scaling = [0.0568 0.0509 1 0.0158 0.25256 0.45819 0.25256 1.221 1.5368];
S.R3M4_kd_est = [56.91208755 35.36250007 0.294555191 1.779555549 0.92358933 0.521582239 0.92358933 0.166345313 0.093940844];

EC_kcat3_scaling = [1,0,0,0,0,0,0,0,0];
PP_H1_kcat3_scaling = [0.5,0,0,0,0,0,0,0,0]; % changed
PP_H2_kcat3_scaling = [0,0,0,0.4,0,0,0,0,0]; % changed

EC_kcat4_scaling = [1,1,1,1,1,1,1,1,1];
PP_1914_kcat4_scaling = [.1,.1,.1,.1,.1,.1,.1,.1,.1]; % changed
PP_2783_kcat4_scaling = [0,0,0,.025,.025,.025,.025,.025,.025]; % changed

ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

%% Figure A AcCoA
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_A = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(3) = 100; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(14) = 1300; % NADH
S.init_cond(16) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
                   0 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta1,Ca1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(1)] = Calc_Function(Ta1,Ca1,S);

[balance_conc_a1, balances_a1, total_conc_a1, carbon_a1] = mass_balance(Ca1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling;  % Using E. coli FabH

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta2,Ca2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(2)] = Calc_Function(Ta2,Ca2,S);

[balance_conc_a2, balances_a2, total_conc_a2, carbon_a2] = mass_balance(Ca2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling; % Using PP FabH1

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta3,Ca3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(3)] = Calc_Function(Ta3,Ca3,S);

[balance_conc_a3, balances_a3, total_conc_a3, carbon_a3] = mass_balance(Ca3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Ta4,Ca4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(4)] = Calc_Function(Ta4,Ca4,S);

[balance_conc_a4, balances_a4, total_conc_a4, carbon_a4] = mass_balance(Ca4,P);

% Plot
figure('Position',[500 200 382 340])
bar(rel_rate_A,'magenta')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
ylim([0 15])
ax = gca;
ax.FontSize = 18; 
text(0.1, 14, 'Acetyl-CoA','FontSize',18)


%% Figure B OcCoA
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_B = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(6) = 100; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(14) = 1300; % NADH
S.init_cond(16) = 500; % malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
                   0 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb1,Cb1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(1)] = Calc_Function(Tb1,Cb1,S);

[balance_conc_b1, balances_b1, total_conc_b1, carbon_b1] = mass_balance(Cb1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling; % Using E. coli FabH

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb2,Cb2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(2)] = Calc_Function(Tb2,Cb2,S);

[balance_conc_b2, balances_b2, total_conc_b2, carbon_b2] = mass_balance(Cb2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling; % Using PP FabH1

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb3,Cb3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(3)] = Calc_Function(Tb3,Cb3,S);

[balance_conc_b3, balances_b3, total_conc_b3, carbon_b3] = mass_balance(Cb3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tb4,Cb4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(4)] = Calc_Function(Tb4,Cb4,S);

[balance_conc_b4, balances_b4, total_conc_b4, carbon_b4] = mass_balance(Cb4,P);

% Plot
figure('Position',[500 200 382 340])
bar(rel_rate_B,'magenta')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
ylim([0 15])
ax = gca;
ax.FontSize = 18; 
text(0.1, 14, 'Octanoyl-CoA','FontSize',18)


%% Figure C No Acyl-CoA
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_C = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(14) = 1300; % NADH
S.init_cond(16) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
                   0 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc1,Cc1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(1)] = Calc_Function(Tc1,Cc1,S);

[balance_conc_c1, balances_c1, total_conc_c1, carbon_c1] = mass_balance(Cc1,P);

% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling; % Using E. coli FabH

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc2,Cc2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(2)] = Calc_Function(Tc2,Cc2,S);

[balance_conc_c2, balances_c2, total_conc_c2, carbon_c2] = mass_balance(Cc2,P);

% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling; % Using PP FabH1

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc3,Cc3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(3)] = Calc_Function(Tc3,Cc3,S);

[balance_conc_c3, balances_c3, total_conc_c3, carbon_c3] = mass_balance(Cc3,P);

% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tc4,Cc4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(4)] = Calc_Function(Tc4,Cc4,S);

[balance_conc_c4, balances_c4, total_conc_c4, carbon_c4] = mass_balance(Cc4,P);

% Plot
figure('Position',[500 200 382 340])
bar(rel_rate_C,'magenta')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabH ';'EC FabH ';'PP FabH1';'PP FabH2'])
ylim([0 15])
ax = gca;
ax.FontSize = 18; 
text(0.1, 14, 'No acyl-CoA','FontSize',18)


%% Testing what params affect the rate

% % PP 4379
% 
% x = linspace(-5,5,11);
% y = 10.^x;
% 
% S.km_table{'k3_1f','parameter_values'} = 40;
% S.km_table{'k3_2f','parameter_values'} = 7.46;
% S.km_table{'k3_3f','parameter_values'} = 5;
% S.param_table{'k3_1f','parameter_values'} = 0.002;
% S.param_table{'k3_1r','parameter_values'} = 0.08;
% S.param_table{'k3_2f','parameter_values'} = 1340;
% S.param_table{'k3_2r','parameter_values'} = 0.01;
% S.param_table{'k3_3f','parameter_values'} = 0.0043;
% S.param_table{'k3_3r','parameter_values'} = 0.0217;
% S.param_table{'k3_4f','parameter_values'} = 0.309;
% S.param_table{'k3_4r','parameter_values'} = 0.01;
% S.param_table{'k3_5f','parameter_values'} = 1.55;
% S.param_table{'k3_5r','parameter_values'} = 0.01;
% S.param_table{'kcat3','parameter_values'} = 3.13;
% S.inhibition_on_rates = [0.3088 1.552];
% S.inhibition_kds = [1.5018    0.0082; 0.2857    0.0104; 0.3352    0.0026; 0.0871    0.0029; 0.0446    0.0009];
% S.ACP_inh(1:2) = [0.0007 0.0916];
% 
% 
% S.enzyme_conc = enz_conc(2,:);
% 
% for i=1:length(y)
% 
%     S.param_table{'kcat3','parameter_values'} = y(i);
% 
%     P = Param_Function(S);
% 
%     parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
% 
%     [T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
% 
%     [~, rel_rate_A(i)] = Calc_Function(Tc3,C,S);
% 
% end
% 
% figure()
% plot(log10(y),rel_rate_A)
