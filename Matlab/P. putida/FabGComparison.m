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

%% Figure D EC_FabH AcCoA
S.kcat_scaling_fabH = EC_kcat3_scaling;  % Using E. coli FabH

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_D = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(3) = 100; % Acetyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(14) = 1300; % NADH
S.init_cond(16) = 500; % Malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 1 0 1 1 10 1 1 1;
                   0 1 1 1 1 1 10 1 1 1]; 

% No FabG
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td1,Cd1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(1)] = Calc_Function(Td1,Cd1,S);

[balance_conc_d1, balances_d1, total_conc_d1, carbon_d1] = mass_balance(Cd1,P);

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling;  % Using E. coli FabG

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td2,Cd2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(2)] = Calc_Function(Td2,Cd2,S);

[balance_conc_d2, balances_d2, total_conc_d2, carbon_d2] = mass_balance(Cd2,P);

% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling; % Using PP 1914 FabG

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td3,Cd3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(3)] = Calc_Function(Td3,Cd3,S);

[balance_conc_d3, balances_d3, total_conc_d3, carbon_d3] = mass_balance(Cd3,P);

% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling; % Using PP 2783 FabG

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Td4,Cd4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(4)] = Calc_Function(Td4,Cd4,S);

[balance_conc_d4, balances_d4, total_conc_d4, carbon_d4] = mass_balance(Cd4,P);

% Plot
figure('Position',[500 200 382 340])
bar(rel_rate_D,'cyan')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabG';'EC FabG';'PP 1914';'PP 2783'])
ylim([0 15])
ax = gca;
ax.FontSize = 18; 
text(0.1, 14, 'Acetyl-CoA','FontSize',18)
text(0.1, 12.5, '1 uM EC FabH','FontSize',18)

%% Figure E PP_FabH2 OcCoA
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;  % Using PP FabH2

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_E = zeros(1,4);

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(6) = 100; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(14) = 1300; % NADH
S.init_cond(16) = 500; % malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 10 0 1 1 10 1 1 1;
                   0 1 10 1 1 1 10 1 1 1]; 

% No FabG
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te1,Ce1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(1)] = Calc_Function(Te1,Ce1,S);

[balance_conc_e1, balances_e1, total_conc_e1, carbon_e1] = mass_balance(Ce1,P);

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using E. coli FabG

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te2,Ce2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(2)] = Calc_Function(Te2,Ce2,S);

[balance_conc_e2, balances_e2, total_conc_e2, carbon_e2] = mass_balance(Ce2,P);

% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling; % Using PP 1914 FabG

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te3,Ce3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(3)] = Calc_Function(Te3,Ce3,S);

[balance_conc_e3, balances_e3, total_conc_e3, carbon_e3] = mass_balance(Ce3,P);

% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling; % Using PP 2783 FabG

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Te4,Ce4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(4)] = Calc_Function(Te4,Ce4,S);

[balance_conc_e4, balances_e4, total_conc_e4, carbon_e4] = mass_balance(Ce4,P);

% Plot
figure('Position',[500 200 382 340])
bar(rel_rate_E,'cyan')
ylabel('Initial Rate (uM C16/m)')
xticklabels(['No FabG';'EC FabG';'PP 1914';'PP 2783'])
ylim([0 15])
ax = gca;
ax.FontSize = 18; 
text(0.1, 14, 'Octanoyl-CoA','FontSize',18)
text(0.1, 12.5, '10 uM PP FabH2','FontSize',18)

%% Figure F PP_FabH2 OcCoA
S.kcat_scaling_fabH = PP_H2_kcat3_scaling; % Using PP FabH2

% changed
S.range = [0 7200]; %2 hours (total production)

% New order from var_name code
S.init_cond = zeros(S.num,1);
S.init_cond(6) = 100; % Octanoyl-CoA
S.init_cond(12) = 10; % holo ACP
S.init_cond(13) = 1300; % NADPH
S.init_cond(14) = 1300; % NADH
S.init_cond(16) = 500; % malonyl-CoA

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 10 1 1 1 10 1 1 1]; 

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling; % Using PP FabH2

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tf1,Cf1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(1,:),~] = Calc_Function(Tf1,Cf1,S);

[balance_conc_f1, balances_f1, total_conc_f1, carbon_f1] = mass_balance(Cf1,P);

%% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling; % Using PP 1914 FabG

S.enzyme_conc = enz_conc;
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tf2,Cf2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(2,:),~] = Calc_Function(Tf2,Cf2,S);

[balance_conc_f2, balances_f2, total_conc_f2, carbon_f2] = mass_balance(Cf2,P);

% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling; % Using PP 2783 FabG

S.enzyme_conc = enz_conc;
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
tic
[Tf3,Cf3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(3,:),~] = Calc_Function(Tf3,Cf3,S);

[balance_conc_f3, balances_f3, total_conc_f3, carbon_f3] = mass_balance(Cf3,P);

F_raw_F_new(:,1:9) = F_raw_F(:,1:9);
for i=10:14
    F_raw_F_new(:,i-5) = F_raw_F(:,i-5)+F_raw_F(:,i);
end
F_raw_F_plot = F_raw_F_new(:,4:8);

% Plot
figure('Position',[500 200 382 340])
stack_labels = {'C10','C12','C14','C16','C18'};
bh = bar(F_raw_F_plot, .9, 'stacked');
xticklabels(['EC FabG';'PP 1914';'PP 2783'])
xtickangle(30)
legend(stack_labels)
ylabel('Production (uM)')
ylim([0 150])
xlim([0.4, 3.6])
set(bh, 'FaceColor', 'Flat')
colors =  [106/255, 173/255, 138/255;...
               238/255, 210/255, 148/255;...
               198/255, 96/255, 93/255;...
               145/255, 145/255, 145/255;...
               5/255, 84/255, 117/255];
colors = mat2cell(colors,ones(5,1),3);
set(bh, {'CData'}, colors)
ax = gca;
ax.FontSize = 18; 

% Plot production vs time
figure()
x = find(endsWith(S.labels, '_FA') | endsWith(S.labels, '_FA_un'));
colors = distinguishable_colors(length(x));
plot(Tf1/60,Cf1(:,x),'LineWidth',1);
lineHandles = findobj(gca, 'Type', 'Line');
set(lineHandles, {'Color'}, num2cell(colors, 2));
legend(cellfun(@(str) strrep(str(3:end), '_', ' '), {S.labels{x}}, 'UniformOutput', false),'Location','bestoutside')
ylabel('Production (uM)')
xlabel('Time (min)')
axis('tight')
ax = gca;
ax.FontSize = 18; 

%% Plotting carbon/reactants and products EC
% figure()
% for i=[90, 326]
%     plot(Tf1(1:1028)/60,Cf1(1:1028,i),"LineWidth",2)
%     hold on
% end
% legend('FabH-AcCoA', 'FabH-OcCoA')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% figure()
% for i=[3 4 8 11 318]
%     plot(Tf1(1:1028)/60,Cf1(1:1028,i),"LineWidth",2)
%     hold on
% end
% legend('AcCoA','ACP','MalCoA','CO2','OcCoA')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% figure()
% for i=[16 21 26 67 77]
%     plot(Tf1(1:1028)/60,Cf1(1:1028,i),"LineWidth",2)
%     hold on
% end
% legend('C4','C6','C8','C18','C20')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% ylim([0 0.5])
% 
% figure()
% for i=[31 37 47 57 42 52 62 72 82]
%     plot(Tf1(1:1028)/60,Cf1(1:1028,i),"LineWidth",2)
%     hold on
% end
% legend('C10','C12','C14','C16','C12 unsat','C14 unsat','C16 unsat','C18 unsat','C20 unsat')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% %% Plotting reactants and products 1914
% figure()
% for i=[3 4 8 11 318]
%     plot(Tf2(1:1071)/60,Cf2(1:1071,i),"LineWidth",2)
%     hold on
% end
% legend('AcCoA','ACP','MalCoA','CO2','OcCoA')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% figure()
% for i=[16 21 26 67 77]
%     plot(Tf2(1:1071)/60,Cf2(1:1071,i),"LineWidth",2)
%     hold on
% end
% legend('C4','C6','C8','C18','C20')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% ylim([0 0.5])
% 
% figure()
% for i=[31 37 47 57 42 52 62 72 82]
%     plot(Tf2(1:1071)/60,Cf2(1:1071,i),"LineWidth",2)
%     hold on
% end
% legend('C10','C12','C14','C16','C12 unsat','C14 unsat','C16 unsat','C18 unsat','C20 unsat')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% %% Plotting reactants and products 2783
% figure()
% for i=[3 4 8 11 318]
%     plot(Tf3(1:925)/60,Cf3(1:925,i),"LineWidth",2)
%     hold on
% end
% legend('AcCoA','ACP','MalCoA','CO2','OcCoA')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% figure()
% for i=[16 21 26]
%     plot(Tf3(1:925)/60,Cf3(1:925,i),"LineWidth",2)
%     hold on
% end
% legend('C4','C6','C8')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% figure()
% for i=[67 77]
%     plot(Tf3(1:925)/60,Cf3(1:925,i),"LineWidth",2)
%     hold on
% end
% legend('C18','C20')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% ylim([0 0.5])
% 
% figure()
% for i=[31 37 47 57 42 52 62 72 82]
%     plot(Tf3(1:925)/60,Cf3(1:925,i),"LineWidth",2)
%     hold on
% end
% legend('C10','C12','C14','C16','C12 unsat','C14 unsat','C16 unsat','C18 unsat','C20 unsat')
% ylabel('Concentration (uM)')
% xlabel('Time (min)')
% xlim([0 120])
% 
% %% Plotting EC FabG Intermediates Separately
% 
% label = {'FabG-NADPH', 'C4 FabG-NADPH-B-ketoacyl-ACPs', 'C6 FabG-NADPH-B-ketoacyl-ACPs',...
% 'C8 FabG-NADPH-B-ketoacyl-ACPs', 'C10 FabG-NADPH-B-ketoacyl-ACPs', 'C12 FabG-NADPH-B-ketoacyl-ACPs',...
% 'C14 FabG-NADPH-B-ketoacyl-ACPs', 'C16 FabG-NADPH-B-ketoacyl-ACPs', 'C18 FabG-NADPH-B-ketoacyl-ACPs',...
% 'C20 FabG-NADPH-B-ketoacyl-ACPs', 'C12 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C14 FabG-NADPH-B-ketoacyl-ACPs Unsat',...
% 'C16 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C18 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C20 FabG-NADPH-B-ketoacyl-ACPs Unsat',...
% FabG-ACP'};
% j=1;
% for i=[93 97 112 127 142 157 172 187 202 217 226 241 256 271 286 294]
%     figure()
%     plot(Tf1,Cf1(:,i))
%     title(label(j))
%     j=j+1;
% end
% 
% %% Plotting EC FabG Intermediates Sat/Unsat
% figure()
% for i=[97 112 127 142 157 172 187 202 217]
%     plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C4','C6','C8','C10','C12','C14','C16','C18','C20')
% title('FabG-NADPH-B-ketoacyl-ACPs')
% 
% figure()
% for i=[226 241 256 271 286]
%     plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C12','C14','C16','C18','C20')
% title('Unsat FabG-NADPH-B-ketoacyl-ACPs')
% 
% %% Plotting EC FabG Intermediates by length
% 
% figure(1)
% for i=[97 112 127]
%     plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C4 TesA=10','C6 TesA=10','C8 TesA=10','C4 TesA=1','C6 TesA=1','C8 TesA=1')
% title('FabG-NADPH-B-ketoacyl-ACPs')
% xlabel('Time (seconds)')
% ylabel('Concentration (uM)')
% axis("padded")
% 
% figure(2)
% for i=[142]
%     plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C10 TesA=10','C10 TesA=1')
% title('FabG-NADPH-B-ketoacyl-ACPs')
% xlabel('Time (seconds)')
% ylabel('Concentration (uM)')
% axis("padded")
% 
% figure(3)
% for i=[157 226]
%         plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C12 TesA=10','C12 Unsat TesA=10','C12 TesA=1','C12 Unsat TesA=1')
% title('FabG-NADPH-B-ketoacyl-ACPs')
% xlabel('Time (seconds)')
% ylabel('Concentration (uM)')
% axis("padded")
% 
% figure(4)
% for i=[172 241 187 256]
%         plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C14 TesA=10','C14 Unsat TesA=10','C16 TesA=10','C16 Unsat TesA=10','C14 TesA=1','C14 Unsat TesA=1','C16 TesA=1','C16 Unsat TesA=1')
% title('FabG-NADPH-B-ketoacyl-ACPs')
% xlabel('Time (seconds)')
% ylabel('Concentration (uM)')
% axis("padded")
% 
% figure(5)
% for i=[202 271 217 286]
%         plot(Tf1,Cf1(:,i))
%     hold on
% end
% legend('C18 TesA=10','C18 Unsat TesA=10','C20 TesA=10','C20 Unsat TesA=10','C18 TesA=1','C18 Unsat TesA=1','C20 TesA=1','C20 Unsat TesA=1')
% title('FabG-NADPH-B-ketoacyl-ACPs')
% xlabel('Time (seconds)')
% ylabel('Concentration (uM)')
% axis("padded")
% 
% %% Testing what params affect the rate
% 
% % PP 1914
% 
% x = linspace(-5,5,11);
% y = 10.^x;
% 
% S.km_table{'k4_1f','parameter_values'} = 10;
% S.km_table{'k4_2f','parameter_values'} = 1; %1
% S.param_table{'k4_1f','parameter_values'} = 0.0079; %No effect
% S.param_table{'k4_1r','parameter_values'} = 0.0793;
% S.param_table{'k4_2f','parameter_values'} = 0.0013; %No effect
% S.param_table{'k4_2r','parameter_values'} = 0.0217;
% S.param_table{'kcat4','parameter_values'} = y(end); %1
% S.ACP_inh(3:4) = [0.0024,0.0916];
% S.kcat_scaling_fabG = [1,1,1,1,1,1,1,1,1];
% 
% S.enzyme_conc = enz_conc(2,:);
% 
% for i=1:length(y)
% 
%     S.kcat_scaling_fabG = [y(i),1,1,1,1,1,1,1,1];
% 
%     P = Param_Function(S);
% 
%     parameterized_ODEs = @(t,c) ODE_Function(t,c,P);
% 
%     [T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
% 
%     rel_rate_D(i) = Calc_Function(T,C,S);
% 
% end
% 
% figure()
% plot(log10(y),rel_rate_D)