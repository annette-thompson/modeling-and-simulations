function S = set_vars()

% Set up structure to store variables
S = struct;

% p_vec = [a1 a2 a3 b1 c2 c3 c1 b2 b3 d1 d2 e f c4 x1 x2 x3 x4];
S.p_vec = [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
    4.645634978 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673 132.8499358...
    2180.050007 0.539756276	0.053673263	34.49718991	11.15058888];

% Potential FA lengths for calculating production profile
% Numbers repeated are saturated and unsaturated
S.FA_dist = [4,6,8,10,12,14,16,18,20,12,14,16,18,20];

S.labels = {'c_ATP', 'c_C1_Bicarbonate', 'c_C2_AcCoA', 'c_C4_SucCoA', 'c_C6_HexCoA', 'c_C8_OcCoA', 'c_C10_DecCoA', 'c_C12_LauCoA', 'c_C14_EthCoA',...
    'c_C16_PalCoA', 'c_C18_OcDecCoA', 'c_ACP', 'c_NADPH', 'c_NADP', 'c_NADH', 'c_NAD', 'c_Fd', 'c_Fd2', 'c_ADP', 'c_C3_MalCoA', 'c_CoA',...
    'c_C3_MalACP', 'c_C1_CO2', 'c_C4_BKeACP', 'c_C6_BKeACP', 'c_C8_BKeACP', 'c_C10_BKeACP', 'c_C12_BKeACP', 'c_C14_BKeACP', 'c_C16_BKeACP',...
    'c_C18_BKeACP', 'c_C20_BKeACP',...
    'c_C4_BHyAcACP', 'c_C6_BHyAcACP', 'c_C8_BHyAcACP', 'c_C10_BHyAcACP', 'c_C12_BHyAcACP', 'c_C14_BHyAcACP', 'c_C16_BHyAcACP',...
    'c_C18_BHyAcACP', 'c_C20_BHyAcACP', 'c_C4_EnAcACP', 'c_C6_EnAcACP', 'c_C8_EnAcACP', 'c_C10_EnAcACP', 'c_C12_EnAcACP', 'c_C14_EnAcACP',...
    'c_C16_EnAcACP', 'c_C18_EnAcACP', 'c_C20_EnAcACP', 'c_C4_AcACP', 'c_C6_AcACP', 'c_C8_AcACP', 'c_C10_AcACP', 'c_C12_AcACP',...
    'c_C14_AcACP', 'c_C16_AcACP', 'c_C18_AcACP', 'c_C20_AcACP', 'c_C12_AcACP_un', 'c_C14_AcACP_un', 'c_C16_AcACP_un', 'c_C18_AcACP_un',...
    'c_C20_AcACP_un', 'c_C4_FA', 'c_C6_FA', 'c_C8_FA', 'c_C10_FA', 'c_C12_FA', 'c_C14_FA', 'c_C16_FA', 'c_C18_FA', 'c_C20_FA', 'c_C12_FA_un',...
    'c_C14_FA_un', 'c_C16_FA_un', 'c_C18_FA_un', 'c_C20_FA_un', 'c_BC_ATP', 'c_C1_BC_ATP_HCO3', 'c_C1_BC_Pi_HCO3',...
    'c_C1_BC_Pi_HCO3_BCCP_Biotin', 'c_C1_BCCP_Biotin_CO2', 'c_C1_CT_BCCP_Biotin_CO2', 'c_C1_CT_Act', 'c_C3_CT_Act_AcCoA',...
    'c_C3_MCMT_MalCoA', 'c_C3_MCMT_Act', 'c_C3_MCMT_Act_ACP', 'c_C2_KASIII_CoA', 'c_C4_KASIII_CoA', 'c_C6_KASIII_CoA', 'c_C8_KASIII_CoA',...
    'c_C10_KASIII_CoA', 'c_C12_KASIII_CoA', 'c_C14_KASIII_CoA', 'c_C16_KASIII_CoA', 'c_C18_KASIII_CoA', 'c_C2_KASIII_Act', 'c_C4_KASIII_Act',...
    'c_C6_KASIII_Act', 'c_C8_KASIII_Act', 'c_C10_KASIII_Act', 'c_C12_KASIII_Act', 'c_C14_KASIII_Act', 'c_C16_KASIII_Act', 'c_C18_KASIII_Act',...
    'c_C5_KASIII_Act_MalACP', 'c_C7_KASIII_Act_MalACP', 'c_C9_KASIII_Act_MalACP', 'c_C11_KASIII_Act_MalACP', 'c_C13_KASIII_Act_MalACP',...
    'c_C15_KASIII_Act_MalACP', 'c_C17_KASIII_Act_MalACP', 'c_C19_KASIII_Act_MalACP', 'c_C21_KASIII_Act_MalACP', 'c_KAR_NADPH',...
    'c_C4_KAR_NADPH_BKeACP', 'c_C6_KAR_NADPH_BKeACP', 'c_C8_KAR_NADPH_BKeACP', 'c_C10_KAR_NADPH_BKeACP',......
    'c_C12_KAR_NADPH_BKeACP', 'c_C14_KAR_NADPH_BKeACP', 'c_C16_KAR_NADPH_BKeACP', 'c_C18_KAR_NADPH_BKeACP',...
    'c_C20_KAR_NADPH_BKeACP', 'c_C4_HAD_BHyAcACP', 'c_C6_HAD_BHyAcACP', 'c_C8_HAD_BHyAcACP', 'c_C10_HAD_BHyAcACP',...
    'c_C12_HAD_BHyAcACP', 'c_C14_HAD_BHyAcACP', 'c_C16_HAD_BHyAcACP', 'c_C18_HAD_BHyAcACP', 'c_C20_HAD_BHyAcACP',...
    'c_C4_HAD_EnAcACP', 'c_C6_HAD_EnAcACP', 'c_C8_HAD_EnAcACP', 'c_C10_HAD_EnAcACP', 'c_C12_HAD_EnAcACP', 'c_C14_HAD_EnAcACP',...
    'c_C16_HAD_EnAcACP', 'c_C18_HAD_EnAcACP', 'c_C20_HAD_EnAcACP', 'c_SAD_Fd', 'c_C12_SAD_Fd_AcACP', 'c_C14_SAD_Fd_AcACP',...
    'c_C16_SAD_Fd_AcACP', 'c_C18_SAD_Fd_AcACP', 'c_C20_SAD_Fd_AcACP', 'c_ER_NADH', 'c_C4_ER_NADH_EnAcACP', 'c_C6_ER_NADH_EnAcACP',...
    'c_C8_ER_NADH_EnAcACP', 'c_C10_ER_NADH_EnAcACP', 'c_C12_ER_NADH_EnAcACP', 'c_C14_ER_NADH_EnAcACP', 'c_C16_ER_NADH_EnAcACP',...
    'c_C18_ER_NADH_EnAcACP', 'c_C20_ER_NADH_EnAcACP', 'c_C4_FatA_AcACP', 'c_C6_FatA_AcACP', 'c_C8_FatA_AcACP', 'c_C10_FatA_AcACP',...
    'c_C12_FatA_AcACP', 'c_C14_FatA_AcACP', 'c_C16_FatA_AcACP', 'c_C18_FatA_AcACP', 'c_C20_FatA_AcACP', 'c_C12_FatA_AcACP_un',...
    'c_C14_FatA_AcACP_un', 'c_C16_FatA_AcACP_un', 'c_C18_FatA_AcACP_un', 'c_C20_FatA_AcACP_un', 'c_C4_KASI_AcACP', 'c_C6_KASI_AcACP',...
    'c_C8_KASI_AcACP', 'c_C10_KASI_AcACP', 'c_C12_KASI_AcACP', 'c_C14_KASI_AcACP', 'c_C16_KASI_AcACP', 'c_C18_KASI_AcACP',...
    'c_C4_KASI_Act', 'c_C6_KASI_Act', 'c_C8_KASI_Act', 'c_C10_KASI_Act', 'c_C12_KASI_Act', 'c_C14_KASI_Act', 'c_C16_KASI_Act', 'c_C18_KASI_Act',...
    'c_C7_KASI_Act_MalACP', 'c_C9_KASI_Act_MalACP', 'c_C11_KASI_Act_MalACP', 'c_C13_KASI_Act_MalACP', 'c_C15_KASI_Act_MalACP',...
    'c_C17_KASI_Act_MalACP', 'c_C19_KASI_Act_MalACP', 'c_C21_KASI_Act_MalACP', 'c_C4_KASII_AcACP', 'c_C6_KASII_AcACP', 'c_C8_KASII_AcACP',...
    'c_C10_KASII_AcACP', 'c_C12_KASII_AcACP', 'c_C14_KASII_AcACP', 'c_C16_KASII_AcACP', 'c_C18_KASII_AcACP', 'c_C4_KASII_Act',...
    'c_C6_KASII_Act', 'c_C8_KASII_Act', 'c_C10_KASII_Act', 'c_C12_KASII_Act', 'c_C14_KASII_Act', 'c_C16_KASII_Act', 'c_C18_KASII_Act',...
    'c_C7_KASII_Act_MalACP', 'c_C9_KASII_Act_MalACP', 'c_C11_KASII_Act_MalACP', 'c_C13_KASII_Act_MalACP', 'c_C15_KASII_Act_MalACP',...
    'c_C17_KASII_Act_MalACP', 'c_C19_KASII_Act_MalACP', 'c_C21_KASII_Act_MalACP', 'c_C4_KASIII_AcACP', 'c_C6_KASIII_AcACP',...
    'c_C8_KASIII_AcACP', 'c_C10_KASIII_AcACP', 'c_C12_KASIII_AcACP', 'c_C14_KASIII_AcACP', 'c_C16_KASIII_AcACP', 'c_C18_KASIII_AcACP',...
    'c_C20_KASIII_AcACP', 'c_C12_KASIII_AcACP_un', 'c_C14_KASIII_AcACP_un', 'c_C16_KASIII_AcACP_un', 'c_C18_KASIII_AcACP_un',...
    'c_C20_KASIII_AcACP_un', 'c_C6_KASIII_Act_AcACP', 'c_C8_KASIII_Act_AcACP', 'c_C10_KASIII_Act_AcACP', 'c_C12_KASIII_Act_AcACP',...
    'c_C14_KASIII_Act_AcACP', 'c_C16_KASIII_Act_AcACP', 'c_C18_KASIII_Act_AcACP', 'c_C20_KASIII_Act_AcACP', 'c_C22_KASIII_Act_AcACP',...
    'c_C14_KASIII_Act_AcACP_un', 'c_C16_KASIII_Act_AcACP_un', 'c_C18_KASIII_Act_AcACP_un', 'c_C20_KASIII_Act_AcACP_un',...
    'c_C22_KASIII_Act_AcACP_un', 'c_FatA_ACP', 'c_KASIII_ACP', 'c_KAR_ACP', 'c_HAD_ACP', 'c_ER_ACP', 'c_KASI_ACP', 'c_SAD_ACP', 'c_KASII_ACP',...
    'c_C2_KASII_AcCoA', 'c_C2_KASII_Act', 'c_C5_KASII_Act_MalACP', 'c_C2_KASI_AcCoA', 'c_C2_KASI_Act', 'c_C5_KASI_Act_MalACP',...
    'c_C3_KASII_MalACP', 'c_C2_AcACP', 'c_C2_KASII_AcACP', 'c_C3_KASI_MalACP', 'c_C2_KASI_AcACP'};

S.num = length(S.labels); %how many diff eqs there are
S.num_elong_steps = 9;%number of elongation steps

% Set thioesterase fitting source
% 'Pf'; Pfleger group measurements
% 'R3M1';
% 'R3M4';
% 'Fox'; Fox group measurements
% 'Non-native'; For alternative thioesterases
S.FatA_fitting_source = 'Pf';

% Load parameter estimates
S.km_table = readtable('km_est.csv','ReadRowNames',true);
S.param_table = readtable('est_param.csv','ReadRowNames',true);

S.acp_bind = S.p_vec(12);%parameter "e"
S.ACP_inh = [S.acp_bind*2.41E-04,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,S.acp_bind*2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02];%ACP_inh lists kon,koff (in this order) for inhibitory ACP binding[KASIII,KAR,HAD,ER,FatA,KASI,SAD,KASII]
S.kd_fits = [S.p_vec(8),S.p_vec(9),S.p_vec(4)]; %[b2,b3,b1] Keq or Kd values that are fit
S.scaling_factor_FabA_unsat = S.p_vec(13);% unused, f, scaled FabA desaturating C10

S.kcat_scaling_FatA = [1,1,1,1,1,1,1,1,1];
S.lin_int = S.p_vec(11);%parameter "d2"
S.lin_slope = S.p_vec(10);%parameter "d1" FatA linear free energy slope and intercept (as a function of chain length)
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

S.scaling_factor_elon = S.p_vec(2);%parameter "a2"
S.scaling_factor_init = S.p_vec(1);%parameter "a1"
S.scaling_factor_kcat = S.p_vec(5);%parameter "c2"
S.scaling_factor_kcat_init = S.p_vec(7);%parameter "c1"
S.scaling_factor_kcat_term = S.p_vec(6);%parameter "c3"
S.scaling_factor_term = S.p_vec(3);%parameter "a3"

S.scaling_factor_kcat8_CO2 = S.p_vec(15);
S.scaling_factor_kcat10_CO2 = S.p_vec(16);
S.scaling_factor_aCoA_8 = S.p_vec(17);
S.scaling_factor_aCoA_10 = S.p_vec(18);

S.kcat_scaling_KASIII = [1,1,1,1,1,1,1,1,1];
S.inhibition_kds = (1/S.acp_bind).*[4335.05,23.67;824.7,30.1;967.5,7.55;251.52,8.484;128.78,2.509];% (Acyl-ACP binding KASIII Kd values in pairs (binding to KASIII and KASIII*), first two values are Kd for 4-12, subsequent values are 14-20)
S.inhibition_on_rates = [0.3088,1.552];

S.kcat_scaling_KAR = [1,1,1,1,1,1,1,1,1];

S.kcat_scaling_HAD = [1,1,1,1,1,1,1,1,1];
S.scaling_factor_HAD_kcat = S.p_vec(14);% scales kcat5, k5_2r (kcat5 reverse)
S.kon_scaling_HAD = [0.469,1,0.296,0.372,0.2,0.0551,0.105,0.105,0.105]; % scales k5_1f, k5_3r

S.kcat_scaling_ER = [1,1,1,1,1,1,1,1,1];

S.kcat_scaling_KASI = [0.914,0.914,0.901,1,0.9,0.289,0.222,0.222,0.222];
S.scaling_factor_KASI = S.p_vec(2);% scales k8_1r/f
S.scaling_factor_KASI_init = 1; % scales k8_4f

S.kcat_scaling_KASII = [0.855,0.855,0.975,0.967,1,0.125,0.0208,0.0208,0.0208];
S.scaling_factor_KASI = S.p_vec(2);% scales k10_1r/f
S.scaling_factor_KASII_init = 1; % scales k10_4f

S.kcat_scaling_SAD = [1,1,1,1,1,1,1,1,1]; % scales kcat9
S.kon_scaling_SAD = [0.0847,0.322,0.717,1,0.751,0.0847,0.0373,0.0373,0.0373]; % scales k9_2