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
S.FA_dist = [4,6,8,10,12,12,14,14,16,16,18,18,20,20];

S.labels = {'s1';'s2';'s3';'s6';'s7';'s8';'p1';'p2';'p3';'p4';'p5';...
    'Q_4';'M_4';'R_4';'T_4';'F_4';'Q_6';'M_6';'R_6';'T_6';'F_6';...
    'Q_8';'M_8';'R_8';'T_8';'F_8';'Q_10';'M_10';'R_10';'T_10';'F_10';'R_10_un';...
    'Q_12';'M_12';'R_12';'T_12';'F_12';'Q_12_un';'M_12_un';'R_12_un';'T_12_un';'F_12_un';...
    'Q_14';'M_14';'R_14';'T_14';'F_14';'Q_14_un';'M_14_un';'R_14_un';'T_14_un';'F_14_un';...
    'Q_16';'M_16';'R_16';'T_16';'F_16';'Q_16_un';'M_16_un';'R_16_un';'T_16_un';'F_16_un';...
    'Q_18';'M_18';'R_18';'T_18';'F_18';'Q_18_un';'M_18_un';'R_18_un';'T_18_un';'F_18_un';...
    'Q_20';'M_20';'R_20';'T_20';'F_20';'Q_20_un';'M_20_un';'R_20_un';'T_20_un';'F_20_un';...
    'e1s1';'e1s1s2';'e1act';'e1acts3';'e2p2';'e2act';'e2acts6';'e3s3';'e3acts3';'e3acts3p4';'e4s7';'e6s8';...
    'e3T_4';'e3actT_4';'e4s7Q_4';'e5M_4';'e5R_4';'e6s8R_4';'e7T_4';'e8T_4';'e8act4';'e8act4p4';'e9M_4';'e9R_4';'e10T_4';'e10act4';'e10act4p4';...
    'e3T_6';'e3actT_6';'e4s7Q_6';'e5M_6';'e5R_6';'e6s8R_6';'e7T_6';'e8T_6';'e8act6';'e8act6p4';'e9M_6';'e9R_6';'e10T_6';'e10act6';'e10act6p4';...
    'e3T_8';'e3actT_8';'e4s7Q_8';'e5M_8';'e5R_8';'e6s8R_8';'e7T_8';'e8T_8';'e8act8';'e8act8p4';'e9M_8';'e9R_8';'e10T_8';'e10act8';'e10act8p4';...
    'e3T_10';'e3actT_10';'e4s7Q_10';'e5M_10';'e5R_10';'e6s8R_10';'e7T_10';'e8T_10';'e8act10';'e8act10p4';'e9M_10';'e9R_10';'e10T_10';'e10act10';'e10act10p4';...
    'e3T_12';'e3actT_12';'e4s7Q_12';'e5M_12';'e5R_12';'e6s8R_12';'e7T_12';'e8T_12';'e8act12';'e8act12p4';'e9M_12';'e9R_12';'e10T_12';'e10act12';'e10act12p4';...
    'e3T_14';'e3actT_14';'e4s7Q_14';'e5M_14';'e5R_14';'e6s8R_14';'e7T_14';'e8T_14';'e8act14';'e8act14p4';'e9M_14';'e9R_14';'e10T_14';'e10act14';'e10act14p4';...
    'e3T_16';'e3actT_16';'e4s7Q_16';'e5M_16';'e5R_16';'e6s8R_16';'e7T_16';'e8T_16';'e8act16';'e8act16p4';'e9M_16';'e9R_16';'e10T_16';'e10act16';'e10act16p4';...
    'e3T_18';'e3actT_18';'e4s7Q_18';'e5M_18';'e5R_18';'e6s8R_18';'e7T_18';'e8T_18';'e8act18';'e8act18p4';'e9M_18';'e9R_18';'e10T_18';'e10act18';'e10act18p4';...
    'e3T_20';'e3actT_20';'e4s7Q_20';'e5M_20';'e5R_20';'e6s8R_20';'e7T_20';'e9M_20';'e9R_20';...
    'e3T_12_un';'e3actT_12_un';'e4s7Q_12_un';'e5M_12_un';'e5R_12_un';'e6s8R_12_un';'e7T_12_un';'e8T_12_un';'e8act12_un';'e8act12_unp4';'e9M_12_un';'e9R_12_un';'e10T_12_un';'e10act12_un';'e10act12_unp4';...
    'e3T_14_un';'e3actT_14_un';'e4s7Q_14_un';'e5M_14_un';'e5R_14_un';'e6s8R_14_un';'e7T_14_un';'e8T_14_un';'e8act14_un';'e8act14_unp4';'e9M_14_un';'e9R_14_un';'e10T_14_un';'e10act14_un';'e10act14_unp4';...
    'e3T_16_un';'e3actT_16_un';'e4s7Q_16_un';'e5M_16_un';'e5R_16_un';'e6s8R_16_un';'e7T_16_un';'e8T_16_un';'e8act16_un';'e8act16_unp4';'e9M_16_un';'e9R_16_un';'e10T_16_un';'e10act16_un';'e10act16_unp4';...
    'e3T_18_un';'e3actT_18_un';'e4s7Q_18_un';'e5M_18_un';'e5R_18_un';'e6s8R_18_un';'e7T_18_un';'e8T_18_un';'e8act18_un';'e8act18_unp4';'e9M_18_un';'e9R_18_un';'e10T_18_un';'e10act18_un';'e10act18_unp4';...
    'e3T_20_un';'e3actT_20_un';'e4s7Q_20_un';'e5M_20_un';'e5R_20_un';'e6s8R_20_un';'e7T_20_un';'e9M_20_un';'e9R_20_un';'e3s6';'e4s6';'e5s6';'e6s6';'e7s6';'e8s6';'e9s6';'e10s6';...
    'e9R_10_un';'e10R_10_un';'e10act10_un';'e10act10_unp4';...
    'e10s3';'e10act2';'e10act2p4';'e8s3';'e8act2';'e8act2p4';...
    'e10p4';'T_2';'e10T_2';'e8p4';'e8T_2';...
    's9';'s10';'s11';'s12';'s13';'s14';'s15';'s16';...
    'e3s9';'e3s10';'e3s11';'e3s12';'e3s13';'e3s14';'e3s15';'e3s16';...
    'e3acts9';'e3acts10';'e3acts11';'e3acts12';'e3acts13';'e3acts14';'e3acts15';'e3acts16';...
    'e3acts9p4';'e3acts10p4';'e3acts11p4';'e3acts12p4';'e3acts13p4';'e3acts14p4';'e3acts15p4';'e3acts16p4';...
};
S.num = length(S.labels); %how many diff eqs there are

% Set thioesterase fitting source
% 'Pf'; Pfleger group measurements
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
S.kcat_scaling_fabF = [0.914,0.914,0.901,1,0.9,0.289,0.0222,0.0222,.0222];
S.kcat_scaling_fabF_unsat = [1,1,1,1,0.9,0.289,0.34,0.34,0.34];%specificity of reaction with unsaturated acyl chains
S.scaling_factor_fabF = S.p_vec(2);%parameter "a2" (option here to modify FabF scaling seperately)
S.kcat_scaling_fabZ = [1,1,1,1,1,1,1,1,1];
S.kon_scaling_fabZ = [0.469,1,0.296,0.372,0.2,0.0551,0.105,0.105,0.105];
S.scaling_factor_fabB_init = 1;
S.scaling_factor_fabF_init = 0.1;

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

%% Figure A AcCoA
EC_kcat3_scaling = [1,1,1,0.5,1,1,1,1,1];
PP_H1_kcat3_scaling = [0.3,1,1,1,1,1,1,1,1];
PP_H2_kcat3_scaling = [0,1,1,0.5,1,1,1,1,1];

load('JpMat.mat','JpMatPrime')
%ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'JPattern',JpMatPrime,'Vectorized','on');
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_A = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(3) = 100; %s3 (Acetyl-CoA)
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1300; %s7 (NADPH)
S.init_cond(6) = 1300; %s8 (NADH) not in figure description
S.init_cond(8) = 500; %p2 (malonyl-CoA)

% FabH parameters for reference
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

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
            0 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Ta1,Ca1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(1)] = Calc_Function(Ta1,Ca1,S);


% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Ta2,Ca2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(2)] = Calc_Function(Ta2,Ca2,S);


% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Ta3,Ca3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(3)] = Calc_Function(Ta3,Ca3,S);


% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Ta4,Ca4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_A(4)] = Calc_Function(Ta4,Ca4,S);


% Plot
figure()
bar(rel_rate_A,'magenta')
ylabel('Initial Rate (uM C16 equiv. per min)')
xticklabels(['No FabH';'EC FabH';'PP 4379';'PP 4545'])
ylim([0 15])

%% Figure B OcCoA

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_B = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(318) = 100; % (Octanoyl-CoA)
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1300; %s7 (NADPH)
S.init_cond(6) = 1300; %s8 (NADH) not in figure description
S.init_cond(8) = 500; %p2 (malonyl-CoA)

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
                   0 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tb1,Cb1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(1)] = Calc_Function(Tb1,Cb1,S);


% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tb2,Cb2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(2)] = Calc_Function(Tb2,Cb2,S);


% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tb3,Cb3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(3)] = Calc_Function(Tb3,Cb3,S);


% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tb4,Cb4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_B(4)] = Calc_Function(Tb4,Cb4,S);

% Plot
figure()
bar(rel_rate_B,'magenta')
ylabel('Initial Rate (uM C16 equiv. per min)')
xticklabels(['No FabH';'EC FabH';'PP 4379';'PP 4545'])
ylim([0 15])

%% Figure C No Acyl-CoA

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_C = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1300; %s7 (NADPH)
S.init_cond(6) = 1300; %s8 (NADH) not in figure description
S.init_cond(8) = 500; %p2 (malonyl-CoA)

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 0 1 1 1 10 1 1 1;
            0 1 1 1 1 1 10 1 1 1]; 

% No FabH
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tc1,Cc1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(1)] = Calc_Function(Tc1,Cc1,S);


% EC FabH
S.kcat_scaling_fabH = EC_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tc2,Cc2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(2)] = Calc_Function(Tc2,Cc2,S);


% PP FabH1
S.kcat_scaling_fabH = PP_H1_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tc3,Cc3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(3)] = Calc_Function(Tc3,Cc3,S);


% PP FabH2
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tc4,Cc4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_C(4)] = Calc_Function(Tc4,Cc4,S);


% Plot
figure()
bar(rel_rate_C,'magenta')
ylabel('Initial Rate (uM C16 equiv. per min)')
xticklabels(['No FabH';'EC FabH';'PP 4379';'PP 4545'])
ylim([0 15])

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
%     parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
% 
%     [T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
% 
%     [~, rel_rate_A(i)] = Calc_Function(Tc3,C,S);
% 
% end
% 
% figure()
% plot(log10(y),rel_rate_A)
