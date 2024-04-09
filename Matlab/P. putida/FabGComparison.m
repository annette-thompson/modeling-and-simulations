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
S.scaling_factor_fabB_init = 0.1;
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

%% Figure D EC_FabH AcCoA
EC_kcat3_scaling = [1,0,0,0,0,0,0,0,0];
PP_H1_kcat3_scaling = [0.5,0,0,0,0,0,0,0,0];
PP_H2_kcat3_scaling = [0,0,0,0.4,0,0,0,0,0];

S.kcat_scaling_fabH = EC_kcat3_scaling;

EC_kcat4_scaling = [1,1,1,1,1,1,1,1,1];
PP_1914_kcat4_scaling = [.1,.1,.1,.1,.1,.1,.1,.1,.1];
PP_2783_kcat4_scaling = [0,0,0,.025,.025,.025,.025,.025,.025];

load('JpMat.mat','JpMatPrime')
%ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'JPattern',JpMatPrime,'Vectorized','on');
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_D = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(3) = 100; %s3 (Acetyl-CoA)
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1300; %s7 (NADPH)
S.init_cond(6) = 1300; %s8 (NADH) not in figure description
S.init_cond(8) = 500; %p2 (malonyl-CoA)

% FabG parameters for reference
% S.km_table{'k4_1f','parameter_values'} = 10;
% S.km_table{'k4_2f','parameter_values'} = 17;
% S.param_table{'k4_1f','parameter_values'} = 0.0079;
% S.param_table{'k4_1r','parameter_values'} = 0.0793;
% S.param_table{'k4_2f','parameter_values'} = 0.0013;
% S.param_table{'k4_2r','parameter_values'} = 0.0217;
% S.param_table{'kcat4','parameter_values'} = 0.59;
% S.ACP_inh(3:4) = [0.0024,0.0916];

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 1 0 1 1 10 1 1 1;
            0 1 1 1 1 1 10 1 1 1]; 

% No FabG
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Td1,Cd1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(1)] = Calc_Function(Td1,Cd1,S);


% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Td2,Cd2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(2)] = Calc_Function(Td2,Cd2,S);


% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling;

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Td3,Cd3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(3)] = Calc_Function(Td3,Cd3,S);


% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling;

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Td4,Cd4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_D(4)] = Calc_Function(Td4,Cd4,S);


%% Plot
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
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;

S.range = [0 150]; %2.5 mins (initial rate)

rel_rate_E = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(318) = 100; % (OcCoA)
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1300; %s7 (NADPH)
S.init_cond(6) = 1300; %s8 (NADH) not in figure description
S.init_cond(8) = 500; %p2 (malonyl-CoA)

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 10 0 1 1 10 1 1 1;
            0 1 10 1 1 1 10 1 1 1]; 

% No FabG
S.enzyme_conc = enz_conc(1,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Te1,Ce1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(1)] = Calc_Function(Te1,Ce1,S);


% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling;

S.enzyme_conc = enz_conc(2,:);

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Te2,Ce2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(2)] = Calc_Function(Te2,Ce2,S);


% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling;

S.enzyme_conc = enz_conc(2,:);
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Te3,Ce3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(3)] = Calc_Function(Te3,Ce3,S);


% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling;

S.enzyme_conc = enz_conc(2,:);
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Te4,Ce4] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[~, rel_rate_E(4)] = Calc_Function(Te4,Ce4,S);

%% Plot
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
S.kcat_scaling_fabH = PP_H2_kcat3_scaling;

S.range = [0 720]; %12 mins (total production)

rel_rate_E = zeros(1,4);

S.init_cond = zeros(S.num,1);
S.init_cond(318) = 100; % (OcCoA)
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1300; %s7 (NADPH)
S.init_cond(6) = 1300; %s8 (NADH) not in figure description
S.init_cond(8) = 500; %p2 (malonyl-CoA)

% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
enz_conc = [0 1 10 1 1 1 10 1 1 1]; 

% EC FabG
S.kcat_scaling_fabG = EC_kcat4_scaling;

S.enzyme_conc = enz_conc;

P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tf1,Cf1] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(1,:),~] = Calc_Function(Tf1,Cf1,S);


% PP 1914
S.kcat_scaling_fabG = PP_1914_kcat4_scaling;

S.enzyme_conc = enz_conc;
    
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tf2,Cf2] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(2,:),~] = Calc_Function(Tf2,Cf2,S);


% PP 2783
S.kcat_scaling_fabG = PP_2783_kcat4_scaling;

S.enzyme_conc = enz_conc;
 
P = Param_Function(S);

parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
tic
[Tf3,Cf3] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
toc
[F_raw_F(3,:),~] = Calc_Function(Tf3,Cf3,S);

F_raw_F_new(:,1:4) = F_raw_F(:,1:4);
j=5;
for i=[5 7 9 11 13]
    F_raw_F_new(:,j) = F_raw_F(:,i)+F_raw_F(:,i+1);
    j=j+1;
end
F_raw_F_plot = F_raw_F_new(:,4:8);

%% Plot
figure('Position',[500 200 382 340])
stack_labels = {'C10','C12','C14','C16','C18'};
bh = bar(F_raw_F_plot,.9,'stacked');
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

%% Plotting EC FabG Intermediates Separately
% 
% label = {'FabG-NADPH', 'C4 FabG-NADPH-B-ketoacyl-ACPs', 'C6 FabG-NADPH-B-ketoacyl-ACPs', 'C8 FabG-NADPH-B-ketoacyl-ACPs', 'C10 FabG-NADPH-B-ketoacyl-ACPs', 'C12 FabG-NADPH-B-ketoacyl-ACPs', 'C14 FabG-NADPH-B-ketoacyl-ACPs', 'C16 FabG-NADPH-B-ketoacyl-ACPs', 'C18 FabG-NADPH-B-ketoacyl-ACPs', 'C20 FabG-NADPH-B-ketoacyl-ACPs', 'C12 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C14 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C16 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C18 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'C20 FabG-NADPH-B-ketoacyl-ACPs Unsat', 'FabG-ACP'};
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

%% Testing what params affect the rate

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
%     parameterized_ODEs = @(t,c) ODE_Function(t,c,P,S.num);
% 
%     [T,C] = ode15s(parameterized_ODEs,S.range,S.init_cond,ODE_options);
% 
%     rel_rate_D(i) = Calc_Function(T,C,S);
% 
% end
% 
% figure()
% plot(log10(y),rel_rate_D)