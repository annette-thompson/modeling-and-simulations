%% Main_Script
% Reorganized model: user only edits options and variables
% clear all
% close all
% clc
% tic


%% Options

% Running the function
run = 'yes';

% Options for the ODE solver
% RelTol: relative error tolerance of solution (minor impact on solve time)
% MaxOrder: Max order of numerical differention equation (default, minor impact on solve time)
% JPattern: Matrix defining sparsity pattern of the Jacobian (1 for all nonzero 
% elements in the Jacobian, 0 otherwise). Included to reduce solve time.
% Vectorized: Specifies that the model is formatted in a vector format
% (signfiicantly impacts solve time)
load('JpMat.mat','JpMatPrime')
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Optimizing p_vec
optimize = 'no';

% What p_vec vars do you want to optimize/add?
p_opt = [1, 2, 3];

% p_vec guess
% p_vec = [a1 a2 a3 b1 c2 c3 c1 b2 b3 d1 d2 e f c4 x1 x2 x3 x4];
p_vec0 = [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
    4.645634978 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673 132.8499358...
    2180.050007 0.539756276	0.053673263	34.49718991	11.15058888];


fun_evals_max = 1000; %max number of function evaluations
max_iter = 180; %max number of optimiztion iterations

opt_options = optimset('MaxFunEvals',fun_evals_max,'Display','iter','MaxIter',max_iter);

% Sensitivity Analysis
sensitivity = 'no';

% Error Analysis
errors = 'no';

%% Variables

% Set up structure to store variables
S = struct;

% Time range (sec)
S.range = [0 720];

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
    's17';'s18'
};
S.num = length(S.labels); %how many diff eqs there are

% Set initial conditions
S.init_cond = zeros(S.num,1);
S.init_cond(1) = 0; %s1 (ATP)
S.init_cond(10) = 0; %p4 (Malonyl-ACP)
S.init_cond(11) = 0; %p5 (CO2)
S.init_cond(2) = 0; %s2 (Bicarbonate)
S.init_cond(3) = 500; %s3 (Acetyl-CoA)
S.init_cond(4) = 10; %s6 (holo ACP)
S.init_cond(5) = 1000; %s7 (NADPH)
S.init_cond(6) = 1000; %s8 (NADH)
S.init_cond(7) = 0; %p1 (ADP)
S.init_cond(8) = 500; %p2 (malonyl-CoA)
S.init_cond(9) = 0; %p3 (CoA)
S.init_cond(318) = 0; %s11 (Octanoyl-CoA)

% Set initial enzyme concentrations
% (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
S.enzyme_conc= [0 1 1 1 1 1 10 1 1 1]; 

% Set thioesterase fitting source
% 'Pf'; Pfleger group measurements
% 'Fox'; Fox group measurements
% 'Non-native'; For alternative thioesterases
% 'R3M1'; Pfleger group measurements for TesA R3M1
% 'R3M4'; Pfleger group measurements for TesA R3M4
S.TesA_fitting_source = 'Pf';

% Load parameter estimates
S.km_table = readtable('km_est.csv','ReadRowNames',true);
S.param_table = readtable('est_param.csv','ReadRowNames',true);

% The rest of the variables will most likely not change

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

% Comparison data with various combinations of FabH, FabF, FabB, and AcCoA 
% removed for optimization
S.opt_init_rate_data = [4.655129552 2.41772584 1.796165966 2.976304272 2.91667542 0.893110994 2.428676471];
% Excel Spreadsheet: 06172021 NADPH Assay No FabH, F, B and AcCoA Testing
S.opt_tot_prod_file = 'Experimental_Dataset.csv'; %Must have column titles Time(min) and Concentration(um) 
% (Yu, X.; Liu, T.; Zhu, F.; Khosla, C. Proceedings of the National Academy
% of Sciences 2011, 108, 18643-18648.)
S.opt_prod_dist_data = [0,0,0.059771046,0.011042448,0.202613717,0.075980144,0.29378989,0.07192787,0.143855739,0.117515956,0.002228751,0.02127444,0,0];
% (Grisewood, M.; Hernández-Lozada, N.; Thoden, J.; Gifford, N.; Mendez-Perez, 
% D.;Lai, R.; Holden, H.; Pfleger, Et. al. ACS Catalysis 2017, 7, 3837-3849.)


%% Optimization
if strcmp(optimize,'yes')
    S.p_vec = p_vec0;
    S.p_vec = Optimization_Function(S,p_opt,opt_options,ODE_options);
end

%% ODES
if strcmp(run,'yes')

    P = Param_Function_new(S);
    
    parameterized_ODEs = @(t,c) ODE_Function_new(t,c,P);
    % Equivalent to param and opt vec

 
    % Intial Conditions
    init_cond = S.init_cond;

    %Time Range
    range = S.range;

    %Numerical solution is found with solver ode15s (ideal for stiff systems)
    %T: Time (sec)
    %C: Matrix of concentration values of each species (columns) at each time point in T (rows)
    [T,C] = ode15s(parameterized_ODEs,range,init_cond,ODE_options);

    check = mass_balance(C,P);
    disp(check)
    % Calcs = Calc_Plot_Function(T,C,run_options); %then can add sections 
    % (options) as I make them and note what figures they are
end

% %% Sensitivity Analysis
% if strcmp(sensitivity,'yes')
%     Analysis = Sensitivity_Analysis_Function(S,sens_options);
% end
% 
% %% Error Calculation
% if strcmp(errors,'yes')
%     Error = Error_Function(S,err_options);
%     % bayesian or other types
% end

% toc
