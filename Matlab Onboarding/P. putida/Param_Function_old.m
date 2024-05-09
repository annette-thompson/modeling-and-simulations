function P = Param_Function_old(S)
%% Parameters
% Includes FabG kcat scaling option

num_elong_steps = S.num_elong_steps; %number of elongation steps

P = struct;

P.labels = S.labels;

% if the parameter is kcat, then scale each estimated kcat value from the 
% table using the appropriate kcat scaling terms

P.kcat7 = S.param_table{'kcat7','parameter_values'}*S.scaling_factor_kcat_term; %termination scaling for TesA

% elongation scaling for FabG,FabI,FabF,FabA,FabB
P.kcat4 = S.param_table{'kcat4','parameter_values'}*S.scaling_factor_kcat;
P.kcat6 = S.param_table{'kcat6','parameter_values'}*S.scaling_factor_kcat;
P.kcat8 = S.param_table{'kcat8','parameter_values'}*S.scaling_factor_kcat;
P.kcat9 = S.param_table{'kcat9','parameter_values'}*S.scaling_factor_kcat;
P.kcat10 = S.param_table{'kcat10','parameter_values'}*S.scaling_factor_kcat;

P.kcat5 = S.param_table{'kcat5','parameter_values'}*S.scaling_factor_fabAZ_kcat;%FabZ (c4) scaling

% initiation scaling for ACC, FabH
P.kcat1_1 = S.param_table{'kcat1_1','parameter_values'}*S.scaling_factor_kcat_init;
P.kcat1_2 = S.param_table{'kcat1_2','parameter_values'}*S.scaling_factor_kcat_init;
P.kcat3 = S.param_table{'kcat3','parameter_values'}*S.scaling_factor_kcat_init;
 
% if the parameter is kon or koff, then scale both using the appropriate 
% scaling factors. Kon is further modified using the estimated Kd values 
% such that kd_est = koff_final/kon_final
 
% For FabD,FabH use "a1" scaling
P.k2_1r = S.param_table{'k2_1r','parameter_values'}*S.scaling_factor_init;
P.k2_3r = S.param_table{'k2_3r','parameter_values'}*S.scaling_factor_init;
P.k3_1r = S.param_table{'k3_1r','parameter_values'}*S.scaling_factor_init;
P.k3_3r = S.param_table{'k3_3r','parameter_values'}*S.scaling_factor_init;
P.k2_1f = P.k2_1r/S.km_table{'k2_1f','parameter_values'};
P.k2_3f = P.k2_3r/S.km_table{'k2_3f','parameter_values'};
P.k3_1f = P.k3_1r/S.km_table{'k3_1f','parameter_values'};
P.k3_3f = P.k3_3r/S.km_table{'k3_3f','parameter_values'};

% For TesA use parameter "a3"
P.k7_1f = S.param_table{'k7_1f','parameter_values'}*S.scaling_factor_term;
P.k7_1r = S.param_table{'k7_1r','parameter_values'}*S.scaling_factor_term;

% For FabF and FabB use paramter "a2" (seperated as scaling_factor_fabF in 
% case different scalings are desired, note scaling_factor_fabF=a2 here)
P.k8_1r = S.param_table{'k8_1r','parameter_values'}*S.scaling_factor_fabF;
P.k10_1r = S.param_table{'k10_1r','parameter_values'}*S.scaling_factor_fabF;
P.k8_1f = P.k8_1r/S.km_table{'k8_1f','parameter_values'};
P.k10_1f = P.k10_1r/S.km_table{'k10_1f','parameter_values'};

% This paramter modification step differs from other steps in thatit is not 
% a binding step, but a reverse catalytic step. For FabZ and FabA reverse 
% use parameter "c4" and "c2." Note that FabZ and FabA are reversible 
% reactions, the forward reaction is denoted by kcat5 and kcat9, the 
% reverse reaction by k5_2r and k9_2r. As both forward and reverse are 
% scaled by the same constant the ratio (Keq) is maintained.
P.k5_2r = S.param_table{'k5_2r','parameter_values'}*S.scaling_factor_fabAZ_kcat;
P.k9_2r = S.param_table{'k9_2r','parameter_values'}*S.scaling_factor_kcat;

 
% For FabD,FabH,FabF,FabB forward and reverse intermediate reaction steps 
% (the intermediate reaction of the ping-pong mechanism) use scaling "b1"
P.k2_2r = S.param_table{'k2_2r','parameter_values'}*S.kd_fits(3);
P.k2_4r = S.param_table{'k2_4r','parameter_values'}*S.kd_fits(3);
P.k3_2r = S.param_table{'k3_2r','parameter_values'}*S.kd_fits(3);
P.k8_2r = S.param_table{'k8_2r','parameter_values'}*S.kd_fits(3);
P.k10_2r = S.param_table{'k10_2r','parameter_values'}*S.kd_fits(3);

% For FabZ use parameter "a2"
P.k5_1r = S.param_table{'k5_1r','parameter_values'}*S.scaling_factor_elon;
P.k5_1f = P.k5_1r/S.km_table{'k5_1f','parameter_values'};

% For FabG,FabI,FabF,FabA,FabB use parameter "a2"
P.k4_1r = S.param_table{'k4_1r','parameter_values'}*S.scaling_factor_elon;
P.k4_2r = S.param_table{'k4_2r','parameter_values'}*S.scaling_factor_elon;
P.k6_1r = S.param_table{'k6_1r','parameter_values'}*S.scaling_factor_elon;
P.k6_2r = S.param_table{'k6_2r','parameter_values'}*S.scaling_factor_elon;
P.k8_3r = S.param_table{'k8_3r','parameter_values'}*S.scaling_factor_elon;
P.k9_1r = S.param_table{'k9_1r','parameter_values'}*S.scaling_factor_elon;
P.k10_3r = S.param_table{'k10_3r','parameter_values'}*S.scaling_factor_elon;

P.k4_1f = P.k4_1r/S.km_table{'k4_1f','parameter_values'};
P.k4_2f = P.k4_2r/S.km_table{'k4_2f','parameter_values'};
P.k6_1f = P.k6_1r/S.km_table{'k6_1f','parameter_values'};
P.k6_2f = P.k6_2r/S.km_table{'k6_2f','parameter_values'};
P.k8_3f = P.k8_3r/S.km_table{'k8_3f','parameter_values'};
P.k9_1f = P.k9_1r/S.km_table{'k9_1f','parameter_values'};
P.k10_3f = P.k10_3r/S.km_table{'k10_3f','parameter_values'};

% ACC and FabH
P.k1_1f = S.param_table{'k1_1f','parameter_values'}*S.scaling_factor_elon;
P.k1_1r = S.param_table{'k1_1r','parameter_values'}*S.scaling_factor_elon;
P.k1_2f = S.param_table{'k1_2f','parameter_values'}*S.scaling_factor_elon;
P.k1_2r = S.param_table{'k1_2r','parameter_values'}*S.scaling_factor_elon;
P.k1_3f = S.param_table{'k1_3f','parameter_values'}*S.scaling_factor_elon;
P.k1_3r = S.param_table{'k1_3r','parameter_values'}*S.scaling_factor_elon;

% After the initial parameter assignment, additional parameters are
% assigned, and modified (for example incorporating substrate specificity)

% For FabH inhibition use appropriate on/off rates (note that inhibition has
% two types, noncompetitive (k3_4*) with respect to acetyl-CoA
% and competitive (k3_5*)with respect to malonyl-ACP

% Noncompetitive inhibition with respect to acetyl-CoA (binding FabH)
P.k3_4f = zeros(1,num_elong_steps);
P.k3_4r = zeros(1,num_elong_steps);

for i = 1:num_elong_steps
    if i <= 5
        P.k3_4f(i) = S.inhibition_on_rates(1); %inhibition for acyl-ACPs of 4-12 has the same value
        P.k3_4r(i) = S.inhibition_kds(1,1)*S.inhibition_on_rates(1);%on rate calculation from Kds are same value for all chain lengths)
    else
        P.k3_4f(i) = S.inhibition_on_rates(1);
        P.k3_4r(i) = S.inhibition_kds(i-4,1)*S.inhibition_on_rates(1);
    end
end
 
% Competitive inhibition with respect to malonyl-ACP (binding FabH*)
P.k3_5f = zeros(1,num_elong_steps);
P.k3_5r = zeros(1,num_elong_steps);

for i = 1:num_elong_steps
    if i <= 5
        P.k3_5f(i) = S.inhibition_on_rates(2);%inhibition for acyl-ACPs of 4-12 has the same value
        P.k3_5r(i) = S.inhibition_kds(1,2)*S.inhibition_on_rates(2);%inhibition for acyl-ACPs of 4-12 has the same value
    else
        P.k3_5f(i) = S.inhibition_on_rates(2);
        P.k3_5r(i) = S.inhibition_kds(i-4,2)*S.inhibition_on_rates(2);
    end
end

% Specify chain length dependence of kon and kcat for TesA
kcat70 = P.kcat7;
P.k7_1f = zeros(1,num_elong_steps);
P.kcat7 = zeros(1,num_elong_steps); 
 
% Specify source of measurments, 'Pf'; Pfleger group measurements, 'Fox'; Fox group measurements
% Implements TesA kcat chain length dependence as relative scaling (of base 
% value which is fit)
if strcmp(S.TesA_fitting_source,'Pf')
    kd_12 = exp(S.lin_slope*(12) + S.lin_int);%kd estimated from linear free energy relationship for chain lengths 12-20
    ratio_val = kd_12/S.Pf_scaling;%ratio used to match kd at chain length 12
    kd_est = (ratio_val).*S.Pf_kd_est_scaling;%Kds for 4,6 and 8-12, estimated Kd is scaled to match linear free energy values
    kcat_scaling = S.Pf_kcat_scaling;
elseif strcmp(S.TesA_fitting_source,'Fox')
    kd_12 = exp(S.lin_slope*(12) + S.lin_int);
    ratio_val = kd_12/S.Fox_scaling;
    kd_est = (ratio_val).*S.Fox_kd_est_scaling;
    kcat_scaling = S.Fox_kcat_scaling;
elseif strcmp(S.TesA_fitting_source,'Non-native')
    kd_est = S.Non_kd_est;
    kcat_scaling = S.Non_kcat_scaling;
elseif strcmp(S.TesA_fitting_source,'R3M1')
    kd_est = S.R3M1_kd_est;
    kcat_scaling = S.R3M1_kcat_scaling;
elseif strcmp(S.TesA_fitting_source,'R3M4')
    kd_est = S.R3M4_kd_est;
    kcat_scaling = S.R3M4_kcat_scaling;
end

for i = 1:num_elong_steps
    kd_long = exp(S.lin_slope*(i*2+2) + S.lin_int);
    P.kcat7(i) = kcat70*kcat_scaling(i);
    if i < 5
        P.k7_1f(i) = P.k7_1r/kd_est(i);%for 4-12 use linear free energy values
    elseif strcmp(S.TesA_fitting_source,'Non-native') || strcmp(S.TesA_fitting_source,'R3M1') || strcmp(S.TesA_fitting_source,'R3M4')
        P.k7_1f(i) = P.k7_1r/kd_est(i);
    else
        P.k7_1f(i) = P.k7_1r/kd_long;%for 12-20 use scaled estimates
    end
end

% FabF and FabB on rates can be further modified by restricting elongation
% here by changing the value of i (i = 5 is chain length 12) and modifying
% the if/else statement
k8_1f0 = P.k8_1f;
k10_1f0 = P.k10_1f;
P.k8_1f = zeros(1,num_elong_steps);
P.k10_1f = zeros(1,num_elong_steps);
for i = 1:num_elong_steps
    if i > 5
        P.k8_1f(i) = k8_1f0;
        P.k10_1f(i) = k10_1f0;
    else
        P.k8_1f(i) = k8_1f0;%change if using restricted elongation mutant
        P.k10_1f(i) = k10_1f0;%change if using restricted elongation mutant
    end
end
 
% Acyl transfer step k_fwd and k_rvs are determined by the fit parameters 
% b2 and b3, which are the Keq values for the transfer step
% b2 is Keq for FabD (first step) only
% b3 is Keq for FabD (second step), FabH, FabF, and FabB

% FabD first transfer step
P.k2_2f = P.k2_2r/S.kd_fits(1);%b2

% FabH transfer step
P.k3_2f = P.k3_2r/S.kd_fits(2);%b3

% FabF transfer step and FabB transfer step
P.k8_2f = zeros(1,num_elong_steps);
P.k10_2f = zeros(1,num_elong_steps);
for i = 1:num_elong_steps
    P.k8_2f(i) = P.k8_2r/S.kd_fits(2);%b3
    P.k10_2f(i) = P.k10_2r/S.kd_fits(2);%b3
end
 
% FabD second transfer step
P.k2_4f = P.k2_4r/S.kd_fits(2);%b3
 
% FabB and FabF parameters for FabH-like activity
P.k10_4f = S.param_table{'k10_4f','parameter_values'}*S.scaling_factor_fabB_init;
P.k10_4r = S.param_table{'k10_4r','parameter_values'};
P.k8_4f = S.param_table{'k8_4f','parameter_values'}*S.scaling_factor_fabF_init;
P.k8_4r = S.param_table{'k8_4r','parameter_values'};
P.k10_9f = S.param_table{'k10_9f','parameter_values'};
P.k10_9r = S.param_table{'k10_9r','parameter_values'};
P.k8_9f = S.param_table{'k8_9f','parameter_values'};
P.k8_9r = S.param_table{'k8_9r','parameter_values'};
P.k10_5f = S.param_table{'k10_5f','parameter_values'};
P.k10_5r = S.param_table{'k10_5r','parameter_values'}*S.scaling_factor_aCoA_10;
P.k8_5f = S.param_table{'k8_5f','parameter_values'};
P.k8_5r = S.param_table{'k8_5r','parameter_values'};
P.k10_6f = S.param_table{'k10_6f','parameter_values'};
P.k10_6r = S.param_table{'k10_6r','parameter_values'};
P.k8_6f = S.param_table{'k8_6f','parameter_values'};
P.k8_6r = S.param_table{'k8_6r','parameter_values'};
P.k10_7f = S.param_table{'k10_7f','parameter_values'};
P.k10_7r = S.param_table{'k10_7r','parameter_values'};
P.k8_7f = S.param_table{'k8_7f','parameter_values'};
P.k8_7r = S.param_table{'k8_7r','parameter_values'};
P.kcat10_H = S.param_table{'kcat10_H','parameter_values'};
P.kcat8_H = S.param_table{'kcat8_H','parameter_values'};
P.kcat10_CO2 = S.param_table{'kcat10_CO2','parameter_values'}*S.scaling_factor_kcat10_CO2;
P.kcat8_CO2 = S.param_table{'kcat8_CO2','parameter_values'}*S.scaling_factor_kcat8_CO2;
P.k10_8f = S.param_table{'k10_8f','parameter_values'};
P.k10_8r = S.param_table{'k10_8r','parameter_values'};
P.k8_8f = S.param_table{'k8_8f','parameter_values'};
P.k8_8r = S.param_table{'k8_8r','parameter_values'}*S.scaling_factor_aCoA_8;

% Chain length specificities for FabH,FabG,FabZ,FabI,FabF,FabA,FabB
kcat30 = P.kcat3;
P.kcat3 = zeros(1,num_elong_steps);
kcat40 = P.kcat4;
P.kcat4 = zeros(1,num_elong_steps);
kcat50 = P.kcat5;
P.kcat5 = zeros(1,num_elong_steps);
kcat60 = P.kcat6;
P.kcat6 = zeros(1,num_elong_steps);
kcat80 = P.kcat8;
P.kcat8 = zeros(1,num_elong_steps);
kcat90 = P.kcat9;
P.kcat9 = zeros(1,num_elong_steps);
kcat100 = P.kcat10;
P.kcat10 = zeros(1,num_elong_steps);
for i = 1:num_elong_steps
    P.kcat3(i) = kcat30*S.kcat_scaling_fabH(i);% FabH chain length kcat scaling
    P.kcat4(i) = kcat40*S.kcat_scaling_fabG(i);% FabG chain length kcat scaling
    P.kcat5(i) = kcat50*S.kcat_scaling_fabZ(i);% FabZ chain length kcat scaling
    P.kcat6(i) = kcat60*S.kcat_scaling_fabI(i);% FabI chain length kcat scaling
    P.kcat8(i) = kcat80*S.kcat_scaling_fabF(i);% FabF chain length kcat scaling
    P.kcat9(i) = kcat90*S.kcat_scaling_fabA(i);% FabA chain length kcat scaling
    P.kcat10(i) = kcat100*S.kcat_scaling_fabB(i);% FabB chain length kcat scaling
end

%For FabZ and FabA chain length specificities of k_rvs (reverse reaction
%rate 'k5_2r' and 'k9_2r') and kon
k5_2r0 = P.k5_2r;
P.k5_2r = zeros(1,num_elong_steps);
k9_2r0 = P.k9_2r;
P.k9_2r = zeros(1,num_elong_steps);
k5_1f0 = P.k5_1f;
P.k5_1f = zeros(1,num_elong_steps);
k9_1f0 = P.k9_1f;
P.k9_1f = zeros(1,num_elong_steps);
for i = 1:num_elong_steps
    P.k5_2r(i) = k5_2r0*S.kcat_scaling_fabZ(i);%FabZ chain length k_rvs scaling
    P.k9_2r(i) = k9_2r0*S.kcat_scaling_fabA(i);%FabA chain length k_rvs scaling
    P.k5_1f(i) = k5_1f0*S.kon_scaling_fabZ(i);%FabZ chain length kon scaling
    P.k9_1f(i) = k9_1f0*S.kon_scaling_fabA(i);%FabA chain length kon scaling
end

% Remaining parameters that need to be vectors for elongation
P.k3_1f = P.k3_1f.*ones(1,num_elong_steps);
P.k3_1r = P.k3_1r.*ones(1,num_elong_steps);
P.k3_2f = P.k3_2f.*ones(1,num_elong_steps);
P.k3_2r = P.k3_2r.*ones(1,num_elong_steps);
P.k3_3f = P.k3_3f.*ones(1,num_elong_steps);
P.k3_3r = P.k3_3r.*ones(1,num_elong_steps);
P.k4_1f = P.k4_1f.*ones(1,num_elong_steps);
P.k4_1r = P.k4_1r.*ones(1,num_elong_steps);
P.k4_2f = P.k4_2f.*ones(1,num_elong_steps);
P.k4_2r = P.k4_2r.*ones(1,num_elong_steps);
P.kcat4 = P.kcat4.*ones(1,num_elong_steps);
P.k5_1r = P.k5_1r.*ones(1,num_elong_steps);
P.k6_1f = P.k6_1f.*ones(1,num_elong_steps);
P.k6_1r = P.k6_1r.*ones(1,num_elong_steps);
P.k6_2f = P.k6_2f.*ones(1,num_elong_steps);
P.k6_2r = P.k6_2r.*ones(1,num_elong_steps);
P.k7_1r = P.k7_1r.*ones(1,num_elong_steps);
P.k8_1r = P.k8_1r.*ones(1,num_elong_steps);
P.k8_2r = P.k8_2r.*ones(1,num_elong_steps);
P.k8_3f = P.k8_3f.*ones(1,num_elong_steps);
P.k8_3r = P.k8_3r.*ones(1,num_elong_steps);
P.k9_1r = P.k9_1r.*ones(1,num_elong_steps);
P.k10_1r = P.k10_1r.*ones(1,num_elong_steps);
P.k10_2r = P.k10_2r.*ones(1,num_elong_steps);
P.k10_3f = P.k10_3f.*ones(1,num_elong_steps);
P.k10_3r = P.k10_3r.*ones(1,num_elong_steps);

% Other params

% ACC (not used)
P.e1tot = S.enzyme_conc(1);

% FabD
P.e2tot = S.enzyme_conc(2);

% FabH
P.k3_inh_f = S.ACP_inh(1);
P.k3_inh_r = S.ACP_inh(2);
P.e3tot = S.enzyme_conc(3);

% FabG
P.k4_inh_f = S.ACP_inh(3);
P.k4_inh_r = S.ACP_inh(4);
P.e4tot = S.enzyme_conc(4);

% FabZ
P.k5_3f = P.k5_1r;
P.k5_3r = P.k5_1f;
P.k5_inh_f = S.ACP_inh(5);
P.k5_inh_r = S.ACP_inh(6);
P.e5tot = S.enzyme_conc(5);

% FabI
P.k6_inh_f = S.ACP_inh(7);
P.k6_inh_r = S.ACP_inh(8);
P.e6tot = S.enzyme_conc(6);

% TesA
P.k7_inh_f = S.ACP_inh(9);
P.k7_inh_r = S.ACP_inh(10);
P.e7tot = S.enzyme_conc(7);

% FabF
P.kcat8_un = P.kcat8(4).*S.kcat_scaling_fabF_unsat;%specificity of reaction with unsaturated acyl chains
P.k8_inh_f = S.ACP_inh(11);
P.k8_inh_r = S.ACP_inh(12);
P.e8tot = S.enzyme_conc(8);

% FabA
P.k9_3f = P.k9_1r;
P.k9_3r = P.k9_1f;
P.k9_1f_un = P.k9_1f.*S.scaling_vector_fabA_unsat;%specificity of reaction with unsaturated acyl chains
P.k9_1r_un = P.k9_1r;
P.kcat9_un = P.kcat9.*S.kcat_scaling_fabA_unsat;
P.k9_2r_un = P.k9_2r;
P.k9_3f_un = P.k9_3f;
P.k9_3r_un = P.k9_3r;
P.k9_inh_f = S.ACP_inh(13);
P.k9_inh_r = S.ACP_inh(14);
P.e9tot = S.enzyme_conc(9);

% FabB
P.kcat10_un = P.kcat10(5).*S.kcat_scaling_fabB_unsat;%specificity of reaction with unsaturated acyl chains
P.k10_inh_f = S.ACP_inh(15);
P.k10_inh_r = S.ACP_inh(16);
P.e10tot = S.enzyme_conc(10);
