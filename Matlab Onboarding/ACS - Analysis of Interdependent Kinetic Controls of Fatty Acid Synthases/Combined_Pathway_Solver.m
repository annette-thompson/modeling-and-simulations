%Combined FAS Pathway Solver
%Initializes the parameters, initial conditions and solver options and
%executes the ode solver ode15s to numerically solve the associated system 
%of equations in Combined_Pathway_Model

%Fit parameters values from optimization
%p_vec = [a1 a2 a3 b1 c2 c3 c1 b2 b3 d1 d2 e]
p_vec = [4498.199121, 15.85799497, 4319.091704, 24.74149287, 87.93965874,...
    0.019016657, 1.650547537, 0.052767764, 0.085412774, -0.289676757,...
    5.443430645, 28.7223893];
acp_bind = p_vec(12);%parameter e

%Initial conditions of species concentrations in init_cond (units of uM)
init_cond = zeros(147,1);
init_cond(1) = 500;%s3 (Acetyl-CoA)
init_cond(2) = 10;%s6 (holo ACP)
init_cond(3) = 1000;%s7 (NADPH)
init_cond(4) = 1000;%s8 (NADH)
init_cond(5) = 500;%p2 (Malonyl-CoA)
init_cond(6) = 0;%p3 (CoA)
init_cond(7) = 0;%p4 (Malonyl-ACP)
init_cond(8) = 0;%p5 (CO2)

%Species labels
labels = {'s3';'s6';'s7';'s8';'p2';'p3';'p4';'p5';'Q_4';'M_4';'R_4';'T_4';...
'F_4';'Q_6';'M_6';'R_6';'T_6';'F_6';'Q_8';'M_8';'R_8';'T_8';'F_8';'Q_10';...
'M_10';'R_10';'T_10';'F_10';'Q_12';'M_12';'R_12';'T_12';'F_12';'Q_14';...
'M_14';'R_14';'T_14';'F_14';'Q_16';'M_16';'R_16';'T_16';'F_16';'Q_18';...
'M_18';'R_18';'T_18';'F_18';'Q_20';'M_20';'R_20';'T_20';'F_20';'e2p2';...
'e2act';'e2acts6';'e3s3';'e3act';'e3actp4';'e4s7';'e6s8';'e3T_4';...
'e3actT_4';'e4s7Q_4';'e5M_4';'e6s8R_4';'e7T_4';'e8T_4';'e8act4';'e8act4p4';...
'e3T_6';'e3actT_6';'e4s7Q_6';'e5M_6';'e6s8R_6';'e7T_6';'e8T_6';'e8act6';...
'e8act6p4';'e3T_8';'e3actT_8';'e4s7Q_8';'e5M_8';'e6s8R_8';'e7T_8';'e8T_8';...
'e8act8';'e8act8p4';'e3T_10';'e3actT_10';'e4s7Q_10';'e5M_10';'e6s8R_10';...
'e7T_10';'e8T_10';'e8act10';'e8act10p4';'e3T_12';'e3actT_12';'e4s7Q_12';...
'e5M_12';'e6s8R_12';'e7T_12';'e8T_12';'e8act12';'e8act12p4';'e3T_14';...
'e3actT_14';'e4s7Q_14';'e5M_14';'e6s8R_14';'e7T_14';'e8T_14';'e8act14';...
'e8act14p4';'e3T_16';'e3actT_16';'e4s7Q_16';'e5M_16';'e6s8R_16';'e7T_16';...
'e8T_16';'e8act16';'e8act16p4';'e3T_18';'e3actT_18';'e4s7Q_18';'e5M_18';...
'e6s8R_18';'e7T_18';'e8T_18';'e8act18';'e8act18p4';'e3T_20';'e3actT_20';...
'e4s7Q_20';'e5M_20';'e6s8R_20';'e7T_20';'e3s6';'e4s6';'e5s6';'e6s6';'e7s6';...
'e8s6';'prodQ_4';'prodQ_6_20'};

%List of kinetic parameter names, order does not matter, except that
%all initiation steps must be listed before all elongation steps (k4 and
%onwards are elongation steps), does not include inhibition parameters
param_names = {{'k2_1f';'k2_1r';'k2_2f';'k2_2r';'k2_3f';'k2_3r';'k2_4f';...
'k2_4r';'k3_1f';'k3_1r';'k3_2f';'k3_2r';'k3_3f';'k3_3r';'k3_4f';'k3_4r';...
'k3_5f';'k3_5r';'kcat3';'k4_1f';'k4_1f';'k4_1r';'k4_2f';'k4_2r';'kcat4';...
'k5_1f';'k5_1r';'kcat5';'k6_1f';'k6_1r';'k6_2f';'k6_2r';'kcat6';'k7_1f';...
'k7_1r';'kcat7';'k8_1f';'k8_1r';'k8_2f';'k8_2r';'k8_3f';'k8_3r';'kcat8';}};


%Sets up parameter options for use in param_func
field1 = 'elong_num'; val1 = 9;%number of elongation steps
field4 = 'dist_opt'; val4 = 'data_defined';%data defined distributions of parameters
field5 = 'opt_name'; val5 = {{'k2_4f','k3_4f','k3_4r','k3_5f','k3_5r',...
    'k7_1f','kcat7','k8_1f','k8_2f'}};%parameters to be modified according to unique rules
field6 = 'param_names'; val6 = param_names;
field10 = 'scaling_factor_init'; val10 = p_vec(1);%parameter a1
field11 = 'scaling_factor_elon'; val11 = p_vec(2);%parameter a2
field12 = 'scaling_factor_kcat'; val12 = p_vec(5);%parameter c2 (elongation kcats)
field13 = 'scaling_factor_term'; val13 = p_vec(3);%parameter a3
field14 = 'scaling_factor_fabf'; val14 = p_vec(2);%same as parameter a2
field15 = 'scaling_factor_kcat_term'; val15 = p_vec(6);%parameter c3

%ACP inhibition on/off rates [FabH on, FabH off, FabG on, FabG off, FabZ
%on, FabZ off, FabI on, FabI off, TesA on, TesA off, FabF on, FabF off]
field16 = 'ACP_inh'; val16 = [acp_bind*2.41E-05,2.17E-02,2.41E-03,2.17E-02,...
    2.41E-03,2.17E-02,2.41E-03,2.17E-02,2.41E-03,2.17E-02,acp_bind*2.41E-03,...
    2.17E-02];
field17 = 'enzyme_conc'; val17 = [1 1 1 1 1 10 1];%enzyme concentrations in uM [FabD FabH FabG FabZ FabI TesA FabF]
field18 = 'inhibition_kds';val18 = (1/acp_bind).*[4335.05,23.67;824.7,30.1;...
    967.5,7.55;251.52,8.484;128.78,2.509];%(chain length 12-20)
field19 = 'inhibition_on_rates';val19 = [0.3088,1.552,0.3088];%inhibition on rates
field20 = 'kd_fits';val20 = [p_vec(8),56,p_vec(9),p_vec(4),p_vec(9)];%parameters [b2,[],b3,b1,b3]
field21 = 'lin_param';val21=[p_vec(10) p_vec(11)];%parameters d1 and d2
field22 = 'scaling_factor_kcat_init';val22 = p_vec(7);%parameter c1

%Parameter options structure for passing options to param_func
opt_struct = struct(field1,val1,field4,val4,field5,val5,field6,val6,...
    field10,val10,field11,val11,field12,val12,field13,val13,field14,val14,...
    field15,val15,field16,val16,field17,val17,field18,val18,field19,val19,...
    field20,val20,field21,val21,field22,val22);

%Calls param_func to return all parameter values given the fit parameters
param_struct = struct;
for i = 1:length(param_names{1})
    param_struct.(param_names{1}{i}) = param_func(param_names{1}{i},i,opt_struct);
end

%Parameterizes the model equations
parameterized_Combined_Pathway_Model = @(t,c) Combined_Pathway_Model(t,c,...
    param_struct,opt_struct);

%Time range (sec)
range = [0 750];

options = odeset('RelTol',1e-4,'MaxOrder',2);%ode options,'MaxOrder' of 
%numerical differentiaion specified to reduce computational cost without 
%signficant decrease in accuracy

%Built-in implicit numerical integrator ode15s returns time and
%concentration values for the specified time range
[T,C] = ode15s(parameterized_Combined_Pathway_Model,range,init_cond,options);

%Calcuates total amounts of ACP intermediate specicies and fatty acids
F_total = zeros(length(T),1);
Q_total = zeros(length(T),1);
M_total = zeros(length(T),1);
R_total = zeros(length(T),1);
T_total = zeros(length(T),1);

for ind = 1:length(labels)
        label_val = char(labels(ind));
        first_char = label_val(1);
        if first_char == char('F')
            F_total = C(:,ind) + F_total;
        end
        if first_char == char('Q')
            Q_total = C(:,ind) + Q_total;
        end
        if first_char == char('M')
            M_total = C(:,ind) + M_total;
        end
        if first_char == char('R')
            R_total = C(:,ind) + R_total;
        end
        if first_char == char('T')
            T_total = C(:,ind) + T_total;
        end
end

%Calculates total fatty acids and palmitic equivalents
F_weighted = zeros(length(T),1);
F_saved = zeros(1,9);
F_raw = zeros(1,9);
weight_vec = [4,6,8,10,12,14,16,18,20]/16;
count = 1;

for ind = 1:length(labels)
    label_val = char(labels(ind));
    first_char = label_val(1);
    if first_char == char('F')
        F_saved(count) = weight_vec(count)*(C(end,ind));
        F_raw(count) = C(end,ind);
        F_weighted = weight_vec(count)*(C(:,ind)) + F_weighted;%Palmitic equivalents
        count = count + 1;
    end
end

total_FA = sum(F_raw);%Total fatty acid
frac_FA = F_raw./total_FA;%Fractional composition of fatty acids

%Estimated experimental distribution of fatty acids
%Source: (Grisewood, M. J. et al. Computational Redesign of Acyl-ACP 
%Thioesterase with Improved Selectivity toward Medium-Chain-Length Fatty Acids. 
%ACS Catal. 2017, 7 (6), 3837–3849.)
norm_pdf_val = [0 0 .059 .0109 .275 .361 .258 .0232 .0129];

%Calculate residuals between gaussian and fatty acid distribution
resid = 100*sum((norm_pdf_val - frac_FA).^2);

%Least squares regresstion against validation data
model_output = F_weighted;
time_val = T/60;
val_file_name = 'Experimental_Dataset.csv';
s = LeastSquaresCalc(time_val,model_output,val_file_name);

%Different objective options

%obj = F_weighted(end); %Palmitic equivalents at last time point
%obj = sum(frac_FA.*weight_vec.*16); %Average fatty acid chain length
%obj = resid*s; %Objective of model used in sensitivity analysis