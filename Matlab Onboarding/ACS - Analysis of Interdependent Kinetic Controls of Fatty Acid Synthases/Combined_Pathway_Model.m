function dc = Combined_Pathway_Model(t,c,param_list,opt_struct)
%Combined_Pathway_Model
%Contains all the differential equations and enzyme balances that define
%the FAS model
%   Input:
%       t: time values (seconds)
%       c: concentration values (uM)
%       param_list: structure containing all kinetic parameters
%   Output:
%       dc: values of all differential equations for given concentrations
%       and the kinetic parameters



%Parameters (assigned based off looking up labels in param_list)

k2_1f = param_list.k2_1f;
k2_1r = param_list.k2_1r;
k2_2f = param_list.k2_2f;
k2_2r = param_list.k2_2r;
k2_3f = param_list.k2_3f;
k2_3r = param_list.k2_3r;
k2_4f = param_list.k2_4f;
k2_4r = param_list.k2_4r;
e2tot = opt_struct.enzyme_conc(1);

k3_1f = param_list.k3_1f;
k3_1r = param_list.k3_1r;
k3_2f = param_list.k3_2f;
k3_2r = param_list.k3_2r;
k3_3f = param_list.k3_3f;
k3_3r = param_list.k3_3r;
k3_4f = param_list.k3_4f;
k3_4r = param_list.k3_4r;
k3_5f = param_list.k3_5f;
k3_5r = param_list.k3_5r;
kcat3 = param_list.kcat3;
k3_inh_f = opt_struct.ACP_inh(1);
k3_inh_r = opt_struct.ACP_inh(2);
e3tot = opt_struct.enzyme_conc(2);

k4_1f = param_list.k4_1f;
k4_1r = param_list.k4_1r;
k4_2f = param_list.k4_2f;
k4_2r = param_list.k4_2r;
kcat4 = param_list.kcat4;
k4_inh_f = opt_struct.ACP_inh(3);
k4_inh_r = opt_struct.ACP_inh(4);
e4tot = opt_struct.enzyme_conc(3);

k5_1f = param_list.k5_1f;
k5_1r = param_list.k5_1r;
kcat5 = param_list.kcat5;
k5_inh_f = opt_struct.ACP_inh(5);
k5_inh_r = opt_struct.ACP_inh(6);
e5tot = opt_struct.enzyme_conc(4);

k6_1f = param_list.k6_1f;
k6_1r = param_list.k6_1r;
k6_2f = param_list.k6_2f;
k6_2r = param_list.k6_2r;
kcat6 = param_list.kcat6;
k6_inh_f = opt_struct.ACP_inh(7);
k6_inh_r = opt_struct.ACP_inh(8);
e6tot = opt_struct.enzyme_conc(5);

k7_1f = param_list.k7_1f;
k7_1r = param_list.k7_1r;
kcat7 = param_list.kcat7;
k7_inh_f = opt_struct.ACP_inh(9);
k7_inh_r = opt_struct.ACP_inh(10);
e7tot = opt_struct.enzyme_conc(6);

k8_1f = param_list.k8_1f;
k8_1r = param_list.k8_1r;
k8_2f = param_list.k8_2f;
k8_2r = param_list.k8_2r;
k8_3f = param_list.k8_3f;
k8_3r = param_list.k8_3r;
kcat8 = param_list.kcat8;
k8_inh_f = opt_struct.ACP_inh(11);
k8_inh_r = opt_struct.ACP_inh(12);
e8tot = opt_struct.enzyme_conc(7);







%List of labels used in differential equations
short_labels = {'s3';'s6';'s7';'s8';'p2';'p3';'p4';'p5';...
'Q_4';'M_4';'R_4';'T_4';'F_4';'Q_6';'M_6';'R_6';'T_6';'F_6';...
'Q_8';'M_8';'R_8';'T_8';'F_8';'Q_10';'M_10';'R_10';'T_10';'F_10';...
'Q_12';'M_12';'R_12';'T_12';'F_12';'Q_14';'M_14';'R_14';'T_14';'F_14';...
'Q_16';'M_16';'R_16';'T_16';'F_16';'Q_18';'M_18';'R_18';'T_18';'F_18';...
'Q_20';'M_20';'R_20';'T_20';'F_20';...
'e2p2';'e2act';'e2acts6';'e3s3';'e3act';'e3actp4';'e4s7';'e6s8';...
'e3T_4';'e3actT_4';'e4s7Q_4';'e5M_4';'e6s8R_4';'e7T_4';'e8T_4';'e8act4';'e8act4p4';...
'e3T_6';'e3actT_6';'e4s7Q_6';'e5M_6';'e6s8R_6';'e7T_6';'e8T_6';'e8act6';'e8act6p4';...
'e3T_8';'e3actT_8';'e4s7Q_8';'e5M_8';'e6s8R_8';'e7T_8';'e8T_8';'e8act8';'e8act8p4';...
'e3T_10';'e3actT_10';'e4s7Q_10';'e5M_10';'e6s8R_10';'e7T_10';'e8T_10';'e8act10';'e8act10p4';...
'e3T_12';'e3actT_12';'e4s7Q_12';'e5M_12';'e6s8R_12';'e7T_12';'e8T_12';'e8act12';'e8act12p4';...
'e3T_14';'e3actT_14';'e4s7Q_14';'e5M_14';'e6s8R_14';'e7T_14';'e8T_14';'e8act14';'e8act14p4';...
'e3T_16';'e3actT_16';'e4s7Q_16';'e5M_16';'e6s8R_16';'e7T_16';'e8T_16';'e8act16';'e8act16p4';...
'e3T_18';'e3actT_18';'e4s7Q_18';'e5M_18';'e6s8R_18';'e7T_18';'e8T_18';'e8act18';'e8act18p4';...
'e3T_20';'e3actT_20';'e4s7Q_20';'e5M_20';'e6s8R_20';'e7T_20';'e3s6';'e4s6';'e5s6';'e6s6';'e7s6';'e8s6';'prodQ_4';'prodQ_6_20'...
};







%Iterates through the list of labels and dynamically assigns c(i)
%value to the variable name in the list (makes it much easier to change
%elements of the pathway)
for i=1:length(short_labels)
    eval(sprintf('%s=%d;\n',short_labels{i},c(i)));
end

dlabels = cell(length(short_labels),1);
for i=1:length(short_labels)
    dlabels{i} = sprintf('%s%s','d',short_labels{i});
end

    

%Enzyme balances
e2 = e2tot - e2p2 - e2act - e2acts6;%FabD
e3 = e3tot - e3s3 - e3act - e3actp4 - e3T_4 - e3T_6 - e3T_8 - e3T_10 - e3T_12 - e3T_14 - e3T_16 - e3T_18 - e3T_20...
    - e3actT_4 - e3actT_6 - e3actT_8 - e3actT_10 - e3actT_12 - e3actT_14 - e3actT_16 - e3actT_18 - e3actT_20 - e3s6;%FabH
e4 = e4tot - e4s7 - e4s7Q_4 - e4s7Q_6 - e4s7Q_8 - e4s7Q_10 - e4s7Q_12 - e4s7Q_14 - e4s7Q_16...
    - e4s7Q_18 - e4s7Q_20 - e4s6;%FabG
e5 = e5tot - e5M_4 - e5M_6 - e5M_8 - e5M_10 - e5M_12 - e5M_14 - e5M_16 - e5M_18 - e5M_20 - e5s6;%FabZ
e6 = e6tot - e6s8 - e6s8R_4 - e6s8R_6 - e6s8R_8 - e6s8R_10 - e6s8R_12 - e6s8R_14 - e6s8R_16...
    - e6s8R_18 - e6s8R_20 - e6s6;%FabI
e7 = e7tot - e7T_4 - e7T_6 - e7T_8 - e7T_10 - e7T_12 - e7T_14 - e7T_16 - e7T_18 - e7T_20 - e7s6;%TesA
e8 = e8tot - e8T_4 - e8act4 - e8act4p4 - e8T_6 - e8act6 - e8act6p4 - e8T_8 - e8act8 - e8act8p4...
     - e8T_10 - e8act10 - e8act10p4 - e8T_12 - e8act12 - e8act12p4 - e8T_14 - e8act14 - e8act14p4...
     - e8T_16 - e8act16 - e8act16p4 - e8T_18 - e8act18 - e8act18p4 - e8s6;%FabF

 
 
 
%Create structure of differential equations
dst = struct;
for i = 1:length(short_labels)
    dst.(dlabels{i}) = 0;
end


%Differential equations of model
dst.ds3 = k3_1r*e3s3 - k3_1f*e3*s3;
dst.ds6 = k2_3r*e2acts6 - k2_3f*e2act*s6...
    + kcat7(1)*e7T_4...
    + kcat7(2)*e7T_6...
    + kcat7(3)*e7T_8...
    + kcat7(4)*e7T_10...
    + kcat7(5)*e7T_12...
    + kcat7(6)*e7T_14...
    + kcat7(7)*e7T_16...
    + kcat7(8)*e7T_18...
    + kcat7(9)*e7T_20...
    + k8_2f(1)*e8T_4 - k8_2r(1)*e8act4*s6...
    + k8_2f(2)*e8T_6 - k8_2r(2)*e8act6*s6...
    + k8_2f(3)*e8T_8 - k8_2r(3)*e8act8*s6...
    + k8_2f(4)*e8T_10- k8_2r(4)*e8act10*s6...
    + k8_2f(5)*e8T_12- k8_2r(5)*e8act12*s6...
    + k8_2f(6)*e8T_14- k8_2r(6)*e8act14*s6...
    + k8_2f(7)*e8T_16- k8_2r(7)*e8act16*s6...
    + k8_2f(8)*e8T_18- k8_2r(8)*e8act18*s6...
    + k3_inh_r*e3s6  - k3_inh_f*e3*s6...
    + k4_inh_r*e4s6  - k4_inh_f*e4*s6...
    + k5_inh_r*e5s6  - k5_inh_f*e5*s6...
    + k6_inh_r*e6s6  - k6_inh_f*e6*s6...
    + k7_inh_r*e7s6  - k7_inh_f*e7*s6...
    + k8_inh_r*e8s6  - k8_inh_f*e8*s6...
    ;
dst.ds7 = k4_1r(1)*e4s7 - k4_1f(1)*e4*s7;
dst.ds8 = k6_1r(1)*e6s8 - k6_1f(1)*e6*s8;

dst.dp2 = k2_1r*e2p2 - k2_1f*e2*p2;
dst.dp3 = k2_2f*e2p2 + k3_2f*e3s3 - k2_2r*e2act*p3 - k3_2r*e3act*p3;
dst.dp4 = k2_4f*e2acts6 + k3_3r*e3actp4 - k2_4r*e2*p4 - k3_3f*e3act*p4...
    + k8_3r(1)*e8act4p4 - k8_3f(1)*e8act4*p4...
    + k8_3r(2)*e8act6p4 - k8_3f(2)*e8act6*p4...
    + k8_3r(3)*e8act8p4 - k8_3f(3)*e8act8*p4...
    + k8_3r(4)*e8act10p4- k8_3f(4)*e8act10*p4...
    + k8_3r(5)*e8act12p4- k8_3f(5)*e8act12*p4...
    + k8_3r(6)*e8act14p4- k8_3f(6)*e8act14*p4...
    + k8_3r(7)*e8act16p4- k8_3f(7)*e8act16*p4...
    + k8_3r(8)*e8act18p4- k8_3f(8)*e8act18*p4...
    ;
dst.dp5 = kcat3*e3actp4 + kcat8(1)*e8act4p4 + kcat8(2)*e8act6p4 + kcat8(3)*e8act8p4 + kcat8(4)*e8act10p4...
    + kcat8(5)*e8act12p4 + kcat8(6)*e8act14p4 + kcat8(7)*e8act16p4 + kcat8(8)*e8act18p4;

dst.dQ_4 = kcat3*e3actp4 + k4_2r(1)*e4s7Q_4 - k4_2f(1)*e4s7*Q_4;
dst.dprodQ_4 = kcat3*e3actp4; %used to track total flux
%into elongation steps if desired
dst.dQ_6 =  kcat8(1)*e8act4p4 + k4_2r(2)*e4s7Q_6 - k4_2f(2)*e4s7*Q_6;
dst.dQ_8 =  kcat8(2)*e8act6p4 + k4_2r(3)*e4s7Q_8 - k4_2f(3)*e4s7*Q_8;
dst.dQ_10 = kcat8(3)*e8act8p4 + k4_2r(4)*e4s7Q_10- k4_2f(4)*e4s7*Q_10;
dst.dQ_12 = kcat8(4)*e8act10p4+ k4_2r(5)*e4s7Q_12- k4_2f(5)*e4s7*Q_12;
dst.dQ_14 = kcat8(5)*e8act12p4+ k4_2r(6)*e4s7Q_14- k4_2f(6)*e4s7*Q_14;
dst.dQ_16 = kcat8(6)*e8act14p4+ k4_2r(7)*e4s7Q_16- k4_2f(7)*e4s7*Q_16;
dst.dQ_18 = kcat8(7)*e8act16p4+ k4_2r(8)*e4s7Q_18- k4_2f(8)*e4s7*Q_18;
dst.dQ_20 = kcat8(8)*e8act18p4+ k4_2r(9)*e4s7Q_20- k4_2f(9)*e4s7*Q_20;
dst.dprodQ_6_20 = kcat8(1)*e8act4p4 + kcat8(2)*e8act6p4 + kcat8(3)*e8act8p4 + kcat8(4)*e8act10p4...
    + kcat8(5)*e8act12p4 + kcat8(6)*e8act14p4 + kcat8(7)*e8act16p4 + kcat8(8)*e8act18p4;%used to track elongation

dst.dM_4 =  kcat4(1)*e4s7Q_4 + k5_1r(1)*e5M_4 - k5_1f(1)*e5*M_4;
dst.dM_6 =  kcat4(2)*e4s7Q_6 + k5_1r(2)*e5M_6 - k5_1f(2)*e5*M_6;
dst.dM_8 =  kcat4(3)*e4s7Q_8 + k5_1r(3)*e5M_8 - k5_1f(3)*e5*M_8;
dst.dM_10 = kcat4(4)*e4s7Q_10+ k5_1r(4)*e5M_10- k5_1f(4)*e5*M_10;
dst.dM_12 = kcat4(5)*e4s7Q_12+ k5_1r(5)*e5M_12- k5_1f(5)*e5*M_12;
dst.dM_14 = kcat4(6)*e4s7Q_14+ k5_1r(6)*e5M_14- k5_1f(6)*e5*M_14;
dst.dM_16 = kcat4(7)*e4s7Q_16+ k5_1r(7)*e5M_16- k5_1f(7)*e5*M_16;
dst.dM_18 = kcat4(8)*e4s7Q_18+ k5_1r(8)*e5M_18- k5_1f(8)*e5*M_18;
dst.dM_20 = kcat4(9)*e4s7Q_20+ k5_1r(9)*e5M_20- k5_1f(9)*e5*M_20;

dst.dR_4 =  kcat5(1)*e5M_4 + k6_2r(1)*e6s8R_4 - k6_2f(1)*e6s8*R_4;
dst.dR_6 =  kcat5(2)*e5M_6 + k6_2r(2)*e6s8R_6 - k6_2f(2)*e6s8*R_6;
dst.dR_8 =  kcat5(3)*e5M_8 + k6_2r(3)*e6s8R_8 - k6_2f(3)*e6s8*R_8;
dst.dR_10 = kcat5(4)*e5M_10+ k6_2r(4)*e6s8R_10- k6_2f(4)*e6s8*R_10;
dst.dR_12 = kcat5(5)*e5M_12+ k6_2r(5)*e6s8R_12- k6_2f(5)*e6s8*R_12;
dst.dR_14 = kcat5(6)*e5M_14+ k6_2r(6)*e6s8R_14- k6_2f(6)*e6s8*R_14;
dst.dR_16 = kcat5(7)*e5M_16+ k6_2r(7)*e6s8R_16- k6_2f(7)*e6s8*R_16;
dst.dR_18 = kcat5(8)*e5M_18+ k6_2r(8)*e6s8R_18- k6_2f(8)*e6s8*R_18;
dst.dR_20 = kcat5(9)*e5M_20+ k6_2r(9)*e6s8R_20- k6_2f(9)*e6s8*R_20;



dst.dT_4 =   kcat6(1)*e6s8R_4 + k7_1r(1)*e7T_4 + k8_1r(1)*e8T_4  - k7_1f(1)*e7*T_4 - k8_1f(1)*e8*T_4 - k3_4f(1)*e3*T_4 + k3_4r(1)*e3T_4 - k3_5f(1)*e3act*T_4 + k3_5r(1)*e3actT_4;
dst.dT_6 =   kcat6(2)*e6s8R_6 + k7_1r(2)*e7T_6 + k8_1r(2)*e8T_6  - k7_1f(2)*e7*T_6 - k8_1f(2)*e8*T_6 - k3_4f(2)*e3*T_6 + k3_4r(2)*e3T_6 - k3_5f(2)*e3act*T_6 + k3_5r(2)*e3actT_6;
dst.dT_8 =   kcat6(3)*e6s8R_8 + k7_1r(3)*e7T_8 + k8_1r(3)*e8T_8  - k7_1f(3)*e7*T_8 - k8_1f(3)*e8*T_8 - k3_4f(3)*e3*T_8 + k3_4r(3)*e3T_8 - k3_5f(3)*e3act*T_8 + k3_5r(3)*e3actT_8;
dst.dT_10 =  kcat6(4)*e6s8R_10+ k7_1r(4)*e7T_10+ k8_1r(4)*e8T_10 - k7_1f(4)*e7*T_10- k8_1f(4)*e8*T_10- k3_4f(4)*e3*T_10+ k3_4r(4)*e3T_10- k3_5f(4)*e3act*T_10+ k3_5r(4)*e3actT_10;
dst.dT_12 =  kcat6(5)*e6s8R_12+ k7_1r(5)*e7T_12+ k8_1r(5)*e8T_12 - k7_1f(5)*e7*T_12- k8_1f(5)*e8*T_12- k3_4f(5)*e3*T_12+ k3_4r(5)*e3T_12- k3_5f(5)*e3act*T_12+ k3_5r(5)*e3actT_12;
dst.dT_14 =  kcat6(6)*e6s8R_14+ k7_1r(6)*e7T_14+ k8_1r(6)*e8T_14 - k7_1f(6)*e7*T_14- k8_1f(6)*e8*T_14- k3_4f(6)*e3*T_14+ k3_4r(6)*e3T_14- k3_5f(6)*e3act*T_14+ k3_5r(6)*e3actT_14;
dst.dT_16 =  kcat6(7)*e6s8R_16+ k7_1r(7)*e7T_16+ k8_1r(7)*e8T_16 - k7_1f(7)*e7*T_16- k8_1f(7)*e8*T_16- k3_4f(7)*e3*T_16+ k3_4r(7)*e3T_16- k3_5f(7)*e3act*T_16+ k3_5r(7)*e3actT_16;
dst.dT_18 =  kcat6(8)*e6s8R_18+ k7_1r(8)*e7T_18+ k8_1r(8)*e8T_18 - k7_1f(8)*e7*T_18- k8_1f(8)*e8*T_18- k3_4f(8)*e3*T_18+ k3_4r(8)*e3T_18- k3_5f(8)*e3act*T_18+ k3_5r(8)*e3actT_18;
dst.dT_20 =  kcat6(9)*e6s8R_20+ k7_1r(9)*e7T_20                  - k7_1f(9)*e7*T_20                  - k3_4f(9)*e3*T_20+ k3_4r(9)*e3T_20- k3_5f(9)*e3act*T_20+ k3_5r(9)*e3actT_20;

dst.dF_4 =  kcat7(1)*e7T_4;
dst.dF_6 =  kcat7(2)*e7T_6;
dst.dF_8 =  kcat7(3)*e7T_8;
dst.dF_10 = kcat7(4)*e7T_10;
dst.dF_12 = kcat7(5)*e7T_12;
dst.dF_14 = kcat7(6)*e7T_14;
dst.dF_16 = kcat7(7)*e7T_16;
dst.dF_18 = kcat7(8)*e7T_18;
dst.dF_20 = kcat7(9)*e7T_20;


dst.de2p2 = k2_1f*e2*p2 + k2_2r*e2act*p3 - k2_1r*e2p2 - k2_2f*e2p2;
dst.de2act = k2_2f*e2p2 + k2_3r*e2acts6 - k2_2r*e2act*p3 - k2_3f*e2act*s6;
dst.de2acts6 = k2_3f*e2act*s6 + k2_4r*e2*p4 - k2_3r*e2acts6 - k2_4f*e2acts6;

dst.de3s3 = k3_1f*e3*s3 + k3_2r*e3act*p3 - k3_1r*e3s3 - k3_2f*e3s3;
dst.de3act = k3_2f*e3s3 + k3_3r*e3actp4 - k3_2r*e3act*p3 - k3_3f*e3act*p4...
    + k3_5r(1)*e3actT_4  - k3_5f(1)*e3act*T_4...
    + k3_5r(2)*e3actT_6  - k3_5f(2)*e3act*T_6...
    + k3_5r(3)*e3actT_8  - k3_5f(3)*e3act*T_8...
    + k3_5r(4)*e3actT_10 - k3_5f(4)*e3act*T_10...
    + k3_5r(5)*e3actT_12 - k3_5f(5)*e3act*T_12...
    + k3_5r(6)*e3actT_14 - k3_5f(6)*e3act*T_14...
    + k3_5r(7)*e3actT_16 - k3_5f(7)*e3act*T_16...
    + k3_5r(8)*e3actT_18 - k3_5f(8)*e3act*T_18...
    + k3_5r(9)*e3actT_20 - k3_5f(9)*e3act*T_20...
    ;
dst.de3actp4 = k3_3f*e3act*p4 - k3_3r*e3actp4 - kcat3*e3actp4;
dst.de4s7 = k4_1f(1)*e4*s7 - k4_1r(1)*e4s7...
    + k4_2r(1)*e4s7Q_4 - k4_2f(1)*e4s7*Q_4...
    + k4_2r(2)*e4s7Q_6 - k4_2f(2)*e4s7*Q_6...
    + k4_2r(3)*e4s7Q_8 - k4_2f(3)*e4s7*Q_8...
    + k4_2r(4)*e4s7Q_10- k4_2f(4)*e4s7*Q_10...
    + k4_2r(5)*e4s7Q_12- k4_2f(5)*e4s7*Q_12...
    + k4_2r(6)*e4s7Q_14- k4_2f(6)*e4s7*Q_14...
    + k4_2r(7)*e4s7Q_16- k4_2f(7)*e4s7*Q_16...
    + k4_2r(8)*e4s7Q_18- k4_2f(8)*e4s7*Q_18...
    + k4_2r(9)*e4s7Q_20- k4_2f(9)*e4s7*Q_20;
dst.de4s7Q_4 =  k4_2f(1)*e4s7*Q_4 - k4_2r(1)*e4s7Q_4 - kcat4(1)*e4s7Q_4;
dst.de4s7Q_6 =  k4_2f(2)*e4s7*Q_6 - k4_2r(2)*e4s7Q_6 - kcat4(2)*e4s7Q_6;
dst.de4s7Q_8 =  k4_2f(3)*e4s7*Q_8 - k4_2r(3)*e4s7Q_8 - kcat4(3)*e4s7Q_8;
dst.de4s7Q_10 = k4_2f(4)*e4s7*Q_10- k4_2r(4)*e4s7Q_10- kcat4(4)*e4s7Q_10;
dst.de4s7Q_12 = k4_2f(5)*e4s7*Q_12- k4_2r(5)*e4s7Q_12- kcat4(5)*e4s7Q_12;
dst.de4s7Q_14 = k4_2f(6)*e4s7*Q_14- k4_2r(6)*e4s7Q_14- kcat4(6)*e4s7Q_14;
dst.de4s7Q_16 = k4_2f(7)*e4s7*Q_16- k4_2r(7)*e4s7Q_16- kcat4(7)*e4s7Q_16;
dst.de4s7Q_18 = k4_2f(8)*e4s7*Q_18- k4_2r(8)*e4s7Q_18- kcat4(8)*e4s7Q_18;
dst.de4s7Q_20 = k4_2f(9)*e4s7*Q_20- k4_2r(9)*e4s7Q_20- kcat4(9)*e4s7Q_20;

dst.de5M_4 =  k5_1f(1)*e5*M_4 - k5_1r(1)*e5M_4 - kcat5(1)*e5M_4;
dst.de5M_6 =  k5_1f(2)*e5*M_6 - k5_1r(2)*e5M_6 - kcat5(2)*e5M_6;
dst.de5M_8 =  k5_1f(3)*e5*M_8 - k5_1r(3)*e5M_8 - kcat5(3)*e5M_8;
dst.de5M_10 = k5_1f(4)*e5*M_10- k5_1r(4)*e5M_10- kcat5(4)*e5M_10;
dst.de5M_12 = k5_1f(5)*e5*M_12- k5_1r(5)*e5M_12- kcat5(5)*e5M_12;
dst.de5M_14 = k5_1f(6)*e5*M_14- k5_1r(6)*e5M_14- kcat5(6)*e5M_14;
dst.de5M_16 = k5_1f(7)*e5*M_16- k5_1r(7)*e5M_16- kcat5(7)*e5M_16;
dst.de5M_18 = k5_1f(8)*e5*M_18- k5_1r(8)*e5M_18- kcat5(8)*e5M_18;
dst.de5M_20 = k5_1f(9)*e5*M_20- k5_1r(9)*e5M_20- kcat5(9)*e5M_20;

dst.de6s8 = k6_1f(1)*e6*s8 - k6_1r(1)*e6s8...
    + k6_2r(1)*e6s8R_4 - k6_2f(1)*e6s8*R_4...
    + k6_2r(2)*e6s8R_6 - k6_2f(2)*e6s8*R_6...
    + k6_2r(3)*e6s8R_8 - k6_2f(3)*e6s8*R_8...
    + k6_2r(4)*e6s8R_10- k6_2f(4)*e6s8*R_10...
    + k6_2r(5)*e6s8R_12- k6_2f(5)*e6s8*R_12...
    + k6_2r(6)*e6s8R_14- k6_2f(6)*e6s8*R_14...
    + k6_2r(7)*e6s8R_16- k6_2f(7)*e6s8*R_16...
    + k6_2r(8)*e6s8R_18- k6_2f(8)*e6s8*R_18...
    + k6_2r(9)*e6s8R_20- k6_2f(9)*e6s8*R_20;
dst.de6s8R_4 =  k6_2f(1)*e6s8*R_4 - k6_2r(1)*e6s8R_4 - kcat6(1)*e6s8R_4;
dst.de6s8R_6 =  k6_2f(2)*e6s8*R_6 - k6_2r(2)*e6s8R_6 - kcat6(2)*e6s8R_6;
dst.de6s8R_8 =  k6_2f(3)*e6s8*R_8 - k6_2r(3)*e6s8R_8 - kcat6(3)*e6s8R_8;
dst.de6s8R_10 = k6_2f(4)*e6s8*R_10- k6_2r(4)*e6s8R_10- kcat6(4)*e6s8R_10;
dst.de6s8R_12 = k6_2f(5)*e6s8*R_12- k6_2r(5)*e6s8R_12- kcat6(5)*e6s8R_12;
dst.de6s8R_14 = k6_2f(6)*e6s8*R_14- k6_2r(6)*e6s8R_14- kcat6(6)*e6s8R_14;
dst.de6s8R_16 = k6_2f(7)*e6s8*R_16- k6_2r(7)*e6s8R_16- kcat6(7)*e6s8R_16;
dst.de6s8R_18 = k6_2f(8)*e6s8*R_18- k6_2r(8)*e6s8R_18- kcat6(8)*e6s8R_18;
dst.de6s8R_20 = k6_2f(9)*e6s8*R_20- k6_2r(9)*e6s8R_20- kcat6(9)*e6s8R_20;

dst.de7T_4 =  k7_1f(1)*e7*T_4 - k7_1r(1)*e7T_4 - kcat7(1)*e7T_4;
dst.de7T_6 =  k7_1f(2)*e7*T_6 - k7_1r(2)*e7T_6 - kcat7(2)*e7T_6;
dst.de7T_8 =  k7_1f(3)*e7*T_8 - k7_1r(3)*e7T_8 - kcat7(3)*e7T_8;
dst.de7T_10 = k7_1f(4)*e7*T_10- k7_1r(4)*e7T_10- kcat7(4)*e7T_10;
dst.de7T_12 = k7_1f(5)*e7*T_12- k7_1r(5)*e7T_12- kcat7(5)*e7T_12;
dst.de7T_14 = k7_1f(6)*e7*T_14- k7_1r(6)*e7T_14- kcat7(6)*e7T_14;
dst.de7T_16 = k7_1f(7)*e7*T_16- k7_1r(7)*e7T_16- kcat7(7)*e7T_16;
dst.de7T_18 = k7_1f(8)*e7*T_18- k7_1r(8)*e7T_18- kcat7(8)*e7T_18;
dst.de7T_20 = k7_1f(9)*e7*T_20- k7_1r(9)*e7T_20- kcat7(9)*e7T_20;

dst.de8T_4 =  k8_1f(1)*e8*T_4 + k8_2r(1)*e8act4*s6 - k8_1r(1)*e8T_4 - k8_2f(1)*e8T_4;
dst.de8T_6 =  k8_1f(2)*e8*T_6 + k8_2r(2)*e8act6*s6 - k8_1r(2)*e8T_6 - k8_2f(2)*e8T_6;
dst.de8T_8 =  k8_1f(3)*e8*T_8 + k8_2r(3)*e8act8*s6 - k8_1r(3)*e8T_8 - k8_2f(3)*e8T_8;
dst.de8T_10 = k8_1f(4)*e8*T_10+ k8_2r(4)*e8act10*s6- k8_1r(4)*e8T_10- k8_2f(4)*e8T_10;
dst.de8T_12 = k8_1f(5)*e8*T_12+ k8_2r(5)*e8act12*s6- k8_1r(5)*e8T_12- k8_2f(5)*e8T_12;
dst.de8T_14 = k8_1f(6)*e8*T_14+ k8_2r(6)*e8act14*s6- k8_1r(6)*e8T_14- k8_2f(6)*e8T_14;
dst.de8T_16 = k8_1f(7)*e8*T_16+ k8_2r(7)*e8act16*s6- k8_1r(7)*e8T_16- k8_2f(7)*e8T_16;
dst.de8T_18 = k8_1f(8)*e8*T_18+ k8_2r(8)*e8act18*s6- k8_1r(8)*e8T_18- k8_2f(8)*e8T_18;


dst.de8act4 =  k8_2f(1)*e8T_4 + k8_3r(1)*e8act4p4 - k8_2r(1)*e8act4*s6 - k8_3f(1)*e8act4*p4;
dst.de8act6 =  k8_2f(2)*e8T_6 + k8_3r(2)*e8act6p4 - k8_2r(2)*e8act6*s6 - k8_3f(2)*e8act6*p4;
dst.de8act8 =  k8_2f(3)*e8T_8 + k8_3r(3)*e8act8p4 - k8_2r(3)*e8act8*s6 - k8_3f(3)*e8act8*p4;
dst.de8act10 = k8_2f(4)*e8T_10+ k8_3r(4)*e8act10p4- k8_2r(4)*e8act10*s6- k8_3f(4)*e8act10*p4;
dst.de8act12 = k8_2f(5)*e8T_12+ k8_3r(5)*e8act12p4- k8_2r(5)*e8act12*s6- k8_3f(5)*e8act12*p4;
dst.de8act14 = k8_2f(6)*e8T_14+ k8_3r(6)*e8act14p4- k8_2r(6)*e8act14*s6- k8_3f(6)*e8act14*p4;
dst.de8act16 = k8_2f(7)*e8T_16+ k8_3r(7)*e8act16p4- k8_2r(7)*e8act16*s6- k8_3f(7)*e8act16*p4;
dst.de8act18 = k8_2f(8)*e8T_18+ k8_3r(8)*e8act18p4- k8_2r(8)*e8act18*s6- k8_3f(8)*e8act18*p4;


dst.de8act4p4 =  k8_3f(1)*e8act4*p4 - k8_3r(1)*e8act4p4 - kcat8(1)*e8act4p4;
dst.de8act6p4 =  k8_3f(2)*e8act6*p4 - k8_3r(2)*e8act6p4 - kcat8(2)*e8act6p4;
dst.de8act8p4 =  k8_3f(3)*e8act8*p4 - k8_3r(3)*e8act8p4 - kcat8(3)*e8act8p4;
dst.de8act10p4 = k8_3f(4)*e8act10*p4- k8_3r(4)*e8act10p4- kcat8(4)*e8act10p4;
dst.de8act12p4 = k8_3f(5)*e8act12*p4- k8_3r(5)*e8act12p4- kcat8(5)*e8act12p4;
dst.de8act14p4 = k8_3f(6)*e8act14*p4- k8_3r(6)*e8act14p4- kcat8(6)*e8act14p4;
dst.de8act16p4 = k8_3f(7)*e8act16*p4- k8_3r(7)*e8act16p4- kcat8(7)*e8act16p4;
dst.de8act18p4 = k8_3f(8)*e8act18*p4- k8_3r(8)*e8act18p4- kcat8(8)*e8act18p4;


dst.de3T_4  = k3_4f(1)*e3*T_4 - k3_4r(1)*e3T_4;
dst.de3T_6  = k3_4f(2)*e3*T_6 - k3_4r(2)*e3T_6;
dst.de3T_8  = k3_4f(3)*e3*T_8 - k3_4r(3)*e3T_8;
dst.de3T_10 = k3_4f(4)*e3*T_10- k3_4r(4)*e3T_10;
dst.de3T_12 = k3_4f(5)*e3*T_12- k3_4r(5)*e3T_12;
dst.de3T_14 = k3_4f(6)*e3*T_14- k3_4r(6)*e3T_14;
dst.de3T_16 = k3_4f(7)*e3*T_16- k3_4r(7)*e3T_16;
dst.de3T_18 = k3_4f(8)*e3*T_18- k3_4r(8)*e3T_18;
dst.de3T_20 = k3_4f(9)*e3*T_20- k3_4r(9)*e3T_20;

dst.de3actT_4  = k3_5f(1)*e3act*T_4  - k3_5r(1)*e3actT_4;
dst.de3actT_6  = k3_5f(2)*e3act*T_6  - k3_5r(2)*e3actT_6;
dst.de3actT_8  = k3_5f(3)*e3act*T_8  - k3_5r(3)*e3actT_8;
dst.de3actT_10 = k3_5f(4)*e3act*T_10 - k3_5r(4)*e3actT_10;
dst.de3actT_12 = k3_5f(5)*e3act*T_12 - k3_5r(5)*e3actT_12;
dst.de3actT_14 = k3_5f(6)*e3act*T_14 - k3_5r(6)*e3actT_14;
dst.de3actT_16 = k3_5f(7)*e3act*T_16 - k3_5r(7)*e3actT_16;
dst.de3actT_18 = k3_5f(8)*e3act*T_18 - k3_5r(8)*e3actT_18;
dst.de3actT_20 = k3_5f(9)*e3act*T_20 - k3_5r(9)*e3actT_20;



%Binding of ACP to TesA
dst.de7s6 = k7_inh_f*e7*s6 - k7_inh_r*e7s6;


%Binding of ACP to FabH FabG FabZ FabI FabF
dst.de3s6 = k3_inh_f*e3*s6 - k3_inh_r*e3s6;
dst.de4s6 = k4_inh_f*e4*s6 - k4_inh_r*e4s6;
dst.de5s6 = k5_inh_f*e5*s6 - k5_inh_r*e5s6;
dst.de6s6 = k6_inh_f*e6*s6 - k6_inh_r*e6s6;
dst.de8s6 = k8_inh_f*e8*s6 - k8_inh_r*e8s6;




%Dynamically returns the calculated values by creating a new list which is
%the same as short list with "d" appended to the start and then assigning
%these to the dc vector
dc = zeros(length(short_labels),1);
for i=1:length(short_labels)
    dc(i) = dst.(dlabels{i});
end



end

