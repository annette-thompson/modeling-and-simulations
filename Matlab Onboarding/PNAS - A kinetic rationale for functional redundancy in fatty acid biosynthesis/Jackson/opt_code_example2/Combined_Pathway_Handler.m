function [total_obj] = Combined_Pathway_Handler(p_vec)


% if p_vec(2) < 100
%     total_obj = 1E8;
%     return
% end


param_names = {{'k1_1f';'k1_1r';'k1_2f';'k1_2r';'kcat1_1';'k1_3f';'k1_3r';...
'kcat1_2';'k2_1f';'k2_1r';'k2_2f';'k2_2r';'k2_3f';'k2_3r';'k2_4f';...
'k2_4r';'k3_1f';'k3_1r';'k3_2f';'k3_2r';'k3_3f';'k3_3r';'k3_4f';'k3_4r';'k3_5f';'k3_5r';'kcat3';...
'k4_1f';'k4_1f';'k4_1r';'k4_2f';'k4_2r';'kcat4';'k5_1f';'k5_1r';...
'kcat5';'k5_2r';'k6_1f';'k6_1r';'k6_2f';'k6_2r';'kcat6';'k7_1f';'k7_1r';...
'kcat7';'k8_1f';'k8_1r';'k8_2f';'k8_2r';'k8_3f';'k8_3r';'kcat8';...
'k9_1f';'k9_1r';'kcat9';'k9_2r';'k10_1f';'k10_1r';'k10_2f';'k10_2r';'k10_3f';'k10_3r';'kcat10'...
}};


%set of parameter values that are varied (or likely subject to change) in
%structure sent to parameterization function param_func.m
acp_bind = p_vec(12);%parameter "e"

field5 = 'opt_name'; val5 = {{'k2_4f','k3_4f','k3_4r','k3_5f','k3_5r','k7_1f','kcat7','k8_1f','k8_2f','k10_1f','k10_2f','kcat5','kcat8','kcat9','kcat10','k5_2r','k9_2r','k5_1f','k9_1f'}};%parameters with unique modifications in param_func
field6 = 'param_names'; val6 = param_names;
field10 = 'scaling_factor_init'; val10 = p_vec(1);%parameter "a1"
field11 = 'scaling_factor_elon'; val11 = p_vec(2);%parameter "a2"
field12 = 'scaling_factor_kcat'; val12 = p_vec(5);%parameter "c2"
field13 = 'scaling_factor_term'; val13 = p_vec(3);%parameter "a3"
field14 = 'scaling_factor_fabf'; val14 = p_vec(2);%parameter "a2" (option here to modify FabF scaling seperately)
field15 = 'scaling_factor_kcat_term'; val15 = p_vec(6);%parameter "c3"
%ACP_inh lists kon,koff (in this order) for inhibitory ACP binding
%[FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB]
field16 = 'ACP_inh'; val16 = [acp_bind*2.41E-04,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,acp_bind*2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02];
%field17 = 'enzyme_conc'; val17 = enz_conc;%Initial enzyme conentrations (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
field18 = 'inhibition_kds';val18 = (1/acp_bind).*[4335.05,23.67;824.7,30.1;967.5,7.55;251.52,8.484;128.78,2.509];%(Acyl-ACP binding FabH Kd values in pairs (binding to FabH and FabH*), first two values are Kd for 4-12, subsequent values are 14-20)
field19 = 'inhibition_on_rates';val19 = [0.3088,1.552];%On rates of acyl-ACP binding FabH (binding to FabH and FabH*)
field20 = 'kd_fits';val20 = [p_vec(8),p_vec(9),p_vec(4),p_vec(9)]; %[b2,b3,b1,b3] %Keq or Kd values that are fit
field21 = 'lin_param';val21=[p_vec(10) p_vec(11)];% [d1 d2]; TesA linear free energy slope and intercept (as a function of chain length)
field22 = 'scaling_factor_kcat_init';val22 = p_vec(7);%parameter "c1"
field23 = 'scaling_factor_FabA_unsat';val23 = p_vec(13);%parameter "f"
field24 = 'scaling_factor_FabAZ_kcat';val24 = p_vec(14);%parameter "c4"
field25 = 'TesA_fitting_source';val25 = 'Pf';%source of thioesterase fitting ('Pf'; Pfleger group measurements, 'Fox'; Fox group measurements', 'Non-native'; For alternative thioesterases)


opt_struct = struct(field5,val5,field6,val6,...
    field10,val10,field11,val11,field12,val12,field13,val13,field14,val14,field15,val15,...
    field16,val16,field18,val18,field19,val19,field20,val20,field21,val21,field22,val22,field23,val23,field24,val24,field25,val25);

%this loads the parameter data in test_data file in a tables stucture
tables = struct;
tables.('param_table') = readtable('est_param.csv','ReadRowNames',true);
tables.('km_table') = readtable('km_est.csv','ReadRowNames',true);
% tables.('kcat_table') = readtable('kcat.csv','ReadRowNames',true);

%generates the parameter values in paramter structure from the parameter
%function
param_struct = struct;
for i = 1:length(param_names{1})
    param_struct.(param_names{1}{i}) = param_func(param_names{1}{i},i,opt_struct,tables);
end

tp_vec = p_vec;
tp_vec(10) = abs(tp_vec(10));
tp_vec(11) = abs(tp_vec(11));
if any(tp_vec(tp_vec<=0))
    total_obj = 1E8;
    return
elseif tp_vec(4) > 6.29E4 || tp_vec(5) > 240 || tp_vec(6) > 1 || tp_vec(7) > 15
    total_obj = 1E8;
    return
elseif tp_vec(1) < .001 || tp_vec(2) < .001 || tp_vec(3) < .001 || tp_vec(5) < 0.01 || tp_vec(6) < 7.41E-05 || tp_vec(7) < 0.001 || tp_vec(8) < 2.46E-7 || tp_vec(9) < 2.46E-7
    total_obj = 1E8;
    return
end

for i = 1:length(param_names{1})
    if ismember(param_names{1}{i},{'k2_1f','k2_3f','k3_1f','k3_3f','k4_1f','k4_2f','k5_1f','k6_1f','k6_2f','k7_1f','k8_1f','k8_3f'})
        if ismember(param_names{1}{i},{'k2_1f','k2_3f','k3_1f'})
            if param_struct.(param_names{1}{i}) > 1650
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            elseif param_struct.(param_names{1}{i+1})/param_struct.(param_names{1}{i}) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(param_names{1}{i},{'k3_3f','k4_2f','k5_1f','k6_2f','k8_1f','k8_3f'})
            if param_struct.(param_names{1}{i}) > 629
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            elseif param_struct.(param_names{1}{i+1})/param_struct.(param_names{1}{i}) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(param_names{1}{i},{'k4_1f','k6_1f'})
            if param_struct.(param_names{1}{i}) > 1650
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            elseif param_struct.(param_names{1}{i+1})/param_struct.(param_names{1}{i}) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(param_names{1}{i},{'k7_1f'})
            if max(param_struct.(param_names{1}{i})) > 629
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            end
        end
    end
end

tic
load('JpMat.mat','JpMatPrime')

%obj_1 = Combined_Pathway_Solver(p_vec,2.878);
%obj_1 = Combined_Pathway_Solver(p_vec,9.712);
%obj_2 = Combined_Pathway_Solver(p_vec,33.81);
%obj_3 = Combined_Pathway_Solver(p_vec,328.4);

enz_conc4 = [0 1 1.092 1 1 1 10 1 1 1];
enz_conc5 = [0 1 31.45 1 1 1 10 1 1 1];
range_init = [0 150];

init_cond_1 = zeros(304,1);
init_cond_1(3) = 200;%s3 (Acetyl-CoA, 1mM)
init_cond_1(4) = 10;%s6 (holo ACP)
init_cond_1(5) = 1000;%s7 (NADPH)
init_cond_1(6) = 1000;%s8 (NADH)
init_cond_1(8) = 500;%p2 (malonyl-CoA)


obj_4 = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc4,range_init,init_cond_1,JpMatPrime);
obj_5 = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc5,range_init,init_cond_1,JpMatPrime);


% obj_sum = obj_4 + obj_5;
% 
% if obj_sum > 1E3
%     total_obj = 1E8;
% else

  
obj_54 = (0.369 - obj_5/obj_4)^2;



enz_conc_uns_1 = [0 1 1 1 10 10 30 1 1 1];
enz_conc_uns_3 = [0 1 1 1 0 10 30 1 10 1];
enz_conc_uns_5 = [0 1 1 1 10 10 30 1 1 10];
enz_conc_uns_6 = [0 1 1 1 10 1 30 1 1 1];

init_cond_unsat = zeros(304,1);
init_cond_unsat(3) = 500;%s3 (Acetyl-CoA, 1mM)
init_cond_unsat(4) = 30;%s6 (holo ACP)
init_cond_unsat(5) = 1000;%s7 (NADPH)
init_cond_unsat(6) = 1000;%s8 (NADH)
init_cond_unsat(8) = 500;%p2 (malonyl-CoA)

range_unsat = [0 30];

obj_unsat1 = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_uns_1,range_unsat,init_cond_unsat,JpMatPrime);
obj_unsat3 = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_uns_3,range_unsat,init_cond_unsat,JpMatPrime);
obj_unsat5 = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_uns_5,range_unsat,init_cond_unsat,JpMatPrime);
obj_unsat6 = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_uns_6,range_unsat,init_cond_unsat,JpMatPrime);

obj_unsat31 = (0.5992 - obj_unsat3/obj_unsat1)^2;
obj_unsat51 = (1.1305 - obj_unsat5/obj_unsat1)^2;
obj_unsat61 = (0.253 - obj_unsat6/obj_unsat1)^2;

obj_sum_rel = (1 + obj_54)*(1 + obj_unsat31 + obj_unsat51 + obj_unsat61);
%obj_sum_rel = 1;

enz_conc6 = [0 1 1 1 1 1 10 1 1 1];
range_long = [0 720];
init_cond_2 = zeros(304,1);
init_cond_2(3) = 500;%s3 (Acetyl-CoA, 1mM)
init_cond_2(4) = 10;%s6 (holo ACP)
init_cond_2(5) = 1000;%s7 (NADPH)
init_cond_2(6) = 1000;%s8 (NADH)
init_cond_2(8) = 500;%p2 (malonyl-CoA)

obj_6 = Combined_Pathway_Solver(p_vec,param_struct,enz_conc6,range_long,init_cond_2,JpMatPrime);
toc
total_obj = (obj_6)*obj_sum_rel;


total_obj
% end

