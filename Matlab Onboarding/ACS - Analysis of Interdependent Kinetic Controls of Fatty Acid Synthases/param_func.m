function param = param_func(name,index,opt_struct)
%param_func Outputs kinetic parameters given the pathway step and specified
%options
%Returns the kinetic parameters for each label in param_names from the
%Combined_Pathway_Solver. The kinetic parameters are specified as a single
%value for the initiation steps, and as a vector for the elongation and
%termination steps.
%   Input:
%       name: the label (a string) being passed
%       index: the index in the list of the label
%       opt_struct: structure containing all of the fit parameter inputs
%       and parameterization options
%   Output:
%       param: A vector or double (depending on the label being passed) of
%       kinetic parameters. Enzymes in elongation steps need a vector as a
%       different kinetic parameters for each subsequent elongation step is
%       an important possibility (right now can only set this for one
%       enzyme specified with opt_name, but ultimately want all elongation
%       enzymes to have different kinetic parameters for different chain
%       lengths)


opt_name = opt_struct.('opt_name');
elong_num = opt_struct.('elong_num');
dist_opt = opt_struct.('dist_opt');
param_names = opt_struct.('param_names');
scaling_factor_init = opt_struct.('scaling_factor_init');
scaling_factor_elon = opt_struct.('scaling_factor_elon');
scaling_factor_term = opt_struct.('scaling_factor_term');
scaling_factor_fabf = opt_struct.('scaling_factor_fabf');
scaling_factor_kcat_term = opt_struct.('scaling_factor_kcat_term');
scaling_factor_kcat = opt_struct.('scaling_factor_kcat');
scaling_factor_kcat_init = opt_struct.('scaling_factor_kcat_init');
inhibition_kds = opt_struct.('inhibition_kds');
inhibition_on_rates = opt_struct.('inhibition_on_rates');
kd_fits = opt_struct.('kd_fits');
lin_param = opt_struct.('lin_param');

lin_slope = lin_param(1);
lin_int = lin_param(2);

 
%Loads the parameter data in csv files as matlab tables
param_table = readtable('est_param.csv','ReadRowNames',true);%kinetic parameter estimates
km_table = readtable('km_est.csv','ReadRowNames',true);%km estimates
kcat_table = readtable('kcat.csv','ReadRowNames',true);%kcat estimates


%Searches for the index in which the first elongation step occurs (this is
%why param_names needs to have all initiation steps listed first)
for i = 1:length(param_names)
    param_val = char(param_names(i));
    first_char = param_val(1:2);
    if first_char == char('k4')
        split_index = i;
        break
    end
end


%Creates structure to hold the data from the param_table, with the table
%variable names assigned as fields
param_struct = struct;
for i = 1:length(param_names)
    if strcmp('kcat',param_names{i}(1:4))
        if strcmp('kcat7',param_names{i})
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat_term;
        elseif strcmp('kcat4',param_names{i}) || strcmp('kcat5',param_names{i}) || strcmp('kcat6',param_names{i})|| strcmp('kcat8',param_names{i})
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat;
        else
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat_init;
        end
        
    elseif ismember(param_names{i},{'k2_1f','k2_1r','k2_3f','k2_3r','k3_1f','k3_1r','k3_3f','k3_3r'})
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_init;
        if ismember(param_names{i},{'k2_1f','k2_3f'})
            kd_est = km_table{{param_names{i}},:}/(1 + kcat_table{{param_names{i}},:}/(scaling_factor_init*param_table{{param_names{i+1}},:}));
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_init;
        elseif ismember(param_names{i},{'k3_1f','k3_3f'})
            kd_est = km_table{{param_names{i}},:}/(1 + scaling_factor_kcat_init*kcat_table{{param_names{i}},:}/(scaling_factor_init*param_table{{param_names{i+1}},:}));
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_init;
        end
    elseif strcmp('k7_1f',param_names{i}) || strcmp('k7_1r',param_names{i})
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_term;
    elseif strcmp('k8_1f',param_names{i}) || strcmp('k8_1r',param_names{i})
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_fabf;
        if ismember(param_names{i},{'k8_1f'})
            kd_est = km_table{{param_names{i}},:}/(1 + scaling_factor_kcat*kcat_table{{param_names{i}},:}/(scaling_factor_fabf*param_table{{param_names{i+1}},:}));
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_fabf;
        end
    elseif strcmp(param_names{i},'k2_2f') || strcmp(param_names{i},'k2_4f') || strcmp(param_names{i},'k3_2f') || strcmp(param_names{i},'k8_2f') || strcmp(param_names{i},'k2_2r') || strcmp(param_names{i},'k2_4r') ||  strcmp(param_names{i},'k3_2r') || strcmp(param_names{i},'k8_2r');
        
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*kd_fits(4);

    else
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_elon;
        if ismember(param_names{i},{'k4_1f','k4_2f','k4_3f','k5_1f','k5_2f','k6_1f','k6_2f','k6_3f','k8_1f','k8_3f','k8_4f'})
            kd_est = km_table{{param_names{i}},:}/(1 + scaling_factor_kcat*kcat_table{{param_names{i}},:}/(scaling_factor_elon*param_table{{param_names{i+1}},:}));
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_elon;
        end
    end
end


%If the parameter being looked for has index lower than split_index then it
%is an initiation step and a double will be returned from the structure
if index < split_index
    param = param_struct.(name);
end

%Number of elongation steps
num_elong_steps = elong_num;

%Parametrization of parameters that differ from above procedure
if strcmp(dist_opt,'data_defined')
    if strcmp(name,'k3_4f') || strcmp(name,'k3_4r')%acyl-ACP binding to FabH
        param = zeros(1,num_elong_steps);
        if strcmp(name,'k3_4r')
            for i = 1:length(param)
                if i <= 5
                    param(i) = inhibition_kds(1,1)*inhibition_on_rates(1);
                else
                    param(i) = inhibition_kds(i-4,1)*inhibition_on_rates(1);
                end
            end
        else
            for i = 1:length(param)
                if i <= 5
                    param(i) = inhibition_on_rates(3);
                else
                    param(i) = inhibition_on_rates(1);
                end
            end
        end
    elseif strcmp(name,'k3_5f') || strcmp(name,'k3_5r')%acyl-ACP binding to activated FabH
        param = zeros(1,num_elong_steps);
        if strcmp(name,'k3_5r')
            for i = 1:length(param)
                if i <= 5
                    param(i) = inhibition_kds(1,2)*inhibition_on_rates(2);
                else
                    param(i) = inhibition_kds(i-4,2)*inhibition_on_rates(2);
                end
            end
        else
            for i = 1:length(param)
                if i <= 5
                    param(i) = inhibition_on_rates(2);
                else
                    param(i) = inhibition_on_rates(2);
                end
            end
        end
        
    elseif strcmp(name,'k7_1f') || strcmp(name,'kcat7')%TesA parameterization
        param = zeros(1,num_elong_steps);
        if strcmp(name,'k7_1f')
            for i = 1:num_elong_steps
                kd_long = exp(lin_slope*(i*2+2) + lin_int);
                kd_12 = exp(lin_slope*(12) + lin_int);
                ratio_val = kd_12/(0.519*14.79);
                kd_est = (ratio_val).*[473 293.9 52.986 14.79];
                if i<5
                    param(i) = param_struct.('k7_1r')/kd_est(i);
                else
                    param(i) = param_struct.('k7_1r')/kd_long;
                end
            end
        elseif strcmp(name,'kcat7')
            kcat_scaling = [0.0568,0.0509,0.1035,0.0158,0.25256,0.45819,1,1.221,1.5368];
            for i = 1:num_elong_steps
                param(i) = param_struct.(name)*kcat_scaling(i);
            end
        end
    elseif strcmp(name,'k8_1f')
        param = zeros(1,num_elong_steps);
        if strcmp(name,'k8_1f')%Useful for specifiying chain lenght specificity of FabF if desired
            for i = 1:num_elong_steps
                if i>5
                    param(i) = param_struct.('k8_1f');
                else
                    param(i) = param_struct.('k8_1f');
                end
            end
        end
    %paramterization of acyl-transfer steps
    elseif strcmp(name,'k2_2f') || strcmp(name,'k3_2f') || strcmp(name,'k8_2f')
        if strcmp(name,'k2_2f') || strcmp(name,'k3_2f')
            if strcmp(name,'k2_2f')
                param = param_struct.('k2_2r')/kd_fits(1);
            else
                param = param_struct.('k3_2r')/kd_fits(3);
            end
        else
            if strcmp(name,'k8_2f')
                for i = 1:num_elong_steps
                    param(i) = param_struct.('k8_2r')/kd_fits(3);
                end
            end
        end
    elseif strcmp(name,'k2_4f')
        param = param_struct.('k2_4r')/kd_fits(5);
    end

end

%Returns a vector of values for kinetic
%parameters in elongation steps (as they should have indicies higher than
%split index)
if index >= split_index && ~any(strcmp(opt_name,name))
    param = zeros(1,num_elong_steps);
    for i = 1:num_elong_steps
        param(i) = param_struct.(name);
    end
end




end

