function s = LeastSquaresCalc(t_vec,mod_res,file_name)
%LEASTSQUARESCALC
%Returns the sum of error squared given model results and a set of
%validation data
    %Input:
        %t_vec: vector of time values from the model at which the model is
            %solved
        %mod_res: vector of concentration values corresponding to the time
            %values from time_vec to be used for validation
        %val_data: name of validatation file written as a string, note the
            %validation file must be csv with only two columns, one with
            %Time(min) as header, and the other with Concentration(um) as
            %header. Must be located in same directory as model
    %Output:
         %s: the sqaure root of the sum of the squared residuals between the
            %model results and the validation data



%Loads data from Experimental_Dataset ((1) Yu, X.; Liu, T.; Zhu, F.;
%Khosla, C. Proceedings of the National Academy of Sciences 2011, 108, 18643-18648.)
warning('off','all')
validation_data = readtable(file_name,'ReadRowNames',false);
warning('on','all')

time_data = validation_data.Time_min_;%Time data in min
conc_data = validation_data.Concentration_um_;%Concentration data in uM

xx = time_data;

%Creates a spline fit of all the model time points and fatty acid
%concentrations
pp = spline(t_vec,mod_res);

yy = ppval(pp,xx);%Calculates the spline fit values at the experimental data points

s = sum((yy-conc_data).^2);%Calculates the sum of squares error between the model and experimental values

end

