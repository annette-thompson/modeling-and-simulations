function [] = plotter(app,results)
%Determines which function to use based on which plot type is requested by
%the user. 

model_scope = app.ModelScopeDropDown.Value;
plot_type = app.PlotTypeDropDown.Value;

% add plot to front page of GUI

if strcmp(plot_type, 'Product Profile')
    length_dist(app.UIAxes,model_scope,results)
else
    time_course(app.UIAxes,model_scope,plot_type,results)
end


