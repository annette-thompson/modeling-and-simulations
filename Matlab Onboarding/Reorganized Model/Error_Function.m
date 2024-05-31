%% Error_Function
function Results = Error_Function(param_struct,err_options)

names = fieldnames(param_struct);
for i = 1:numel(names)
    assignin('caller', names{i}, param_struct.(names{i}));
end

end