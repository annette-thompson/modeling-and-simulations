%Optimizes model parameters toward matching model results (least squares) to data or
%specified outputs. Utilizes multi-parameter approaches (constrained
%or unconstrained)
function [p_vec_sol] = Combined_Pathway_Optimizer(p_vec0)


fun_evals_max = 1000;%max number of function evaluations
max_iter = 30000;%max number of optimiztion iterations

%set up optimization call

fitfunc = @(x) Combined_Pathway_Handler(x);

options = optimset('MaxFunEvals',fun_evals_max,'Display','iter','MaxIter',max_iter);
[p_vec_sol,~,~,~] = fminsearch(fitfunc,p_vec0,options);
