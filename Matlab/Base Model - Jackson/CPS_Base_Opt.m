%Optimizes model parameters toward matching model results (least squares) to data or
%specified outputs. Utilizes multi-parameter approaches (constrained
%or unconstrained)
function [p_vec_sol] = CPS_Base_Opt(p_vec)

max_iter = 90;
fun_evals_max = 1000;

%set up optimization call
%fitfunc = @(x) CPS_Base_handler([106494.8474 x(1) 4.617637215 62895.56962 87.77959928 0.005962119 x(2:6) 9.912128884 132.8499358 2164.845252 0.539756276 0.053673263 34.49718991 11.15058888 0.17946003 0.123650181]);
%p_vec0 = [7726.392489 4.288620753 0.005308707 0.295368266 -0.296472666 3.658257772];

fitfunc = @(x) CPS_Base_handler(x);
%p_vec0 = [106494.8474 7726.392489 4.617637215 62895.56962 87.77959928 0.005962119 4.288620753 0.005308707 0.295368266 -0.296472666 3.658257772 9.912128884 132.8499358 2164.845252 0.539756276	0.053673263	34.49718991	11.15058888 0.17946003 0.123650181];

options = optimset('MaxFunEvals',fun_evals_max,'Display','iter','MaxIter',max_iter);
[p_vec_sol,~,~,~] = fminsearch(fitfunc,p_vec,options);


%Just put a for loop in here to run fminsearch? have to redefine fit func
%as p_vec changes. 