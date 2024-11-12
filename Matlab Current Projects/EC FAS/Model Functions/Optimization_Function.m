%% Optimization_Function
function [x_sol, fval] = Optimization_Function()
% Runs parameter optimization
%   Output:
%       x_sol: parameter solution
%       fval: total objective value being minimized

%clear all;close all;clc

% Give access to all necessary folders
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))

% Set ODE solver options
ODE_options = odeset('RelTol',1e-6,'MaxOrder',5,'Vectorized','on');

% Set optimization options
format shortG % show significant figures in shortest way
fun_evals_max = 1000; %max number of function evaluations
max_iter = 300; %max number of optimiztion iterations
opt_options = optimset('MaxFunEvals',fun_evals_max,'MaxIter',max_iter,'Display','iter','OutputFcn',@optout);

% p_vec = [a1 a2 a3 b1 c2 c3...
%           c1 b2 b3 d1 d2 e...
%           f c4 x1 x2 x3 x4];
% p_vec = [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
%           4.645634978 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673...
%           132.8499358 2180.050007 0.539756276 0.053673263 34.49718991 11.15058888];

% Initial guess for parameters being optimized
% x0 = [161.8306e-003 1.6737e+000];
x0 = [4.645634978 0.1];

% p_vec with x(1),x(2),... for parameters being optimized
fitfunc = @(x) Optimization_Function_Handler(...
    [142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274...
    abs(x(1)) 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673...
    132.8499358 2180.050007 0.539756276 0.053673263 34.49718991 11.15058888 abs(x(2))],...
    ODE_options);

% Search for parameter values that make the simulation closest to
% experiments
[x_sol, fval] = fminsearch(fitfunc,x0,opt_options);

end