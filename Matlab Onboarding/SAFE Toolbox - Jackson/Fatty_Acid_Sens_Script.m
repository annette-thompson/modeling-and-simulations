clear all

%Script for performing sensitivity analysis on fatty acid model
M = 9; %Number of model parameters

%Define parameter ranges for sampling
%(FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
xmin = [0.1 0.1 0.1 0 0.1 0.1 0 0 0];
xmax = [10 10 10 10 10 10 10 10 10];

%Creat parameter distributions
DistrFun = 'unif';%uniform distribution in defined range

%create cell array with distribution parameters (min and max values for 
%uniform distribution)
DistrPar = cell(M,1);                                                          
for i=1:M
    DistrPar{i} = [ xmin(i) xmax(i) ] ; 
end

%define parameter label names
X_labels = {'FabD','FabH','FabG','FabZ','FabA','FabF','FabB','FabI','TesA'};
%define model function
myfun = 'Combined_Pathway_Solver_SA';


%number of elementary effects
r = 1000;


SampStrategy = 'lhs';%sampling of latin hypercube
design_type = 'radial';%radial method of sampling

%generate samples
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);

%evaluate samples
[total_FA,avg_chain,unsat_frac] = model_evaluation(myfun,X);


%calculate elementary effects
[EE_mean_total_FA,EE_dev_total_FA] = EET_indices(r,xmin,xmax,X,total_FA,design_type);
[EE_mean_avg_chain,EE_dev_avg_chain] = EET_indices(r,xmin,xmax,X,avg_chain,design_type);
[EE_mean_unsat_frac,EE_dev_unsat_frac] = EET_indices(r,xmin,xmax,X,unsat_frac,design_type);
norm_EE_mean_total_FA = EE_mean_total_FA./max(EE_mean_total_FA);
norm_EE_mean_avg_chain = EE_mean_avg_chain./max(EE_mean_avg_chain);
norm_EE_mean_unsat_frac = EE_mean_unsat_frac./max(EE_mean_unsat_frac);

% %plot elementary effects and std_dev
% EET_plot(mean_course,std_dev_course,X_labels)
% 
% %bootstrapping for confidence interval estimation
% Nboot=100;                                                                     
% [mi,sigma,EE,mi_sd,sigma_sd,mi_lb,sigma_lb,mi_ub,sigma_ub] = EET_indices(r,xmin,xmax,X,Y,design_type,Nboot); 
% 
% %plot bootstrapped results
% EET_plot(mi,sigma,X_labels,mi_lb,mi_ub,sigma_lb,sigma_ub)
% 
% 
% %perform convergence analysis with reduced sample(1/5th of original sample)
% rr = [ r/5:r/5:r ] ;                                                           
% m_r = EET_convergence(EE,rr);
% 
% figure; plot_convergence(m_r,rr*(M+1),[],[],[],...                             
% 'no of model evaluations','mean of EEs',X_labels)    
% 
% %bootstrap for confidence estimation of convergence analysis
% Nboot = 100;                                                                   
% rr = [ r/5:r/5:r ] ;                                                           
% [m_r,s_r,m_lb_r,m_ub_r] = EET_convergence(EE,rr,Nboot);                        
% % Plot the sensitivity measure (mean of elementary effects) as a function      
% % of model evaluations:                                                        
% figure; plot_convergence(m_r,rr*(M+1),m_lb_r,m_ub_r,[],...                     
% 'no of model evaluations','mean of EEs',X_labels)   
% 
