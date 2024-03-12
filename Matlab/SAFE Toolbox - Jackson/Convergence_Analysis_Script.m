Script for performing sensitivity analysis on fatty acid model


M = 15; %Number of model parameters

%Define parameter ranges for sampling
xmin = [1E2 1E0 1E2 1E0 1E0 1E0 1E0 1E0 1E-3 1E-1 1E-3 1E-3 -.5 0 .1];
xmax = [1E4 1E2 1E4 1E2 2.4E2 2.4E2 2.4E2 2.4E2 1E-1 1.5E1 1E-1 1E-1 0 10 200];



%Creat parameter distributions
DistrFun = 'unif';%uniform distribution in defined range

%create cell array with distribution parameters (min and max values for 
%uniform distribution)
DistrPar = cell(M,1);                                                          
for i=1:M
    DistrPar{i} = [ xmin(i) xmax(i) ] ; 
end

%define parameter label names
X_labels = {'InitSc','ElongSc','TermSc','AcylSc','Kc','TermKc','InitKc','Kd1','Kd2','slope','int','acp'};

%define model function
myfun = 'Combined_Pathway_Solver_function';


%number of elementary effects
r = 100;


SampStrategy = 'lhs';%sampling of latin hypercube
design_type = 'radial';%radial method of sampling

%generate samples
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);

%evaluate samples
Y = model_evaluation(myfun,X);


%calculate elementary effects
[EE_mean,EE_dev] = EET_indices(r,xmin,xmax,X,Y(:,1),design_type);

%plot elementary effects and std_dev
% EET_plot(mean_course,std_dev_course,X_labels)

%bootstrapping for confidence interval estimation
% Nboot=1000;                                                                     
% [mi,sigma,EE,mi_sd,sigma_sd,mi_lb,sigma_lb,mi_ub,sigma_ub] = EET_indices(r,xmin,xmax,X,Y(:,1),design_type,Nboot); 

%plot bootstrapped results
% EET_plot(mi,sigma,X_labels,mi_lb,mi_ub,sigma_lb,sigma_ub)


%perform convergence analysis with reduced sample(1/5th of original sample)
% rr = [ r/5:r/5:r ] ;                                                           
% m_r = EET_convergence(EE,rr);

% figure; plot_convergence(m_r,rr*(M+1),[],[],[],...                             
% 'no of model evaluations','mean of EEs',X_labels)    

%bootstrap for confidence estimation of convergence analysis
% Nboot = 1000;                                                                   
% rr = [ r/5:r/5:r ] ;                                                           
% [m_r,s_r,m_lb_r,m_ub_r] = EET_convergence(EE,rr,Nboot);                        
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
% figure; plot_convergence(m_r,rr*(M+1),m_lb_r,m_ub_r,[],...                     
% 'no of model evaluations','mean of EEs',X_labels)   

