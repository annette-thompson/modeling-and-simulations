%% Step 1 (add paths)                                                          
my_dir = '/Users/Annette/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Matlab/PNAS - A kinetic rationale for functional redundancy in fatty acid biosynthesis' ;                      
                                                                               
% Set current directory to 'my_dir' and add path to sub-folders:               
cd(my_dir)                                                              
addpath(genpath(my_dir))
                                                                               
%% Step 2 (setup the Hymod model)                                              

% Number of uncertain parameters subject to SA:                                
M = 9 ;                                                                     
% Parameter ranges (from literature):                                          
xmin = [0.1 0.1 0.1 0 0.1 0.1 0 0 0];                                                     
xmax = [10 10 10 10 10 10 10 10 10];                                                     
% Parameter distributions:                                                     
DistrFun  = 'unif'  ;                                                          
DistrPar = cell(M,1);                                                          
for i=1:M; DistrPar{i} = [xmin(i) xmax(i)] ; end                             
% Name of parameters (will be used to costumize plots):                        
X_labels = {'D','H','G','Z','I','T','F','A','B'} ;                                    
                                                                             
% Define output:                                                               
myfun = 'Solver_func';

%% Step 3 (sample inputs space)                                                
                                                                               
r = 1000 ; % Number of Elementary Effects                                       
% [The final number of model evaluations will be equal to r*(M+1)]                                                                                              
                                                                                                              
SampStrategy = 'lhs' ; % Latin Hypercube                                       
design_type = 'radial';   

X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);              
                                                                               
%% Step 4 (run the model) 

[Y1,Y2,Y3] = model_execution(myfun,X); % size (r*(M+1),1)   
%I edited this function from the original

%% Step 5 (Computation of the Elementary effects)                              
                                                  
[miS, sigmaS] = EET_indices(r,xmin,xmax,X,Y3,design_type);  
[miA, sigmaA] = EET_indices(r,xmin,xmax,X,Y2,design_type); 
[miP, sigmaP] = EET_indices(r,xmin,xmax,X,Y1,design_type);                   
