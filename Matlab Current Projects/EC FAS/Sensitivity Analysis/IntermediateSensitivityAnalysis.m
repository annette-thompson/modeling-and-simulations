%% Step 1 (add paths)                                                          
% Set current directory to 'my_dir' and add path to sub-folders:               
my_dir ='/Users/annettethompson/Library/CloudStorage/OneDrive-UCB-O365/Annie Thompson/Git Repository/Matlab Current Projects/EC FAS';
cd(my_dir)
addpath(genpath(my_dir))
                                                                               
%% Step 2 (setup the Hymod model)                                              

% Number of uncertain parameters subject to SA:                                
%M = 9 ;   
M = 4; 
% Parameter ranges (from literature, 0 for redundant, pnas said up to 100 yielded same effects):                                          
%xmin = [0.1 0.1 0.1 0 0.1 0.1 0 0 0];                                                     
%xmax = [10 10 10 10 10 10 10 10 10];    
xmin = [0 0 0 1];
xmax = [1000 1000 1000 720];
% Parameter distributions:                                                     
DistrFun  = 'unif'  ;                                                          
DistrPar = cell(M,1);                                                          
for i=1:M; DistrPar{i} = [xmin(i) xmax(i)] ; end                             
% Name of parameters (will be used to customize plots):                        
%X_labels = {'D','H','G','Z','I','T','F','A','B'};  
X_labels = {'Bicarbonate','Acetyl-CoA','Malonyl-CoA','Time'}; 
                                                                             
% Define output:                                                               
myfun = 'Intermediate_Solver_func';

%% Step 3 (sample inputs space)                                                
                                                                               
r = 100; % Number of Elementary Effects                                       
% [The final number of model evaluations will be equal to r*(M+1)]                                                                                              
                                                                                                              
SampStrategy = 'lhs' ; % Latin Hypercube                                       
design_type = 'radial';   

Y = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);              
                                                                               
%% Step 4 (run the model) 

[Y1,Y2,Y3,Y4] = model_execution(myfun,X); % size (r*(M+1),1)   
%I edited this function from the original - made it parallel

%% Step 5 (compute elementary effects)                              
                                                  
[miA, sigmaA] = EET_indices(r,xmin,xmax,X,Y1,design_type);  
[miC, sigmaC] = EET_indices(r,xmin,xmax,X,Y2,design_type); 
[miR, sigmaR] = EET_indices(r,xmin,xmax,X,Y3,design_type);        
[miU, sigmaU] = EET_indices(r,xmin,xmax,X,Y4,design_type);        

save('Int_SA_results')