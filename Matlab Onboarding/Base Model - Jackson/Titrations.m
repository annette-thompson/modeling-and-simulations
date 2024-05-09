clear all
clc

%vector of initial conditions (units of uM)
init_cond = zeros(336,1);

init_cond(1) = 0;%s1 (ATP, not used)
init_cond(2) = 0;%s2 (Bicarbonate, not used)
init_cond(3) = 500;%s3 (Acetyl-CoA) %500 uM
init_cond(4) = 10;%s6 (holo ACP) %10 uM
init_cond(5) = 1000;%s7 (NADPH) %1000 uM
init_cond(6) = 1000;%s8 (NADH) %1000 uM

init_cond(7) = 0;%p1 (ADP, not generated)
init_cond(8) = 500;%p2 (malonyl-CoA) %500 uM
init_cond(9) = 0;%p3 (CoA)
init_cond(10) = 0;%p4 (malonyl-ACP)
init_cond(11) = 0;%p5 (CO2)

time_range = [0 720];


%% FabB Titration
%Check intitial conditions concentrations. 
avg_chain= zeros(101,1);
final_palm_equiv = zeros(101,1);
total_FA = zeros(101,1);
weight_vec = [4,6,8,10,12,12,14,14,16,16,18,18,20,20]/16;


for i = 1:36
    %(ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
    clear F_weighted
    clear frac_FA
    enz_conc = [0 1 10 1 1 1 10 i-1 1 1];
    [total_FA(i,1),F_weighted,frac_FA] = CPS_Base_solv(init_cond,enz_conc,time_range);
    avg_chain(i,1) = sum(frac_FA.*weight_vec.*16);
    final_palm_equiv(i,1) = F_weighted(end);

end
excel_copy_paste = [avg_chain, final_palm_equiv, total_FA]; 