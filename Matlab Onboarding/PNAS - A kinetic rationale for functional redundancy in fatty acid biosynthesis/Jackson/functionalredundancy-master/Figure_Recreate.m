%vector of initial conditions (units of uM)
init_cond = zeros(304,1);

init_cond(1) = 0;%s1 (ATP, not used)
init_cond(2) = 0;%s2 (Bicarbonate, not used)
init_cond(3) = 500;%s3 (Acetyl-CoA) %500 uM
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1000;%s7 (NADPH)
init_cond(6) = 1000;%s8 (NADH)

init_cond(7) = 0;%p1 (ADP, not generated)
init_cond(8) = 500;%p2 (malonyl-CoA)
init_cond(9) = 0;%p3 (CoA)
init_cond(10) = 0;%p4 (malonyl-ACP)
init_cond(11) = 0;%p5 (CO2)

time_range = [0 720];

%% FigS15
% avg_chain= zeros(101,1);
% final_palm_equiv = zeros(101,1);
% total_FA = zeros(101,1);
% for i = 1:101
%     enz_conc = [0 1 i-1 1 1 1 10 1 1 4];
%     [avg_chain(i,1),final_palm_equiv(i,1),total_FA(i,1),~] = Combined_Pathway_Solver_figrecreate(init_cond,enz_conc,time_range);
% end

%% FigS19
% avg_chain = zeros(5,3);
% final_palm_equiv = zeros(5,3);
% total_FA = zeros(5,3);
% for i = 1:3
%     for j = 1:4
%         enz_conc = [0 1 1 1 1 1 10^(i-1) 10^(j-1) 1 1];
%         [avg_chain(j,i),final_palm_equiv(j,i),total_FA(j,i),~] = Combined_Pathway_Solver_figrecreate(init_cond,enz_conc,time_range);
%         %enz_conc = [0 1 1 1 1 1 10^(i-1) 0 1 1];
%         %[avg_chain(j,i),final_palm_equiv(j,i),total_FA(j,i),~] = Combined_Pathway_Solver_figrecreate(init_cond,enz_conc,time_range);
%     
%     end
% end

%% FigS20
avg_chain = zeros(11,1);
final_palm_equiv = zeros(11,1);
total_FA = zeros(11,1);
frac_FA = zeros(11,14);
for i = 1:11
    enz_conc = [0 1 (i/10-0.1) 1 1 1 0.1 10 1 1];
    [avg_chain(i,1),final_palm_equiv(i,1),total_FA(i,1),frac_FA(i,:)] = Combined_Pathway_Solver_figrecreate(init_cond,enz_conc,time_range);
end

