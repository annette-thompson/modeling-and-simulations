load('r1000.mat');
%order of norm_EE_mean (FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
bar_labels = {'D','H','G','Z','A','F','B','I','T'};
reorder_norm_EE_mean_total_FA = [norm_EE_mean_total_FA(1:4) norm_EE_mean_total_FA(8) norm_EE_mean_total_FA(7) norm_EE_mean_total_FA(9) norm_EE_mean_total_FA(5:6)];
reorder_norm_EE_mean_unsat_frac = [norm_EE_mean_unsat_frac(1:4) norm_EE_mean_unsat_frac(8) norm_EE_mean_unsat_frac(7) norm_EE_mean_unsat_frac(9) norm_EE_mean_unsat_frac(5:6)];
reorder_norm_EE_mean_avg_chain = [norm_EE_mean_avg_chain(1:4) norm_EE_mean_avg_chain(8) norm_EE_mean_avg_chain(7) norm_EE_mean_avg_chain(9) norm_EE_mean_avg_chain(5:6)];
figure()
color = bar([reorder_norm_EE_mean_total_FA;reorder_norm_EE_mean_unsat_frac;reorder_norm_EE_mean_avg_chain]');
color(1).FaceColor = 'k';
color(2).FaceColor = 'r';
color(3).FaceColor = 'c';
set(gca,'xticklabel',bar_labels)
legend('Production','Unsaturation','Chain Length','Location','best')
ylim([0 1.1])
ylabel('Normalized EE','FontSize',18)
title('Fig 3E Sensitivity Analysis','FontSize',18)  