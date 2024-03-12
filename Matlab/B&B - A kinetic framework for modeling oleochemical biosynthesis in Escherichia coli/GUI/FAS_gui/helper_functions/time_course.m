function [] = time_course(myAxes,model_scope,plot_type,results)
%Plots time course data for each expanded model
%Enables visualization of intermediates within each expanded model. 

cla(myAxes,'reset')

switch model_scope
    case "Fatty Acids"
        switch plot_type
            case "Acyl-ACPs"
               T = results.time;
               Acyl_ACPs = results.acyl_ACP_time;
               plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
               T = results.time;
               FFAs = results.FFA_time;
               plot(myAxes,T,FFAs)
        end
    
    case "Alcohols"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "Fatty Aldehydes"
                T = results.time;
                Ald = results.aldehyde_time;
                plot(myAxes,T,Ald)
            case "Fatty Alcohols"
                T = results.time;
                Alc = results.alcohol_time;
                plot(myAxes,T,Alc)                
        end 
     case "Alcohols - ACR1"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "Acyl-CoAs"
                T = results.time;
                Acyl_CoA = results.aCoA_time;
                plot(myAxes,T,Acyl_CoA)
            case "Fatty Aldehydes"
                T = results.time;
                Ald = results.aldehyde_time;
                plot(myAxes,T,Ald)
            case "Fatty Alcohols"
                T = results.time;
                Alc = results.alcohol_time;
                plot(myAxes,T,Alc)                
        end 
     case "Alcohols - ACR2"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "Acyl-CoAs"
                T = results.time;
                Acyl_CoA = results.aCoA_time;
                plot(myAxes,T,Acyl_CoA)
            case "Fatty Alcohols"
                T = results.time;
                Alc = results.alcohol_time;
                plot(myAxes,T,Alc)                
        end 
     case "Alcohols - ATR"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Alcohols"
                T = results.time;
                Alc = results.alcohol_time;
                plot(myAxes,T,Alc)                
        end 

    case "Alkanes"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Aldehydes"
                T = results.time;
                Ald = results.aldehyde_time;
                plot(myAxes,T,Ald)
            case "Alkanes"
                T = results.time;
                Alk = results.alkane_time;
                plot(myAxes,T,Alk)                
        end 

    case "Methyl Ketones"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "Acyl-CoAs"
                T = results.time;
                Acyl_CoAs = results.acyl_CoA_time;
                plot(myAxes,T,Acyl_CoAs)
            case "Methyl Ketones"
                T = results.time;
                MK = results.MK_time;
                plot(myAxes,T,MK)                
        end 
    case "Methyl Ketones - Thiolase"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "Acyl-CoAs"
                T = results.time;
                Acyl_CoAs = results.acyl_CoA_time;
                plot(myAxes,T,Acyl_CoAs)
            case "Methyl Ketones"
                T = results.time;
                MK = results.MK_time;
                plot(myAxes,T,MK)                
        end 

    case "FAMEs"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "FAMEs"
                T = results.time;
                FAME = results.FAME_time;
                plot(myAxes,T,FAME)                
        end 

    case "FAEEs"
        switch plot_type
            case "Acyl-ACPs"
                T = results.time;
                Acyl_ACPs = results.acyl_ACP_time;
                plot(myAxes,T,Acyl_ACPs)
            case "Fatty Acids"
                T = results.time;
                FFA = results.FFA_time;
                plot(myAxes,T,FFA)
            case "Acyl-CoAs"
                T = results.time;
                Acyl_CoAs = results.acyl_CoA_time;
                plot(myAxes,T,Acyl_CoAs)
            case "FAEEs"
                T = results.time;
                FAEE = results.FAEE_time;
                plot(myAxes,T,FAEE)                
        end 
end
                xlabel(myAxes, 'Time (s)')
                ylabel(myAxes, 'Concentration (\muM)')
                title(myAxes,'Product Profile')
end