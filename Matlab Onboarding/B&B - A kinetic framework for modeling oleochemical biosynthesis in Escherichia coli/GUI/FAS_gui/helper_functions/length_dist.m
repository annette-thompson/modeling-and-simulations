function [] = length_dist(myAxes,model_scope,results)
%Plots product profiles for each expanded model, and sets the appopriate
%axis tiles. 
switch model_scope
    case "Fatty Acids"
        FFA = results.FFA_dist;
        barnames = {'C4','C6','C8','C10','C12','C12:1','C14','C14:1','C16','C16:1','C18','C18:1'};
        bar(myAxes,FFA)
        
    case "Alcohols"
        Alc = results.alcohol_dist;
        barnames = {'C4','C6','C8','C10','C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
        bar(myAxes,Alc)
    case "Alcohols - ACR1"
        Alc = results.alcohol_dist;
        barnames = {'C6','C8','C10','C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
        bar(myAxes,Alc)
    case "Alcohols - ACR2"
        Alc = results.alcohol_dist;
        barnames = {'C6','C8','C10','C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
        bar(myAxes,Alc)
    case "Alcohols - ATR"
        Alc = results.alcohol_dist;
        barnames = {'C6','C8','C10','C12','C14','C16','C12:1','C14:1','C16:1'};
        bar(myAxes,Alc)

    case "Alkanes"
        Alk = results.alkane_dist;
        barnames = {'C15','C17','C15:1','C17:1'};
        bar(myAxes,Alk)

    case "Methyl Ketones"
        MK = results.MK_dist;
        barnames = {'C5','C7','C9','C11','C13','C15','C17','C11:1','C13:1','C15:1','C17:1'};
        bar(myAxes,MK)
    case "Methyl Ketones - Thiolase"
        MK = results.MK_dist;
        barnames = {'C5','C7','C9','C11','C13','C15','C17','C19','C11:1','C13:1','C15:1','C17:1','C19:1'};
        bar(myAxes,MK)

    case "FAMEs"
        FAME = results.FAME_dist;
        barnames = {'C12','C14','C16','C12:1','C14:1','C16:1',};
        bar(myAxes,FAME)
        
    case "FAEEs"
        FAEE = results.FAEE_dist;
        barnames = {'C12','C14','C16','C18','C12:1','C14:1','C16:1','C18:1'};
        bar(myAxes,FAEE)
end
        set(myAxes,'xticklabel',barnames)
        xlabel(myAxes, 'Chain Length')
        ylabel(myAxes, 'Concentration (\muM)')
        title(myAxes,'Product Profile')
end

