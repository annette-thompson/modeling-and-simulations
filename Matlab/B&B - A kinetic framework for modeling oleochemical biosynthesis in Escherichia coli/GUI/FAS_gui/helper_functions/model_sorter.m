function results = model_sorter(user)
%Identifies where each expaced model is and runs the appropriate solver
%file. Finds the folder, enters it, runs the model, and then backs out of
%the folder. 

switch user.model_scope
    case 'Fatty Acids'
        fpath = 'Base Model';
        addpath(fpath)
        results = CPS_Base_solv_GUI(user);
        rmpath(fpath);
    case 'Alcohols'
        fpath = 'Alcohol';
        addpath(fpath)
        results = CPS_Alc_solv_GUI(user);
        rmpath(fpath);
    case 'Alcohols - ACR1'
        fpath = 'Alcohol - ACR1';
        addpath(fpath)
        results = CPS_Alc_solv_GUI(user);
        rmpath(fpath);
    case 'Alcohols - ACR2'
        fpath = 'Alcohol - ACR2';
        addpath(fpath)
        results = CPS_Alc_solv_GUI(user);
        rmpath(fpath);
    case 'Alcohols - ATR'
        fpath = 'Alcohol - ATR';
        addpath(fpath)
        results = CPS_Alc_solv_GUI(user);
        rmpath(fpath);
    case 'Alkanes'
        fpath = 'Alkane';
        addpath(fpath)
        results = CPS_Alk_solv_GUI(user);
        rmpath(fpath);
    case 'Methyl Ketones'
        fpath = 'Methyl Ketone';
        addpath(fpath)
        results = CPS_Ket_solv_GUI(user);
        rmpath(fpath);
    case 'Methyl Ketones - Thiolase'
        fpath = 'Methyl Ketone - Thiolase';
        addpath(fpath)
        results = CPS_Thiol_solv_GUI(user);
        rmpath(fpath);
    case 'FAMEs'
        fpath = 'FAME';
        addpath(fpath)
        results = CPS_FAME_solv_GUI(user);
        rmpath(fpath);
    case 'FAEEs'
        fpath = 'FAEE';
        addpath(fpath)
        results = CPS_FAEE_solv_GUI(user);
        rmpath(fpath);
    otherwise
        error('No model detected...')
end
end 

