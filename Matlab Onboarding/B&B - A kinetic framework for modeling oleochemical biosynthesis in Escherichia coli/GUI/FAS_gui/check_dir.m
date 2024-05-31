function [error_code] = check_dir()
%CHECK_DIR is used to automatically navigate to the correct workspace
%directory, and if the proper directory cannot be found, 
error_code = 0;
thisDir = cd;
if ismac
    allDirs = split(thisDir,'/');
elseif ispc 
    allDirs = split(thisDir,'\');
end 

if ~strcmp(allDirs(end),'FAS_gui')
    try
        cd 'FAS_gui'
    catch
        try
            if ismac
                cd ../../'FAS_gui'
            elseif ispc
                cd ..\..\'FAS_gui'
            end
        catch
            error_code = 1;
        end
    end
end
end

