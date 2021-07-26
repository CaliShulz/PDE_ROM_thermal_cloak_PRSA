%SETPATH add the subfolder FEM_library and RB_library to the path 

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'RB_library')));
addpath(genpath(strcat(pwd,sslash,'FEM_library')));
rmpath(genpath(strcat(pwd,sslash,'FEM_library',sslash,'LinearSolver',sslash,'Mumps',sslash,'Libraries')));

fprintf('\n------------------------------------------')
fprintf('\n      *** Welcome to redbKIT ***')
fprintf('\n------------------------------------------\n\n')
