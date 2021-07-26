function [] = setup_path(main,parent)
%add necessary paths


if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

path_str = strcat(main,sslash,parent);
addpath(path_str)

end

