function [sslash] = path_setup()
% setup redbkit and library path
% return sslash

    % useful to save files
    if ispc
        sslash = '\';
    elseif isunix
        sslash = '/';
    end

    % setup redbKit paths
    eval(' cd redbKit-2.2 ');
    eval(' setPath '       );
    eval(' cd .. '         );

    addpath('useful_functions');
    
end

