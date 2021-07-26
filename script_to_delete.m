clear all
clc
close all

sslash = path_setup() ; % setup path 

mesh_name = 'mesh_circular_cloak';
load(strcat('archive_data',sslash,'FOM_setup_',mesh_name));

h_min   = min(FOM.MESH.h);


% save basis vector and matrices for affine evaluation
FOM.A_d_0       = FOM.A_d;
FOM.A_d_ocp_0   = FOM.A_d_ocp;
FOM.A_d_dir_0   = FOM.A_d_dir;

FOM.F_0     = FOM.F;
FOM.F_ocp_0 = FOM.F_ocp;
% Transient Parameters

param.max_iter  = 2;



[N_z,~]   = size(FOM.A_d);
[N_q,~] = size(FOM.B);

mu = 5;
mu_test = [mu 0 100 1e-10] ;

%dt_max  = (h_min / mu * 2/3)/4;

param.dt = 0.01;
param.T       = 5;
param.dimt    = param.T/param.dt+1;

[FOM] = evaluate_theta_terms(mu_test,FOM);
[FOM] = assemble_ato_SS(FOM);

FOM.q0      = zeros(N_q,1);

FOM.M   = FOM.M_ocp;
FOM.A   = FOM.A_d_ocp+FOM.A_d_robin_ocp;
FOM.F   = FOM.F_ocp +  FOM.F_dir;

[N_q,N_u] = size(FOM.B);
u         = 0*ones(N_u,param.dimt);

tic
[q]       = integrate_state(u,param,FOM,1);
toc

MESH         = FOM.MESH;
vertices     = MESH.vertices;
elements     = MESH.elements;


q_opt_FOM    = FOM.T_dir*ones(length(vertices),param.dimt);  % Cambiare!

q_opt_FOM(FOM.nodes_ocp_in,:) = q;

[state_elements,state_boundaries] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
FOM_state_data.reduced.vertices     = MESH.vertices(:,FOM.nodes_ocp);
FOM_state_data.reduced.elements     = state_elements; 
FOM_state_data.reduced.indexes      = FOM.nodes_ocp;

FOM_state_data.name       = "state_FOM";
FOM_state_data.y          = q_opt_FOM;

FOM_state_plot_data.dt    = param.dt;
FOM_state_plot_data.limits = [min([q_opt_FOM(:)]) max([q_opt_FOM(:)])];
FOM_state_plot_data.title = "State (FOM)";


% Setup fonts for plots
font_label = 18;
font_title = 19;
font_legend = 10;


fonts_data.font_title = font_title;
fonts_data.font_label = font_label;

video_controlled_FOM     = setup_video("video_q_fom");
fig = gobjects(0);
set(0,'DefaultFigureVisible','off');

for tt=1:5:param.dimt
        
    % Controlled System
    fig(length(fig)+1)     = figure;
    FOM_state_plot_data.tt = tt;
    [fig] = plot_field(fig,FOM_state_data,FOM_state_plot_data,fonts_data);
    frame = getframe(gcf);
    writeVideo(video_controlled_FOM,frame);
    
    
    
    
end

close(video_controlled_FOM);




function [video_object] = setup_video(name_video)

video_object = VideoWriter(name_video);
video_object.FrameRate = 4;
open(video_object);

end
