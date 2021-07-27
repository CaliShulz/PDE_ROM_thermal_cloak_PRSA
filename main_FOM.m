function [] = main_FOM(mesh_name)

%% Solve FOM OCP transient and SS and compare results

% Modified OtD method
% Final condition on adjoint p_N = p_SS


sslash = path_setup() ; % setup path 

load(strcat('archive_data',sslash,'FOM_setup_',mesh_name));

% save basis vector and matrices for affine evaluation
FOM.A_d_0       = FOM.A_d;
FOM.A_d_ocp_0   = FOM.A_d_ocp;
FOM.A_d_dir_0   = FOM.A_d_dir;

FOM.F_0     = FOM.F;
FOM.F_ocp_0 = FOM.F_ocp;

% Transient Parameters
param.dt      = 0.05;
param.T       = 5;
param.dimt    = param.T/param.dt + 1;

dt = param.dt;

FOM.beta   = 1e-07;
FOM.alfa_T = 0;                                     % weight on terminal cost
FOM.alfa_R = 1;                                     % weight on running cost
param.max_iter  = 25;
param.tol       = 1e-06;

% System dimension
[N_z,~]   = size(FOM.A_d);
[N_q,N_u] = size(FOM.B);

% Solve test FOM - Steady-state
mu_test = [3.5 10000 0 0];

tic;
[FOM] = evaluate_theta_terms(mu_test,FOM);
[FOM] = assemble_ato_SS(FOM);
% Solve for snapshot AtO
x_opt = FOM.A_big \ FOM.F_big;
t_FOM_SS = toc;
z_SS = x_opt(1:N_z,1);
q_SS = x_opt(N_z         + (1:N_q),1);
p_SS = x_opt((N_z+N_q)   + (1:N_q), 1);
u_SS = x_opt((N_z+2*N_q) + (1:N_u),1);
delta_q_SS  = q_SS - FOM.E*z_SS;
J_SS = 0.5*(FOM.beta * ( transpose(u_SS) * FOM.M_u * u_SS ) + FOM.beta_g * ( transpose(u_SS) * FOM.A_u * u_SS )+ ...
                +   transpose(delta_q_SS)* FOM.M_obs * delta_q_SS ) ;
            
FOM.p_SS = p_SS;
FOM.q_SS = q_SS;
FOM.u_SS = u_SS;
FOM.z_SS = z_SS;


delta_q_SS_norm = transpose(delta_q_SS)  * FOM.M_obs * delta_q_SS;

J_SS    = 0.5* ( FOM.beta   * u_SS' * FOM.M_u * u_SS + ...
               mu_test(4)   * u_SS' * FOM.A_u * u_SS + ...
               (q_SS- FOM.E*z_SS)' * FOM.M_obs * (q_SS - FOM.E*z_SS));
        


% Solve FOM - Transient
tic;
[z_T , q_T , p_T , u_T , J_T , hist]     = solve_HF_OCP(mu_test,FOM,param);
tFOM_T = toc;

T    = param.T;
dimt = param.dimt;
delta_q     = q_T   - FOM.E*z_T;   % Tracking difference


J_T = 0.5*(trapz(diag(  FOM.beta    * ( transpose(u_T)    * FOM.M_u * u_T ) + FOM.beta_g * ( transpose(u_T) * FOM.A_u * u_T ) + ...
                   FOM.alfa_R  * ( transpose(delta_q)* FOM.M_obs * delta_q))) *dt                                        + ...
                +  FOM.alfa_T  * ( transpose(delta_q(:,end)) * FOM.M_obs * delta_q(:,end))) ;

J_T_time =  0.5*diag(  FOM.beta    * ( transpose(u_T)    * FOM.M_u * u_T ) + FOM.beta_g * ( transpose(u_T) * FOM.A_u * u_T ) + ...
                   FOM.alfa_R  * ( transpose(delta_q)* FOM.M_obs * delta_q))           ; 
            
%% Figures
close all

% Setup fonts for plots
font_label = 18;
font_title = 20;
font_legend =15;


fonts_data.font_title = font_title;
fonts_data.font_label = font_label;

% Setup figures
fig = gobjects(0);
set(0,'DefaultFigureVisible','on');


            
time_vect       = linspace(0,T,dimt);       
delta_q_norm    = diag(transpose(delta_q)* FOM.M_obs * delta_q);
delta_q_SS_norm = transpose(delta_q_SS)  * FOM.M_obs * delta_q_SS;



p_T_norm        = diag( transpose(p_T) * FOM.M_ocp * p_T );
p_SS_norm       =       transpose(p_SS)* FOM.M_ocp * p_SS;


fig(length(fig)+1) = figure;
fig(length(fig)).Name = "adjoint_norm_T_SS";
plot(time_vect,p_T_norm,'-','LineWidth',1.5)
hold on
plot(time_vect,p_SS_norm*time_vect./time_vect,'--r','LineWidth',1.5)
legend_temp = legend({'Transient','Steady-State'},'interpreter','latex','Location','southeast');
set(legend_temp,'FontSize',font_legend);
ylim([0,p_SS_norm*1.1]);
grid('minor');
label_temp = xlabel('Time [s]','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$ ||\mathbf{p}_h ||_{\tilde{M}}  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
title_temp = title('Adjoint norm','interpreter','latex');
set(title_temp,'FontSize',font_title);
axis square
grid on

yh = get(gca,'ylabel'); % handle to the label object
adjoint_label_pos = get(yh,'position'); % get the current position property



fig(length(fig)+1) = figure;
fig(length(fig)).Name = "track_error_T_SS";
semilogy(time_vect,delta_q_norm,'-','LineWidth',1.5)
hold on
semilogy(time_vect,delta_q_SS_norm*time_vect./time_vect,'--r','LineWidth',1.5)
grid('minor');
axis square
ylim([min(delta_q_norm)*1.1 max(delta_q_norm)*1.2])
label_temp = xlabel('Time [s]','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$ ||\mathbf{q}_h - E \mathbf{z}_h||_{M_{o}}  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
legend_temp = legend({'Transient','Steady-State'},'interpreter','latex','Location','southeast');
set(legend_temp,'FontSize',font_legend);
title_temp = title('Tracking error','interpreter','latex');
set(title_temp,'FontSize',font_title);

yh = get(gca,'ylabel'); % handle to the label object
current_pos = get(yh,'position'); % get the current position property
current_pos(1) = adjoint_label_pos(1);
set(yh,'position',current_pos)   % set the new position



grid on

q_T_norm        = diag( transpose(q_T) * FOM.M_ocp * q_T );
q_SS_norm       =       transpose(q_SS)* FOM.M_ocp * q_SS;

fig(length(fig)+1) = figure;
fig(length(fig)).Name = "state_norm_T_SS";
plot(time_vect,q_T_norm,'-','LineWidth',1.5)
hold on
plot(time_vect,q_SS_norm*time_vect./time_vect,'--r','LineWidth',1.5)
legend_temp = legend({'Transient','Steady-State'},'interpreter','latex','Location','southeast');
set(legend_temp,'FontSize',font_legend);
grid('minor');
label_temp = xlabel('Time [s]','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$ ||\mathbf{q}_h ||_{\tilde{M}}  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
title_temp = title('State norm','interpreter','latex');
set(title_temp,'FontSize',font_title);
axis square
grid on
yh = get(gca,'ylabel'); % handle to the label object
current_pos = get(yh,'position') ; % get the current position property
current_pos(1) = adjoint_label_pos(1);
set(yh,'position',current_pos)   % set the new position


fig(length(fig)+1) = figure;
fig(length(fig)).Name = "cost_T_SS";
semilogy(time_vect,J_T_time,'-','LineWidth',1.5)
hold on
semilogy(time_vect,J_SS*time_vect./time_vect,'--r','LineWidth',1.5)
grid('minor');
label_temp = xlabel('Time [s]','interpreter','latex');
set(label_temp,'FontSize',font_label);
legend_temp = legend({'Transient','Steady-State'},'interpreter','latex','Location','southeast');
set(legend_temp,'FontSize',font_legend);
label_temp = ylabel('$ J(t)  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
title_temp = title('Running cost','interpreter','latex');
set(title_temp,'FontSize',font_title);
axis square
grid on
yh = get(gca,'ylabel'); % handle to the label object
current_pos = get(yh,'position'); % get the current position property
current_pos(1) = adjoint_label_pos(1);
set(yh,'position',current_pos)   % set the new position

% u_T norm
u_T_norm        = diag( transpose(u_T) * FOM.M_u * u_T );
u_SS_norm       =       transpose(u_SS)* FOM.M_u * u_SS;

fig(length(fig)+1) = figure;
fig(length(fig)).Name = "control_norm_T_SS";
plot(time_vect,u_T_norm,'-','LineWidth',1.5)
hold on
plot(time_vect,u_SS_norm*time_vect./time_vect,'--r','LineWidth',1.5)
grid('minor');
legend_temp = legend({'Transient','Steady-State'},'interpreter','latex','Location','southeast');
set(legend_temp,'FontSize',font_legend);
ylim([0,u_SS_norm*1.1]);
label_temp = xlabel('Time [s]','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$ ||\mathbf{u}_h ||_{M_{u}}  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
title_temp = title('Control norm','interpreter','latex');
set(title_temp,'FontSize',font_title);
axis square
grid on
yh = get(gca,'ylabel'); % handle to the label object
current_pos = get(yh,'position'); % get the current position property
current_pos(1) = adjoint_label_pos(1);
set(yh,'position',current_pos)   % set the new position




z_T_norm        = diag( transpose(z_T) * FOM.M * z_T );
z_SS_norm       =       transpose(z_SS)* FOM.M * z_SS;

fig(length(fig)+1) = figure;
fig(length(fig)).Name = "reference_norm_T_SS";
plot(time_vect,z_T_norm,'-','LineWidth',1.5)
hold on
plot(time_vect,z_SS_norm*time_vect./time_vect,'--r','LineWidth',1.5)
grid('minor');
legend_temp = legend({'Transient','Steady-State'},'interpreter','latex','Location','southeast');
set(legend_temp,'FontSize',font_legend);

label_temp = xlabel('Time [s]','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$ ||\mathbf{z}_h ||_{M}  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
title_temp = title('Reference norm','interpreter','latex');
set(title_temp,'FontSize',font_title);
axis square
yh = get(gca,'ylabel'); % handle to the label object
current_pos = get(yh,'position'); % get the current position property
current_pos(1) = adjoint_label_pos(1);
set(yh,'position',current_pos)   % set the new position


% Save figures
name_sim = strcat('archive_sim',sslash,'FOM_only_T_SS');

folder_name = name_sim;
mkdir(folder_name);    
for jj=1:length(fig)

    name_tmp_jj = strcat(pwd,sslash,folder_name,sslash,fig(jj).Name,'.png');
    exportgraphics(fig(jj),name_tmp_jj);
    
end 


% % Save txt file with parameters of simulation
% 
% %generate report txt file
file_id = fopen(strcat(pwd,sslash,folder_name,sslash,'report','.txt'),'w');

fprintf(file_id,'################################### \n');
fprintf(file_id,'\n');
fprintf(file_id,'PARAMETERS MU TEST \n');
fprintf(file_id,'mu1 = %d, mu2 =%d , mu3 = %d , mu4 =%d \n ',mu_test);
fprintf(file_id,'TIME PARAMETERS \n');
fprintf(file_id,'dt = %d ,  T = %d , dimt = %d \n',[param.dt param.T param.dimt]);



fprintf(file_id,'################################### \n');
fprintf(file_id,'\n');
fprintf(file_id,'FOM \n');
fprintf(file_id,'################################### \n');
fprintf(file_id,evalc('disp(FOM)'));
fprintf(file_id,'################################### \n');



MESH         = FOM.MESH;
vertices     = MESH.vertices;
elements     = MESH.elements;


q_opt_FOM    = FOM.T_dir*ones(length(vertices),param.dimt);  % Cambiare!

q_opt_FOM(FOM.nodes_ocp_in,:) = q_T;

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


end

function [video_object] = setup_video(name_video)

video_object = VideoWriter(name_video);
video_object.FrameRate = 4;
open(video_object);

end





