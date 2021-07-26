%% Assemble_problem.m 

% Assemble data structures for selected mesh data
% Plot on screen relevant computational mesh figures
% Generate m files in archive_data folder (that is generated if
% nonexistent)

clear all
clc
close all

sslash = path_setup() ; % setup path 
set(0,'DefaultFigureVisible','on');     % plot figures on screen

% Load mesh according to test case to assemble
% change this to "mesh_circular_cloak" , "mesh_circular_cloak_coarse 
% "mesh_disconnected_cloak" or "mesh_boar_cloak" ( note that mesh_boar_cloak_fine is very fine and need a lot of RAM)

mesh_name = "mesh_circular_cloak";  
mesh_pet  = load(strcat(mesh_name,".mat"));
vertices         = mesh_pet.to_save.p;
boundaries       = mesh_pet.to_save.e;
elements         = mesh_pet.to_save.t;
control_dom_labels = mesh_pet.to_save.control_dom_labels;
obs_dom_labels     = mesh_pet.to_save.obs_dom_labels;
ocp_dom_labels     = mesh_pet.to_save.ocp_dom_labels;
boundary_ocp_labels= mesh_pet.to_save.boundary_ocp_labels;
outer_boundary_labels = mesh_pet.to_save.outer_boundary_labels;

% Visualize mesh
figure
pdeplot(vertices,boundaries,elements);


%% Assemble FOM matrices

dim_problem  = 2;
fem = 'P1';
quad_order = 4;

DATA       = read_DataFile('datafile_empty', dim_problem);
DATA.param = [];
DATA.flag_robin   = outer_boundary_labels;                      % outer boundaries are of Robin-type
DATA.bcRob_alpha    = @(x,y,t,param)(1+ 0.*x.*y);

mu_diff        = 1;
DATA.diffusion = @(x,y,t,param)(mu_diff + 0.*x.*y);

[ MESH ] = buildMESH( dim_problem, elements, vertices, boundaries, fem, quad_order, DATA ); % build full reference mesh

% get elements and boundaries of the ocp problem
elements_ocp  = [];
for dom_label = ocp_dom_labels
   
    elements_ocp = unique([elements_ocp find(elements(4,:)==dom_label)]);
    
end

nodes_ocp     = unique( [ elements(1,elements_ocp) elements(2,elements_ocp) elements(3,elements_ocp) ] ) ;

boundary_dof_ocp = 0;
for bound_label = unique(boundary_ocp_labels)
    boundary_dof_ocp = boundary_dof_ocp + (boundaries(5,:)==bound_label);
end
boundary_dof_ocp = find(boundary_dof_ocp==1);
boundary_dof_ocp = unique([boundaries(1,boundary_dof_ocp) boundaries(2,boundary_dof_ocp)]);


% Remove internal boundary nodes from nodes to keep
nodes_ocp_in  = setdiff(nodes_ocp,boundary_dof_ocp);


vertices_ocp = vertices(:,nodes_ocp);

[elements_ocp,boundaries_ocp] = get_reduced_mesh(MESH,nodes_ocp);

MESH_OCP.vertices   = vertices_ocp;
MESH_OCP.elements   = elements_ocp;
MESH_OCP.boundaries = boundaries_ocp;


figure
pdeplot(vertices_ocp,boundaries_ocp,elements_ocp)  % Visualize OCP mesh
hold on
scatter(vertices(1,boundary_dof_ocp),vertices(2,boundary_dof_ocp)) % Visualize boundaries of OCP mesh


% Extraction matrix, maps full nodes to ocp nodes
E = nodes_ocp_in' ./ [1:length(vertices)];
E(E~=1) = 0;
E       = sparse(E);

% Extraction matrix, maps full nodes to dirichlet nodes
E_dir           = boundary_dof_ocp' ./ [1:length(vertices)];
E_dir(E_dir~=1) = 0;
E_dir           = sparse(E_dir);


[ FE_SPACE ]   = buildFESpace( MESH, fem, 1, quad_order );
param.FE_SPACE = FE_SPACE;


% Assemble source right-hand side located in a small circular domain with
% constant intensity
DATA.force = @(x,y,t,param)( 1 + 0.*x.*y );
[~, F_dom, ~] = ADR_Assembler(MESH, DATA, FE_SPACE, [], [], [], [], 0);

x0 =  0.6;
y0 =  0;
r  = 0.025;
source_strength = 1;
id_source = @(x,y) ( ((x-x0).^2 + (y-y0).^2 - r) <=0 ) * source_strength;
DATA.force = @(x,y,t,param)( id_source(x,y) );
[~, F, ~] = ADR_Assembler(MESH, DATA, FE_SPACE, [], [], [], [], 0);


% Assemble mass and diffusion matrices for reference MESH
[M,~,A_d,~] = compute_state_input_matrices_simple(MESH, DATA, FE_SPACE); % assemble matrices
[A_d_robin, F, ~] =  ADR_ApplyBC(0*A_d, F, FE_SPACE, MESH, DATA, 0);  % Should save A_d (affine in mu) + A_robin (not mu dependent)

% Extraction matrix E, E_dir to get M,A and F for ocp domain
A_d_dir       = E * A_d       * transpose(E_dir);
A_d_robin_dir = E * A_d_robin * transpose(E_dir);
A_d_OCP       = E * A_d       * transpose(E);
A_d_robin_OCP = E * A_d_robin * transpose(E);
F_OCP         = E * F;
M_OCP         = E * M         * transpose(E);


% Set T boundary dirichlet
T_dir      = 1;
T_dir_vect = T_dir * sparse(ones(length(boundary_dof_ocp),1)); 

% Assemble control matrix B ( N_q x N_c mass matrix ) 

% get elements in which control is defined ( Omega_c in the paper)
elements_control  = [];
for dom_label = control_dom_labels
   
    elements_control = unique([elements_control find(elements(4,:)==dom_label)]);
    
end

control_basis_index = unique( [ elements(1,elements_control) elements(2,elements_control) elements(3,elements_control) ] ) ;

B     = M(nodes_ocp_in,control_basis_index);          % control matrix
M_u   = M(control_basis_index,control_basis_index);   % control mass matrix
A_u   = A_d(control_basis_index,control_basis_index); % control diffusion matrix ( to weight the gradient of the control in Omega_c )

% Visualize control nodes
%plot(vertices(1,control_basis_index),vertices(2,control_basis_index),'ok','MarkerFaceColor','c') 




vertices_control = vertices(:,control_basis_index);
[elements_control,~] = get_reduced_mesh(MESH,control_basis_index);
figure
pdeplot(vertices_control,elements_control(1:3,:));

% Assemble observation mass matrix M_obs 

% get observation elements
elements_obs  = [];
for dom_label = obs_dom_labels
   
    elements_obs = unique([elements_obs find(elements(4,:)==dom_label)]);
    
end

observation_basis_index = unique( [ elements(1,elements_obs) elements(2,elements_obs) elements(3,elements_obs) ] ) ;


% Observation matrix construction

% Extraction matrix, maps ocp_nodes to observation nodes
E_obs           = observation_basis_index' ./ nodes_ocp_in;
E_obs(E_obs~=1) = 0;
E_obs           = sparse(E_obs);


M_obs = sparse(size(M,1),size(M,2));
M_obs(observation_basis_index,observation_basis_index) = M(observation_basis_index,observation_basis_index);
M_obs = M_obs(nodes_ocp_in,nodes_ocp_in);

vertices_obs = vertices(:,observation_basis_index);
[elements_obs,~] = get_reduced_mesh(MESH,observation_basis_index);
figure
pdeplot(vertices_obs,elements_obs(1:3,:));

[N_z,~]   = size(A_d);
[N_q,N_u] = size(B);


% Plot different domains with different colors, replicate Figure 
figure
pdeplot(vertices,elements(1:3,:));
hold on
pdeplot(vertices_ocp,elements_ocp(1:3,:),'edgecolor','k');
pdeplot(vertices_obs,elements_obs(1:3,:),'edgecolor','g');
pdeplot(vertices_control,elements_control(1:3,:),'edgecolor','r');
scatter(vertices(1,boundary_dof_ocp),vertices(2,boundary_dof_ocp),'filled','MarkerFaceColor','b') 

beta   = 1E-07; % weights control norm
beta_g = 1E-06; % weights gradient control norm

A_big = [   A_d+A_d_robin     sparse(N_z,N_q)             sparse(N_z,N_q)             sparse(N_z,N_u) ; ...
          sparse(N_q,N_z)     A_d_OCP+A_d_robin_OCP       sparse(N_q,N_q)                  -B             ; ...
          M_obs*E                -M_obs                  A_d_OCP+A_d_robin_OCP        sparse(N_q,N_u) ; ...
          sparse(N_u,N_z)    sparse(N_u,N_q)              transpose(B)                beta*M_u+beta_g*A_u;        ];

F_dir =  -( A_d_dir + A_d_robin_dir) * T_dir_vect ;
F_big = [ F ; F_OCP+F_dir; sparse(N_q,1) ; sparse(N_u,1) ];

q_free_in = A_d_OCP \ (F_OCP - A_d_dir * T_dir_vect);


% Save necessary matrices for model order reduction
FOM.beta            = beta;
FOM.beta_g          = beta_g;
FOM.DATA            = DATA;
FOM.FE_SPACE        = FE_SPACE;
FOM.MESH            = MESH;
FOM.A_d_robin       = A_d_robin;
FOM.A_d_robin_dir   = A_d_robin_dir;
FOM.A_d             = A_d;
FOM.A_d_dir         = A_d_dir;

FOM.A_d_ocp         = A_d_OCP;
FOM.A_d_robin_ocp   = A_d_robin_OCP;



FOM.E               = E;
FOM.E_obs           = E_obs;
FOM.E_dir           = E_dir;
FOM.F               = F;
FOM.F_dom           = F_dom;
FOM.F_ocp           = F_OCP;
FOM.M_ref           = M_obs*E;
FOM.M               = M;
FOM.B               = B;
FOM.M_u             = M_u;
FOM.A_u             = A_u;

FOM.A_big           = A_big; 
FOM.M_ocp           = M_OCP;
FOM.M_obs           = M_obs;


FOM.control_basis_index     = control_basis_index; 
FOM.observation_basis_index = observation_basis_index;
FOM.nodes_ocp               = nodes_ocp;
FOM.nodes_ocp_in            = nodes_ocp_in;

FOM.name = mesh_pet.to_save.name;

if isfield(mesh_pet.to_save,'shape')
   
   FOM.shape = mesh_pet.to_save.shape;
    
end

mkdir('archive_data')
data_set_name = strcat('archive_data',sslash,'FOM_setup_',mesh_name);
save(data_set_name,'FOM','-v7.3');


 






