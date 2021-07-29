
function [] = T_POD_training(mesh_name)

%% T_POD_training.m

%  Generate snapshot dataset for transient dynamics for the selected mesh
%  depending on the number of snapshots takes (METTERE) and the dataset is
%  rather big (). The dataset generated is saved in archive_data

% To change input parameters range change FOM.mu_min, FOM.mu_max
% To change dimension of the set change mu_train_Dimension



sslash = path_setup() ; % setup path 
load(strcat('archive_data',sslash,'FOM_setup_',mesh_name));

% use three parameters [ diffusivity source strength T_dir beta_g]
FOM.P      = 4;
FOM.mu_max = [5 15000   100 1e-10];        
FOM.mu_min = [1 500     0   1e-10];

mu_train_Dimension = 25;
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));
FOM.mu_train       = mu_train;

% save basis vector and matrices for affine evaluation
FOM.A_d_0       = FOM.A_d;
FOM.A_d_ocp_0   = FOM.A_d_ocp;
FOM.A_d_dir_0   = FOM.A_d_dir;

FOM.F_0     = FOM.F;
FOM.F_ocp_0 = FOM.F_ocp;

% Transient Parameters
param.dt        = 0.025;
param.T         = 5;
param.dimt      = param.T/param.dt + 1;
param.max_iter  = 25;
param.tol       = 1e-06;


FOM.beta   = 1e-10;
FOM.alfa_T = 0;                      % weight on terminal cost
FOM.alfa_R = 1;                      % weight on running cost


Z_ref     = [];
Q_opt     = [];
P_opt     = [];
U_opt     = [];

parfor (jj = 1:size(mu_train,1),10)    
%for jj = 1:size(mu_train,1)
    mu_test      = mu_train(jj,:);
    FOM_temp     = FOM;            % create a copy of FOM to parallelize
    
    % Solve steady-state OCP for adjoint final condition
    [~,~,~,~,~,FOM_temp] = solve_HF_OCP_SS(mu_test,FOM_temp);
    
    % Solve transient OCP
    [z , q_opt_in , p_opt_in , u_opt_in , J] = solve_HF_OCP(mu_test,FOM_temp,param);

    
    Z_ref = [Z_ref z       ];
    Q_opt = [Q_opt q_opt_in];
    P_opt = [P_opt p_opt_in];
    U_opt = [U_opt u_opt_in];
    
    J_cells{jj} = J;
    
end

% save dataset
dataset.Z_ref   = Z_ref;
dataset.Q_opt   = Q_opt;
dataset.P_opt   = P_opt;
dataset.U_opt   = U_opt;


% Solve test FOM
mu_test = [3.2 10000 50 1e-10];


tFOM_Start = tic();
% Solve steady-state problem for adjoint final condition
[z_SS,q_SS,p_SS,u_SS,J_SS,FOM] = solve_HF_OCP_SS(mu_test,FOM);
[z , q_opt_in , p_opt_in , u_opt_in , J] = solve_HF_OCP(mu_test,FOM,param);
tFOM = toc(tFOM_Start);

% save test case
test_case.mu_test = mu_test;
test_case.FOM.z   = z;
test_case.FOM.q   = q_opt_in;
test_case.FOM.p   = p_opt_in;
test_case.FOM.u   = u_opt_in;
test_case.FOM.J   = J;
test_case.FOM.t   = tFOM;


to_save.test_case = test_case;
to_save.dataset   = dataset;
to_save.FOM       = FOM;
to_save.param     = param;

data_set_name = strcat('archive_data',sslash,'T_dataset_',mesh_name);
save(data_set_name,'to_save','-v7.3');

end

