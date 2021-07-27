function [eta_plot] = beta_comparison(mesh_name,mu_test)

sslash = path_setup() ; % setup path 
load(strcat('archive_data',sslash,'FOM_setup_',mesh_name));

% save basis vector and matrices for affine evaluation
FOM.A_d_0       = FOM.A_d;
FOM.A_d_ocp_0   = FOM.A_d_ocp;
FOM.A_d_dir_0   = FOM.A_d_dir;

FOM.F_0     = FOM.F;
FOM.F_ocp_0 = FOM.F_ocp;

[N_z,~]   = size(FOM.A_d);
[N_q,N_u] = size(FOM.B);

betas = 10 .^ [ -1:-0.25:-12 ];

jj = 1;

for beta = betas
    
    FOM.beta   = beta;

    [FOM] = evaluate_theta_terms(mu_test,FOM);
    [FOM] = assemble_ato_SS(FOM);
    x_opt = FOM.A_big \ FOM.F_big;


    z = x_opt(1:N_z,1);
    q = x_opt(N_z         + (1:N_q),1);
    p = x_opt((N_z+N_q)   + (1:N_q), 1);
    u = x_opt((N_z+2*N_q) + (1:N_u),1);


    q_free = (FOM.A_d_ocp+FOM.A_d_robin_ocp) \ (FOM.F_ocp+FOM.F_dir);

    area_omega_obs    = sum(FOM.E_obs*FOM.E*FOM.F_dom);
    track_error       = transpose(q - FOM.E*z)      * FOM.M_obs * (q - FOM.E*z);
    track_error_free  = transpose(q_free - FOM.E*z) * FOM.M_obs * (q_free - FOM.E*z);

    MTE_unc        = sqrt(track_error_free/area_omega_obs);
    MTE_opt        = sqrt(track_error     /area_omega_obs);

    eta_cloak(jj) = (abs(MTE_unc-MTE_opt)) / MTE_unc; 
    jj = jj+1;
    
end

eta_plot = [eta_cloak ; betas ];

end
