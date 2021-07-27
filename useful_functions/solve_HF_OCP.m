function [z , q_opt , p_opt , u_opt , J , hist] = solve_HF_OCP(mu_test,FOM,param)

% Solve high fidelity time-dependent OCP for parameter mu_test

[N_z,~]   = size(FOM.A_d);
[N_q,~] = size(FOM.B);

[FOM] = evaluate_theta_terms(mu_test,FOM);
[FOM] = assemble_ato_SS(FOM);

% Solve reference dynamics
[FOM_ref] = assemble_time_dep_model(FOM.M,FOM.A_d_robin+FOM.A_d,zeros(N_z,1),FOM.F,zeros(N_z,1));
[z] = integrate_state(zeros(1,param.dimt),param,FOM_ref,1);

% Project to get reference state in the OCP mesh
FOM.z0_num  = FOM.E * z(:,end);
FOM.z_r_num = FOM.E * z; 
FOM.q0      = zeros(N_q,1);

FOM.M   = FOM.M_ocp;
FOM.A   = FOM.A_d_ocp+FOM.A_d_robin_ocp;
FOM.F   = FOM.F_ocp + FOM.F_dir;


% initial guess
u_0 = repmat(FOM.u_SS,1,param.dimt); 
[q_opt,p_opt,u_opt,J,hist_N]  = modifiedNewton(u_0,param,FOM,param.max_iter,param.tol);

hist.Newton = hist_N;


end

