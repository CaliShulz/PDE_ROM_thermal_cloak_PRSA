function [z_SS,q_SS,p_SS,u_SS,J_SS,FOM] = solve_HF_OCP_SS(mu_test,FOM)
%Solve steady-state OCP with AtO Method


[FOM] = evaluate_theta_terms(mu_test,FOM)
[FOM] = assemble_ato_SS(FOM)

[N_z,~]   = size(FOM.A_d);
[N_q,N_u] = size(FOM.B);


% Solve for snapshot AtO
x_opt = FOM.A_big \ FOM.F_big;


z_SS = x_opt(1:N_z,1);
q_SS = x_opt(N_z         + (1:N_q),1);
p_SS = x_opt((N_z+N_q)   + (1:N_q), 1);
u_SS = x_opt((N_z+2*N_q) + (1:N_u),1);

% Salva anche in FOM
FOM.z_SS = z_SS;
FOM.q_SS = q_SS;
FOM.p_SS = p_SS;
FOM.u_SS = u_SS;


delta_q_SS  = q_SS - FOM.E*z_SS;
J_SS = 0.5*(FOM.beta * ( transpose(u_SS) * FOM.M_u * u_SS ) + FOM.beta_g * ( transpose(u_SS) * FOM.A_u * u_SS )+ ...
                +   transpose(delta_q_SS)* FOM.M_obs * delta_q_SS ) ;
      

end

