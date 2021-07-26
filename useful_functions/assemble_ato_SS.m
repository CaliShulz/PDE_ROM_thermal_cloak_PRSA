function [model] = assemble_ato_SS(model)

%Assemble optimality system + reference A_big, F_big for steady-state
%problem


[N_q,N_u]   = size(model.B);
[N_z,~]     = size(model.A_d);
[~,N_b]     = size(model.A_d_dir);


% control matrix and weighting
B           = model.B;
beta        = model.beta;
beta_g      = model.beta_g;

T_dir_vect  = model.T_dir * ones(N_b,1);

% forcing rhs
F           = model.F;
F_OCP       = model.F_ocp;

% diffusion matrices
A_d         = model.A_d;
A_d_OCP     = model.A_d_ocp;
A_d_dir     = model.A_d_dir;


% mass matrices
M_ref       = model.M_ref;
M_obs       = model.M_obs;
M_u         = model.M_u;
A_u         = model.A_u;

% robin boundary matrices
A_d_robin     = model.A_d_robin;
A_d_robin_OCP = model.A_d_robin_ocp;
A_d_robin_dir = model.A_d_robin_dir;


% Assemble matrices for all-at-once
A_big = [   A_d+A_d_robin     sparse(N_z,N_q)             sparse(N_z,N_q)             sparse(N_z,N_u) ; ...
          sparse(N_q,N_z)     A_d_OCP+A_d_robin_OCP       sparse(N_q,N_q)                  -B         ; ...
          M_ref                -M_obs                  A_d_OCP+A_d_robin_OCP           sparse(N_q,N_u); ...
          sparse(N_u,N_z)    sparse(N_u,N_q)              transpose(B)                beta*M_u+beta_g*A_u   ];

F_dir =  -( A_d_dir + A_d_robin_dir) * T_dir_vect ;
F_big = [ F ; F_OCP+F_dir; sparse(N_q,1) ; sparse(N_u,1) ];

model.A_big = A_big;
model.F_big = F_big;


model.F_dir = F_dir;




end

