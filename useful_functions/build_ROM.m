function [ROM] = build_ROM(N_ref,N_pq_basis,N_u_basis,FOM,dataset,N_t,flag)

Q_opt = dataset.Q_opt;
P_opt = dataset.P_opt;
U_opt = dataset.U_opt;
Z_ref = dataset.Z_ref;

%Given N_basis,
[N_h,N_all] = size(Q_opt);
N_snapshots = N_all/N_t;

if flag == "sequential"
    [V_ref,Sigma_ref]  = my_POD([Z_ref(:,1:N_t)],N_ref) ;
    [V_pq, Sigma_pq ]  = my_POD([Q_opt(:,1:N_t) P_opt(:,1:N_t) ],N_pq_basis);
    [V_u , Sigma_u  ]  = my_POD([U_opt(:,1:N_t)],N_u_basis);

    for jj = 2:N_snapshots

        [V_ref,Sigma_ref]  = my_POD([V_ref, Z_ref(:, (jj-1)*N_t + (1:N_t) )],N_ref)    ;


        [V_pq, Sigma_pq]   = my_POD([V_pq, Q_opt(:, (jj-1)*N_t + (1:N_t) ) P_opt(:, (jj-1)*N_t + (1:N_t) )],N_pq_basis);

        [V_u , Sigma_u]    = my_POD([V_u ,  U_opt(:, (jj-1)*N_t + (1:N_t) )],N_u_basis);


    end
end

if flag == "ato"
    
    [V_ref,Sigma_ref]  = my_POD([Z_ref],N_ref) ;
    [V_pq, Sigma_pq]  = my_POD([Q_opt P_opt],N_pq_basis); % try adding Z_ref to state basis
    [V_u , Sigma_u]  = my_POD([U_opt],N_u_basis);
    
end

% Steady-state OCP , all-at-once

% A            x          x      x         = F               z
% x          A_ocp        x     -B         = E F             q
% M_ocp*E   -M_ocp      A_ocp    x         = 0               p 
% x            x         B^T    aM_u       = 0               u


% OCP reduction
ROM.V_pq              = V_pq;
ROM.V_u               = V_u;
ROM.V_ref             = V_ref;


ROM.beta     = FOM.beta;
ROM.MESH     = FOM.MESH;
ROM.FE_SPACE = FOM.FE_SPACE;
ROM.DATA     = FOM.DATA;

if isfield(FOM,'alfa_T')
    
    ROM.alfa_T = FOM.alfa_T;
    ROM.alfa_R = FOM.alfa_R;
end

end

