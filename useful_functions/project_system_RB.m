function [ROM] = project_system_RB(V_N,ROM,FOM)

%Project system onto lower dimensional space


V_u   = V_N.V_u;
V_pq  = V_N.V_pq;
V_ref = V_N.V_ref;

if isfield(FOM,'rhs')
    ROM.rhs               = cellfun(@(B) transpose(V_ref)*B,FOM.rhs,'UniformOutput',false);
    ROM.HyRED             = FOM.HyRED;
    ROM.Qf                = FOM.Qf;
end

ROM.B                 = transpose(V_pq) *FOM.B*V_u;

% mass matrices
ROM.M_u               = transpose(V_u)  *FOM.M_u  *V_u;
ROM.A_u               = transpose(V_u)  *FOM.A_u  *V_u;
ROM.M_obs             = transpose(V_pq) *FOM.M_obs*V_pq;
ROM.M_ref             = transpose(V_pq) *FOM.M_ref*V_ref;
ROM.M                 = transpose(V_ref)*FOM.M    *V_ref;
ROM.M_ocp             = transpose(V_pq) *FOM.M_ocp*V_pq;

% robin boundary matrices
ROM.A_d_robin     = transpose(V_ref)*FOM.A_d_robin    *V_ref;
ROM.A_d_robin_ocp = transpose(V_pq) *FOM.A_d_robin_ocp*V_pq;
ROM.A_d_robin_dir = transpose(V_pq) *FOM.A_d_robin_dir;

% Reduction of normalized diffusivity terms
ROM.A_d_0     = transpose(V_ref)*FOM.A_d_0*V_ref;
ROM.A_d_ocp_0 = transpose(V_pq)*FOM.A_d_ocp_0*V_pq;
ROM.A_d_dir_0 = transpose(V_pq)*FOM.A_d_dir_0;


% Reduction of normalized right-hand sides
ROM.F_0     = transpose(V_ref)*FOM.F_0;
ROM.F_ocp_0 = transpose(V_pq)*FOM.F_ocp_0;

% Reduction of restriction matrix (may be not needed)
ROM.E       = transpose(V_pq)*FOM.E*V_ref;

end

