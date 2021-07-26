function [M_full,M_ff,A_d_ff,A_d,B] = compute_state_input_matrices_simple(MESH, DATA, FE_SPACE)

%just get M and A_d for linear parabolic problem

[M,A_d,~,~,~]  =  ADR_Assembler_mod(MESH, DATA, FE_SPACE, [], [], [], [], 0);

% get reduced matrix for homogenous dirichlet problem

M_full = M;

if isempty(MESH.Dirichlet_dof)
    dirichlet_dof = 0;
else
    dirichlet_dof = MESH.Dirichlet_dof(end);
end

M_ff   = M(MESH.internal_dof,MESH.internal_dof); 
A_d_ff = A_d(MESH.internal_dof,MESH.internal_dof);

 
B = M_ff;

end

