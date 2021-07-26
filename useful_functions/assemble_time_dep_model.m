function [model] = assemble_time_dep_model(M,A,B,F,q0)
%generate model structure to solve with integrate_state

model.M  = M;
model.A  = A;
model.B  = B;
model.F  = F;
model.q0 = q0;

end

