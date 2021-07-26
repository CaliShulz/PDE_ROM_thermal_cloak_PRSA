function [model] = evaluate_theta_terms(mu,model)

mu_diff         = mu(1);
source_strength = mu(2);
T_dir           = mu(3);
beta_g          = mu(4);

% Evaluate matrices that linearly depends on diffusivity (TODO guardare RB
% EIM GAUSSIAN per come sono trattate le valutazioni)

model.A_d       = mu_diff * model.A_d_0;
model.A_d_ocp   = mu_diff * model.A_d_ocp_0;
model.A_d_dir   = mu_diff * model.A_d_dir_0;

% Evaluate RHS parametrize with source strength
model.F         = source_strength * model.F_0;
model.F_ocp     = source_strength * model.F_ocp_0;

% Evaluate RHS with Dirichlet Temperature
model.T_dir     = T_dir;

% Beta_g parameter
model.beta_g    = beta_g;




end

