function [p] = integrate_adjoint(q,param,model, theta)

% Integrate adjoint, discretization of Optimize-then-Discretize

% Adjoint equation is M p_dot + A p = 

if nargin < 4
    theta = 0.5;
end

alfa_R = model.alfa_R;

M     = model.M;
A     = model.A;
M_obs = model.M_obs;

z_r_num = model.z_r_num;

dt = param.dt;


A_minus = (M/dt - (1-theta)*A);
A_plus  = (M/dt + theta*A);

dimt = param.dimt;

p(:,dimt) = model.p_SS;
%p(:,dimt) = model.p_SS*0;   % uncomment to see the effect of homogenous adjoint boundaries

for ii = dimt:-1:2
    
    b = A_minus*p(:,ii) + theta     * ( alfa_R* (M_obs*q(:,ii-1)- M_obs*z_r_num(:,ii-1)) ) ... 
                        + (1-theta) * ( alfa_R* (M_obs*q(:,ii)  - M_obs*z_r_num(:,ii  )) );
                   
    p(:,ii-1) = A_plus \ b;
    
    
end




end

