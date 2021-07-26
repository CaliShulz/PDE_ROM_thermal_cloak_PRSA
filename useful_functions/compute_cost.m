function [J] = compute_cost(u,param,model)


    % parameters independent on FOM/ROM

    dimt = param.dimt;
    dt = param.dt;
    
    beta   = model.beta;
    alfa_R = model.alfa_R;
    alfa_T = model.alfa_T;
    
    %% Given a set of initial parameters, return the cost for a given u;
    M_u      = model.M_u;
    z0_num   = model.z0_num;
    z_r_num  = model.z_r_num;
    M_obs    = model.M_obs;
    
    [dim_c,~] = size(M_u);
    
    u = reshape(u,dim_c,dimt);
        
    % Solve state equation   
    q = integrate_state(u,param,model);
    
    
    % compute cost functional
    % u has dimension dim_control x dim_t
        
    control_cost_vect = diag(transpose(u)*M_u*u);
    state_running_cost_vect = diag(transpose(q-z_r_num)*M_obs*(q-z_r_num));
    
    J = alfa_T/2*transpose(q(:,end)-z0_num)*M_obs*(q(:,end)-z0_num) + beta/2*trapz(control_cost_vect)*dt + alfa_R/2*trapz(state_running_cost_vect)*dt;

end

