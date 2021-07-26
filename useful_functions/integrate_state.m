function [q] = integrate_state(u,param,model,theta)

    %% given control and parameters solve forward problem for q using Crank-Nicolson
    %% for linear parabolic heat quation Mq_dot +Aq = B u + F
    
    % CN produces spurious oscillation
    % Backward Euler does not
    
    if nargin < 4
        theta = 0.5;
    end
    
    % added input forcing F

    M      = model.M;
    A      = model.A;
    B      = model.B;
    F      = model.F;
    q0     = model.q0;
    
    
    dt = param.dt;
    T = param.T;

    num_steps = T/dt;
    
    q = q0;
    A_minus = (M/dt - (1-theta)*A);
    A_plus  = (M/dt + theta*A);
    
    
    for tt=1:num_steps
        
        n = tt; % current time instant
        np1 = tt+1; % time instant + 1

        u_np1 = u(:,np1);
        u_n = u(:,n);
        
        b = A_minus*q(:,n) + (1-theta)*B*u_n +theta*B*u_np1   + F;
        q(:,np1) = A_plus \ b;

        
        
    end


end
