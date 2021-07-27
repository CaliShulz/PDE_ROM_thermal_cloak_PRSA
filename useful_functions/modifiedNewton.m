function [q,p,u,J,hist] = modifiedNewton(u,param,model,max_iter,tol)

J = [];
dimt   = param.dimt;
beta   = model.beta;
beta_g = model.beta_g;
A_u    = model.A_u;
M_u    = model.M_u;
B      = model.B;

% convergence data
hist.t = [];
hist.f = [];

[~,N_u] = size(B);

fprintf( ' Modified Newton Method \n ' ) ;
fprintf( ' Current Cost  || Iteration  || Grad norm || Tau \n ' )
tau = 1;

for it=1:max_iter

        tic 
        
        J= [ J compute_cost(u,param,model) ] ; 
        [q]     = integrate_state(u,param,model,1);
        [p]     = integrate_adjoint(q,param,model,1);
        
        grad_J  = (beta*M_u+beta_g*A_u)*u + transpose(B)*p;
        grad_J  = reshape(grad_J,N_u,dimt);
        J_test = J(end) +1 ;

        [N_u,dimt] = size(grad_J);
                          
        fprintf( '     %5.2f     ||     %d/%d  ||     %d  ||   %d  \n ' , J(end) , it , max_iter , norm(grad_J(:)) , tau);

        
        if beta_g == 0
            delta_u = M_u \ (-1/beta*grad_J);
        else
            delta_u = (beta*M_u+beta_g*A_u) \ (-grad_J);
        end

        tau = 1;
        while J_test > J(end) 


            u_test = u + tau*delta_u;                  % modified newton step
            J_test = compute_cost(u_test,param,model);
            tau    = 0.5*tau;

            if tau < 1e-08
                fprintf(" Cost is not decreasing along d direction ");
                return
            end

        end

        u = u + tau*delta_u;
        
        
        %TODO exit cycle once tolerance is reached
        time_iter = toc;
        hist.t = [hist.t time_iter];
        hist.f = [hist.f J(end)];
        
        if norm(grad_J(:)) < tol
            fprintf(" Convergence reached with tolerance %d ", tol);
            return
            
        end
        
end

end

