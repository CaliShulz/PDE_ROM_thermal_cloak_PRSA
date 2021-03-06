function [I_f] = RBF_evaluate(x, RBF_data)
%RBF_EVALUATE evaluates RBF interpolant in a given set of points
%
%   [I_F] = RBF_EVALUATE(X, RBF_data) given a matrix X of size dimX
%   x nPoints and the struct RBF_data generated by RBF_SETUP, returns a row
%   vector I_F of size 1 x nPoints containing the values of the RBF
%   interpolant in the points X.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

interp_points     = RBF_data.x;

[dimI,  ~]      = size(interp_points);
[dimX,  ~] = size(x);

if (dimI~=dimX)
  error('x should have the same number of rows as RBF_data.x');
end

I_f = RBF_evaluate_Fast(RBF_data.RBF_function_type, interp_points, x,  RBF_data.constant, RBF_data.coeff);

% I_f = zeros(1, nPoints);       % Unrecognized function or variable 'nPoints'.


% for i = 1 : nPoints
%     
%     r =  (x(:,i)*ones(1,nI)) - interp_points;
%     r = sqrt(sum(r.*r, 1));
%     
%     s = RBF_data.coeff(nI+1) + sum(RBF_data.coeff(1:nI)'.*RBF_data.RBF_function(r, RBF_data.constant));
%     
%     for k = 1 : dimX
%         s = s + RBF_data.coeff(k+nI+1)*x(k,i);     
%     end
%     
%     I_f(i) = s;
% end

end