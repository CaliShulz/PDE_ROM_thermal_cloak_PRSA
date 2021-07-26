function [V,Sigma] = my_POD(S,N_tol)
%compute POD of S

% Form correlation matrix
K =  transpose(S)*S;
[PSI, Sigma]  = svd(full(K));
Sigma         = diag(Sigma);     % squared singular values


if N_tol < 1
      
      Sigma_N  = cumsum(Sigma);
  
      N        = find(N_tol.^2>=1-Sigma_N/sum(Sigma),1,'first');
else
      
      N        = min(N_tol,size(PSI,2));
    
end

PSI   =  PSI(:,1:N);

% Form POD basis
V   =  S*PSI;

for cc = 1 : N
    V(:,cc) = 1/(sqrt(Sigma(cc))) * V(:,cc);
end



end

