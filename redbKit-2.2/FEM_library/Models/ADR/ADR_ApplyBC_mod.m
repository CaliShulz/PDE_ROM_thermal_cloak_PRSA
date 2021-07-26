function [L] =  ADR_ApplyBC_mod(FE_SPACE, MESH, DATA, t)
% modified version of ADR_ApplyBC

% compute contribution of boundary basis control
%
%


param = DATA.param;


[csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
[phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
eta            =  1 - csi;
nqn            =  length(csi);


nof         = length(MESH.Robin_side);
nbn         = MESH.numBoundaryDof;

Arows       = zeros(nbn*nbn*nof,1);
Acols       = Arows;
Acoef       = Arows;
Rrows       = zeros(nbn*nof,1);
Rcoef       = Rrows;
[rows,cols] = meshgrid(1:nbn,1:nbn);
rows        = rows(:);
cols        = cols(:);

xlt = zeros(nof,nqn); ylt = xlt;
coord_ref = [eta; csi];
for j = 1 : 2
    dof = MESH.boundaries(j,MESH.Robin_side);
    vtemp = MESH.vertices(1,dof);
    xlt = xlt + vtemp'*coord_ref(j,:);
    vtemp = MESH.vertices(2,dof);
    ylt = ylt + vtemp'*coord_ref(j,:);
end

u_Robin = DATA.bcRob_fun(xlt,ylt,t,param);
alphaR  = DATA.bcRob_alpha(xlt,ylt,t,param);

one       = ones(nof,nqn);
u_Robin   = u_Robin.*one;
alphaR    = alphaR.*one;

x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Robin_side));
y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Robin_side));

side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);

for l = 1 : nof
    face = MESH.Robin_side(l);

    u_Robin_loc  = u_Robin(l,:).*wi;
    u_Robin_loc  = u_Robin_loc(1,:).';

    alphaR_loc   = alphaR(l,:).*wi;
    alphaR_loc   = alphaR_loc(ones(nbn, 1),:);

    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
    Rcoef(1+(l-1)*nbn:l*nbn)    = side_length(l)*phi*u_Robin_loc;

    Arows(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  MESH.boundaries(rows,face);
    Acols(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  MESH.boundaries(cols,face);
    Acoef(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  side_length(l)*(alphaR_loc.*phi)*phi';

end

L = sparse(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);

    
end

