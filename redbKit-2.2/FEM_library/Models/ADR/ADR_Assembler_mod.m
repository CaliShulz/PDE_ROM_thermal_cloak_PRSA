function [M,A_d,B_x,B_y,C_r] = ADR_Assembler_mod(MESH, DATA, FE_SPACE, OPERATOR, TC_d, TC_t, subdomain, t, stabilization, dt)
%
% modified version of ADR_ASSEMBLER, 
%
% return diffusion, transport x,y and reaction matrices
%
%


if nargin < 4 || isempty(OPERATOR)
    OPERATOR = 'all';
end

if nargin < 5 || isempty(TC_d)
    TC_d = [10 10];% diagonal components of the diffusion operator by default
end

if nargin < 6 || isempty(TC_t)
    TC_t = 10; % all components of the transport operator by default
end

if nargin < 7
    subdomain = [];
end

if nargin < 8
    t = [];
end

if nargin < 9
    stabilization = [];
end

if nargin < 10
    dt = 1;
end

if ~isempty(subdomain)    
    index_subd = [];
    for q = 1 : length(subdomain)
        index_subd = [index_subd find(MESH.elements(FE_SPACE.numElemDof+1,:) == subdomain(q))];
    end
    MESH.elements = MESH.elements(:,index_subd);
    MESH.numElem  = size(MESH.elements,2);
else
    index_subd = [1:MESH.numElem];
end


%% Computations of all quadrature nodes in the elements
coord_ref = MESH.chi;

switch MESH.dim
    case 2
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x;

        % riempe x ed y in un colpo senza for
        for j = 1 : 3   %triangoli
            i = MESH.elements(j,:);   % indici dei nodi j-th di tutti gli elementi
            vtemp = MESH.vertices(1,i); % coordinata x del nodo j-th di tutti gli elementi
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
        end
        
        %% Evaluation of coefficients in the quadrature nodes
        mu  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.diffusion, index_subd);
        si  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.reaction, index_subd);
        f   = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.force, index_subd);
        
        bx  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.transport{1}, index_subd);
        by  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.transport{2}, index_subd);
        
        b = [bx by];
        
    case 3
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x; z = x;

        for j = 1:4     %tetraedri
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,i);
            z = z + vtemp'*coord_ref(j,:);
        end
        
        %% Evaluation of coefficients in the quadrature nodes
        mu  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.diffusion, index_subd);
        si  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.reaction, index_subd);
        f   = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.force, index_subd);
  
        bx  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.transport{1}, index_subd);
        by  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.transport{2}, index_subd);
        bz  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.transport{3}, index_subd);
        
        b = [bx by bz];
end

%% Assembly
if isempty( stabilization )
    
    % C assembly, returns matrices in sparse vector format
    [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef,B_xcoef,B_ycoef,A_dcoef,C_coef] =  ADR_assembler_C_omp_mod(MESH.dim, OPERATOR, TC_d, TC_t, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,...
        FE_SPACE.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);
    
        
    % Build sparse matrices and rhs
    A    = GlobalAssemble(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);
    A_d  = GlobalAssemble(Arows,Acols,A_dcoef,MESH.numNodes,MESH.numNodes);
    B_x  = GlobalAssemble(Arows,Acols,B_xcoef,MESH.numNodes,MESH.numNodes);
    B_y  = GlobalAssemble(Arows,Acols,B_ycoef,MESH.numNodes,MESH.numNodes);
    C_r  = GlobalAssemble(Arows,Acols,C_coef,MESH.numNodes,MESH.numNodes);
    M    = GlobalAssemble(Arows,Acols,Mcoef,MESH.numNodes,MESH.numNodes);
    F    = GlobalAssemble(Rrows,1,Rcoef,MESH.numNodes,1);
    
else
        
    switch stabilization
       
        case 'SUPG'
            
            % C assembly, returns matrices in sparse vector format
            [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef] = ADR_SUPGassembler_C_omp(MESH.dim, stabilization, dt, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,...
                FE_SPACE.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);
            
            M = [];
            
        case 'SUPGt'
            
            % C assembly, returns matrices in sparse vector format
            [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef] = ADR_SUPGassembler_C_omp(MESH.dim, stabilization, dt, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,...
                FE_SPACE.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);
        
            M    = GlobalAssemble(Arows,Acols,Mcoef,MESH.numNodes,MESH.numNodes);
    end
    
    % Build sparse matrices and rhs
    A    = GlobalAssemble(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);
    F    = GlobalAssemble(Rrows,1,Rcoef,MESH.numNodes,1);
    
end

return