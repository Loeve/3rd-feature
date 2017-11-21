function [F,logZ] = betheMRF2(nodePot, edgePot,nodeBel,edgeBel,V,E)
% nodePot(n,k) potential on state k for node n
% edgePot(k,k',e) potential for states (k,k') on edge e
% nodeBel(n,k) belief on state k for node n
% edgeBel(k,k',e) belief for states (k,k') on edge e
% [V,E] star edge representation
%
% free energy F = energy - entropy = -ln Z = -loglik
% Bethe free energy approximates the entropy term.
% Exact Energy E= E2 + E1
%  E2 = -sum_{ij} sum_{x_i,x_j) b_ij(xi,xj) ln pot(xi,xj)
%  E1 = -sum_{i} sum_{x_i) b_i(xi) ln pot(xi)
% Approximate Entropy S = H2 + H1
%  H2 = -sum_{ij} sum_{x_i,x_j) b_ij(xi,xj) ln b_ij(xi,xj)
%  H1 = sum_{i} (qi-1) * sum_{x_i) b_i(xi) ln b_i(xi) %negative entropy
% where qi = num nbrs of node n

E1 = 0; E2 = 0;
H1 = 0; H2 = 0;
[nnodes nstates] = size(nodePot);


nodePot=nodePot+eps;
edgePot=edgePot+eps;
for i=1:nnodes
    [nbrs,nbrsEdges] = StarEdge_FindNeighbors(V,E,i);
    nnbrs = length(nbrs);
    b = nodeBel(i,:);
    b(b<.000001)=1;
    H1 = H1 + (nnbrs-1)*(sum(b .* log(b)));
    b = nodeBel(i,:);
    E1 = E1 - sum(b .* log(nodePot(i,:)));
    for n=1:length(nbrs)
        j=nbrs(n);
        if i > j
            e=nbrsEdges(n);
            pot_ij = edgePot(:,:,e);
            b = edgeBel(:,:,e);
            b(b<.000001)=1;
            H2 = H2 - sum(b(:) .* log(b(:)));
            b = edgeBel(:,:,e);
            E2 = E2 - sum(b(:) .* log(pot_ij(:)));
        end
    end
end

F = (E1+E2) - (H1+H2);
logZ = -F;
%fprintf('E1 = %.3f, E2 = %.3f, H1 = %.3f, H2 = %.3f\n',E1,E2,H1,H2);
