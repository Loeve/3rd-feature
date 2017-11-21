function edgeBel = computeEdgeBel(nodeBel,V,E)
% nodeBel(n,k) belief on state k for node n
% edgePot(k,k',e) potential for states (k,k') on edge e
% msg(e,k) message for state k on edge e
% [V,E] star edge representation
%
% Output:
% edgeBel(k,k',e)
%
% Compute pairwise beliefs
% For Mean Field, this is the product of beliefs of nodes
% (star edge representation version)

[nnodes nstates] = size(nodeBel);
nedges = length(E);
edgeBel = zeros(nstates,nstates,nedges);

for i=1:nnodes
    [nbrs,nbrsEdges] = StarEdge_FindNeighbors(V,E,i);
    for n=1:length(nbrs)
        j=nbrs(n);
        e1 = nbrsEdges(n); %E(i,j)
        for k_i = 1:nstates
            for k_j = 1:nstates
                edgeBel(k_i,k_j,e1) = nodeBel(i,k_i)*nodeBel(j,k_j);
            end
        end
    end
end