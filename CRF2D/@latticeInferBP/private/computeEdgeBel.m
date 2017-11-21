function edgeBel = computeEdgeBel(nodeBel,nodePot,edgePot,msg,V,E)
% nodeBel(n,k) belief on state k for node n
% edgePot(k,k',e) potential for states (k,k') on edge e
% msg(e,k) message for state k on edge e
% [V,E] star edge representation
%
% Output:
% belE(k,k',e)
%
% Compute belief on edges given belief on marginal nodes and messages
% (star edge representation version)

[nnodes nstates] = size(nodeBel);
nedges = length(E);
edgeBel = zeros(nstates,nstates,nedges);

for i = 1:nnodes
    [nbrs,nbrsEdges] = StarEdge_FindNeighbors(V,E,i);
    for n=1:length(nbrs)
        j=nbrs(n);


        e = nbrsEdges(n); %E(i,j)

        for s_i = 1:nstates
            for s_j = 1:nstates

                % Product of Local Potentials with Edge Potential
                b = edgePot(s_i,s_j,e)*nodePot(i,s_i)*nodePot(j,s_j);

                % Times messages from neighbors of i except j
                [nbrs_i,nbrsEdges_i] = StarEdge_FindNeighbors(V,E,i);
                for n_i=1:length(nbrs_i)
                    k=nbrs_i(n_i);
                    if k == j
                        continue;
                    end
                    b = b*msg(StarEdge_EdgeNum(V,E,k,i),s_i);
                end

                % Times messages from neighbors of j except i
                [nbrs_j,nbrsEdges_j] = StarEdge_FindNeighbors(V,E,j);
                for n_j=1:length(nbrs_j)
                    l=nbrs_j(n_j);
                    if l == i
                        continue;
                    end
                    b = b*msg(StarEdge_EdgeNum(V,E,l,j),s_j);
                end

                edgeBel(s_i,s_j,e) = b;
            end

        end
        
        % Force the pairwise bel to be a probability distribution
        edgeBel(:,:,e)=edgeBel(:,:,e)/sum(sum(edgeBel(:,:,e)));
        
    end
end
