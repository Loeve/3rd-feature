function [nodeBel,niter] = MF_General(edgePot,nodePot,maxIter,optTol,V,E)
% edgePot(k,k',e) potential for states (k,k') on edge e
% nodePot(n,k) - Potential at node n for state k
% maxIter - Maximum number of iterations
% opTol - Optimality Tolerance
% starEdge_V - Vertex vector in star edge representation
% starEdge_E - Edge vector in star edge representation
%
% Output:
% nodeBel(n,k) - Belief at node i for state k
% niter - Number of iterations

[nnodes nstates] = size(nodePot);

% Initialize Beliefs to node Potentials
nodeBel = nodePot;

% Normalise beliefs
nodeBel = nodeBel ./ repmat(sum(nodeBel.').',[1 nstates]);

% DEBUG only:
% FreeEnergy = MFGibbsFreeEnergy(nodePot,edgePot,nodeBel,V,E)

for niter = 1:maxIter
    old_nodeBel = nodeBel;
    for i = 1:nnodes
        [nbrs,nbrsEdges] = StarEdge_FindNeighbors(V,E,i);
        for k_i = 1:nstates
            b = 0;
            % For all Neighbors
            nnbrs = length(nbrs);
            for n=1:length(nbrs)
                j=nbrs(n);
                e=nbrsEdges(n);
                  for k_j=1:nstates
                      b = b+(nodeBel(j,k_j)*log(edgePot(k_i,k_j,e)));
                  end
            end
            nb(i,k_i) = nodePot(i,k_i)*exp(b);
        end
        nodeBel(i,:) = nb(i,:)/sum(nb(i,:));
        nodeBel(i,:);
    end
    
    % DEBUG only:
    %FreeEnergy = MFGibbsFreeEnergy(nodePot,edgePot,nodeBel,V,E)
    
    % Test convergence
    if sum(sum(abs(nodeBel-old_nodeBel))) < optTol
        break;
    end
    
end


if niter < maxIter
    fprintf('Converged in %d iterations\n',niter);
else
    fprintf('Did not converge after max = %d iterations\n',maxIter);
end
