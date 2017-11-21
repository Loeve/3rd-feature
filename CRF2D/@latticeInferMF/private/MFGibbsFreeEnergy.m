function [F,logZ] = MFGibbsfreeEnergy(nodePot,edgePot,nodeBel,V,E)
% Compute mean-fields Gibbs free energy as in:
%   "Understanding Belief Propagation and its Generalizations"
%       Yedidia, Freeman, Weiss
%
% nodePot(n,k) potential on state k for node n
% edgePot(k,k',e) potential for states (k,k') on edge e
% nodeBel(n,k) belief on state k for node n
% [V,E] star edge representation

[nnodes nstates] = size(nodePot);

U1 = 0;
U2 = 0;
S1 = 0;

nodePot=nodePot+eps;
edgePot=edgePot+eps;

for i =1:nnodes
 
  % Local Mean-Field Average Energy Term
  b = nodeBel(i,:);
  U1 = U1 + sum(b .* log(nodePot(i,:)));
  
  % Mean-Field Entropy Term
  b = nodeBel(i,:);
  b(b<.000001)=1;
  S1 = S1 + sum(b .* log(b));
  
  % For all Neighbors
  [nbrs,nbrsEdges] = StarEdge_FindNeighbors(V,E,i);
  nnbrs = length(nbrs);
  for n=1:length(nbrs)
      j=nbrs(n);
      if i > j
      e=nbrsEdges(n);
      
      b_i = repmat(nodeBel(i,:),[nstates 1]);
      b_j = repmat(nodeBel(j,:).',[1 nstates]);
      pot_ij = edgePot(:,:,e);
      
    % Pairwise Mean-Field Average Engery Term
      U2 = U2 + sum(b_i(:).*b_j(:).*log(pot_ij(:)));
      end
  end
end

F = -U2 - U1 + S1;
logZ = -F;
