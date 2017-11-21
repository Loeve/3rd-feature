function  [belRC, MMRC, niter, edgeBel, logZ] = infer(infEng, nodePot, edgePot)
% INPUT
% nodePot(r,c,k) where k indexes states
% edgePot(k,k',e) where e is computed by assign_edge_nums_lattice
%   Hence each undirected edge has 2 numbers, depending on direction.
%   WE ASSUME YOU HAVE CREATED A POTENTIAL FOR EACH DIRECTION
%    (need be NOT be transposes of each other)
%
% OUTPUT
% belRC(r,c,k) 
% MMRC(r,c) = max marginal at node r,c = max_q P(X(r,c)=q|y)
% niter = number of iterations used
%
% The following outputs are only computed if you request them
% egdeBel(k,k',e)
% logZ


[nr,nc,nstates]  = size(nodePot);
[K,K2,nedges] = size(edgePot);

% Create general graph starFormat nodes and edge numberings
% (assign_edge_lattice.m has already puts the edges in the right order)
% So make you use this to make your edge order!
nodePot = reshape(nodePot,[nr*nc nstates]);
[V,E] = StarEdge_MakeEdgeNums_Lattice2(nr,nc);

% C code only handles doubles
nodePot = double(reshape(nodePot,[nr*nc nstates]));
edgePot = double(edgePot);

if(isreal(nodePot)&&isreal(edgePot))
    [nodeBel,niter,msgs,edgeBel] = BP_General_C(...
        edgePot,nodePot,infEng.maximize,infEng.maxIter,...
        infEng.tol,int32(V),int32(E));
else
    nodePot = complex(real(nodePot), imag(nodePot));
    edgePot = complex(real(edgePot), imag(edgePot));
    [nodeBel,niter,msgs,edgeBel] = BP_GeneralComplex_C(...
        edgePot,nodePot,infEng.maximize,infEng.maxIter,...
        infEng.tol,int32(V),int32(E));
end

belRC = reshape(nodeBel,[nr nc nstates]);
[junk, MMRC] = max(belRC, [], 3);
MMRC = -MMRC*2+3;

if nargout > 4
  [F,logZ] = betheMRF2log0(nodePot,edgePot,nodeBel,edgeBel,V,E);
end
