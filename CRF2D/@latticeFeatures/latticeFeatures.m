function obj = latticeFeatures(expandNode, expandEdge)

% Construct feature vectors for 2D lattice CRF as used in
% "Discriminative Random Fields: A Discriminative Framework for Contextual
%  Interaction in Classification", S. Kumar and M. Hebert", ICCV03
% We use the interaction potential described in
% "Discriminative Fields for Modeling Spatial Dependencies in Natural Images",
% Kumar and Herber, NIPS04
%
% Input:
% expandNode = 1 if we do a quadratic expansion in mkNodeFeatures, otherwise just prepend 1
% expandEdge = 1 if we do a quadratic expansion in mkEdgeFeatures, otherwise just prepend 1
% Note:  size of graph is determined during enter_evidence (since it may change)
%
% We assume the hidden states xi in {-1,+1}
% Node potentials are
% nodepot(xi, i) = exp[ xi w^T h(:,i) ]
% where h(:,i) is (derived from) the features at node i
%
% Edge potentials (for edge e=i-j) are
% edgepot(xi,xj,e) = exp[ xi xj v^T mu(:,e) ]
% 
% In the code, h(:,i) = nodeFeatures(:,i)
% mu(:,e) = edgeFeatures(:,e)

if nargin==0 % Used when objects are loaded from disk
  obj = init_fields;
  obj = class(obj, 'latticeFeatures');
  return;
end
firstArg = expandNode;
if isa(firstArg, 'latticeFeatures') %  used when objects are passed as arguments
  obj = firstArg;
  return;
end

obj = init_fields;
obj.nstates = 2;
obj.expandNode = expandNode;
obj.expandEdge = expandEdge;

obj = class(obj, 'latticeFeatures'); 

%%%%%%%

function obj = init_fields();

% make an empty object

obj.nstates = [];
obj.nnodes = [];
obj.nedges  = [];
obj.nrows = [];
obj.ncols = [];
obj.Dnode = [];
obj.Dedge = [];
obj.nodeFeatures = [];
obj.edgeFeatures = [];
obj.nodeLabels = [];
obj.edgeDirNum = [];
obj.outNbr = [];
obj.edgeEndsIJ = [];
obj.edgeEndsRC = [];
