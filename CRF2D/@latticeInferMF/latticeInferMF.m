function obj = latticeInferMF(nstates, varargin)

% Mean field inference on a 2D lattice with pairwise potentials.
% We assume all nodes have the same number of states.
% However, the edge potentials can vary.

if nargin==0 % Used when objects are loaded from disk
  obj = init_fields;
  obj = class(obj, 'latticeInferMF');
  return;
end
firstArg = nstates;
if isa(firstArg, 'latticeInferMF') %  used when objects are passed as arguments
  obj = firstArg;
  return;
end

[maxIter, tol, maximize] = process_options(...
    varargin, 'maxIter', 200, 'tol', 1e-3, 'maximize', 0);

obj.nstates = nstates;
obj.maxIter = maxIter;
obj.tol = tol;
obj.maximize = maximize;

obj = class(obj, 'latticeInferMF'); 

%%%%%%%%%

function obj = init_fields()

obj.nstates = [];
obj.maxIter = [];
obj.tol = [];
obj.maximize = [];
