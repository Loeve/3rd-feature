function  [nodePot, edgePot] = mkPotentials(featureEng, weights)

% Make data-dependent node and edge potentials
% in a format suitable for BP_edgeDir/infer
% nodePot(r,c,q)
% edgePot(q,q,ed) for e=1:num DIRECTED edges
%
% You must call enter_evidence first.
% We assume node/ edge weights are tied across nodes/ edges.

[nodeWeights, edgeWeights] = splitWeights(featureEng, weights);

tmp = nodeWeights(:).' * featureEng.nodeFeatures; % tmp(i) = sum_d w(d) f(d,i)
nodePot = [exp(tmp); exp(-tmp)]; % nodePot(q,i)
nodePot = reshape(nodePot.', [featureEng.nrows, featureEng.ncols, featureEng.nstates]);

tmp = edgeWeights(:).' * featureEng.edgeFeatures; % tmp(i) = sum_d w(d) f(d,i)
edgePot(1,1,:) = exp(tmp); % +1 * +1
edgePot(1,2,:) = exp(-tmp); % +1 * -1
edgePot(2,1,:) = exp(-tmp); % -1 * +1
edgePot(2,2,:) = exp(tmp); % -1 * -1
