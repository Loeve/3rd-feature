function  ONF = mkObsNodeFeatures(featureEng);
% NF(d,i) = xi * h(d,i) where xi = +/-1 is true label of node i

L = featureEng.nodeLabels; % L(i)
Dnode = featureEng.Dnode;
W = repmat(L(:).', Dnode, 1); % W(d,i) = -1 or +1
ONF = featureEng.nodeFeatures .* W;
