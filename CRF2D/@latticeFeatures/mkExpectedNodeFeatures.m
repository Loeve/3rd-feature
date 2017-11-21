function ENF = mkExpectedNodeFeatures(featureEng, belRC)
% belRC(r,c,q), q=1 or 2
% EF(d,i) = <xi> * h(d,i), where <xi> = P(xi=1)*1 + P(xi=2)*-1

[nr nc nstates] = size(belRC);
nodesBel = reshape(belRC, nr*nc, nstates).'; % nodeBel(q,i)

EL = nodesBel(1,:) - nodesBel(2,:);
Dnode = featureEng.Dnode;
W = repmat(EL(:).', Dnode, 1); % W(d,i)
ENF = featureEng.nodeFeatures .* W;
