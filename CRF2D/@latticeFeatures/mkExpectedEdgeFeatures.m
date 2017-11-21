function  EEF = mkExpectedEdgeFeatures(featureEng, edgeBel);
% edgeBel(q,q,eu) for eu=1:numUndirEdges
% EEF(d,e) = <xixj> * mu(d,eu) 


% EEL(e) = <xi*xj> = P(xi=1,xj=1)*1 + P(xi=2,xj=2)*1 -P(xi=1,xj=2) -P(xi=2,xj=1)
%  = expected label of edge e

numUndirEdges = size(edgeBel, 3); 

EEL = edgeBel(1,1,:) + edgeBel(2,2,:) - edgeBel(1,2,:) - edgeBel(2,1,:);

if 0
EEL2 = zeros(numUndirEdges, 1); 
for e=1:numUndirEdges
  EEL2(e) = edgeBel(1,1,e) + edgeBel(2,2,e) -edgeBel(1,2,e) - edgeBel(2,1,e);
end
assert(approxeq(EEL,EEL2))
end

Dedge = featureEng.Dedge;
W = repmat(EEL(:).', Dedge, 1); % W(d,e) =expcted value of edge e
ndx = 1:numUndirEdges;
EEF = featureEng.edgeFeatures(:,ndx) .* W;
