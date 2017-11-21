function  OEF = mkObsEdgeFeatures(featureEng);
% OEF(d,eu) = xixj * mu(d,eu)
% where xixi = xi * xj = +/-1 is true label of edge e
% and eu = 1:numUndirEdges

L = featureEng.nodeLabels; % L(i)
%OEL(e) = xi*xj = +/-1 = observed edge label

numDirEdges = size(featureEng.edgeEndsIJ,1);
numUndirEdges = numDirEdges;%/2;
ndx = 1:numUndirEdges;
OEL = L(featureEng.edgeEndsIJ(ndx,1)) .* L(featureEng.edgeEndsIJ(ndx,2));

if 0
E = size(featureEng.edgeEndsIJ,1);
OEL2 = zeros(E, 1); 
for e=1:E
  i = featureEng.edgeEndsIJ(e,1);
  j = featureEng.edgeEndsIJ(e,2);
  %if i>j, continue; end 
  OEL2(e) = L(i)*L(j);
end
assert(approxeq(OEL,OEL2))
end

Dedge = featureEng.Dedge;
W = repmat(OEL(:).', Dedge, 1); % W(d,i) = -1 or +1
OEF = featureEng.edgeFeatures(:,ndx) .* W;
