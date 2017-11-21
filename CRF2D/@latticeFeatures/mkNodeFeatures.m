function  nodeFeatures = mkNodeFeatures(featureEng, raw);

% raw(:,r,c,s)
% NF(:,i,s) = [1, raw(:,r,c,s))] if expandNode=0
% NF(:,i,s) = [1, quadratice(raw(:,r,c,s))] if expandNode=1

[Draw nr nc ncases] = size(raw);
if featureEng.expandNode
  D = length(expandFeatureVec(raw(:,1,1,1)));
else
  D = Draw+1;
end

nodeFeatures = zeros(D, nr, nc, ncases);
%size(nodeFeatures)
%size(raw)
%subdata = raw(:, 1, 1, 1)
%anydata = [1, subdata'];
%size(anydata)

for s=1:ncases
  for r=1:nr
    for c=1:nc
      if featureEng.expandNode
	h = expandFeatureVec(raw(:,r,c,s));
	nodeFeatures(:,r,c,s) = h/norm(h);
      else
    % BUG: Change made here: Vicente Ordonez (2009).
	nodeFeatures(:,r,c,s) = [1, raw(:,r,c,s)'];
      end
    end
  end
end
nodeFeatures = reshape(nodeFeatures, [D nr*nc ncases]); 
