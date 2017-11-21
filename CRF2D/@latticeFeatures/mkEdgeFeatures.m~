function  edgeFeatures = mkEdgeFeatures(featureEng, raw);

% raw(:,r,c,s)
% EF(:,e,s) = [1, |raw(:,r,c) - raw(:,r',c')| ] if expandEdge = 0
% EF(:,e,s) = [1, Quadratic(|raw(:,r,c) - raw(:,r',c')|) ] if expandEdge = 1

[Draw nr nc ncases] = size(raw);
expand = featureEng.expandEdge;
D = length(getMu(raw(:,1,1,1), raw(:,1,1,1), expand));
    
[out_edge, outNbr, edgeEndsIJ, edgeEndsRC, in_edge, nedges] = ...
    assign_edge_nums_lattice(nr, nc);

nedgesDir = size(edgeEndsIJ,1);
edgeFeatures = zeros(D, nedgesDir, ncases);
for s=1:ncases
  for e=1:nedgesDir 
    r1 = edgeEndsRC(e,1); c1 = edgeEndsRC(e,2); r2 = edgeEndsRC(e,3); c2 = edgeEndsRC(e,4);
    mu = getMu(raw(:,r1,c1,s), raw(:,r2,c2,s), expand);
    edgeFeatures(:,e,s) = mu;
  end
end
