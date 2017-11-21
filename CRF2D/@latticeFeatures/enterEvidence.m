function featureEng = enterEvidence(featureEng, data, s)
% data.nodeFeatures(d,i,s)
% data.edgeFeatures(d,e,s)
% data.nodeLabels(r,c,s) : infer size of image from this

[nr nc ncases] = size(data.nodeLabels);

featureEng.nodeLabels = reshape(data.nodeLabels(:,:,s), nr*nc, 1);
featureEng.nodeFeatures = data.nodeFeatures(:,:,s);
featureEng.edgeFeatures = data.edgeFeatures(:,:,s);

[featureEng.edgeDirNum, featureEng.outNbr, featureEng.edgeEndsIJ, featureEng.edgeEndsRC]  = ...
    assign_edge_nums_lattice(nr, nc);

featureEng.nrows = nr;
featureEng.ncols = nc;
featureEng.nnodes = nr*nc;
%featureEng.nedges = size(featureEng.edgeEndsRC,1);
featureEng.Dnode = size(data.nodeFeatures,1);
featureEng.Dedge = size(data.edgeFeatures,1);
