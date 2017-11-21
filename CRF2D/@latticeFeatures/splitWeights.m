function [nodeWeights, edgeWeights] = splitWeights(featureEng, weights);

ndx = featureEng.Dnode;
nodeWeights = weights(1:ndx);
edgeWeights = weights(ndx+1:end);
