function [err, grad] = pseudoLikGradient(weights, ndx, featureEng, data, lambda)
% Input
% weights - interpreted by featureEng
% ndx - specifies which set of data to use (for mini-batch training)
% featureEngine
% data - interpreted by featureEng
% lambda - quadratic regularizer on weights
%
% Output
% err = *negative* log likelihood
% grad(i) = sum_{s in ndx} d err(s)/d w(i)


NcasesTotal = data.ncases;
NcasesBatch = length(ndx); 
Nfeatures = length(weights);
gradient = zeros(Nfeatures, NcasesBatch);

for s=1:NcasesBatch
  casenum = ndx(s);
  featureEng = enterEvidence(featureEng, data, casenum);
  [ll(s), gradient(:,s)] = plGradient(featureEng, weights);
end

% compute final values + regularizer
% err = -log-likelihood + lambda*(batch_size/train_set_size)*sum(w(i)^2)
% grad(i) = -dll/dw(i) - 2*lambda*(batch_size/train_set_size*w(i)

err = -sum(ll) + ((lambda/2)*sum(weights.^2))*NcasesBatch/NcasesTotal;
grad = -sum(gradient,2) + (lambda*weights)*NcasesBatch/NcasesTotal;
