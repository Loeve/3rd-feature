function [errRate, nerr] = classifPerformance(weights, ndx, featureEng, infEng, data)
% compute the fraction of misclassified pixesl

Ncases = length(ndx);
for s=1:Ncases
  casenum = ndx(s);
  featureEng = enterEvidence(featureEng, data, casenum);
  [nodePot, edgePot] = mkPotentials(featureEng, weights);
  [nodeBel, MAPlabels] = infer(infEng, nodePot, edgePot);
  npixels(s) = length(MAPlabels(:));
  nerr(s) = compareLabels(featureEng, MAPlabels);
end
errRate = sum(nerr)/sum(npixels);

