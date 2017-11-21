function err = compareLabels(featureEng, MAPlabels)
% Compare MAPlabel estimates with the true labels for current case (after enterEvidence)
% We just sum the number of misclassified nodes.

err = sum(MAPlabels(:) ~= featureEng.nodeLabels(:));
