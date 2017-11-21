function w = initWeights(featureEng, Dnode, Dedge)

NW = 0.1*randn(Dnode, 1);
EW = 0.1*randn(Dedge, 1);
w = [NW(:); EW(:)];
