function demoMulti()

% Compare batch and stochastic training methods
% using BP as inference engine

% Load Training/Testing Data
trainNdx = 1:10;
testNdx = 1:5;
%trainNdx = 1:50;
%testNdx = 1:10;
label = sign(double(imread('X.png'))-1);
label=label(:,:,1);
nstates = 2;

for i = trainNdx
    train(:,:,i) = label+randn(32,32);
end
for i = testNdx
    test(:,:,i) = label+randn(32,32);
end

figure;
for i = testNdx
    subplot(2,5,i);
    imshow(test(:,:,i))
end
a = test(:,:, 1)
whos a
suptitle('test data')

% Make Features and Feature Engine
featureEng = latticeFeatures(0,0);

trainFeatures = permute(train,[4 1 2 3]);
whos trainFeatures
traindata.nodeFeatures = mkNodeFeatures(featureEng,trainFeatures);
traindata.edgeFeatures = mkEdgeFeatures(featureEng,trainFeatures);
training_labels = repmat(label,[1 1 length(trainNdx)]);
whos training_labels
traindata.nodeLabels = training_labels;
traindata.ncases = length(trainNdx);
trainNdx = 1:traindata.ncases;
nNodeFeatures = size(traindata.nodeFeatures,1);
nEdgeFeatures = size(traindata.edgeFeatures,1);
winit = initWeights(featureEng,nNodeFeatures,nEdgeFeatures);

testFeatures = permute(test,[4 1 2 3]);
testdata.nodeFeatures = mkNodeFeatures(featureEng,testFeatures);
testdata.edgeFeatures = mkEdgeFeatures(featureEng,testFeatures);
testdata.nodeLabels = repmat(label,[1 1 length(trainNdx)]);
testdata.ncases = length(testNdx);
testNdx = 1:testdata.ncases;

%%%%%%%% Random params

infEng = latticeInferBP(nstates);
showResults(winit, testNdx, featureEng, infEng, testdata, 'rnd params');

%%%%%%%%% BFGS training with BP

reg = 1;
maxIter = 3;
options = optimset('Display','iter','Diagnostics','off','GradObj','on',...
    'LargeScale','off','MaxFunEval',maxIter);

gradFunc = @scrfGradient;
gradArgs = {featureEng, infEng, traindata, reg};

weights = fminunc(gradFunc,winit,options,trainNdx,gradArgs{:});

showResults(weights, testNdx, featureEng, infEng, testdata, 'BFGS+BP');


%%%%%%%%%%% SG trainign with BP

reg = 1;
maxIter = 1;
eta = 0.0001;
batch_size = 1;
anneal = 0;
tau = 0;

weights = stochgrad(gradFunc,winit,trainNdx,'gradArgs',gradArgs,...
    'maxIter',maxIter,'eta',eta,'batch_size',batch_size);

showResults(weights, testNdx, featureEng, infEng, testdata, 'SG+BP');


%%%%%%%%%%%% SMD with BP

reg = 1;
maxIter = 1;
eta0 = 0.0001;
mu = 0.001;
lambda = 0.9;
batch_size = 1;

weights = smd(gradFunc,winit,trainNdx,'gradArgs',gradArgs,...
    'maxIter',maxIter,'eta0',eta0,'batch_size',batch_size,...
    'lambda', lambda, 'mu', mu);

showResults(weights, testNdx, featureEng, infEng, testdata, 'SMD+BP');


%%%%%%%%%%%%

function showResults(weights, testNdx, featureEng, infEng, testdata, ttl)

figure;
for i = testNdx
    subplot(2,5,i);
    featureEng = enterEvidence(featureEng, testdata, i);
    [nodePot, edgePot] = mkPotentials(featureEng, weights);
    [nodeBel, MAPlabels] = infer(infEng, nodePot, edgePot);
    imshow(MAPlabels);
    title(sprintf('%d',i));
end
testErrorRate = classifPerformance(weights, testNdx, featureEng, infEng, testdata);
suptitle(sprintf('%s, error rate = %5.3f', ttl, testErrorRate))
drawnow
