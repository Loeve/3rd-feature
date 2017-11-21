function [w, err, errTrace, weightsTrace,errTrace_iter1,weightsTrace_iter1] = stochgrad(gradient, winit, trainNdx, varargin)
% Stochastic gradient descent
%
% gradient is a fn of the form [f,g]=gradient(w, ndx, gradArgs{:})
%   where ndx specifies the indices in the current mini batch
% winit is the initial weight vector
% trainNdx is a vector of integers which specify the complete
%  training set  (usually 1:N); this is used for deriving mini-batches
%  (otherwise we don't know how big the dataset is)
% varargin supports various optional arguments, in particular:
%
% eta - (initial) step size
% anneal - set to 1 to turn on annealing
% tau - annealing parameter
% maxIter
% batch_size
% gradArgs - a cell array passed to gradient (default {})
% displayFn - a fn of the form displayFn(iter, weights, f, g, displayArgs{:})
% displayArgs - a cell array passed to displayFn (default {})
%
% EXAMPLE
%[weights, err, nllSG_train, SGweights] = stochgrad(@scrfGradient, winit, trainNdx, 'gradArgs', gradArgs, ...
%    'maxIter', T+1, 'eta', eta0, 'batch_size', batch_size, 'displayFn', @displayIter);
% where
% gradArgs = {featureEng, infEng, traindata, reg};
% [err, grad] = scrfGradient(weights, ndx, featureEng, infEng, data, lambda)
%
% ***** NOTE: Returning errTrace is expensive, so this code will run much
% faster if you don't need it.

[eta0, maxIter, batch_size, gradArgs,  displayFn, displayArgs,anneal,tau] = process_options(...
    varargin, 'eta', 0.0005, ...
    'maxIter', [], 'batch_size', 5, 'gradArgs', {}, ...
    'displayFn', [], 'displayArgs', {}, 'anneal',0,'tau',0);

maxIter = maxIter+1;
batchIndices = mkMinibatches(length(trainNdx), batch_size);
Nbatches = length(batchIndices);
p = length(winit);
T = maxIter;
w = winit; % Initialize weights

if nargout >= 3
    % let's see how good the initial parameter guess is...
    [err] = feval(gradient, w, trainNdx, gradArgs{:});
    errTrace(1) = err;
else
    errTrace = [];
end

if nargout >= 4
    weightsTrace(:,1) = w;
else
    weightsTrace = [];
end

if nargout >= 5
    errTrace_iter1(1) = err;
else
    errTrace_iter1 = [];
end

if nargout >= 6
    weightsTrace_iter1(:,1) = w;
else
    weightsTrace_iter1 = [];
end

for i = 2:T
    for b=1:Nbatches
        
        % Update learning rate
        if anneal
            eta = eta0/(tau+i-1);
        else
            eta = eta0;
        end
        
        % Compute and take the step
        batchNdx = batchIndices{b};
        [f(b),g] = feval(gradient, w, batchNdx, gradArgs{:});
        w = w - eta*g;
        
        fprintf('Iter = %d, Batch = %d, f(b) = %d\n',i-1,b,f(b));

        if ~isempty(errTrace_iter1) && i == 2
            errTrace_iter1(b+1)=feval(gradient, w, trainNdx, gradArgs{:});
            if ~isempty(weightsTrace_iter1)
                weightsTrace_iter1(:,b+1)=w;
            end
        end
    end
    if ~isempty(displayFn) || ~isempty(errTrace)
        err = feval(gradient, w, trainNdx, gradArgs{:});
    end
    if ~isempty(displayFn), feval(displayFn, i, b, w, err, g, displayArgs{:}); end
    if ~isempty(weightsTrace)
        weightsTrace(:,i) = w;
    end
    if ~isempty(errTrace)
        errTrace(i) = err;
    end
end

if nargout >= 2
    [err] = feval(gradient, w, trainNdx, gradArgs{:});
end