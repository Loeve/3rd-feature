function [w, err, errTrace,  weightsTrace, errTrace_iter1, weightsTrace_iter1,etaTrace, vTrace] = smd(gradient, winit, trainNdx, varargin)
% Stochastic meta descent
%
% Arguments are almost identical to stochgrad, except we can also specify
%
% lambda
% mu

[lambda, mu, eta0, maxIter, batch_size, gradArgs,  ...
    displayFn, displayArgs] = process_options(...
    varargin, 'lambda', 0.99, 'mu', 0.02, 'eta0', 0.0005, ...
    'maxIter', [], 'batch_size', 5, 'gradArgs', {}, ...
    'displayFn', [], 'displayArgs', {});

maxIter=maxIter+1;
batchIndices = mkMinibatches(length(trainNdx), batch_size);
Nbatches = length(batchIndices);

ii = 1e-150*sqrt(-1);
p = length(winit);
T = maxIter;
v = zeros(p,1); % auxiliary vector
eta = eta0*ones(p,1); % learning rate per parameter
w = winit;

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
if nargout >= 7
    etaTrace(:,1) = eta;
else
    etaTrace = [];
end
if nargout >= 8
    vTrace(:,1) = v;
else
    vTrace = [];
end

for i = 2:T
    for b=1:Nbatches
        batchNdx = batchIndices{b};
        % Nic's code - uses complex number trick
        [f(b),g] = feval(gradient, w + ii*v, batchNdx, gradArgs{:});
        eta = eta.*max(1/2,1+mu*v.*real(g));
        w = w - eta.*real(g);
        v = lambda*v+eta.*(real(g)-lambda*imag(g)*1e150);
        
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
    if ~isempty(displayFn), feval(displayFn, i, b, w, real(err), g, displayArgs{:}); end
    if ~isempty(weightsTrace)
        weightsTrace(:,i) = w;
    end
    if ~isempty(errTrace)
        errTrace(i) = err;
    end
    if ~isempty(etaTrace)
        etaTrace(:,i) = eta;
    end
    if ~isempty(vTrace)
        vTrace(:,i) = v;
    end
end
if nargout >= 2
    [err] = feval(gradient, w, trainNdx, gradArgs{:});
end
