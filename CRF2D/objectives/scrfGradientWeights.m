function [err, grad] = scrfGradientWeights(weights, ndx, featureEng, infEng, data, lambda)

% This is the same as scrfGradient except it stores the weights and err at each iteration

[err, grad] = scrfGradient(weights, ndx, featureEng, infEng, data, lambda);

global FMINweightsTrace FMINerrTrace
i = size(FMINweightsTrace,2)+1;
fprintf('iter %d, err %5.3f, norm(grad) = %5.3f\n', i, err, norm(grad,inf));
if err < min(FMINerrTrace(:))
    FMINweightsTrace(:,i) = weights;
    FMINerrTrace(i) = err;
else
    FMINweightsTrace(:,i) = FMINweightsTrace(:,i-1);
    FMINerrTrace(i) = FMINerrTrace(i-1);
end
