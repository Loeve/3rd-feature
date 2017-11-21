function batchIndices = mkMinibatches(Ntrain, batchsize)
% batchIndices{b} has size batchsize, unless it is the last one

% This could probably be done more efficiently!
i = 1;
b = 1;
done = 0;
while ~done
  last = min(Ntrain, i+batchsize-1);
  batchIndices{b} = i:last;
  b = b + 1;
  i = i + batchsize;
  if i>Ntrain, done=1; end
end
