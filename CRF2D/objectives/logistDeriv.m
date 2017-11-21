function [f,g,H] = logistDeriv(w,X,y,lambda)
% Value, Gradient and Hessian of logisitc fucntion sigma(w^T x)
% w = weight vector
% X(i,:)
% y(i) = 0 or 1

if nargin < 4, lambda = 0; end

[n k] = size(X);

p = 1 ./ (1 + exp(-X*w));

f = -sum( (y.*log(p+eps) + (1-y).*log(1-p+eps)) );

g = -X'*((y-p));

wt = p .* (1-p);
H = zeros(k,k);
for i = 1:k,
  wxi = wt .* X(:,i);
  for j = i:k,
    hij = wxi' * X(:,j);
    H(i,j) = -hij;
    H(j,i) = -hij;
  end
end
H=-H;

% penalized ll = ll - lambda/2 ||w||^2 (since p(w) propto exp(-lam/2 ||w||^2)
% penalized f = f + lambda/2 ||w||^2
% penalized g = g + lambda w
% penalized H = H + diag(lambda)
f = f + lambda/2*sum(w.^2);
g = g + lambda*w;
H = H + eye(k)*lambda;
