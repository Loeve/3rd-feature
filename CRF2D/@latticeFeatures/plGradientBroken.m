function [ll, g] = plGradientBroken(featureEng, weights)

% Compute log likelihood and gradient of pseudo likelihood
% for this training case
% We assume that x(i) is the label {-1,+1} 
% Node potentoials
% phi(xi,i) = exp [w'* NF(:,i) xi],  NF(:,i) are node features
% edge potentals
% psi(xi,xj,e) = exp [v'* EF(:,e) xi xj]  EF(:,e) are edge features
%
% Define
% sn(:,i) = sum_{j in nbrs i} x(j) EF(:,ij) = sum_neighbors
% snx(:,i) = sum_{j in nbrs i} x(i) x(j) EF(:,ij) = sum_neighbors_x_i
% snv(i) = sum_{j in nbrs i} x(j) v' * EF(:,ij) = sum_neighbors_v
%
% A(i) = w'*NF(:,i) + snv(i)
% esumA(i) = e^{Ai} + e^{-Ai}
% edifA(i) = e^{Ai} - e^{-Ai}}
%
% Then 
% ll(i) = log(exp [xi(w^T NF(:,i) + sum_j v^T xj EF(:,ij))] /
%               " with xi=1 + " with xi=-1 )
% = xi Ai - log(esumA(i))
%
% d ll(:,i)/dw = NF(:,i)(xi- edifA(i)/esumA(i))
% d ll(:,i)/dv = snx(:,i) - sn(:,i)*edifA(i)/esumA(i)

[w,v] = splitWeights(featureEng, weights);
nnodes = featureEng.nnodes;
edgeDirNum = featureEng.edgeDirNum;
outNbr = featureEng.outNbr;
Dnode = featureEng.Dnode;
Dedge = featureEng.Dedge;
nr = featureEng.nrows;
nc = featureEng.ncols;
X = reshape(featureEng.nodeLabels, nr, nc);
EF = featureEng.edgeFeatures;
NF = featureEng.nodeFeatures;

ll = 0;
grad_w = zeros(length(w),1);
grad_v = zeros(length(v),1);

sn = zeros(Dedge, nnodes);
snv = zeros(1, nnodes);
snx = zeros(Dedge, nnodes);
for r=1:nr
  for c=1:nc
    i = sub2ind([nr nc], r, c);
    for dir=1:4
      e  = edgeDirNum(r,c,dir);
      if e==0, continue; end
      r2 = outNbr(r,c,dir,1);
      c2 = outNbr(r,c,dir,2);
      sn(:,i) = sn(:,i) + X(r2,c2)*EF(:,e);
      snx(:,i) = snx(:,i) + X(r,c)*X(r2,c2)*EF(:,e);
      snv(i) = snv(i) + X(r2,c2)*v.'*EF(:,e);
    end
  end
  A(i) = w(:).'*NF(:,i) + snv(i);
  esumA(i) = exp(A(i)) + exp(-A(i));
  ediffA(i) = exp(A(i)) - exp(-A(i));
  
  ll = ll + A(i)*X(i) - log(esumA(i));
  grad_w = grad_w + X(i)*NF(:,i) - NF(:,i)*ediffA(i)/esumA(i);
  grad_v = grad_v + snx(:,i) - sn(:,i)*ediffA(i)/esumA(i);
end

g = [grad_w;grad_v];
