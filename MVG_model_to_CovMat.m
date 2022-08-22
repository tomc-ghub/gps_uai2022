function  Sigma = MVG_model_to_CovMat(B,V)
% Convert MVG model X = B*X + E (+ C), with E ~ N(.|0,V), and const. C (=0)
% into covariance matrix Sigma of the corresponding MVG pdf as:
% => Sigma = (I-B)^(-1)*V*(I-B)^(-T), with V a diagonal matrix
% 
% input:
% - B     = matrix of weights in MVG model x_i = B(i,:)*X + v_i
% - V     = independent variance (noise term) as vector or diag.cov.matrix
%           
% output:
% - Sigma = covariance matrix corresponding to MVG model 
% =========================================================================

  % 1 - initialize
  N = size(B,1); 
  % make sure V is diagonal matrix
  if ~isdiag(V), V = diag(V); end;

  % compute covariance matrix
  invI_B = (eye(N) - B) \ eye(N);  % more stable than inv(I-B)
  Sigma = invI_B*V*invI_B';       

end  % function mvg_model_to_covmat
