function [Sigma,Mu,B,V] = ag_to_random_MVG_PDF(M,mode)
% ag_to_random_MVG generates a random linear Gaussian model corresponding
% to the ancestral graph M with accompanying covariance matrix Sigma and
% zero-mean vector Mu. 
% Model: X = B*X + E with noise E ~ N(.|0,V), B(j,i)~=0 iff i->j in G
% Generate corresponding pdf (cov.matrix Sigma) via canonical DAG, where each
% bidirected edge x<->y becomes an (unobserved) confounder x<-z->y, and
% each undirected edge x---y becomes an (unobserved) conditioned x->z<-y
% 
% input: 
% - M     = mixed graph representing (M)AG structure 
% - mode
%   .b_range  = [lower,upper] bound on (nonzero) edge parameters B(i,j)
%   .s_range  = [lower,upper] bound on std.dev. of intrinsic noise term e_i
%   .minabs_b = minimum abs.size non-zero edge strengths in B (default 0)
%           
% output:
% - Sigma = random covariance matrix corresponding to ancestal graph M (?)
% - Mu    = (column) of means (default just zero)
% - B     = matrix of edge weights in MVG model of G: x_i = B_G(i,:) * X + e_i
% - V     = idem, vector of variance of independent noise terms e_i ~ N(.|0,v_i)
% =========================================================================
  
  % 1 - initialize  
  N = size(M,1);
  % get parameter properties
  if (nargin < 2), mode.default = 1; end;
  % get range of edge parameter properties
  if isfield(mode,'b_range'), 
    min_b = mode.b_range(1); max_b = mode.b_range(2); 
  else
    min_b = -1.0; max_b = +1.0; % default edge range [-1.0,+1.0]
  end;
  if isfield(mode,'minabs_b'), 
    minabs_b = mode.minabs_b; 
  else
    minabs_b = 0;   % default no min.par. |b| >= 0
  end;
  % range of std.dev. of intrinsic noise per variable in (canonical) DAG
  if isfield(mode,'s_range'), 
    min_s = mode.s_range(1); max_s = mode.s_range(2); 
  else
    min_s = 0.01; max_s = 0.2; % default range of std.dev. of intrinsic noise [0.01,0.2]
  end;
      
  % ===========================================
  % 2: convert (M)AG to canonical DAG, incl. 'padded' latent/selection variables
  [G,Lat,Sel] = mxg_to_dag(M); 
  
  % 3: add parameters to canonical DAG
  S = G';   % transposed to match parameter encoding
  S(S > 1) = 0;  % S(i,j) = 1 <-> edge j--*i in canonical DAG G
  nG = size(G,1);
  % get edge parameter matrix (G)B for canonical DAG G, over non-zero
  % adjacencies (parents) in S)
  % sampled uniform from interval {[min_b,-minabs_b],[+minabs_b,max_b]}
  B = S.*( (max_b - min_b - 2*minabs_b)*rand(nG) + min_b + minabs_b );
  B = B + sign(B)*minabs_b;
  % sample std.dev. of intrinsic noise terms on variables in can.DAG G
  V = (max_s - min_s)*rand(nG,1) + min_s;
  V = V.^2; % convert to variance 
  
  % NOTE: parameters for bidirected/selection edges could take sqrt() to
  % match strength of other edges (product of 2 edges), but not now.
  
  % convert to covariance matrix over can. DAG G (ignore mean)
  C = zeros(nG,1);      % dummy means
  Sigma_G = MVG_model_to_CovMat(B,V);
  % impose selection bias (if any)
  if ~isempty(Sel),
    % implicit conditioning
    order = [1:N, Lat, Sel];
    nObsLat = N + length(Lat);
    % reorder covariance matrix for partitioning ([Obs,Lat,Sel]
    Sigma_G = Sigma_G(order,order);     % should not be needed?
    % get precision matrix
    Lambda_G = Sigma_G \ eye(nG);  % better than Lambda_G  = inv(Sigma_G);   
    % marginalize selection variables and convert back to covariance 
    Lambda_OL = Lambda_G(1:nObsLat,1:nObsLat);  % 
    Sigma_G = Lambda_OL \ eye(nObsLat); % = inv(Lambda_OL);
  end;
  
  % restrict to observed variables (N in MAG, not nG in can.DAG)
  Mu    = zeros(N,1);     % could generate random means, but ignored anyway
  Sigma = Sigma_G([1:N],[1:N]);
  % no explicit positive definite check
  
end   % function ag_to_random_MVG_PDF

