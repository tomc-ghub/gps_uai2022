function C = reachability_graph(G,self)
  % REACHABILITY_GRAPH C(i,j) = 1 iff there is a path from i to j in DAG G
  % input:
  % - G: MAG/PAG G (circle marks are treated as potential tails)
  % - self: (= 0/1) treat node x as reachable to itself (default: no)
  % output:
  % - C (binary) reachability graph, 
  %   C(i,j) = 1 iff (potential) directed path from i to j in G
  if (nargin < 2), self = 0; end;   % default exclude self
  % based on BayesNetToolbox\FullBNT-1.0.4\graph\reachability_graph.m
  % tipping point around n = 14
  n = length(G);
  % for MAG/PAG: remove arrowheads (== 2) and set circle (==3) to tail
  G(G == 2) = 0;
  G(G == 3) = 1;
  % check possible directed paths
  if (n > 14)
    % expm(G) = I + G + G^2 / 2! + G^3 / 3! + ...
    if (self == 0),
      M = expm(double(full(G))) - eye(n); % exclude self (only proper ancestors)
    else
      M = expm(double(full(G))); % do not exclude self : x \in An(x)
    end;
    C = (M>0);
  else
    % This computes C = G + G^2 + ... + G^{n-1} but include self for acyclic
    A = G;
    if (self == 0),
      C = zeros(n);     % exclude self (default)
    else
      C = eye(n);       % include self
    end;
    for i=1:n
      C = C + A;
      A = A * G;
    end
    C = (C > 0);
  end
end % reachability_graph
