function [G,lat_set,sel_set] = mxg_to_dag(M)
% MXG_TO_DAG Convert mixed graph (MAG/ADMG) into canonical dag representation
% Input:
% M = mixed graph encoded as 
%   [M(i,j), M(j,i)] = [0,0]  : not adjacent    i     j
%                    = [1,1]  : undirected edge i --- j
%                    = [1,2]  : arrow           i --> j
%                    = [2,1]  : arrow           i <-- j
%                    = [2,2]  : bidirected edge i <-> j
% Output:
% G = dag representation of M (mixed graph with only arcs)
% lat_set = set of latent nodes added to G      (i <-> j into i <-- k --> j)
% sel_set = set of selection nodes added to G   (i --- j into i --> s <-- j)
%
% Description:
% The canonical DAG representation is an ancestral graph that entails the same 
% independence relations between its nodes as G, but contains only directed 
% edges. As a result a number of implicit (hidden) nodes are added that act
% as latent or selection nodes (Rich/Spir,2002):
% - add selection child for every undirected edge to remove
% - add latent parent for every bidirected edge to remove
% Note: assumes M is valid acylic latent projection/ancestral graph; 
% all added nodes have indices > length(G) so that index in M = index in G.

% initialise
n = length(M);
cn = n;
lat_set = [];
sel_set = [];
edgeM = M + M'; % symmetrical matrix to characterise the edges

% create output matrix for G from M, including zero placeholders for possible extra nodes
n_add = sum(sum( triu(((edgeM == 2) + (edgeM == 4)),1) ));
if (n_add == 0)
  G = M;
  return;
else
  G = [M zeros(n,n_add); zeros(n_add, (n+n_add))];
end;

% loop over upper triangular part of M
for i = 1:(n-1)
  for j = (i+1):n
    if (edgeM(i,j) == 2)
      % undirected edge between i and j
      % add selection child node cn with arrows from i and j
      cn = cn+1;
      G([i j], cn) = 1;         % tails from i and j
      G(cn, [i j]) = 2;         % arrowheads at new node
      G(i,j) = 0; G(j,i) = 0;   % remove link i-j
      sel_set = [sel_set cn];
    elseif (edgeM(i,j) == 4)
      % bi-directed edge between i and j
      % add latent parent node cn with arrows to i and j
      cn = cn+1;
      G([i j], cn) = 2;         % arrowheads at i and j
      G(cn, [i j]) = 1;         % tails from new node
      G(i,j) = 0; G(j,i) = 0;   % remove link i-j
      lat_set = [lat_set cn];
    end;  % if 
  end;  % for j
end;  % for i

% done: return canonical dag G

end  % function mxg_to_dag
