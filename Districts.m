function D = Districts(G)
% compute matrix with districts / c-components (+ parents) of G 
% input:
% - G     : NxN graph (assumes edges symmetric nonzero entries)
% output:
% - D     : NxnD matrix with D(i,j)==1 : node i present in district j
%                            D(i,j)==2 : node i parent of district j
% 
% example (R4b for discr.path x->z<->w<->y for w)
% .    x,z,w,y
% G = [0,0,1,0;
% .    0,0,2,2;
% .    2,2,0,1;
% .    0,2,2,0];
% then 2 districts: 1={x}, 2={z,w,y}, with parent x
% D = [1,2;
% .    0,1;
% .    0,1;
% .    0,1];
DEBUG = ~true;

  % initialize
  N  = size(G,1);           % nr. of variables
  D  = zeros(N,N);          % output D(node,district), max.N components preallocated
  nD = 0;                   % nr. of districts
  visited = zeros(1,N);     % visit x before?
  root = 1;
  % start at next root node, expand district, and repeat 
  while ~isempty(root)
    % next district
    nD = nD + 1;
    district = zeros(1,N);
    district(root) = 1;
    % get all neighbours connected to 'root' with bidirected edge in G
    new_neighbours = find(G(:,root) == 2 & G(root,:)' == 2);
    while ~isempty(new_neighbours)
      % expand district
      district(new_neighbours) = 1;
      % find new neighbours (of nodes in 'new_neighbours') not already in district
      % note: no need to check for '~visited' in previous districts
      new_neighbours = find(sum(G(:,new_neighbours) == 2 & G(new_neighbours,:)' == 2 & ~district',2));
    end; % while (new_neighbours)
    
    % update D + visited
    nodes = find(district);
    D(nodes,nD) = 1;
    visited(nodes) = 1;
    
    % collect parents of district in G (but not part of the district itself)
    parents = find(sum(G(:,nodes) == 1 & G(nodes,:)' == 2 & ~district',2));
    % store in D as well
    D(parents,nD) = 2;
    
    % find next root node to process (if any)
    root = find(~visited,1);
  end; % while (root)
  
  % clean up and return
  D(:,[nD+1:end]) = [];
end  % function Districts
