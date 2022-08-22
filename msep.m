function sep = msep(X,Y,Z,M,method)
% MSEP Is X m-separated (independent) from Y given Z in ancestral graph M?
% sep = msep(X, Y, Z, G)  
%
% input:
% - x,y   : nodes to separate
% - Z
% - M     : true ancestral graph (not necessarily maximal)
% - method  =  (in tests Bayes Ball seemed faster)
%   - 1 = (exploding) Bayes Ball, 
%   - 2 = reachable in moralized ancestral graph via canonical DAG
% output:
% - G     : causal graph (partially oriented skeleton / PAG)
% - nCI   : array of nr. of CI tests over i variables
% - CIL   : list of established (minimal) independencies

% Note on graph encoding:
% - M = NxN ancestral graph encoded as M(i,j) =
%       0      : not adjacent    i     j
%       1      : tail            i --* j    (with * = {-,>,o}
%       2      : arrowhead       i <-* j

  % Initialize
  if nargin < 5, method = 1; end;     % default exploding Bayes-ball?
  N = length(M);
  % no need to search for separating sets of adjacent nodes in M
  if find(M(X,Y) ~= 0,1), sep = false; return; end;
  % temporary fix for using a cPAG in call to msep
  % => if not sufficient then use 'cpag_to_mag' first!
  M(M == 3) = 1;

  % assume separated until connecting path is found
  sep = true;
  if (method == 1)
    % (exploding) Bayes Ball: like regular Bayes Ball, except at each step
    % all possible next steps are followed simultaneously. Start at X and walk all
    % untraversed edge(halves) from nodes reached, crossing nodes not in Z via non-collider 
    % paths and nodes in Z via collider paths; terminate at a node if no available 
    % untraversed edges available. If Y is reached msep is false (so X and Y
    % dependent given Z), if no more available edges Y cannot be reached and
    % msep = true (so X and Y independent given Z).

    % construct the graph to update
    G = M;
    % all outgoing arcs with no arrowhead at nodes in Z are eliminated from G
    for z = Z
      Vtmp = find(M(z,:) == 1);
      G(z,Vtmp) = 0;
      G(Vtmp,z) = 0;
    end;  % for s

    % start with all X as noncollider-endpoint
    Vhead = [];
    Vtail = X;
    % while new found (but not yet y) keep trying
    while ~isempty([Vhead,Vtail])
      % reset collector variables for next nodes
      Vnewh = [];
      Vnewt = [];
      % loop over into-nodes 
      for v = Vhead
        % find all nodes (still) reachable from v in G
        % process depending on whether v is selection
        if ~isempty(find(Z == v,1)), 
          % v is selection node: (only arrowheads at v in G)
          v_to = find(G(v,:) ~= 0);
          % check M for edge marks
          v_head = v_to(M(v_to,v) == 2);
          v_tail = v_to(M(v_to,v) == 1);  % or complement
        else
          % v is not selection node .. only tails from v
          v_to = find(G(v,:) == 1);
          % check M for edge marks
          v_head = v_to(M(v_to,v) == 2);
          v_tail = v_to(M(v_to,v) == 1);  % or complement
        end;
        % set traversed edge-half to zero in G and add to nodes found
        G(v,v_to) = 0;
        Vnewh = [Vnewh,v_head];
        Vnewt = [Vnewt,v_tail];
      end;  % for Vhead

      % loop over out-of-nodes 
      for v = Vtail
        % find all nodes (still) reachable from v in G
        if ~isempty(find(Z == v,1)), 
          % non-collider always blocked by selection node
          continue;
        else
          % v is not selection node .. all edge-halves traversable
          v_to = find(G(v,:) ~= 0);
          % check M for edge marks
          v_head = v_to(M(v_to,v) == 2);
          v_tail = v_to(M(v_to,v) == 1);  % or complement
        end;
        % set traversed edge-half to zero in G and add to nodes found
        G(v,v_to) = 0;
        Vnewh = [Vnewh,v_head];
        Vnewt = [Vnewt,v_tail];
      end;  % for Vtail

      % gather next nodes reached at heads/tails
      Vhead = myunique2(Vnewh);
      Vtail = myunique2(Vnewt);    

      % check if we reached y (or set Y) in this step
      if ~isempty(myintersect2([Vhead,Vtail], Y)),
        sep = false;
        return;
      end;
    end; % while
    % if here then not connected
    sep = true;
    return;

  elseif (method == 2)
    % Convert to canonical DAG, take subgraph over nodes that can reach X, Y 
    % and/or Z (incl. self), convert into moralized undirected graph and see if 
    % removing all connections to/from Z leaves X and Y unconnected.
    %
    % compute canonical ancestral graph from G (DAG)
    [R,lat_set,sel_set] = ag_to_cag(M);
    % convert to adjacency notation (0/1) by removing all arrowheads (2 => 0)
    R = ag_to_adj(R);
    % get ancestral part of dag
    C = reachability_graph(R);
    Z = myunion(Z, sel_set);      % all nodes in sel_set are implicitly selected
    G = myunion([X,Y], Z);
    [A,junk] = find(C(:, G));     % nodes that can reach nodes M are relevant
    A = unique(A);                % (but once is enough)
    A = myunion(A, G);            % in combination with the nodes from M itself
    % now construct the moralized ancestral graph
    GM = moralize(R(A,A));
    % and see whether X and Y are separated by Z in this undirected graph
    sep = graph_separated(GM, find_equiv_posns(X,A), find_equiv_posns(Y,A), find_equiv_posns(Z,A));
    return;
  end;  % if method==2

end  % function msep


% =========================================================================
% local functions to optimize speed 
% TODO: check for further improvmenets.
% =========================================================================
function B = myunique2(A)
% MYUNIQUE Sorts unique elements in set (array) A of positive integers (much faster than built-in unique)
% B = myunique(A)
  if isempty(A)
    B = [];
  else
    bits = zeros(1, max(A));
    bits(A) = 1;
    B = find(bits);
  end
end % myunique2

function C = myintersect2(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = myintersect(A,B)
  if isempty(A) || isempty(B), 
    C = []; 
  else
    bits = zeros(1, max([A,B]));
    bits(A) = 1;
    C = B(logical(bits(B)));  
  end;
end  % myintersect2



