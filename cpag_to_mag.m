function M = cpag_to_mag(G,Mode)
  % constructs a particuler Maximal Ancestral Graph (MAG) M from the equivalence class 
  % represented by a given completed Partial Ancestral Graph (cPAG) G
  %
  % Input:
  % G = completed partial ancestral graph encoded in matrix form with 
  %   G(i,j) = 0  : no mark at i (no edge i-j) i     j
  %          = 1  : tail mark at i (to j)      i --* j    
  %          = 2  : arrowhead at i (to j)      i <-* j    (i not ancestor of j)
  %          = 3  : circle mark at i (to j)    i o-* j    (non-invariant)  
  % Mode = (optional) {1,2} select between (default) arc-augmentation (=1)
  %        and tail-augmentation (=2) procedure (see below).
  % Output:
  % M = (maximal) ancestral graph encoded in the form  
  %    (Gij, Gji) = (0,0)  : not adjacent    i     j
  %               = (1,1)  : undirected edge i --- j
  %               = (1,2)  : arrow           i --> j
  %               = (2,1)  : arrow           i <-- j
  %               = (2,2)  : bidirected edge i <-> j
  %
  % Description:
  % CPAGs encode the equivalence class of a (M)AG with every invariant edge
  % mark indicated by either arrowhead '>' (definite non-ancestor) or tail '-' 
  % (definite ancestor = causal link). Non-invariant edges are indicated by
  % a circle mark 'o' (= maybe ancestor ... maybe not) (see Zh/Sp,2005)
  % The default conversion is based on the 'arc augmentation' procedure
  % - orient all o-> edges into -->
  % - orient all --o edges into -->
  % - orient remaining subgraph over o-o edges into a DAG with no
  %   unshielded colliders (is chordal, so always possible), by:
  %   - choose node with highest number of processed parents and orient all edges as outward arcs
  %   - repeat with next unhandled, highest #parents until all have been done
  % As a result, no undirected or bi-directed edges are added.
  % See def.12, p.41 in (Zh/Sp,2005), and (Meek,95) for DAG orientation part
  
  % Alternative mode: 'tail-augmented' procedure (see Zhang's PhD thesis)
  % - orient all o-> edges into -->
  % - orient all --o edges into ---
  % - orient all o-o edges with no arrowheads attached into ---
  % - orient remaining subgraph over o-o edges into a DAG with no unsh.col

  % Initialise
  global DEBUG;
  if isempty(DEBUG), debug = 0; else debug = DEBUG; end;
  ARC  = 1;
  TAIL = 2;
  if (nargin < 2), Mode = ARC; end; % default arc augmentation
  n = length(G);
  
  % start with cpag
  M = G;
  % 1: orient all o-> edges into -->
  M(myintersect(find(G == 3), find(G' == 2) )) = 1;     % 
  
  % 2: orient all --o edges into -->
  if (Mode == ARC)
    M(myintersect(find(G == 3), find(G' == 1) )) = 2;     % 
  elseif (Mode == TAIL)
    M(myintersect(find(G == 3), find(G' == 1) )) = 1;     % 
    % 2b: all o-o edges with no arrowheads as ---
    AH = find(sum(G == 2,2) > 0);   % all nodes with 
    C = zeros(n);
    C(M == 3 & M' == 3) = 1;
    C(AH,AH) = 0;       % should be good as they are disjoint
    M(C == 1) = 1;
  end;
  
  % 3a: find circle component(s)
  C = zeros(n);
  C(myunion(find(M == 3), find(M' == 3) )) = 3;     % 
  Pc = C;               % (store copy of original circle component)
  [Cnodes,junk] = ind2sub(size(G),find(C));
  Cnodes = unique(Cnodes);
  Cc = C(Cnodes,Cnodes);               % (only circle component)
  Cunprocessed = Cnodes;
  Cprocessed = [];
  
  % 3b: loop over nodes in C
  while not(isempty(Cunprocessed))
    % find node with highest number of processed (parents) in C
    if isempty(Cprocessed)
      % first one: pick (first) node with highest degree in Cc
      node = Cnodes(find(sum(Cc,1) == max(sum(Cc,1)),1));
    else
      % find one with highest number of parents in Cprocessed
      % parents is M(i,j) == 1 equals i is parent of j 
      maxnumpar = max(sum(C(Cprocessed,Cunprocessed),1));
      node = Cunprocessed(find(sum(C(Cprocessed,Cunprocessed),1) == maxnumpar,1));
    end;
    
    % direct all undirected edges away from this node (inC)
    nodes_to = find(C(node,:) == 3);
    C(node,nodes_to) = 1;
    C(nodes_to,node) = 2;
    
    % and record as done
    Cprocessed   = [Cprocessed; node];
    Cunprocessed(find(Cunprocessed(:) == node)) = [];
  end;  % while
  
  % 4: transfer C to M and return
  M(find(Pc)) = C(find(Pc));      % or simply: M = M + C;  %!

  % done: return
  
end  % function cpag_to_mag