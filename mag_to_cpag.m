function [P,nRules] = mag_to_cpag(G)
  % calculates the completed Partial Ancestral Graph representation for a given 
  % maximal ancestral graph.
  % input:
  % - G = (maximal) ancestral graph encoded in the form  
  %    (Gij, Gji) = (0,0)  : not adjacent    i     j
  %               = (1,1)  : undirected edge i --- j
  %               = (1,2)  : arrow           i --> j
  %               = (2,1)  : arrow           i <-- j
  %               = (2,2)  : bidirected edge i <-> j
  % output:
  % - P = completed ancestral graph encoded in matrix form with 
  %   P(i,j) = 0  : no mark at i (no edge i-j) i     j
  %          = 1  : tail mark at i (to j)      i --* j    (causal link i=>j? NO)
  %          = 2  : arrowhead at i (to j)      i <-* j    (i not ancestor of j)
  %          = 3  : circle mark at i (to j)    i o-* j    (non-invariant)  
  % - nRules = array with counters of orientation rules triggered (FCI numbering, see mag_to_cpag)
  %
  % Description:
  % Follows the augmented FCI orientation strategy from [Zhang,2008]. Identical to 
  % mec_to_cpag except R0b (v-structures) and R4 (discr.path).
  % Note the implementation below also includes the more efficient way to
  % check for discriminated nodes (without explicit reconstruction of the
  % discr.path, see R4' in the GPS article), in combination with a more
  % efficient implementation for the selection bias rules that now combines
  % R5-R7 per node in one step. Same approach is also used to improve speed
  % for checking R9(/R10).

  % =======================================================================
  % 1 - Initialise
  % =======================================================================
  global DEBUG;
  if isempty(DEBUG), debug = 0; else debug = DEBUG;  end;
  % note: set debug>1 for detailed rule tracing in output during orientation
  % debug = 2;
  
  % start from P as skeleton with o-o edges
  P = 3*(G > 0);
  nRules = zeros(1,14);  % function global counter of rules R0b-R10 that triggered
  
  % =======================================================================
  % 2 - Run blocks of orientation rules
  % =======================================================================
  % Flow: repeat each step until no more change
  % G -> R0b -> (R1-R4b) -> R5 -> (R6-R7) -> (R8a-R10) -> P
  if (debug > 1), figure(3); clf; draw_cpmag(P);  title(sprintf('mag_to_cpag: PAG - stage 1')); end;
  P = FCI_Orientation_R0b(P,G);       % add v-structures
  
  if (debug > 1), figure(3); clf; draw_cpmag(P);  title(sprintf('mag_to_cpag: PAG - stage 2')); end;
  P = FCI_Orientation_R1_R4b(P,G);    % add remaining arrowheads (+ some tails)
  
  if (debug > 1), figure(3); clf; draw_cpmag(P);  title(sprintf('mag_to_cpag: PAG - stage 3')); end;
  P = FCI_Orientation_R5_R7(P);         % do selection bias (undirected edges)
  
  if (debug > 1), figure(3); clf; draw_cpmag(P);  title(sprintf('mag_to_cpag: PAG - stage 4')); end;
  P = FCI_Orientation_R8a_R10(P);       % add remaining tails

  % PAG P completed
  if (debug > 1), figure(3); clf; draw_cpmag(P);  title(sprintf('mag_to_cpag: PAG - stage 5')); end;
  return;
  
% =======================================================================
% Subfunctions
% =======================================================================
% =======================================================================
% rule R0b
% =======================================================================
function P = FCI_Orientation_R0b(P,G)
%function P = FCI_Orientation_R0b(P)
  % 1: initialize
  N = size(P,1); 
  % loop over all nodes z
  for z = 1:N
    S = find(G(z,:) == 2);       % all nodes with arrowheads to z in G
    nS = length(S);
    for i = 1:nS-1
      x = S(i);
      for j = (i+1):nS
        y = S(j);
        if (G(x,y) == 0),       % for each nonadjacent pair {x,y} ...
          % implies z not in sepset(x,y), so orient v-structure x *-> z <-* y
          P(z,x) = 2;
          P(z,y) = 2;
          if (debug > 1), fprintf('FCI: R0b orients  %d *-> %d <-* %d \n',x,z,y); end;
          nRules(1) = nRules(1) + 1;  % R0b = 1
        end; 
      end;  % for y
    end;  % for x
  end;  % for z
end  % function FCI_Orientation_R0b
  
% =======================================================================
% Orientation rules R1-R4b
% =======================================================================
function P = FCI_Orientation_R1_R4b(P,G)
  % 1: initialize
  N = size(P,1);  
  % start from all nodes in P with arrowheads (from R0)
  check_nodes = find(sum(P == 2,2)' > 0);

  % 2: double loop until no more changes
  while ~isempty(check_nodes)
    % 2a: repeatedly process rules 1-3:
    % add arrowheads (and some tails) until no more can be found
    while ~isempty(check_nodes)    
      % keep track of new arrownodes 
      new_arrow = zeros(1,N);

      % ============================
      % R1: (near) y-structure
      % if (x *-> z o-* y) and x,y not adjacent then orient x *-> z --> y
      for z = check_nodes(:)'
        % loop over nodes with arrowheads at z (at least check_nodes)
        X = find(P(z,:) == 2);
        for x = X(:)'
          % find nodes with circle marks at z and no edge to x
          Y = find(P(z,:) == 3 & P(x,:) == 0);
          if ~isempty(Y),
            P(z,Y) = 1;
            P(Y,z) = 2;
            % add to new_arrow
            new_arrow(Y) = 1;
            if (debug > 1), for y = Y(:)', fprintf('FCI: R1 orients  (%d *->) %d --> %d \n',x,z,y); end; end; 
            nRules(2) = nRules(2) + length(Y);  % R1 = 2
          end;
        end;  % for x
      end;  % for z

      % ============================
      % R2: definite non-ancestor 
      % if z --> x *-> y (R2a) or z *-> x --> y (R2b) and z *-o y then orient z *-> y
      for x = 1:N
        % R2a: get nodes with z --> x 
        Z = find(P(x,:) == 2 & P(:,x)' == 1);
        for z = Z(:)'   % forces Z==[] as 'empty matrix 1x0' (otherwise runs for z=[]) 
          % find x *-> y + z *-o y
          Y = find( P(:,x)' == 2 & P(:,z)' == 3 );
          if ~isempty(Y),
            P(Y,z) = 2;
            % add to new_arrow
            new_arrow(Y) = 1;
            if (debug > 1), for y = Y, fprintf('FCI: R2a orients %d *-> %d (via %d --> %d *-> %d) \n',z,y,z,x,y); end; end; 
            nRules(3) = nRules(3) + length(Y);  % R2a = 3
          end;
        end;  
        % R2b: get nodes with x --> y 
        Y = find(P(x,:) == 1 & P(:,x)' == 2);
        for y = Y(:)'
          % find z *-> x + z *-o y
          Z = find( P(x,:) == 2 & P(y,:) == 3 );
          if ~isempty(Z),
            P(y,Z) = 2;
            % add to new_arrow
            new_arrow(y) = 1;
            if (debug > 1), for z = Z(:)', fprintf('FCI: R2b orients %d *-> %d (via %d *-> %d --> %d) \n',z,y,z,x,y); end; end; 
            nRules(4) = nRules(4) + length(Z);  % R2b = 4
          end;
        end;      
      end;  % for x

      % ============================
      % R3: conditional dependence 
      % if (x *-> z <-* y), (x *-o w o-* y), {x,y} not adjacent and (w *-o z) then orient (w *-> z)
      for z = check_nodes(:)'
        % loop over nodes with arrowheads at z
        S = find(P(z,:) == 2);
        nS = length(S);
        % loop over all pairs (combinations)
        for i = 1:(nS-1)
          for j = (i+1):nS
            x = S(i); y = S(j);
            if (P(x,y) > 0), continue; end;   % nonadjacent
            W = find( P(z,:) == 3 & P(:,x)' == 3 & P(:,y)' == 3 );
            if ~isempty(W),
              P(z,W) = 2;
              % add to new_arrow
              new_arrow(z) = 1;
              if (debug > 1), for w = W(:)', fprintf('FCI: R3 orients  %d *-> %d from %d *-> %d <-* %d \n',w,z,x,z,y); end; end; 
              nRules(5) = nRules(5) + length(W);  % R3 = 5
            end;
          end; 
        end;
      end;  % for z

      % see if we're done
      check_nodes = find(new_arrow > 0);  % only nodes with (new) arrowheads need to be checked
    end;  % while (check_nodes)

    % ===================================
    % 2b: process R4 discriminating paths
    % ===================================
    % if (z o-* y) and discr.path=(x,...,w,z,y), then orient z --> y in P iff z --* y in G, 
    % otherwise orient w <-> z <-> y 
    % NOTE: New version without explicit construction of the discr.path: 
    %   per y: collect all parents W, group in bidirected components, per
    %   node w (in component Ww), if there is a valid x *-> (Ww) not adjacent to y,
    %   then for all w <-* z o-* y orient z --> y. (see also R4c in mec_to_cpag)
    % ERRATUM: in proof of R4' in the supplement to the GPS article, this was erroneously
    % described as 'intersection of district and parents', but it should be
    % as described in the main article+code below 'district between parents'.
    new_arrow = zeros(1,N);
    % get bidirected subgraph of P
    B = zeros(N,N);
    B(P == 2 & P' == 2) = 1;   % arrowhead marks (=2) on both sides
    % loop over all nodes y
    for y = 1:N
      % get all identified parents W of y in P
      W = find(P(y,:) == 2 & P(:,y)' == 1);   % W --> y most limiting feature first
      if isempty(W), continue; end;
      % partition W in bidirected clusters of w 
      BW = B(W,W);  % get bidirected subgraph over W
      BW = double(expm(BW) > 0);    % expm ensures BW(w) = 1
      % loop over W
      for w = W(:)'
        % get all W in the same bidirected component (in BW) as w
        Ww = W(find(BW(W == w,:) > 0));
        % check for some x *-> w' in Ww, that is not adjacent to y
        x = find(sum(P(Ww,:) == 2,1) > 0 & P(y,:) == 0,1);
        if ~isempty(x)
          % collect all discriminated nodes not yet oriented in P
          Z = find(P(w,:) == 2 & P(:,y)' == 3);  % w <-* Z o-* y
          for z = Z(:)'
            % discr.path x *-> (W <->) w <-* z o-* y  
            % check edge mark z to y in G (equivalent to 'z in sepset(x,y)?')
            if (G(z,y) == 1),   % tail mark at z --* y in G, so R4a (noncollider)
              % put arc at z --> y
              P(z,y) = 1;
              P(y,z) = 2;
              % add to new_arrow    NOTE: may need to check R4 again with z
              new_arrow(y) = 1;
              if (debug > 1), fprintf('FCI: R4a orients %d --> %d \n',z,y); end; 
              nRules(6) = nRules(6) + 1;  % R4a = 6
            else  % arrowhead mark at z <-* y in G, so R4b (collider)
              % put w <-> z <-> y
              P(y,z) = 2;
              P(z,y) = 2;
              P(z,w) = 2;
              % add to new_arrow
              new_arrow([z,y]) = 1;
              if (debug > 1), fprintf('FCI: R4b orients %d <-> %d <-> %d \n',w,z,y); end; 
              nRules(7) = nRules(7) + 1;  % R4b = 7
            end;          
          end;  % for z
        end;  % if (x)
      end;  % for w
    end;  % for y   
    % === end of Rule 4 ====    
    
    % see if we're done
    check_nodes = find(new_arrow > 0);  % only nodes with (new) arrowheads need to be checked    
  end;  % while (check_nodes)
  
end  % function FCI_Orientation_R1_R4b

% =======================================================================
% Orientation rules R5-R7 
% NOTE: new version covers all 3 rules at once per node (much faster)
% Approach: do not try to find an uncovered noncollider path >4 and then
% execute R5 (and afterwards R6+R7), but start from a node and see if
% propagating (hypothetical) tail marks along uncovered noncolliders
% reaches back. If so then that node has definite selection bias (R5) and
% all hypothetical tail marks correspond to invariant tail marks; if not
% then that node has no sel. bias and the hypothetical tails should remain circles.  
% =======================================================================
function P = FCI_Orientation_R5_R7(P)    
  % identify undirected edges from uncovered circle paths (R5, sel.bias),
  % and propagate invariant tails (R6+R7)
  % 1: initialise
  N = size(P,1); 
  % do temp orientations on copy of P (not needed, but just to make sure)
  U = P; % all edges
  % get nodes without arrowheads but with circle edges
  PosSelNodes = find(sum(P == 2,2) == 0 & sum((P+P' == 6),2) > 0);
  % track identifiable sel.bias status of nodes in Usel (array over N nodes) 
  % Usel(i) = -1 => definitely no sel.bias (skip), 
  %         = +1 => definitely sel.bias, propagate invariant tails (R5/6/7)
  %         = +2 => definitely sel.bias, fully processed already (skip)
  %         =  0 => not yet known, so propagate hypothetical tail marks to
  %         see if they reach back: if so then def. sel.bias, if not then not.
  Usel = -1*ones(1,N);      % initialise as definitely no sel.bias ...
  Usel(PosSelNodes) = 0;    % ... except for the ones we don't know yet
    
  % find next node to process: first check for unprocessed def. sel.bias
  % (Usel = 1), or else unknown (Usel = 0) 
  i = find(Usel == 1,1);  % def. sel.bias, not fully processed yet
  if isempty(i), i = find(Usel == 0,1); end; % otherwise first unknown
  % reserve block to process edges 
  check_Edges = zeros(N*5,2); % maybe adjust to make sure it is enough
  
  % 2: loop over nodes to process (2a) or check (2b)
  % NOTE: new invariant tails in U are only set in 2a loop (def.sel.bias node),
  % the 2b loop only identifies def.sel.bias (or not) for unknown nodes
  while ~isempty(i)    
    % determine if (unprocessed) def.sel bias or unknown
    if (Usel(i) == 1),
      % 2a: definite sel.bias on i, but not fully processed yet, so propagate invariant tails 
      % first all nodes adjacent to i get definite tail edge (R5+R6)
      Y = find(U(i,:) == 3); % all adjacent to i, but still circle mark
      nEdges = length(Y);
      % add def. tail marks at i between each nonadjacent pair in X
      for j = 1:nEdges
        y = Y(j);
        U(i,y) = 1;
        check_Edges(j,:) = [i,y]; % to propagate via R7 in block below
        % see if we created another undirected edge y-z
        if U(y,i) == 1
          % yes! signal sel.bias.on both (if not already done, but should not be)
          if Usel(y) == 0, Usel(y) = 1; end;
        end;
      end; % for j
      if (debug > 1), for y = Y(:)', fprintf('R5/R6 orients  %d --* %d \n',i,y); end; end; 
      nRules(8) = nRules(8) + length(Y);  % R5 = 8 (no distinction R5/R6)

      % loop over all tail marks on edges to propagate from i (R7)
      idxEdge = 1; 
      while (idxEdge <= nEdges)
        % get next edge 
        x = check_Edges(idxEdge,1);
        z = check_Edges(idxEdge,2);
        idxEdge = idxEdge + 1;
        % find nodes with a circle OR HYPO mark at Z not adjacent to x, but not x
        idxY = ((U(z,:) == 3 | U(z,:) < 0 ) & U(x,:) == 0);
        idxY(x) = 0;
        Y = find(idxY);
        for y = Y(:)'
          % check if we reached the another (or starting?) hypotail
          U(z,y) = 1; % add new def.tail 
          nEdges = nEdges + 1;
          check_Edges(nEdges,:) = [z,y];
          
          % see if we created another undirected edge y-z
          if U(y,z) == 1
            % yes! signal sel.bias.on both (if not already done, but should not be)
            if Usel(y) == 0, Usel(y) = 1; end;
            if Usel(z) == 0, Usel(z) = 1; end;
          end;
        end;  % for y
        if (debug > 1), for y = Y(:)', fprintf('R7 orients (%d --o) %d --* %d \n',x,z,y); end; end; 
        nRules(10) = nRules(10) + length(Y);  % R7 = 10 
          
      end; % while idxEdge

      % signal def.selnode 'i' processed
      Usel(i) = 2;
        
    else
      % 2b: Usel(i) = 0 => unknown node: propagate hypothetical tails along
      % each edge, and see if we reach back or not. 
      
      % add hypo(thetical) tail at each until reaching back
      % note: could already have a def.tail! .. so proc only circles
      J = find(U(i,:) == 3); % find gives row vector
      Usel_found = false;
      hypo_mark = 0;    % use negative marks for hypothetical tails (that can 
                        % easily be removed again if they do not reach back)

      % add hypo tail marks at i between each nonadjacent pair in X
      % process each in turn and propagte to see if it reaches back
      for j = J(:)'
        % first check exist k: k-i-j, k and j not adjacent (but not j itself)
        K = find( U(i,:) ~= 0 & U(j,:) == 0);
        if isempty(find(K ~= j,1)), continue; end;
        
        % ok: start and process hypo tail 
        hypo_mark = hypo_mark - 1;  % avoid trying edges multiple times via different paths
        U(i,j) = hypo_mark;   % note: was circle mark in U
        
        check_Edges(1,:) = [i,j];
        idxEdge = 1; nEdges = 1;
        Usel_found = false;
        
        while (idxEdge <= nEdges) && ~Usel_found
          % pop edge to process
          x = check_Edges(idxEdge,1);
          z = check_Edges(idxEdge,2);
          idxEdge = idxEdge + 1;
          % find nodes with a circle OR (current) hypo mark at Z not adjacent to x, but not x
          idxY = ((U(z,:) == hypo_mark | U(z,:) == 3 ) & U(x,:) == 0);
          idxY(x) = 0;
          Y = find(idxY);
          for y = Y(:)'
            % check if we reached another 'hypotail'
            if (U(z,y) == hypo_mark)  
              if (z == i) && (y == j)  % reached starting mark?
                % yes! signal sel.bias.
                Usel(z) = 1;
                Usel_found = true;
                break;  % break (twice) out of J loop (then process as def.sel.bias node in 2a loop) 
              end;
            % propagate new hypo on circle amrk
            else
              % so U(z,y) == 1
              U(z,y) = hypo_mark; % add new hypo edge
              nEdges = nEdges + 1;
              check_Edges(nEdges,:) = [z,y];
            end;
          end;  % for y
          
        end; % while idxEdge
        
        if (Usel_found), break; end; % break out of J loop (then process as def.sel.bias node in 2a loop) 
      end;  % for j        
      
      % reset circle marks back for all hypo tails in U
      % NOTE: could be adapted so if Usel found then hypo tails set are
      % kept as tails, but for simplicity we only set invariant tails in
      % '2a loop'
      U(U < 0) = 3;
      % processed all possible paths for i: either sel.node or not
      if ~(Usel_found),
        % if not here, then never .. signal no sel. bias
        Usel(i) = -1;
      end;
    end;  % proc node i

    % see if we're done: try to get the next i
    i = find(Usel == 1,1);  % def. sel.bias, not fully processed yet
    if isempty(i), i = find(Usel == 0,1); end; % otherwise first unknown
  
  end;  % while ~isempty(i)

  % all done: transfer back to P as undirected graph component
  P(U == 1) = 1;

end;  % subfunction FCI_Orientation_R5_R7


% =======================================================================
% Orientation rules R8-R10 
% =======================================================================
function P = FCI_Orientation_R8a_R10(P)
  % 1: initialize
  N = size(P,1);
  
  % start from all nodes with arrowheads
  check_nodes = find(sum(P == 2,2)' > 0);
  
  % 2: repeat Rule 8-10, adding tails until no more can be found
  while ~isempty(check_nodes)    
    % keep track of new arrownodes 
    new_tail = zeros(1,N);

    % ============================
    % R8: transitivity
    % if z --> x --> y (R8a) or z --o x --> y (R8b) and z o-> y then orient z --> y
    for x = 1:N
      % get x --> Y
      Y = find( P(x,:) == 1 & P(:,x)' == 2 );
      for y = Y(:)'
        % R8a: get nodes with z --> x and z o-> y 
        Z = find(P(x,:) == 2 & P(:,x)' == 1 & P(y,:) == 2 & P(:,y)' == 3);
        if ~isempty(Z),
          P(Z,y) = 1;
          new_tail(Z) = 1;
          if (debug > 1), for z = Z(:)', fprintf('R8a orients %d --> %d via %d --> %d --> %d \n',z,y,z,x,y); end; end; 
          nRules(11) = nRules(11) + length(Z);  % R8a = 11
        end;
        % R8b: get nodes with z --o x and z o-> y 
        Z = find(P(x,:) == 3 & P(:,x)' == 1 & P(y,:) == 2 & P(:,y)' == 3);
        if ~isempty(Z),
          P(Z,y) = 1;
          new_tail(Z) = 1;
          if (debug > 1), for z = Z(:)', fprintf('R8b orients %d --> %d via %d --o %d --> %d \n',z,y,z,x,y); end; end; 
          nRules(12) = nRules(12) + length(Z);  % R8b = 12
        end;        
      end;  % for y
    end;  % for x
    
    % ============================
    % R9: transitive path (uncovered p.d. path, new faster version like R5 edges)
    % if x o-> y, u=(x,z,v,..,y) is u.p.d path, z,y not adjacent then orient x --> y    
    % preallocate block to process edges 
    check_Edges = zeros(N*5,2); % maybe adjust to make sure it is enough
    % get all x o-> y
    [X,Y] = find(P == 3 & P' == 2);
    nXY = size(X,1);
    for i = 1:nXY
      % process next arc
      x = X(i); y = Y(i);
      % keep track of edges processed
      U = zeros(size(P));     
      % note: not needed for proper PAGs, but in case of invalid MECs it
      % could lead to directed cycles while constructing, leading to infinite loops

      % find first nodes z along possibly u.p.d. 
      idxZ = (P(x,:) == 1 | P(x,:) == 3) & P(y,:) == 0;
      idxZ(y) = 0;    % but not y itself obviously ..
      Z = find(idxZ);
      % fill starting edges for x - Z 
      nEdges = length(Z);
      % add edges to follow
      for j = 1:nEdges
        z = Z(j);
        check_Edges(j,:) = [x,z]; 
        U(x,z) = 1;
      end;

      % now process
      idxEdge = 1; reached_y = false;
      while (idxEdge <= nEdges) && ~reached_y
        % pop edge
        z1 = check_Edges(idxEdge,1);
        z2 = check_Edges(idxEdge,2);
        idxEdge = idxEdge + 1;
        % find nodes with a circle or tail mark at z2 not adjacent to z1, but not z1
        idxZ3 = ((P(z2,:) == 3 | P(z2,:) == 1 ) & P(z1,:) == 0);
        idxZ3(z1) = 0;
        Z3 = find(idxZ3);
        for z3 = Z3(:)'
          % check if we reached the target node y
          if (z3 == y)
            % yes! signal and break
            reached_y = true;
            break;
          elseif (U(z2,z3) == 0)
            % not yet and new: add and propagate
            nEdges = nEdges + 1;
            check_Edges(nEdges,:) = [z2,z3];
            U(z2,z3) = 1;
            % track of edges already processed to avoid infinite loops for invalid MECs
          end;
        end;  % for z3
      end; % while idxEdge

      % update if reached
      if (reached_y)
        % set tail at x to y
        P(x,y) = 1;
        new_tail(x) = 1;
        if (debug > 1), fprintf('5.R9 orients %d --> %d\n',x,y); end;
        nRules(13) = nRules(13) + 1;  % R9 = 13
      end;

    end;  % for i in [X,Y]

  
    % ============================
    % R10: double uncovered pdp
    % if x o-> y, q --> y <-- z, u1=(x,(s,..,)z) and u2=(x,(q,..,)r) u.p.d. paths 
    % then if s (possibly z) and q (possibly r) not adjacent, then orient x --> y
    % note: the arcs z --> y <-- v result from R9
    % TODO: alter R10 like R5/R9 (but rarely triggered anyway, so leave for now)
    for y = 1:N
      % get all nodes with z --> y
      Z = find( P(y,:) == 2 & P(:,y)' == 1 );
      nZ = length(Z);
      if (nZ < 2), continue; end;
      % get all x o-> y
      X = find( P(y,:) == 2 & P(:,y)' == 3);
      for x = X(:)'
        % get all nodes with no tail/circle adjacent to z
        % note: no arrowhead, as e.g. s *-> x would imply x --> p (by R1),
        % and so on towards y, and then back to x via R8a, and so x --> y
        % note: actually also adjacent to y, otherwise R9 would have triggered
        T = find( (P(x,:) == 1 | P(x,:) == 3) & (P(:,x)' == 1 | P(:,x)' == 3) ...
                & (P(:,y)'== 1 | P(:,y)'== 3) );
        nT = length(T);
        if (nT < 2), continue; end;
        % find first uncovered pdp, if found try second
        % try from start S, stop if a node in Z is reached? no. just try in
        % order (not that many anyway)
        for i = 1:(nT-1)
          t = T(i);         
          % find all nodes in S not adjacent to s (used in up2)
          Q = mysetdiff(T(P(t,T) == 0),t);
          nQ = length(Q);
          % ok: try first upd path x - S(i) .. Z(j)
          for j = 1:nZ      % note: not nZ-1
            z = Z(j);
            % find uncovered pdp x -o s - .. z
            % note: we cannot look for just 'any' path (as in R9), as it is
            % about a specific combination of two disjoint paths 
            % store beginning of possible uncovered possibly directed paths [x (-/o)-(o/>) s (..) y]
            up1{1} = [x,t,z];

            % process edge
            while ~isempty(up1)
              u1 = up1{1}(end-2);   % yes: handle s==z (-> go to up2
              v1 = up1{1}(end-1);
              % get possible extension nodes: u o-> v o-> w, with w not adjacent to u)
              idxW1 = (P(v1,:) == 1 | P(v1,:) == 3 ) & (P(:,v1)' == 2 | P(:,v1)' == 3 ) & (P(u1,:) == 0);
              % check if we're there: reach target z (no additional checks needed)
              if (idxW1(z) > 0) || (t == z), 
                % yes: uncovered possibly directed path: try second path!
                % get all other nodes in Z
                R = mysetdiff(Z,z);
                % now loop for second path
                for k = 1:nQ            % note: not nQ-1
                  q = Q(k);
                  for l = 1:(nZ-1)      % nR = nZ - 1
                    r = R(l);
                    up2{1} = [x,q,r];   % note: up2{one}
                    
                    % process edge (identical to first path
                    while ~isempty(up2)
                      u2 = up2{1}(end-2);    % avoid overwriting {u1,v1}
                      v2 = up2{1}(end-1);
                      % get possible extension nodes: u o-> v o-> w, with w not adjacent to u)
                      idxW2 = (P(v2,:) == 1 | P(v2,:) == 3 ) & (P(:,v2)' == 2 | P(:,v2)' == 3 ) & (P(u2,:) == 0);
                      % check if we're there: reach target z (no additional checks needed)
                      if (idxW2(r) > 0) || (q == r), 
                        % yes: second uncovered possibly directed path!
                        P(x,y) = 1;
                        new_tail(x) = 1;
                        if (debug > 1), fprintf('5.R10 orients %d --> %d via upds [%d,..,%d] and [%d,..,%d] \n',...
                           x,y,up1{1}(2),up1{1}(end-1),up2{1}(2),up2{1}(end-1)); end;   % include up2
                        nRules(14) = nRules(14) + 1;  % R10 = 14
                        break;    % go to next unchecked edge (outer while loop)
                      else
                        % no: exclude nodes already on the path (and/or the
                        % first path? ... yes: otherwise other rule should
                        % have triggered already
                        idxW2(up2{1}) = 0;
                        W2 = find(idxW2 > 0);
                        % extend path with 
                        for w2 = W2(:)',
                          up2{end+1} = [up2{1}(1:end-1),w2,up2{1}(end)];
                        end
                        % try next path
                        up2 = up2(2:end);
                      end;
                    end;  % while up2
                    if (P(x,y) == 1), break; end;
                  end;  % for l
                  if (P(x,y) == 1), break; end;
                end;  % for k
              else  
                % z not reached in idxW1
                % exclude nodes already on the path 
                idxW1(up1{1}) = 0;
                W1 = find(idxW1 > 0);
                % extend path with 
                for w1 = W1(:)',
                  up1{end+1} = [up1{1}(1:end-1),w1,up1{1}(end)];
                end
              end;  % if (z in idxW1)
              % try next path
              up1 = up1(2:end);
              if (P(x,y) == 1), break; end;
            end;  % while up1
            if (P(x,y) == 1), break; end;
          end;  % for j
          if (P(x,y) == 1), break; end;
        end;  % for i 
        
      end;  % for x
    end;  % for y
    
    % see if we're done
    check_nodes = find(new_tail > 0);  % only nodes with (new) arrowheads need to be checked
  end;  % while 

end;  % subfunction FCI_Orientation_R8a_R10
  
end  % mag_to_cpag