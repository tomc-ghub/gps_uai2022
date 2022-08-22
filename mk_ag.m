function [M,Reach,Vt,Sel] = mk_ag(N,maxK,avgK,p_bi,p_sel)
% MK_AG Create random ancestral graph over N variables with max. node
% degree maxK, and probability on bi-directed edge = p_bi.
% M = mk_ag(N, maxK, p_bi)
% input:
% - N       : nr.of nodes
% - maxK    : max.degree
% - avgK    : average node degree (2 < avgK < maxK)
% - p_bi    : probability of edge being bidirected
% - p_sel   : probability of selection bias per node (default 1/N)
% output:
% - M       : ancestral (causal) graph 
% - Reach   : reachability/ancestor matrix of directed paths in M
% - Vt      : consistent time ordering of nodes in V
% - Sel     : nodes with (direct) selection bias
% Note: on laptop 'mk_ag(500,8,3.5,0.3);' takes about 1 sec. 

% Approach: 
% 1: generate undirected skeleton over N variables 
%    - start from empty graph, 
%    - add random edges (ensure other properties?)
% 2: pick edge, choose random orientation
global DEBUG;
if isempty(DEBUG), debug = 0; else debug = DEBUG; end;

  % 1 - Initialize  
  % check/process input parameters
  if (nargin < 2), maxK = 5; end; 
  % and convert N, maxK into min/max/avg edge limits
  posEdge = N*(N-1)/2;
  minEdge = (N-1); 
  maxEdge = N*maxK/2;
  
  % check edge probabilities, if necessary adjust p_Edge to satisfy
  % connected but not too dense graphs
  if (nargin < 3), avgK = (1 + maxK/2); 
  elseif (avgK < 2), avgK = 2;
  elseif (avgK > maxK), avgK = maxK; end;
  p_Edge = avgK/(N-1);
  % sample number of edges until nEdge falls within valid interval
  % NOTE: not guaranteed exactly nEdge edges in M. If close to maximum for
  % given maxK, then in some configurations no more edges can be added.
  % This is a different type of problem, and will not be tackled here.
  nEdge = 0;
  while (nEdge < minEdge) || (nEdge > maxEdge)
    nEdge = binornd(posEdge,p_Edge);  % 
  end;  % while
  if (debug > 1), fprintf('%i bi-directed edges\n', nEdge); end;
  
  % check probability of bi-directed edges
  if (nargin < 4), p_bi   = 0.3; end; 
  % check probability of selection bias per node
  if (nargin < 5), p_sel  = 1/N; end; 
  % large N: 36.8% 0 nodes, 36.8% 1 node, 18.4% 2 nodes, 6.1% 3 nodes ... etc.
  
  % initialize variables
  M = zeros(N);       % empty graph
  D = zeros(1,N);     % zero degree for all nodes
  V0 = [1:N];         % array of nodes not yet connected (prio)
  Vk = [];            % array of nodes to choose from (prio 2)
  VK = [];            % array of nodes with degree = maxK (no more candidate)
  
  % =======================================================================
  % 2 - Skeleton
  % created connected skeletom over nEdge edges
  % start with random node in Vk (and remove from V0)
  i = randi(N); Vk = [i]; V0(i) = [];
  % first connect all nodes (perhaps weighted by degree, order?)
  for nE = 1:(N-1)
    % start with one from Vk then one from V0 or Vk
    i = randi(length(Vk));
    x = Vk(i);    
    % new not-yet-connected node
    j = randi(length(V0));
    y = V0(j);
    % add undirected edge to M
    M(x,y) = 1; M(y,x) = 1;
    % update arrays/degrees
    V0(j) = [];
    Vk    = [Vk,y];
    % increase degree of both
    D([x,y]) = D([x,y])+1;
    % move nodes at max.degree to VK (maxK > 1, so never Vj)
    if (D(x) >= maxK), 
      VK    = [VK,x];
      Vk(i) = [];
    end;
  end;  % for all unconnected
  if (~isempty(V0)), disp('Error mk_ag.1'); end;
  
  % now add all other edges
  for nE = N:nEdge
    % draw node from Vk  (weighted by degree, order?)
    i = randi(length(Vk));
    x = Vk(i);
    % draw second node via temp array of valid candidates
    Vx = Vk; 
    for j = length(Vk):-1:1
      % remove invalid candidates: no self-loops or already existing edges
      y = Vk(j);
      if (x == y) || (M(x,y) > 0), Vx(j) = []; end;
    end;
    % it is possible that Vx becomes empty, so then just continue
    if isempty(Vx), continue; end;
    % get node and get j in terms of index in Vk
    j = randi(length(Vx));      % with rand('state',3) Vx gets to be empty ...?
    y = Vx(j);
    j = find(Vk == y,1);    
    % add undirected edge to M
    M(x,y) = 1; M(y,x) = 1;
    % increase degree of both
    D([x,y]) = D([x,y])+1;
    % check degree, first y / j
    if (D(y) >= maxK), 
      VK    = [VK,y];
      Vk(j) = [];
    end;
    % check degree, then x / i
    if (D(x) >= maxK), 
      VK    = [VK,x];
      Vk(Vk == x) = [];     % in case j < i
    end;
  end;  % for i=N:nEdge
  
  % =======================================================================
  % 3 - orient edges
  % bookkeeping: make M = circle skeleton of M 
  M = 3*M;
  % set random time-order (to ensure consistent orientations for arcs)
  Vt = randperm(N);
  if (debug > 1), disp(Vt); end;
  % get inverse permutation?
  B = sortrows([1:N;Vt]',2)';
  Vinvt = B(1,:);
  % ancestor graph A(x,y) = {0,1,2} = {2 = x not ancestor of y, 1 = is}
  A = tril(2*ones(N) - eye(N));
  A = A(Vinvt,Vinvt);
  % set counters
  nArc = 0;         % nr.of arcs allocated in M
  nBi  = 0;         % nr.of bi-directed arcs in M
  % per node: orientations to do is degree unoriented 
  %Vx = [1:N];       % nodes with edges to orient
  Dx = D;           % nr.of edges at node to orient
  % loop until all are done
  while (sum(Dx) > 0)
    % select node with unoriented edge
    Vx = find(Dx > 0);
    i = randi(length(Vx));
    x = Vx(i);
    % select one from set of nodes on unoriented edge to x
    Vy = find(M(x,:) == 3);
    j = randi(length(Vy));
    y = Vy(j);
    % ensure x->y consistent with Vt
    if (Vinvt(x) > Vinvt(y)), 
      % swap direction
      tmp = x; x = y; y = tmp;
    end;
    
    % 1. update nr.of unoriented edges
    Dx([x,y]) = Dx([x,y]) - 1;
    add_bi = (rand < p_bi);
    implied = false;
    % determine arc/bi
    if (add_bi), 
      AddBi(x,y); 
    else
      % add arc x->y
      AddArc(x,y);
    end;  % add arc/bi
    
    % process
    implied = true;
    if (add_bi)
      % if x -> y assigned as x <-> y, then a.d.cycle with x -> .. -> u -> v -> .. -> y
      % if added bi-dir x <-> y, with x < y in Vt, then
      % - get De(x)/An(y) and add all unoriented edges De(x) <-> An(y)
      % check implied
      DeX = myunion(x,find(A(x,:) == 1) );
      AnY = myunion(y,find(A(:,y) == 1)');
      % first check x << y     (superfluous)
      if (Vinvt(x) > Vinvt(y)),
        disp('ERROR - add bi-dir.2'); 
      end;
      % find and process all unoriented edges De(x)-An(y)
      [UnDx,UnDy] = find(M(DeX,AnY) == 3);
      for i = 1:length(UnDx)
        u = DeX(UnDx(i));   % convert indices in DeX to actual nodes
        v = AnY(UnDy(i));   % convert indices in AnY to actual nodes
        % if arc would be u -> v, then a.d.cycle
        if (Vinvt(u) < Vinvt(v)),
          % add bi-directed edge 
          AddBi(u,v);
          % update edges-to-orient counter
          Dx([u,v]) = Dx([u,v]) - 1;
        end;
      end;  % for i:UnDx
      
    else
      % add arc
        % if added arc x -> y then 
        % - get An(x)/De(y) and add all unoriented edges An(x) -> De(y)
        % - get all bi-dirs anx_i <-> z to An(x), and then get An(z), 
        %   if unoriented edge anz_j - dey_k, and dey_k < anz_j in Vt,
        %   then add anz_j <-> dey_k
        % - get all bi-dirs dey_i <-> z to De(y), and then get De(z),
        %   if unoriented edge dez_j - anx_k, and dez_j < anx_k in Vt,
        %   then add anz_j <-> dey_k
        
      % check implied
      AnX = myunion(x,find(A(:,x) == 1)');
      DeY = myunion(y,find(A(y,:) == 1) );
      % 1 - find and process all unoriented edges An(x)-De(y) (similar to add_bi)
      [UnDx,UnDy] = find(M(AnX,DeY) == 3);
      for i = 1:length(UnDx)
        u = AnX(UnDx(i));   % convert indices in AnX to actual nodes
        v = DeY(UnDy(i));   % convert indices in DeY to actual nodes
        % add arc u -> v (instead of bi-directed edge)
        AddArc(u,v);
        % update edges-to-orient counter
        Dx([u,v]) = Dx([u,v]) - 1;
      end;  % for i:UnDx
      
      % 2 - get all bi-dirs anx_i <-> z to An(x), and then get An(z), 
      %   if unoriented edge anz_j - dey_k, and dey_k < anz_j in Vt,
      %   then add anz_j <-> dey_k
      BiZidx = myintersect(find(M(AnX,:)==2), find(M(:,AnX)'==2));
      [BiX,BiZ] = ind2sub([length(AnX),N],BiZidx);
      % loop over all these bi-directed edges
      for i = 1:length(BiX)
        bix = AnX(BiX(i));  % convert indices in AnX to actual nodes
        z   = BiZ(i);       % indices in BiZ are already actual nodes
        AnZ = find(A(:,z)' == 1);   % ancestors of z (incl.z)
        
        % now find unoriented edges AnZ - DeY
        [UnDz,UnDy] = find(M(AnZ,DeY) == 3);
        for j = 1:length(UnDz)
          u = AnZ(UnDz(j));   % convert indices in AnZ to actual nodes
          v = DeY(UnDy(j));   % convert indices in DeY to actual nodes
          % check v << u in Vt
          if (Vinvt(u) > Vinvt(v)),     % so v << u in Vt
            % add bi-dir u <-> v
            AddBi(u,v);
            % update edges-to-orient counter
            Dx([u,v]) = Dx([u,v]) - 1;
          end;
        end;  % for j
      end;  % for i
      
      % 3 - get all bi-dirs dey_i <-> z to De(y), and then get De(z),
      %   if unoriented edge dez_j - anx_k, and dez_j < anx_k in Vt,
      %   then add anz_j <-> dey_k
      BiZidx = myintersect(find(M(DeY,:)==2), find(M(:,DeY)'==2));
      [BiY,BiZ] = ind2sub([length(DeY),N],BiZidx);
      % loop over all these bi-directed edges
      for i = 1:length(BiY)
        biy = DeY(BiY(i));  % convert indices in AnX to actual nodes
        z   = BiZ(i);       % indices in BiZ are already actual nodes
        DeZ = find(A(z,:) == 1);   % descendants of z (incl.z)
        
        % now find unoriented edges DeZ - AnX
        [UnDz,UnDx] = find(M(DeZ,AnX) == 3);
        for j = 1:length(UnDz)
          u = DeZ(UnDz(j));   % convert indices in AnZ to actual nodes
          v = AnX(UnDx(j));   % convert indices in DeY to actual nodes
          % check u << v in Vt
          if (Vinvt(u) < Vinvt(v)),     % so u << v in Vt
            % add bi-dir u <-> v
            AddBi(u,v);
            % update edges-to-orient counter
            Dx([u,v]) = Dx([u,v]) - 1;
          end;
        end;  % for j
      end;  % for i

    end;  % if add_bi else
        
  end;  % while edges to do
  
  % compute reachability graph (ancestor in M)
  B = zeros(size(M));
  B(M == 1) = 1;            % adjacency graph (only keep tails from M)
  Reach = reachability_graph(B);

  
  % finally add selection nodes (leading to undirected edges) if necessary
  if (p_sel > 0) && (p_sel <= 1),
    % sample selection nodes
    Sel = find(binornd(1,p_sel,1,N) > 0);
    % find all nodes with directed paths to [Sel]
    AllSel = myunion(find(sum(Reach(:,[Sel]),2) > 0),Sel);
    M2 = zeros(size(M));
    M2([AllSel],:) = M([AllSel],:);
    % set all 'from' edge marks to tails
    M(M2 > 1) = 1;
  else 
    % no selection
    Sel = [];
  end;


  % =======================================================================
  % verify graph M ?
  if (debug > 1)
    % no unprocessed edges
    if ~isempty(find(M > 2,1)),
      disp('ERROR - unprocessed edges');
    end;

    % get reachability graph from M
    B  = M;
    B(M(:,:) == 2) = 0;   % get adjacency graph (remove arrowheads)
    % R = reachability_graph(B);
    Bn = B;
    R2  = zeros(N);
    for i = 1:N
      R2  = R2 + Bn;
      Bn = Bn * B;
    end;
    R2(R2(:,:) > 0) = 1;
    % check self-loops
    if ~isempty(find(diag(R2) > 0, 1)),
      disp('ERROR - self loop in M'); 
    end;

    % now loop over all bi-directed edges
    B = triu(M + M');
    [X,Y] = find(B == 4);
    % loop over all these bi-directed edges
    for i = 1:length(X)
      % get next bidirected edge to check
      x = X(i); y = Y(i);
      if (R2(x,y) > 0) || (R2(y,x) > 0),
        fprintf('ERROR - almost directed cycle %i <-> %i\n',x,y);
      end;
    end;
  end;  % if verify
  
  
%   % finally: show if needed (only useful if N < 20)
%   if DEBUG, figure(6); clf; draw_cpmag(M); end;

  
  % ======================================================================
  function AddBi(x,y)
  % add bi-directed edge x <-> y to M, verify valid, and update A, nBi 
    M(x,y) = 2; M(y,x) = 2;
    if (debug > 1), 
      if (implied),
        fprintf('Add implied %i <-> %i \n',x,y); 
      else
        fprintf('Add %i <-> %i \n',x,y); 
      end;
    end;
    nBi = nBi + 1;
    % update A(ncestral matrix) (but verify first)
    if (A(x,y) == 1) || (A(y,x) == 1), 
      fprintf('ERROR - AddBi: %i <-> %i\n',x,y); 
    end;
    A(x,y) = 2; A(y,x) = 2;  
  end % local function AddBi
  
  function AddArc(x,y)
  % add arc x --> y to M, verify valid, and update A, nArc
    M(x,y) = 1; M(y,x) = 2;
    if (debug > 1), 
      if (implied),
        fprintf('Add implied %i --> %i \n',x,y); 
      else
        fprintf('Add %i --> %i \n',x,y); 
      end;
    end;
    nArc = nArc + 1;
    % update A(ncestral matrix) for all AnX -> DeY (local? ... )
    anX = myunion(x,find(A(:,x) == 1)');
    deY = myunion(y,find(A(y,:) == 1) );
    if ~isempty(find(A(anX,deY) == 2,1)) || ~isempty(find(A(deY,anX) == 1,1)), 
      fprintf('ERROR - AddArc: %i --> %i\n',x,y); 
    end;
    A(anX,deY) = 1; A(deY,anX) = 2;
  end % local function AddArc

  function C = reachability_graph(G,self)
    % REACHABILITY_GRAPH C(i,j) = 1 iff there is a path from i to j in DAG G
    % C = reachability_graph(G)
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
        expG = expm(double(full(G))) - eye(n); % exclude self (only proper ancestors)
      else
        expG = expm(double(full(G))); % do not exclude self : x \in An(x)
      end;
      C = (expG>0);
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
  end % local function reachability_graph

  function C = myintersect(A,B)
    % MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
    % C = myintersect(A,B)
    A = A(:)'; B = B(:)';

    if isempty(A)
      ma = 0;
    else
      ma = max(A);
    end

    if isempty(B)
      mb = 0;
    else
      mb = max(B);
    end

    if ma==0 | mb==0
      C = [];
    else
      %bits = sparse(1, max(ma,mb));
      bits = zeros(1, max(ma,mb));
      bits(A) = 1;
      C = B(logical(bits(B)));  
    end
  end % local function myintersect

  function C = myunion(A,B)
    % MYUNION Union of two sets of positive integers (much faster than built-in union)
    % C = myunion(A,B)
    if isempty(A)
      ma = 0;
    else
      ma = max(A);
    end

    if isempty(B)
      mb = 0;
    else
      mb = max(B);
    end

    if ma==0 & mb==0
      C = [];
    % elseif ma==0 & mb>0  % always enforce merge+order
    %   C = B;
    % elseif ma>0 & mb==0
    %   C = A;
    else
      %bits = sparse(1, max(ma,mb));
      bits = zeros(1, max(ma,mb));
      bits(A) = 1;
      bits(B) = 1;
      C = find(bits);
    end
  end  % local function myunion

end  % function mk_ag








