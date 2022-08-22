function MEC = mag_to_mec(G)
% Get Markov equivalence class (MEC = {S,C,D}) representation from 
% (maximal/partial) ancestral graph G over N variables. 
% Handles MAGs and (core)PAGs in the same way. If G is a modified PAG then
% undetermined triples with order are treated as default 'noncollider'. 
% Example: suppose G is obtained by adding edge z o-o y to 
%  PAG P = x o-> w <-o z + w --> y. Then [w,z,y] is a triple with order 1
%  (corresponding to discriminating path <x,w,z,y>), but from G it is not
%  unambiguously clear whether [w,z,y] is a collider or noncollider (as
%  edge mark z o-* y should be invariant either way). In that case triple
%  [w,z,y] is taken to be a *noncollider*. In other words: non-invariant
%  circle marks in G encountered when building the MEC are treated as if
%  they were tail marks.
% Note: another option is to consider a collection of MECs for a graph G
% with ambiguous triples, but here we opt for a single default MEC.
%
% Input:
% G = (maximal) ancestral graph encoded in the form  
%    (Gij, Gji) = (0,0)  : not adjacent    i     j
%               = (1,1)  : undirected edge i --- j
%               = (1,2)  : arrow           i --> j
%               = (2,1)  : arrow           i <-- j
%               = (2,2)  : bidirected edge i <-> j
% Note: starting from PAG, circle marks i o-* j encoded as Gij = 3
% Output:
% - MEC : MEC of G
%  MEC.S = skeleton as NxN matrix, 0=no edge, 1=edge
%  MEC.C = collider triples with order [k,x,z,y,q], where
%     k=order, [x,z,y]= triple, q=shared node triggering triples (order>0) 
%     combo collider [q,x,z] + noncollider [q,x,y] (not needed, tracking only)
%  MEC.D = noncollider triples with order (idem)
%
% Flow: 
% 1: initialise, extract skeleton
% 2: get unshielded triples, initialise process list
% 3: repeat (process next triple) until list is empty
% Note: in practice for random graphs k>2 is rarely encountered. For
% efficiency an index list IndC/D is used to avoid searching through the
% entire C/D lists for matching triples, so we only need to scan max.'d' 
% entries for sparse graphs with max. node degree = d.
% =========================================================================
global DEBUG;
if ~isempty(DEBUG), debug = DEBUG; else debug = 0; end; % local debug lvl

% 1: initialize
% check/process input parameters
if (nargin < 1), return; end; 
% initialize variables
N = size(G,1);
% get skeleton
S = zeros(N,N);
S(G ~= 0) = 1;
% get max. degree (note: N^2 complexity unless given or edge-based skeleton)
d = max(sum(S,1));
% intitialize triple lists + index
BlockSize = N*d;    % size to add to C/D lists on exceeding entries
C = zeros(BlockSize,5);     % preallocate BlockSize entries
D = zeros(BlockSize,5);
nC = 0;
nD = 0;
IndC = zeros(N,N,d);   % index for quick access (actually (d-1) should suffice)
IndD = zeros(N,N,d);   % 
% initialise process list (no index needed)
T = zeros(BlockSize,5);
nT = 0;

% 2: unshielded triples
% loop over all nodes z
for z = 1:N
  % find set of all nodes with edges to z 
  Adj_z = find(G(z,:) > 0);
  nAdj_z = length(Adj_z);
  % loop over all pairs {x,y} adjacent to z
  for i = 1:(nAdj_z-1)
    x = Adj_z(i); 
    for j = (i+1):nAdj_z
      y = Adj_z(j);
      % process if nonadjacent
      if (G(x,y) == 0)
        % x-z-y is unshielded triple, 
        isCol = (G(z,x) == 2) && (G(z,y) == 2);
        id = Add_TripleWithOrder(0,x,z,y,-1,isCol);      % updates nC/nD counter + logs
      end; % if (unshielded triple)
    end; % for j
  end; % for i      
end; % for z

% initialise process Triples list 
% start from colliders 
for i = 1:nC
  % get entry
  %[k,x,z,y,q] = C(i,:);
  k = C(i,1);
  x = C(i,2);
  z = C(i,3);
  y = C(i,4);
  q = C(i,5);
  
  % search for [x,z,q] or [q,z,x] or [q,z,y] or [y,z,q] in D
  % (note this could be faster with a 2-D index array, but ok for now)

  % check triple order 1 for y (allow for circle mark at z in case of modified PAG)
  Tri = find( ((G(z,:) == 1) | (G(z,:) == 3)) & (G(x,:) == 0) & (G(y,:) > 0) );
  for i = Tri
    % add T[k,a,b,c,q]
    id = Add_ProcessTriple(k+1,z,y,i,x);  % updates nT
  end;
  % and reverse triple order 1 for x
  Tri = find( ((G(z,:) == 1) | (G(z,:) == 3)) & (G(x,:) > 0) & (G(y,:) == 0) );
  for i = Tri
    % add T[k,a,b,c,q]
    id = Add_ProcessTriple(k+1,z,x,i,y);  % updates nT
  end;
  
end;  % for i


% 3: process list T until empty (recursively)
idxT = 1;
while (idxT <= nT)
  % pop top entry (actually keep counter to avoid modifying list)
  % [k,x,z,y,q] = T(idxT,:);
  k = T(idxT,1);
  x = T(idxT,2);
  z = T(idxT,3);
  y = T(idxT,4);
  q = T(idxT,5);
  idxT = idxT + 1;

  % get type (collider/noncollider) from G
  isCol = (G(z,y) == 2);  % effectively treats undetermined triples as noncollider
  if (debug > 1), 
    if (G(z,y) >= 3) 
      fprintf('WARNING:  Undetermined higher order triple [%d,%d,%d,%d,%d]: %d; treated as noncollider.\n',k,x,z,y,q,G(z,y)); 
    end;
  end;
  
  % add triple with order to C/D (updates index/counter as wel)
  [id,isNew] = Add_TripleWithOrder(k,x,z,y,q,isCol);

  % check for new triples to process
  if ~isNew
    % already found: no need to check for other triples 
    continue;
  elseif (isCol)
    % collider x*->z<-*y (+ x-->y ?) ... matches noncollider(s) in D
    % - x*->z-->q, if y-q in G  => triple [z,y,q] not already in C or D
    % - y*->z-->q, if x-q in G  => triple [z,x,q] not already in C or D

    % look for possible matching noncollider triples and add
    % match on [x,z,q] in D with g-y in G
    for id = IndD(x,z,:)
      if (id == 0), break; end;    % for sparse graphs max.'d' entries
      % get entry 
      % [m,a,b,c,u] = D(id,:);  % a=x,b=z
      m = D(id,1);
      a = D(id,2);
      b = D(id,3);
      c = D(id,4);
      u = D(id,5);
      if (debug > 1), 
        % (verify a==x,b==z .. or if possibly c==x?)
        if (b ~= z) || ((a ~= x) && (c ~= x)),
          fprintf('ERROR 3.Col non-matching triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] \n',id,m,a,b,c,u,k,x,z,y,q); 
        else
          fprintf('*** 3.Col processing triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] \n',id,m,a,b,c,u,k,x,z,y,q); 
        end;
      end;
      % check edge c-y
      if (G(c,y) > 0) && (a == x) && (b == z),
        % add process triple [l,z,y,c,x]
        l = min(k,m) + 1;   % order of triple
        id = Add_ProcessTriple(l,z,y,c,x);  % also updates nT
      end;
      % check reverse for order 0 noncollider and edge a-y
      if (m == 0) && (G(a,y) > 0) && (c == x) && (b == z),
        % add process triple [l,z,y,a,x]
        l = min(k,m) + 1;
        id = Add_ProcessTriple(l,z,y,a,x);  % also updates nT
      end;
    end;  % for id
    
    % and reverse match on [y,z,q] in D with g-x in G
    for id = IndD(y,z,:)
      if (id == 0), break; end;
      % get entry 
      % [m,a,b,c,u] = D(id,:);  % a=x,b=z
      m = D(id,1);
      a = D(id,2);
      b = D(id,3);
      c = D(id,4);
      u = D(id,5);
      if (debug > 1), 
        % (verify a==y,b==z .. or if possibly c==y?
        if (b ~= z) || ((a ~= y) && (c ~= y)),
          fprintf('ERROR 3.Col non-matching triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] ',id,m,a,b,c,u,k,x,z,y,q); 
        else
          fprintf('*** 3.Col processing triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] ',id,m,a,b,c,u,k,x,z,y,q); 
        end;
      end;
      % check edge c-x
      if (G(c,x) > 0) && (a == y) && (b == z),
        % add process triple [l,z,y,c,x]
        l = min(k,m) + 1;   % order of triple
        id = Add_ProcessTriple(l,z,x,c,y);  % also updates nT
      end;
      % check reverse for order 0 noncollider and edge a-x
      if (m == 0) && (G(a,x) > 0) && (c == y) && (b == z),
        % add process triple [l,z,x,a,y]
        l = min(k,m) + 1;
        id = Add_ProcessTriple(l,z,x,a,y);  % also updates nT
      end;
    end;  % for id
    
  else
    % process noncollider   [k,x,z,y,q] = T(idxT,:);
    % noncollider x*->z-->y (or x<--z*-y for k==0) ... matches collider(s) in C
    % - x*->z<--q, if y-q in G  => triple [z,q,y] (not already in C or D)
    % or for k==0:
    % - y*->z<-*q, if x-q in G  => triple [z,q,x] (not already in C or D

    % match on [x,z,q] in C with q-y in G
    for id = IndC(x,z,:)
      if (id == 0), break; end;
      % get entry 
      % [m,a,b,c,u] = C(id,:);  % a=x,b=z
      m = C(id,1);
      a = C(id,2);
      b = C(id,3);
      c = C(id,4);
      u = C(id,5);
      % check match 
      if (a == x) && (b == z), 
        if (G(c,y) > 0),
          % add process triple [l,z,y,c,x]
          l = min(k,m) + 1;   % order of triple
          id = Add_ProcessTriple(l,z,c,y,x);  % also updates nT
        end;
      elseif (c == x) && (b == z),
        if (G(a,y) > 0),
          % add process triple [l,z,y,c,x]
          l = min(k,m) + 1;   % order of triple
          id = Add_ProcessTriple(l,z,a,y,x);  % also updates nT
        end;
      else
        fprintf('ERROR 3.Noncol non-matching triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] ',id,m,a,b,c,u,k,x,z,y,q); 
      end;
      if (debug > 1), 
        fprintf('*** 3.Noncol processing triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] ',id,m,a,b,c,u,k,x,z,y,q); 
      end;         
    end;  % for id    
    
    % and reverse for k==0
    if (k == 0),
      % match on [y,z,q] in C with q-x in G
      for id = IndC(y,z,:)
        if (id == 0), break; end;
        % get entry 
        % [m,a,b,c,u] = C(id,:);  % a=x/y,b=z
        m = C(id,1);
        a = C(id,2);
        b = C(id,3);
        c = C(id,4);
        u = C(id,5);
        % check match 
        if (a == y) && (b == z), 
          if (G(c,x) > 0),
            % add process triple [l,z,y,c,x]
            l = min(k,m) + 1;   % order of triple
            id = Add_ProcessTriple(l,z,c,x,y);  % also updates nT
          end;
        elseif (c == y) && (b == z),
          if (G(a,x) > 0),
            % add process triple [l,z,y,c,x]
            l = min(k,m) + 1;   % order of triple
            id = Add_ProcessTriple(l,z,a,x,y);  % also updates nT
          end;
        else
          fprintf('ERROR 3.Noncol non-matching reverse triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] ',id,m,a,b,c,u,k,x,z,y,q); 
        end;
        if (debug > 1), 
          fprintf('*** 3.Noncol processing reverse triple id=%d [%d,%d,%d,%d,%d] on proc. [%d,%d,%d,%d,%d] ',id,m,a,b,c,u,k,x,z,y,q); 
        end;         
      end;  % for id    
    end;  % if (k==0)
  end;  % if iscol .. else
  
end;  % while

% 4: finalise and return
% clean up superfluous entries from C / D
C(nC+1:end,:) = [];
D(nD+1:end,:) = [];

% assign output and return
MEC.S = S;
MEC.C = C;
MEC.D = D;
if (debug > 2), fprintf('mag_to_mec finished: nC = %d, nD = %d',nC,nD); end;
return;

% =========================================================================
% helper routines for list management (update list/index/counters/logging)
function [id,isNew] = Add_TripleWithOrder(k,a,b,c,q,isCol)
  id = 0;
  % check if already present (in index list, avoids searching through C/D)
  isNew = true;
  if isCol
    % check C
    IDs = IndC(a,b,:); 
    for i = IDs
      if (i == 0), break; end;  % zero indicates no more filled ids
      if (C(i,4) == c) || (C(i,4) == a) && (C(i,2) == c),  % usually false, but not always
        isNew = false;
        id = i;
        % (no need to check order)
      end;
    end;  % for i
  else
    % isnoncol
    % check D
    IDs = IndD(a,b,:); 
    for i = IDs
      if (i == 0), break; end;  % zero indicates no more filled ids
      if (D(i,4) == c) || (D(i,4) == a) && (D(i,2) == c),  % usually false, but not always
        isNew = false;
        id = i;
        % (no need to check order)
      end;
    end;  % for i
  end;
  if ~(isNew), return; end;
  
  % process new entry: add to C/D list and corresponding index list
  if (isCol),    
    nC = nC + 1;
    id = nC; % return id of entry 
    % extend triple list? (preallocate next block of triples)
    if (nC > size(C,1)),
      C(end+BlockSize,end) = 0;
    end;
    % add triple to collider list
    C(nC,:) = [k,a,b,c,q];
    % add to index (3D matrix)
    % check if we need to extend IndC array
    if (IndC(a,b,end)>0) || ((k==0) && (IndC(c,b,end)>0) )
      % too long? extend array with zeros
      fprintf('WARNING AddTriple: IndC(%d,%d(,%d)) order %k full!? .. (double size)\n',a,b,c,k); 
      IndC(a,b,2*size(IndC,3)) = 0;    % extends for all entries ... but should not be
    end;
    % add entry
    n = find(IndC(a,b,:) == 0,1);   % find next unused entry 
    IndC(a,b,n) = nC;
    % also add reverse for zero order
    if (k == 0),
      n = find(IndC(c,b,:) == 0,1); 
      IndC(c,b,n) = nC;
    end;
    
    % log message
    if (debug > 1), fprintf('Add COLLIDER order %d: %d *-> %d <-* %d (via %d)\n',k,a,b,c,q); end; 
  else
    % add noncollider
    nD = nD + 1;
    id = nD; % return index of entry ?
    % extend triple list? (preallocate next block of triples)
    if (nD > size(D,1)),
      D(end+BlockSize,end) = 0;
    end;
    % add triple to noncollider list
    D(nD,:) = [k,a,b,c,q];

    % add to index (3D matrix)
    % check if we need to extend IndC array
    if (IndD(a,b,end)>0) || ((k==0) && (IndD(c,b,end)>0) )
      % too long? extend array with zeros
      fprintf('WARNING AddTriple: IndD(%d,%d(,%d)) order %k full!? .. (double size)\n',a,b,c,k); 
      IndD(a,b,2*size(IndD,3)) = 0;    % extends for all entries ... but should not be
    end;
    % add entry
    n = find(IndD(a,b,:) == 0,1);
    IndD(a,b,n) = nD;
    % also add reverse for zero order (NOT for all noncolliders)
    if (k == 0),
      n = find(IndD(c,b,:) == 0,1);
      IndD(c,b,n) = nD;
    end;
    if (debug > 1), fprintf('Add NONCOLLIDER order %d: %d *-- %d --* %d (via %d)\n',k,a,b,c,q); end; 
  end;  % if (iscol) else ..
end % add_triplewithorder

function id = Add_ProcessTriple(k,a,b,c,q)
  nT = nT + 1;
  id = nT;
  % extend process list T? (preallocate next block of triples)
  if (nT > size(T,1)),
    T(end+BlockSize,end) = 0;
  end;
  % add triple to process list
  T(nT,:) = [k,a,b,c,q];  
  % log message
  if (debug > 1), fprintf('Add process triple order %d: %d *-> %d <-* %d (via %d)\n',k,a,b,c,q); end; 
end  % function add_processtriple

end % function mag_to_mec

