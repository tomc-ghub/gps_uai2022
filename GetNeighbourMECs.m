function [MECS,nCounts] = GetNeighbourMECs(M,Operators,MAGFlags,P)
% Generates collection of neighbouring MECs of MEC M/PAG P, using either
% basic (=1) or extended (=2) add/delete/mknoncol/mkcol operators.
% Input:
% - M = MEC encoded as struct with fields
%   .S = skeleton as NxN 0/1 adjacency matrix 
%   .C = collider triples with order (i,:) = [k,x,z,y,q]
%     k: order 
%     x,z,y: triple for collider z in vstructure or on discr.path to y
%     q: triggering triple node (given [q,x,z] in C, [q,x,y] in D)
%   .D = noncollider triples with order (idem)
% - Operators = {0:3} flags for which neighbours to consider, as
%              [edge add, edge del, col -> noncol, noncol -> col]   
%              where 1=basic (one candidate),2=extended
% - MAGFlags  = 0/1 flags for allowed invariant edge types in MAGS 
%              [allow bidirected, allow undirected]   
%   note: used to avoid scoring undirected components in GSMAG score 
% Output:
% - MECS = array of structs of neighbouring MECS with extra fields
%   .S/C/D = skeleton/coll./noncoll with order
%   .R = reachability matrix
%   .P = corresponding PAG
%   .G = arc augmented MAG instance
%   .action = [operator,k,x,z,y,q]
%    - operator: 1=add edge, 2=del edge, 3=make noncol, 4=make col 
%    - [k,x,z,y,q]: edge x-y, or (output) triple with order 
% - nCounts = struct with counters on 
%   .nOper[1,4]  : nr.oper[add,delete,mkNoncol,mkCol] tried, 
%   .nMECs[9,4]  : nr.MECs[add,delete,mkNoncol,mkCol] tried/valid per operator
%      (row) 1=operator tried, 2=cand. valid, 3=cand. invalid, 4=mag valid, 
%            5=mag invalid, 6=pag dupl, 7=pag add, 8=limit cand, 9=max.cand
% ======================================================================
% Set local debug level, either via global, or overrule here:
% => 0 = no logging, 1 = basic, 2 = extended/detailed
global DEBUG
if ~isempty(DEBUG), debug = DEBUG; else debug = 0; end;

  % 1: Initialise
  if (debug > 1), disp('Function GetneighbourMECS()'); end; 
  if (nargin < 1), return; end;
  if (nargin < 2), Operators = [1,1,1,1]; end;  % default all neighbours
  % output
  MECS = [];
  nCounts.nMECs = zeros(9,4);   

  % define constants for readability later
  ARC = 1; TAIL = 2;    % for cpag_to_mag
  % operator index (also into nCounts.nMECs tracker)
  opADD = 1;
  opDEL = 2;
  opMNC = 3;
  opMCL = 4;
  MAX_CAND = 128;   % maximum nr. of candidates to consider per operator
  % collect operator levels (0=off, 1=baseline, 2=extended (all versions)
  ADD       = Operators(opADD);
  DELETE    = Operators(opDEL);
  MK_NONCOL = Operators(opMNC);
  MK_COL    = Operators(opMCL);

  % get/set variables used
  nM = 0; 
  S  = M.S; C = M.C; D = M.D;    % extract {S,C,D} from M
  N  = size(S,1);
  nC = size(C,1);
  nD = size(D,1);

  % preprocessing 
  taMAG = cpag_to_mag(P,TAIL);    % tail augmented MAG for Add/DeleteEdge
  % ancestors in taMAG/PAG
  R_taMAG = reachability_graph(taMAG,1);  % ancestor in MAG (include self?)
  rP = P; rP(P == 3) = 0;     % blank circles in PAG (only directed paths)
  R_P = reachability_graph(rP,1);  % ancestor in PAG, include self

  % 2: generate & collect core PAG of all candidate MEC versions (cMECS)
  % preallocate N^2 elements (for starters)  
  cMECS(N*N,1).P = [];  % struct('P',[],'G',[],'action',zeros(1,5),'version',0);  
  nCM = 0;  % index into (still empty) array of candidate MECs 
  
  % 2a: Operator 1 - add edges
  if (ADD > 0), 
    if (debug > 1), disp('Operator 1: edge add'); end;
    % loop over all pairs of nodes (alt: use [X,Y] = find(S == 0))
    for i = 1:(N-1)
      for j = (i+1):N
        % check if pair {i,j} needs to be processed
        if (S(i,j) == 0)
          % update nr. operators tried (1,:)
          nCounts.nMECs(1,opADD) = nCounts.nMECs(1,opADD) + 1; 
          % for one(basic) / four (extended) possible edge types in turn
          if (ADD == 1),
            % 2a.1: type 1, add i o-o j as two noncolliders 
            nCM = nCM + 1; version = 1;
            cMECS(nCM).G = taMAG; % start from tail-augmented MAG
            % count invariant arrowheads at i, j in P
            nrAH_i = sum(P(i,:) == 2);
            nrAH_j = sum(P(j,:) == 2);
            if (nrAH_i == 0) && (nrAH_j == 0)
              cMECS(nCM).G(i,j) = 1; cMECS(nCM).G(j,i) = 1; 
              cMECS(nCM).action = [1,0,i,0,j,1];  % added edge i-j, type 1
              cMECS(nCM).version = version;      
            elseif (R_taMAG(i,j) > 0) || ((nrAH_i == 0) && (nrAH_j > 0))
              % path from i to j or no arh at i, but j <-*
              cMECS(nCM).G(i,j) = 1; cMECS(nCM).G(j,i) = 2; 
              cMECS(nCM).action = [1,0,i,0,j,2];  % added edge i->j, type 2
              cMECS(nCM).version = version;      
            elseif (R_taMAG(j,i) > 0) || ((nrAH_i > 0) && (nrAH_j == 0))
              cMECS(nCM).G(i,j) = 2; cMECS(nCM).G(j,i) = 1; 
              cMECS(nCM).action = [1,0,i,0,j,3];  % added edge i<-j, type 3
              cMECS(nCM).version = version; 
            elseif (R_taMAG(i,j) == 0) &&  (R_taMAG(j,i) == 0)
              cMECS(nCM).G(i,j) = 2; cMECS(nCM).G(j,i) = 2; 
              cMECS(nCM).action = [1,0,i,0,j,4];  % added edge i<->j, type 4
              cMECS(nCM).version = version; 
            else
              disp('ERROR: invalid ADD??');     % should not occur
            end;
            % update nr. candidates valid (2,:)
            nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 

          else % if ADD > 1
            % extended: multiple edge versions + alternatives for circle marks
            version = 0;
            % 2a.1: type 2, add i --- j if no arrowheads
            valid = isempty(find(P(i,:) == 2 | P(j,:) == 2,1));
            if (valid)
              % start from original core PAG for M with added edge i --- j
              G = P; 
              G(i,j) = 1; G(j,i) = 1;
              % add base version for this edge to candidate PAGs
              nCM = nCM + 1; version = version + 1;
              cMECS(nCM).G = G;
              cMECS(nCM).action = [1,0,i,0,j,1];  % added edge i-j, type 1
              cMECS(nCM).version = version;
              % update nr. candidates valid (2,:)
              nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 
            else
              % invalid: update nr. candidates invalid (3,:) and skip
              nCounts.nMECs(3,opADD) = nCounts.nMECs(3,opADD) + 1;               
            end;
            
            % 2a.2: type 2, add i <-- j as arc
            % start from original core PAG for M with added edge i <-- j
            G = P; 
            % validate: no undirected edge at i and no almost directed cycle?
            valid = isempty(find(P(i,:) == 1 & P(:,i)' == 1,1));
            if valid
              % check no almost directed cycle in P (Lemma 7)
              % (actually: no protected ancestor, but ok)
              Z = find(R_P(:,j)' > 0);  % all ancestors of j in P
              W = find(R_P(i,:)  > 0);  % all descendants of i in P
              % note: no edge i - j in P, so no need to exclude  via copy as in step 2c (make-noncol)
              if ~isempty(find(P(Z,W) == 2,1))
                valid = false;
              end;
            end;
            if (valid)
              % if ok then process: add
              G(i,j) = 2; G(j,i) = 1;
              % add base version for this edge to candidate PAGs
              nCM = nCM + 1; version = version + 1;
              cMECS(nCM).G = G;
              cMECS(nCM).action = [1,0,i,0,j,2];  % added edge i-j, type 2
              cMECS(nCM).version = version;
              % update nr. candidates valid (2,:)
              nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 

              % consider alternatives for all circle marks at i in G
              % note: none are on an undirected edge in P, so no need to check!
              U = find( G(i,:) == 3);
              nU = length(U);
              nCand = 2^nU;
              % track max nr. of candidates and limit if needed
              if (nCand > nCounts.nMECs(9,opADD)), nCounts.nMECs(9,opADD) = nCand; end;
              if (nCand > MAX_CAND),
                if (debug > 1),fprintf('Limit too many candidates for ADD-U: %d > %d\n',nCand, MAX_CAND); end;
                nCand = MAX_CAND;
                % update nr. of times nr. of candidates were limited (8,:)
                nCounts.nMECs(8,opADD) = nCounts.nMECs(8,opADD) + 1; 
              end;
              % edgemarks = ind2subv(2*ones(1,nU), 1:2^nU);   % generate all combis for U nodes
              edgemarks = ind2subv(2*ones(1,nU), 1:nCand);   % generate (maxed) combis for U nodes
              % add all combinations as separate candidate PAG
              for k = 2:(nCand)      % note: k==1 (all tails) == base version above
                nCM = nCM + 1; version = version + 1;
                cMECS(nCM).G = G;
                % modify in PAG 
                cMECS(nCM).G(i,U) = edgemarks(k,:);   % new combi of tails/arrowheads
                cMECS(nCM).action = [1,0,i,0,j,3];  % type 3
                cMECS(nCM).version = version;
                % update nr. candidates valid (2,:)
                nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 
              end;  % for k
            else  
              % invalid: update nr. candidates invalid (3,:) and skip
              nCounts.nMECs(3,opADD) = nCounts.nMECs(3,opADD) + 1;               
            end;  % if valid

            % 2a.3: type 3, add i --> j as arc
            % start from original PAG for M with added edge i --> j
            G = P; 
            % validate: no undirected edge at j and no almost directed cycle?
            valid = isempty(find(P(j,:) == 1 & P(:,j)' == 1,1));
            if valid
              % check no almost directed cycle in P (Lemma 7)
              % (actually: no protected ancestor, but ok)
              Z = find(R_P(:,i)' > 0);  % all ancestors of i in P
              W = find(R_P(j,:)  > 0);  % all descendants of j in P
              if ~isempty(find(P(Z,W) == 2,1))
                valid = false;
              end;
            end;
            if (valid)
              G(i,j) = 1; G(j,i) = 2;
              % add base version for this edge to candidate PAGs
              nCM = nCM + 1; version = version + 1;
              cMECS(nCM).G = G;
              cMECS(nCM).action = [1,0,i,0,j,3];  % added edge i-j, type 3
              cMECS(nCM).version = version;   
              % update nr. candidates valid (2,:)
              nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 
              % consider alternatives for all circle marks at j 
              % now find all circle marks at j in 
              V = find( G(j,:) == 3);
              nV = length(V);
              nCand = 2^nV;
              % track max nr. of candidates and limit if needed
              if (nCand > nCounts.nMECs(9,opADD)), nCounts.nMECs(9,opADD) = nCand; end;
              if (nCand > MAX_CAND),
                if (debug > 1),fprintf('Limit too many candidates for ADD-V: %d > %d\n',nCand, MAX_CAND); end;
                nCand = MAX_CAND;
                % update nr. of times nr. of candidates were limited (8,:)
                nCounts.nMECs(8,opADD) = nCounts.nMECs(8,opADD) + 1; 
              end;
              %edgemarks = ind2subv(2*ones(1,nV), 1:2^nV);   % generate all combis for V nodes
              edgemarks = ind2subv(2*ones(1,nV), 1:nCand);   % generate all combis for V nodes
              % add all combinations as separate candidate PAG
              for k = 2:nCand      % note: k==1 (all tails) == base version above
                nCM = nCM + 1; version = version + 1;
                cMECS(nCM).G = G;
                % modify marks 
                cMECS(nCM).G(j,V) = edgemarks(k,:);   % new combi of tails/arrowheads
                cMECS(nCM).action = [1,0,i,0,j,3];  % type 3
                cMECS(nCM).version = version;
                % update nr. candidates valid (2,:)
                nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 
              end;  % for k
            else  
              % invalid: update nr. candidates invalid (3,:) and skip
              nCounts.nMECs(3,opADD) = nCounts.nMECs(3,opADD) + 1;               
            end;  % if valid

            % 2a.4: type 4, add i <-> j as bidirected arc
            % start from original PAG for M with added edge i <-> j
            G = P; 
            % validate: no undirected edge at i or j and no directed path between them?
            valid = isempty(find( (P(i,:) == 1 & P(:,i)' == 1) | (P(j,:) == 1 & P(:,j)' == 1),1));
            if valid
              % check no directed path in P (Lemma 7)
              % (actually: no protected ancestor, but ok)
              if (R_P(i,j) > 0) || (R_P(j,i) > 0)
                valid = false;
              end;
            end;
            if (valid)
              G(i,j) = 2; G(j,i) = 2;
              % add base version for this edge to candidate PAGs
              nCM = nCM + 1; version = version + 1;
              cMECS(nCM).G = G;
              cMECS(nCM).action = [1,0,i,0,j,4];  % added edge i-j, type 4
              cMECS(nCM).version = version;   
              % update nr. candidates valid (2,:)
              nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 
              % consider alternatives for all circle marks at i AND j 
              % now do the same for all combinations of U and V 
              % WARNING: potentially expensive! (ok in practice)
              U = find( G(i,:) == 3);   % recompute, as stage 2 may not have run invalid)
              V = find( G(j,:) == 3);   % idem stage 3
              nU = length(U);
              nV = length(V);
              W = [U,V];
              nW = nU + nV;
              nCand = 2^nW;
              % track max nr. of candidates and limit if needed
              if (nCand > nCounts.nMECs(9,opADD)), nCounts.nMECs(9,opADD) = nCand; end;
              if (nCand > MAX_CAND),
                if (debug > 1),fprintf('Limit too many candidates for ADD-W: %d > %d\n',nCand, MAX_CAND); end;
                nCand = MAX_CAND;
                % update nr. of times nr. of candidates were limited (8,:)
                nCounts.nMECs(8,opADD) = nCounts.nMECs(8,opADD) + 1; 
              end;            
              % edgemarks = ind2subv(2*ones(1,nW), 1:2^nW);   % generate all combis for W nodes
              edgemarks = ind2subv(2*ones(1,nW), 1:nCand);   % generate all combis for W nodes
              % add all combinations as separate candidate PAG
              % NOTE: this can be a killer ... avoid by limit
              for k = 2:nCand      % note: k==1 (all tails) == base version above
                nCM = nCM + 1; version = version + 1;
                cMECS(nCM).G = G;
                % modify in core PAG 
                cMECS(nCM).G(i,U) = edgemarks(k,[1:nU]);   % new combi of tails/arrowheads
                cMECS(nCM).G(j,V) = edgemarks(k,[nU+1:nW]);   % new combi of tails/arrowheads
                cMECS(nCM).action = [1,0,i,0,j,4];  % type 3
                cMECS(nCM).version = version;
                % update nr. candidates valid (2,:)
                nCounts.nMECs(2,opADD) = nCounts.nMECs(2,opADD) + 1; 
              end;  % for k
            else  
              % invalid: update nr. candidates invalid (3,:) and skip
              nCounts.nMECs(3,opADD) = nCounts.nMECs(3,opADD) + 1;               
            end;  % if valid
              
          end;  % if ADD==1, else ..
        end;  % if process {i,j}
      end; % for j
    end; % for i
  end;  % operator 1 - edge add

  % 2b: Operator 2 - delete edges
  if (DELETE > 0), 
    if (debug > 1), disp('Operator 2: edge delete'); end; 
    % loop over all pairs of nodes (alt: use [X,Y] = find(S ~= 0))
    for i = 1:(N-1)
      for j = (i+1):N
        % check if pair {i,j} needs to be processed
        if (S(i,j) ~= 0) 
          % update nr. operators tried (1,:)
          nCounts.nMECs(1,opDEL) = nCounts.nMECs(1,opDEL) + 1; 
          % delete edge i *-* j 
          if (DELETE == 1), 
            nCM = nCM + 1;
            cMECS(nCM).G = taMAG; % start from tail-augmented MAG
            % delete edge i *-* j 
            cMECS(nCM).G(i,j) = 0; 
            cMECS(nCM).G(j,i) = 0; 
            cMECS(nCM).action = [2,0,i,0,j,0];  % delete edge i-j, 
            cMECS(nCM).version = 1;  
            % update nr. candidates valid (2,:)
            nCounts.nMECs(2,opDEL) = nCounts.nMECs(2,opDEL) + 1; 
          
          else  % DELETE > 1
            % extended: consider alternatives for all circle marks at i or j 
            G = P;      % extended starts from PAG, NOT taMAG!!

            % NOTE: there is no a priori straightforward sanity check for
            % validity, as all circle marks can by definition appear in a
            % valid MAG as either tail or arrowhead, and in that MAG we can
            % always remove any edge to still have a valid MAG. So for now
            % => NO SEPARATE VALIDITY CHECK FOR THE DELETE OPERATOR
            
            G(i,j) = 0; G(j,i) = 0; % delete edge
            % find all potentially ambiguous zero order triples
            % meaning either circle+arrowhead or circle+circle
            U = find( (G(:,i) == 3 & G(:,j) >= 2) | (G(:,i) >= 2 & G(:,j) == 3) );
            nU = length(U);
            % find all combinations of possible colliders 
            nCand = 2^nU;
            % track max nr. of candidates and limit if needed
            if (nCand > nCounts.nMECs(9,opDEL)), nCounts.nMECs(9,opDEL) = nCand; end;
            if (nCand > MAX_CAND),
              if (debug > 1),fprintf('Limit too many candidates for ADD-U: %d > %d\n',nCand, MAX_CAND); end;
              nCand = MAX_CAND;
              % update nr. of times nr. of candidates were limited (8,:)
              nCounts.nMECs(8,opDEL) = nCounts.nMECs(8,opDEL) + 1; 
            end;
            edgemarks = ind2subv(2*ones(1,nU), 1:nCand);   % generate (maxed) combis for U nodes
            for k = 1:nCand
              V = U(edgemarks(k,:) == 2);   % get new collider nodes
              % potentially collider: add explicit collider instance
              nCM = nCM + 1;
              cMECS(nCM).G = G;
              % add colliders i *-> {V} <-* j for all nodes in V
              if ~isempty(V),
                cMECS(nCM).G(V,[i,j]) = 2; 
              end;
              cMECS(nCM).action = [2,0,i,0,j,0];  % delete edge i-j, 
              cMECS(nCM).version = k;      
              % update nr. candidates valid (2,:)
              nCounts.nMECs(2,opDEL) = nCounts.nMECs(2,opDEL) + 1; 
            end;
          end;  % if delete=1, else          
        end;  % if process {i,j}
      end; % for j
    end; % for i
  end;  % operator 2 - edge delete  
  
  % 2c: Operator 3: make noncollider
  if (MK_NONCOL > 0), 
    if (debug > 1), disp('Operators 3: make noncol'); end; 
    P2 = P;     % create copy of P for validation purposes (Lemma 7)
    % count all entries in lists C/D to process
    nC = size(C,1); 
    % update nr. operators tried (1,:)
    nCounts.nMECs(1,opMNC) = nCounts.nMECs(1,opMNC) + nC; 
    % loop over all triples in C
    for i = 1:nC
      % get collider triple C(i) to turn into noncollider
      k = C(i,1);
      x = C(i,2);
      z = C(i,3);
      y = C(i,4);
      q = C(i,5);
      version = 0;
      
      % try z --* y (always, i.e. for all orders 'k')
      % validate: do no created undirected edge z --- y (given arrowhead x *-> z)   
      %valid = (length(find( P(z,:) == 2 )) > 2) && (P(y,z) ~= 1);
      valid = (P(y,z) ~= 1);    % NOTE: corrected superfluous check
      if valid
        % check y <-- z not (almost) directed cycle in (new) P (Lemma 7)
        % (actually: no protected ancestor, but ok)
        Z = find(R_P(:,z)' > 0);  % all ancestors of z in P
        W = find(R_P(y,:)  > 0);  % all descendants of y in P
        % exclude mark z <-* y in P (modify copy P2, check, and restore)
        P2(z,y) = 1;
        if ~isempty(find(P2(Z,W) == 2,1))
          valid = false;
        end;
        P2(z,y) = 2; % restore
      end;
      if (valid)
        % start again from original PAG for M
        nCM = nCM + 1; version = version + 1;
        cMECS(nCM).G = P;
        % modify to noncollider in PAG (note: for k=0 three options!
        cMECS(nCM).G(z,y) = 1; % NOTE: corrected from 'P(z,y)' in original GPS code
        cMECS(nCM).action = [3,k,x,z,y,q];  % try mkNoncol 
        cMECS(nCM).version = version;      
        % update nr. candidates valid (2,:)
        nCounts.nMECs(2,opMNC) = nCounts.nMECs(2,opMNC) + 1; 
      else  
        % invalid: update nr. candidates invalid (3,:) and skip
        nCounts.nMECs(3,opMNC) = nCounts.nMECs(3,opMNC) + 1;               
      end;  % if valid
        
      % try x *-- z (only for k == 0, for higher order triples z --> y is necessary)
      % validate: do no created undirected edge x --- z (given arrowhead z <-* y)   
      %valid = (length(find( P(z,:) == 2 )) > 2) && (P(x,z) ~= 1);
      valid = ((k == 0) && (P(x,z) ~= 1));
      if valid
        % check x <-- z not (almost) directed cycle in (new) P (Lemma 7)
        % (actually: no protected ancestor, but ok)
        Z = find(R_P(:,z)' > 0);  % all ancestors of z in P
        W = find(R_P(x,:)  > 0);  % all descendants of x in P
        % exclude mark z <-* x in P (modify copy P2, check, and restore)
        P2(z,x) = 1;
        if ~isempty(find(P2(Z,W) == 2,1))
          valid = false;
        end;
        P2(z,x) = 2; % restore
      end;
      if (valid)
        nCM = nCM + 1; version = version + 1;
        cMECS(nCM).G = P;
        cMECS(nCM).G(z,x) = 1; 
        cMECS(nCM).action = [3,k,x,z,y,q];  % try extra mkNoncol 
        cMECS(nCM).version = version;       % could still be identical to above
        % update nr. candidates valid (2,:)
        nCounts.nMECs(2,opMNC) = nCounts.nMECs(2,opMNC) + 1; 
      else  
        % invalid: update nr. candidates invalid (3,:) and skip
        nCounts.nMECs(3,opMNC) = nCounts.nMECs(3,opMNC) + 1;               
      end;  % if valid
        
      % and try x *-- z --* y 
      % multiple validation steps
      valid = true;
      if (P(x,z) == 1),
        % will result in x --- z, so check no other arrowheads at x or z
        if ~isempty(find(P(x,:) == 2,1)) || (length(find(P(z,:) == 2)) > 2),
          valid = false;
        end;
      end;
      if valid && (P(y,z) == 1),
        % will result in y --- z, so check no other arrowheads at y or z
        if ~isempty(find(P(y,:) == 2,1)) || (length(find(P(z,:) == 2)) > 2),
          valid = false;
        end;
      end;
      if valid && ( P(x,z) == 2 || P(y,z) == 2 ) % if directed edge created ..
        % .. then check for (almost) directed paths
        rG = P;  % extra copy for reachability graph (now two edges have changed)
        rG(z,[x,y]) = 1;
        rG(rG == 3) = 0;     % blank circles in modified reachability graph rG
        R_G = reachability_graph(rG,1);  % def.ancestors in modified rG
        if (P(x,z) == 2)  % only check if needed (directed edge x <-- z created)
          % check x <-- z not (almost) directed cycle in (new) G (Lemma 7)
          % (actually: no protected ancestor, but ok)
          Z = find(R_G(:,z)' > 0);  % all def.ancestors of z in P
          W = find(R_G(x,:)  > 0);  % all def.descendants of x in P
          % exclude marks z <-* x/y in P (modify copy P2, check, and restore)
          P2(z,[x,y]) = 1;
          if ~isempty(find(P2(Z,W) == 2,1))
            valid = false;
          end;
          P2(z,[x,y]) = 2; % restore
        end;
        if valid && (P(y,z) == 2)
          % check y <-- z not (almost) directed cycle in (new) G (Lemma 7)
          % (actually: no protected ancestor, but ok)
          Z = find(R_G(:,z)' > 0);  % all ancestors of z in P
          W = find(R_G(y,:)  > 0);  % all descendants of y in P
          % exclude marks z <-* x/y in P (modify copy P2, check, and restore)
          P2(z,[x,y]) = 1;
          if ~isempty(find(P2(Z,W) == 2,1))
            valid = false;
          end;
          P2(z,[x,y]) = 2; % restore
        end;
      end;  % validation complete
      % add candidate if valid
      if (valid)
        nCM = nCM + 1; version = version + 1;
        cMECS(nCM).G = P;
        cMECS(nCM).G(z,[x,y]) = 1; 
        cMECS(nCM).action = [3,k,x,z,y,q];  % try extra mkNoncol 
        cMECS(nCM).version = version;       % could still be identical to above
        % update nr. candidates valid (2,:)
        nCounts.nMECs(2,opMNC) = nCounts.nMECs(2,opMNC) + 1; 
      else  
        % invalid: update nr. candidates invalid (3,:) and skip
        nCounts.nMECs(3,opMNC) = nCounts.nMECs(3,opMNC) + 1;               
      end;  % if valid
      
      % extended version? not needed, as no new colliders can be created!
    end;  % for i=1:nC
  end;  % operators 3 - make noncollider

  % 2d: Operator 4: make collider
  if (MK_COL > 0), 
    if (debug > 1), disp('Operator 4: make collider'); end; 
    % count all entries in list D to process
    nD = size(D,1);
    % update nr. operators tried (1,:)
    nCounts.nMECs(1,opMCL) = nCounts.nMECs(1,opMCL) + nD; 
    % further preallocate? ... not really necessary, but maybe    
    % loop over all triples in D
    for i = 1:nD
      % get noncollider of triple D(i)
      k = D(i,1);
      x = D(i,2);
      z = D(i,3);
      y = D(i,4);
      q = D(i,5);
      version = 0;
      % validate
      % no other node with undirected edge to z
      U = (P(z,:) == 1) & (P(:,z)' == 1);
      U([x,y]) = 0;     % except for {x,y}
      valid = isempty(find(U > 0,1));
      if valid && ( P(z,x) == 1 || P(z,y) == 1 )
        % actually: check z not protected ancestor of x or y in P
        % implies z--*x in P , and directed path z to x in G'
        rG = P;
        rG(z,[x,y]) = 2; % add arrowheads x->z<-y
        % find reachability graph
        rG(rG == 3) = 0;
        R_G = reachability_graph(rG);  % exclude self (does not matter here)
        if (R_G(z,x) > 0) || (R_G(z,y) > 0)
          valid = false;
        end;
      end;  % validation
      
      % add candidate if valid
      if (valid)
        nCM = nCM + 1; version = version + 1;
        % start again from original PAG for M
        cMECS(nCM).G = P;
        % modify to collider 
        cMECS(nCM).G(z,[x,y]) = 2; % add arrowheads x->z<-y
        cMECS(nCM).action = [4,k,x,z,y,q];  % try mkCol 
        cMECS(nCM).version = version;
        % update nr. candidates valid (2,:)
        nCounts.nMECs(2,opMCL) = nCounts.nMECs(2,opMCL) + 1; 
      else  
        % invalid: update nr. candidates invalid (3,:) and skip
        nCounts.nMECs(3,opMCL) = nCounts.nMECs(3,opMCL) + 1;               
      end;  % if valid
      
      if (MK_COL == 1) || ~(valid), continue; end;
      
      % extended: check all other (noncol) triples with a shared edge at z-[x,y]
      % NOTE: no added validity checks at extended stage
      G = cMECS(nCM).G;
      % find all circle marks at z
      U = find( G(z,:) == 3);
      nU = length(U);
      nCand = 2^nU;
      % track max nr. of candidates and limit if needed
      if (nCand > nCounts.nMECs(9,opMCL)), nCounts.nMECs(9,opMCL) = nCand; end;
      if (nCand > MAX_CAND),
        if (debug > 1),fprintf('Limit too many candidates for MCL-U: %d > %d\n',nCand, MAX_CAND); end;
        nCand = MAX_CAND;
        % update nr. of times nr. of candidates were limited (8,:)
        nCounts.nMECs(8,opMCL) = nCounts.nMECs(8,opMCL) + 1; 
      end;
      % edgemarks = ind2subv(2*ones(1,nU), 1:2^nU);   % generate all combis for U nodes
      edgemarks = ind2subv(2*ones(1,nU), 1:nCand);   % generate all combis for U nodes
      % add all combinations as separate candidate PAG (up to 2^|maxK-2|)
      for j = 2:nCand       % first = all tails == circle marks, so skip
        nCM = nCM + 1; version = version + 1;
        cMECS(nCM).G = P;
        % modify in G
        cMECS(nCM).G(z,U) = edgemarks(j,:);   % new combi of tails/arrowheads
        cMECS(nCM).action = [4,k,x,z,y,q];  % try mkCol 
        cMECS(nCM).version = version;               
        % update nr. candidates valid (2,:)
        nCounts.nMECs(2,opMCL) = nCounts.nMECs(2,opMCL) + 1; 
      end;  % for j
    end;  % for i
  end;  % operators 4 - make collider
  % ======================================================================  
  
  % 3: process candidates and add to neighbour list if valid and new
  MECS(nCM).M = [];         % preallocate for current nr.candidates
  nM = 0;
  for k = 1:nCM
    if (debug > 1), fprintf('Validate candidate MEC %d\n',k); end;
    % get candidate info
    G0      = cMECS(k).G;               % starting modified graph
    version = cMECS(k).version;         % check for multiple PAG instances per action
    oper    = cMECS(k).action(1);       % 1=ADD,2=DEL,3=mkNcol,4=mkCol
    % reset PAG comparison if new version = 1 encountered (not necessary valid)
    if (version == 1), nM_v1 = 0; end;
    
    % convert to PAG(s): G0 -> MEC -> PAG -> MAG-> validate->duplicates
    M = mag_to_mec(G0);         % handles mag/pag_to_mec (undetermined triples become noncollider)
    [P,coreP] = mec_to_cpag(M); % get (core) PAG from MEC (check for duplicates if valid)
    G = cpag_to_mag(P,ARC);   % arc-augmented MAG for validation and scoring
    
    % check validity
    if mag_valid(G,MAGFlags)
      % yes: update counts and store new neighbour entry if not duplicate
      nCounts.nMECs(4,oper) = nCounts.nMECs(4,oper) + 1;    % (4,:) = mag valid
      
      % check duplicate PAGs? 
      if (nM_v1 == 0)
        % new (first) valid version, just store index into MEC collection
        nM_v1 = nM + 1;
      else
        % higher version, so compare to all previous versions
        duplicate = false;
        for m = nM_v1:nM
          % check against previous PAGs
          if isempty(find(MECS(m).P ~= P,1))
            % duplicate PAG, skip
            duplicate = true;
            nCounts.nMECs(6,oper) = nCounts.nMECs(6,oper) + 1;    % (6,:) = pag duplicates            
            break;  % 
          end;        
        end; % for m
        if (duplicate), continue; end; % move to next candidate MEC k+1
      end; % if (version=1), check duplicates
      
      % all good: add to neighbouring MEC collection
      nM = nM + 1;
      MECS(nM).M      = M;                  % MEC
      MECS(nM).G      = G;                  % arc-augmented MAG
      MECS(nM).P      = P;                  % CPAG 
      MECS(nM).coreP  = coreP;              % core P
      MECS(nM).action = cMECS(k).action;    % [1,0,i,0,j,k]; = added edge i-j
      nCounts.nMECs(7,oper) = nCounts.nMECs(7,oper) + 1;    % (7,:) = pag add
      
    else
      % mag invalid: track number of invalid candidates per action
      nCounts.nMECs(5,oper) = nCounts.nMECs(5,oper) + 1;    % (5,:) = mag invalid
      if (debug > 1), fprintf('WARNING: skip invalid candidate MEC %d\n',k); end;
      continue;
    end;  % if (mag_valid)

  end;  % for k

  
  % 4: finalise 
  % remove unused entries (if any)
  if nM < size(MECS,2),
    MECS(nM+1:end) = [];
  end;
  
  % return collection of neighbouring MECS
  if (debug > 1), fprintf('GetNeighbourMECs: return %d new MECs\n',nM); end; 
  return;
    
end % main function getneighbourmecs

