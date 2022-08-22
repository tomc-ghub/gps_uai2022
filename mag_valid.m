function ok = mag_valid(G,MAGFlags)
% check if MAG G is valid, which implies:
% - no arrowheads at nodes with undirected edges
% - no (almost) directed cycles
% optional check MAGFlags[allow_bidrected,allow_undirected]:
% - 0 = not allowed, 1 = allowed

global DEBUG;
if isempty(DEBUG), debug = 0; else debug = DEBUG;  end;

% 1: initialize
ok = false;      % guilty until proven innocent 
if (nargin < 1), return; end; 
if (nargin < 2), MAGFlags = [1,1]; end; % default allow all
N = size(G,1);
ALLOW_BIDIRECTED = (MAGFlags(1) > 0);
ALLOW_UNDIRECTED = (MAGFlags(2) > 0);

% 2: no arrowheads at nodes with undirected edge
% loop over all nodes
for i = 1:N
  % find if it is on an undirected edge
  if ~isempty(find( (G(i,:) == 1) & (G(:,i)' == 1),1))
    % check if it has arrowheads
    if ~ALLOW_UNDIRECTED || ~isempty(find( (G(i,:) == 2),1))
      if (debug > 1), fprintf('MAG invalid: und.edge(+arrowhead) at %d\n',i); end; 
      return;
    end;
  end;
end;  % for i

% 3: no (almost) directed cycles
% remove undirected edges
P = G;
P(G+G' == 2) = 0;

% check bidirected allowed (rarely used)
if ~ALLOW_BIDIRECTED
  [row,col] = find(G == 2 & G' == 2,1);
  if ~isempty(row)
    if (debug > 1), fprintf('WARNING: invalid MAG, bidirected edge %d-%d not allowed\n',row,col); end; 
    return;
  end;
end;

% compute reachability matrix
R = reachability_graph(P);
% loop over all nodes
for i = 1:N
  % check cycle
  if R(i,i) > 0
    if (debug > 1), fprintf('WARNING: invalid MAG, cycle at %d\n',i); end; 
    return;
  end;
  % check almost directed cycle
  % find reachable node connected with bidirected arc
  j = find( (G(i,:) == 2) & (G(:,i)' == 2) & (R(i,:) > 0),1);
  if ~isempty(j)
    if (debug > 1), fprintf('WARNING: invalid MAG, almost directed cycle %d -> %d\n',i,j); end; 
    return;
  end;
end;  % for i

% 4: finish: if we reach here then it was ok
ok = true;
if (debug > 2), disp('mag_valid = ok'); end; 

end % function mag_valid


