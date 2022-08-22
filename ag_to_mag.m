function [M, R, nIndEdge] = ag_to_mag(G)
% ag_to_mag: Make ancestral graph maximal.
% input:
% - G       : ancestral graph (AG)
% output:
% - M       : maximal ancestral graph (MAG)
% - R       : anterior reachability graph
% - nIndEdge: nr. of induced edges (in M, not in G)

global DEBUG
if ~isempty(DEBUG), debug = DEBUG; else debug = 1; end;

  % 1 - Initialize  
  N   = length(G);
  M   = G;
  R   = reachability_graph(G);  % R(i,j) > 0 iff i anterior to j in G
  Sel = []; 
  nIndEdge = 0;
  
  % 2 - loop over all combinations and test for m-separation given ancestors
  for x = 1:(N-1)
    for y = (x+1):N
      % skip nodes already adjacent in G
      if G(x,y) > 0, continue; end;
      % find all ancestors of x and/or y and/or S
      AxyS = mysetdiff(find(sum(R(:,[x,y,Sel]) > 0, 2)' > 0), [x,y]);
      % test for m-separation given all possible ancestors (anterior nodes)
      sep = msep(x,y,[AxyS,Sel],G);   % include S explicitly
      if ~(sep)
        % inducing path from x to y: not separable, so add 
        if (debug > 1), fprintf('Induced edge (%i *-* %i) in ancestral graph\n',x,y); end;
        if R(x,y) > 0, M(x,y) = 1; else M(x,y) = 2; end; % should always be x<->y in M
        if R(y,x) > 0, M(y,x) = 1; else M(y,x) = 2; end;
        nIndEdge = nIndEdge + 1;
      end;
    end;  % for y
  end;  % for x
  
  % return M  
end  % function mk_mag

  
  
  