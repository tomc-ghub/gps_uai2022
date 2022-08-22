function [SHD, wSHD] = GetPAGSHD(truePAG,P,WeightedUnknown)
% compute Structural Hamming Distance between truePAG and PAG P
% SHD = nr.Diff.Marks + nr.Diff.Edge
% - WeightedUnknown: treat 'unknown' for tail/arrowhead (or v.v.) as less
%   severe error than 'wrong'; but not vs. no-edge
%   change: (tail <-1/2-> circle <-1/2-> arrowhead) <-> no-edge
% output:
% - SHD = Structural Hamming Distance 
% - wSHD (optional) = SHD with weighted unknown

  % 1: initialise
  if (nargin < 3), WeightedUnknown = false; end;  % default just raw SHD
  SHD = 0; wSHD = 0;
  
  % 2: compute SHD
  % nr. of different skeleton edges
  trueS = (truePAG ~= 0);
  S = (P ~= 0);
  diff_edge = sum(sum(trueS ~= S)) / 2; % each edge has two marks

  % penalty for different marks
  diff_mark = sum(sum(truePAG ~= P));
  % sum to SHD
  SHD = diff_edge + diff_mark;
    
  % 3: compute weighted version (if needed)
  if (WeightedUnknown)
    % yes: tail/arrow vs. unknown is less bad than tail vs. arrow
    WrongUnknown = sum(sum( ((truePAG == 3) & ( (P==1) | (P==2))) | ...
                            ((P == 3) & ( (truePAG==1) | (truePAG==2))) ...
                          ));
    % different, but only penalized half, so subtract from SHD
    wSHD = SHD - WrongUnknown/2;    
  end;
  
  % 4: done
  return;
end  % function GetPAGSHD
