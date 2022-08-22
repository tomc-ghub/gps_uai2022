function [Score,Comps,CompScores] = GetMAGScore(G,CovMat,nSamples,tol,refG,refComps,refCompScores)
% Compute Gaussian MAG score, adapted from GSMAG (Triantafillou, 2016), and
% greedyBAPs (Nowzohour, 2017)
% Version only recomputing changed c-components/districts
% steps:
% 1: compute districts / subgraphs / subCov
% 2: loop over components
%   - compute score with RICFfit per component
% 3: collect scores and sum
% ======================================================================

  % recompute when no reference graph
  FULL_RECOMP = (nargin < 7);
  tol = max(tol, 1e-8);     % use reasonable threshold for convergence

  % 1: initialize components/scores
  Score = 0;
  Comps = Districts(G);      % NxnD matrix 
  [N,nD] = size(Comps);
  CompScores = zeros(1,nD);

  % 2: loop over all components/districts
  for i = 1:nD
    % get component
    CNodes = find(Comps(:,i) > 0);  % district members+parents
    DNodes = find(Comps(:,i) == 1);    % district members only
    PaNodes = find(Comps(:,i) == 2);   % just parents of district
    % get component MAG
    compMAG = G;                        
    compMAG(PaNodes,PaNodes) = 0;       % blank all interparent links
    compMAG = compMAG(CNodes,CNodes);   % keep subMAG over CNodes
    
    % check if we really need to recompute
    if ~FULL_RECOMP
      % try to find matching district in ref 
      % every node appears exactly once in a district
      % find matching ref comp 
      ref_j = find(refComps(DNodes(1),:) == 1,1);
      if isempty(find(Comps(:,i) ~= refComps(:,ref_j),1)),
        % ok: districts match, now check graph details
        refCNodes = find(refComps(:,ref_j) > 0); % yeah bit of overlap ..
        refPaNodes = find(refComps(:,ref_j) == 2);
        refMAG = refG;
        refMAG(refPaNodes,refPaNodes) = 0; % blank all interparent links
        refMAG = refMAG(refCNodes,refCNodes); % keep subMAG over CNodes
        % compare
        if isempty(find(compMAG ~= refMAG,1))
          % identical district: copy score and go to next one
          CompScores(i) = refCompScores(ref_j);
          continue;     % next district (break to 'for i' loop)
        end;
      end; % if matching comp
    end; % if not(fullrecomp)
    
    % when here: no available match so compute score for this component
    compCov = CovMat(CNodes,CNodes);
    % parents in component (rest is district)
    idxParents = zeros(1,N); 
    idxParents(PaNodes) = 1; 
    idxParents = idxParents(CNodes);  % index of parents in component
    % get score (includes RICF step)
    CompScores(i) = MAGComp_MVGML_Score(compMAG,idxParents,compCov,nSamples,tol);
  end; % for i
  
  % compute overall MAG score (negative so we can minimize the score)
  MLScore = (-nSamples/2)*sum(CompScores);  % neg.log-ML
  % BIC penalty
  nEdges = length(find(G > 0));   
  BICp = bicPenalty(nSamples, N, nEdges); 
  % overall score (alpha = tunable parameter to induce sparsity (rel. ML weight)
  alpha = 2;    % GSMAG = 2, GSBAP = 1 ... lower should make it sparser
  Score = (-alpha*MLScore + BICp)/nSamples;  % scaled, and negative so minimize like penalty
end % function getmagscore

function bp = bicPenalty(nSamples, nVars, nEdges) 
  % note: BIC penalty choice bit arbitrary, Nowz BAP uses (nodes+edges)
  bp = log(nSamples)*(2*nVars+nEdges);    % GSMAG version
  %  bp = log(nSamples)*(nVars+nEdges);    % GSBAP version
end

