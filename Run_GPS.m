function [P,Info] = Run_GPS(CovMat,nSamples,mode,P0,truePAG)
% Run GPS algorithm for Gaussian data (using MAG score) or true graph
% (using SHD score), with configurable run mode / parameters. 
% Input:
% - CovMat   = covariance matrix obtained from data / model
% - nSamples = (equivalent) nr. of samples behind CovMat
% - mode     = struct with running mode for GPS
%   .ScoreType  : {'MML' (Gaussian MAG score, default), 'SHD' (structural Hamming distance)}
%   .RICF_tol   : tolerance threshold for convergence in RICF step (MML, default=1e-4)
%   .MAGFlags   : [allow bidirected edges (0/1), allow undirected edges (0/1)] (default=[1,0])
%      note: existing Gaussian MAG score cannot handle undirected edges, so
%            default [1,0] avoids generating PAGs with invariant undirected edges
%   .Stages     : which operators active when & how, encoded as [stages x operators] matrix 
%       GPS-flow: start at first stage, if no improvement found progress to next stage 
%                 until no more stages available, then return top-scoring PAG
%       Stages(i,:) = [AddEdge, DeleteEdge, MakeNoncol, MakeCol, nextStage], with:
%       - value per operator = {0,1,2} with
%          0 = off (no candidate PAGs for this operator), 
%          1 = basic (1 or 2 candidate PAGs per operator, default 'noncollider' for new triples with order), 
%          2 = extended (candidate PAGs for all non/collider combinations for new triples with order per operator)
%       - nextStage = row index of next Stage (in Stages matrix) to go to after success
%       => Examples:
%          Baseline: Stages = [1,1,1,1,1];
%          Extended: Stages = [2,2,2,2,1];
%          Hybrid  : Stages = [1,1,1,1,1;
%                              2,2,2,2,1];   note: final '1' reverts back to baseline (stage at row 1) 
%                                                  after 'extended' finds new top score
% - P0       = (optional) starting PAG (default empty graph)
% - truePAG  = (optional) ground truth PAG (for use with SHD, so not feasible in practice)
% Output:
% - P    = completed partial ancestral graph encoded in matrix form with 
%   P(i,j) = 0  : no mark at i (no edge i-j) i     j
%          = 1  : tail mark at i (to j)      i --* j    
%          = 2  : arrowhead at i (to j)      i <-* j    
%          = 3  : circle mark at i (to j)    i o-* j    (non-invariant)  
% - Info = (optional) struct with counts/intermediate results/tracking etc. for analysis
%          => ignore (or see subroutines for details)
% =========================================================================
% Set local debug level, either via global, or overrule here:
% => 0 = no logging, 1 = basic, 2 = extended/detailed
global DEBUG
if ~isempty(DEBUG), debug = DEBUG; else debug = 0; end;

  % 1: Initialise 
  if (debug > 1), disp('Function Run_GPS()'); end; 
  % process input (no sanity checks included)
  if (nargin < 3), 
    % default mode: MML
    ScoreType   = 'MML';        % MAG likelihood score
    RICF_tol    = 1e-4;         % default tolerance lvl for convergence
    Stages      = [1,1,1,1,1];  % baseline: run all basic operators
                                % [edge add, edge delete, col -> noncol, noncol -> col,next-stage-id]
    MAGFlags    = [1,0];        % [allow bidirected, allow undirected edges
  else
    % use mode info
    ScoreType = mode.ScoreType;
    if (ScoreType == 'MML'), RICF_tol  = mode.RICF_tol; else RICF_tol  = 0; end;
    Stages    = mode.Stages;    % arrat [nrStages x 4 operators] array
    MAGFlags  = mode.MAGFlags;
  end;
  if (nargin < 4), P0 = []; end;
  
  % init output
  Info = [];
  if ~isempty(P0), 
    P = P0;          % e.g. when starting from FCI output PAG
  else
    P = zeros(D,D);  % default start from empty graph
  end;
  % note: counts/trackers/etc. can be ignored (strictly for analysis/verification purposes)
  TotCounts.nMECs = zeros(9,4);  % row 1=operator tried, 2=cand. valid, 3=cand. invalid, 
             % 4=mag valid, 5=mag invalid, 6=pag dupl, 7=pag add, 8=limit cand, 9=max.cand
  
  % set/extract relevant parameters/variables
  % graph size/dimension DxD 
  if ~isempty(CovMat), D = size(CovMat,1); else D = size(truePAG,1); end;
  nrStages = size(Stages,1);
  loopcounter = zeros(1,nrStages);  % track nr. of loops per stage
    
  % 2: start run GPS
  % initialise search stage with starting best/top PAG/MAG/MEC 
  topMEC.P = P;
  topMEC.G = cpag_to_mag(topMEC.P,1);       % 1 = arc-augmented MAG (default)
  topMEC.M = mag_to_mec(topMEC.G);          % MEC = {S,C,D}
  % get initial score (best score so far, minimized in loop below)
  switch(ScoreType) 
    case 'MML'      % MAG maximum likelihood score (on RICF fit) + BIC penalty
      [topMEC.Score,topMEC.Comps,topMEC.CompScores] = GetMAGScore(topMEC.G,CovMat,nSamples,RICF_tol);
    case 'SHD'      % SHD score relative to true PAG
      topMEC.Score = GetPAGSHD(truePAG,topMEC.P);   
    otherwise
      disp('ERROR - Run_GPS: Unknown score %s',ScoreType);
      return;
  end; % switch(score)
  topScore = topMEC.Score; 
  
  % 3: outer loop over all stages (until no more improvements)
  idxStage = 1;
  while (idxStage <= nrStages)
    % get active operators
    Operators = Stages(idxStage,[1:4]);
    done = false;
    % run inner stage loop until no more improvement
    while ~(done)
      loopcounter(idxStage) = loopcounter(idxStage) + 1;
      if (debug >= 1), fprintf('Stage %d, Loop %d: topscore = %d\n',idxStage,loopcounter(idxStage),topScore); end;

      % generate collection of neighbouring MECs/PAGs/MAGs from M
      [MECS,nCounts] = GetNeighbourMECs(topMEC.M,Operators,MAGFlags,topMEC.P);
      % update counts
      TotCounts.nMECs  = TotCounts.nMECs  + nCounts.nMECs;

      % score neighbouring MECs/PAGs for this run (track best)
      nM = size(MECS,2);
      minScore = +Inf; minId = 0;
      for i = 1:nM
        % get score on PAG (SHD) or MAG (likelihood)
        switch(ScoreType)
          case 'MML'      % MAG maximum likelihood score (on RICF fit)
          % (only recompute changed c-components/districts)
            [MECS(i).Score, MECS(i).Comps, MECS(i).CompScores] = ...
                GetMAGScore(MECS(i).G,CovMat,nSamples,RICF_tol,topMEC.G,topMEC.Comps,topMEC.CompScores);
          case 'SHD'      % SHD to true PAG
            MECS(i).Score = GetPAGSHD(truePAG,MECS(i).P);  
        end; % switch(score)
        
        % keep track of best so far (smaller, lower error)
        if MECS(i).Score < minScore,
          minScore = MECS(i).Score;
          minId = i;
        end;
      end;  % for i (loop over all neighbour MECs)
      
      % if improve over previous: set as new best and try again
      if minScore < topScore        % score is 'error', so lower = better
        % set as new starting MEC
        topScore = minScore;
        topMEC   = MECS(minId);

        % show progress (if needed)
        if (debug >= 1), figure(4); clf; draw_cpmag(topMEC.P);  title(sprintf('topPAG (iter=%d)',sum(loopcounter))); end;
        % check stage to return to on success (allows for hybrid versions
        % to escape local optima before returning to simpler strat again
        if (size(Stages,2) > 4) && (Stages(idxStage,5) ~= idxStage),
          % revert back to other (simpler) Stage
          idxStage = Stages(idxStage,5);    % stage(:,5) = next stage on gain
          Operators = Stages(idxStage,[1:4]); % simpler operator version
        end;        
      else
        % no more improvement: exit inner loop
        done = true;
      end;
    end;  % while (stage loop)

    % stage complete, try next one?
    idxStage = idxStage + 1;
  end;  % while (idxstage (in case of different search stages)
   
  
  % 4: collect final output 
  P = topMEC.P;
  if (nargout > 1),
    Info.topMEC.M       = topMEC.M; % save MEC skeleton/coll/noncol lists
    Info.topScore       = topScore;
    Info.loopcounter    = loopcounter;
    Info.ScoreType      = ScoreType;        % MML/SHD
    if (ScoreType == 'MML'), Info.RICF_tol  = RICF_tol; end;
    Info.Stages         = Stages;           % single stage, try everything
    Info.MAGFlags       = MAGFlags;        % [allow bidirected, allow undirected edges
    Info.TotCounts      = TotCounts;
  end

  % done.
  return;
      
end  % function Run_GPS
      