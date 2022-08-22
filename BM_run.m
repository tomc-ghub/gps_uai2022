function [BR,resBR,resEvalBR,resCounts] = BM_run(BR,BMdir,Resdir)
% or SCRIPT to run a single benchmark ?
% wrapper to run <BRmode> for GPS tests
% input:
% - BR.ID       = '1';             % benchmark run ID 
% - BR.BMset    = 'GPS_A';         % benchmark collection ID (idem)
% - BR.Alg      = 'GPS_MML';       % algorithm/score to use
% - BR.MM       = [1:100];         % 100 models
% - BR.DD       = [6,10,20];       % 3 sizes
% - BR.KK       = [2.5,3.5];       % 2 densities
% - BR.NN       = [];              % nr.of data points 
% - BMdir    e.g. './BM/data';  % parent location of BM data sets for GPS tests
% - Resdir   e.g. './BM';  % output directory for result file

% what do we get from a PC/BCCD/ACI/GES?
% - single PAG + possibly: distribution over each mark {-,>,o,' '}

% evaluation per model:
% - confusion matrix for PAG - truePAG (PAG <-> trueAG ?)
%
% NOT:
% - causal matrix MC - truePAGreach and MC-trueMC
%  => precision = TP/(TP+FP), recall = TP/(TP+FN), accuracy = (TP+TN)/(P+N)
%  => FDR = FP/(TP+FP), Type-I err = FP, Type-II err = FN
% set based quantities:
% - AUC: (order predicitons of causal/noncausal/edgemarks, plot Pr vs. Re
% - overall precision,recall,accuracy,FDR
% weighted measures:
% - for n predictions with k possible outcomes, divide by proportion of k
% =========================================================================
global DEBUG
if ~isempty(DEBUG), debug = DEBUG; else debug = 0; end; % local debug lvl

% =========================================================================
% 1: intialise
if (nargin < 4), BMdir  = './BM/data'; end; 	% default output directory => <dir>/BM<ID>/<BM_...>
if (nargin < 5), Resdir = './BM/result'; end;

% 1b: BM set (separate file)
BR.BM_Info = BM_getInfo(BR.BMset);
% load BM data file (variables: BM_Info, BM); ignore BM_Info (= BR.BM_Info)
load([BMdir,'/',BR.BMset,'/',BR.BM_Info.BMfile]);  % fills variable 'BM'

% 1c: algorithm/score used in benchmark run
BR.Mode = getAlgorithmMode(BR.Alg);   % local function below

% 1d: extract/collect batchrun info from BR
MM      = BR.MM;    % idxs of models to run e.g. [51:100]
DD      = BR.DD;    % nr. of variables per model [5:5:25]
KK      = BR.KK;    % set of densities to consider [2.5,3.5]
NN      = BR.NN;    % nr. of data points [200,1000,5000] 
if isempty(NN), NN = [0]; end;  % SHD uses PAG, not data
% get nr. of tests per dimension
nMM = length(MM);
nDD = length(DD);
nKK = length(KK);
nNN = length(NN);
nr_test = nMM*nDD*nKK*nNN;

% 1c: prepare output
BRoutfile = [Resdir,'/',sprintf('BR%s_%s_M%d_D%d-%d_K%.1f-%.1f_N%d-%d.mat', ...
  BR.ID,BR.Alg,nMM, DD(1),DD(end), KK(1),KK(end), NN(1),NN(end))];

% create output targets to fill (models inferred)
resBR = cell(nMM,nDD,nKK,nNN);  % model-id x dimension x density x data size
% idem result of evaluation per model (as single table)
%            1   2 3 4 5   6    7    8   9-24
% Result = [(idx,M,D,K,N,) accP,accS,SHD,confP(:)'];
resEvalBR = zeros(nr_test,24);    

% notify
fprintf('\nRunning BM %s: for %s\n', BR.ID, BR.Alg);
fprintf(' - %d models\n', length(MM));
fprintf(' - graphs size %d-%d\n', DD(1),DD(end));
fprintf(' - avg.density %d-%d\n', KK(1),KK(end));
fprintf(' - data size %d-%d\n',NN(1),NN(end));
fprintf('result BRoutfile = %s \n',BRoutfile);


% 2: process  models
idx = 0;
% outer loop over models of increasing size d
for idx_d = 1:nDD
  d = DD(idx_d);
  % loop over densities k
  for idx_k = 1:nKK
    k = KK(idx_k);
    fprintf('Running %d models, size D=%d, density K=%.1f \n',nMM, d,k);
    % loop over individual models
    for idx_m = 1:nMM
      m = MM(idx_m);      % note: BR can do subrange of models (batch)
      % get (true) model M (contains MAG, PAG, etc.)
      bm_d = find(BM_Info.DD == d,1);   % find mathcing entry in BM for D
      bm_k = find(BM_Info.KK == k,1);   % idem for K
      BMmdk = BM{m,bm_d,bm_k};    % imported benchmark model data 

      % choose number of data points to use
      for idx_n = 1:nNN
        % get relevant input
        n = NN(idx_n);
        % for GPS-MVG datasize N is just parameter with cov, not actual data set
  %       datafile = tBMmn.datafile;
  %       D = load([input_dir,datafile],'-ascii');
  %       D  = D(1:nD,:);      % sample data set (or read from file)
        % BK = background knowledge? No.

        time = clock;
        fprintf('[%d:%d:%2.1f] === Running model %d (D=%d,K=%d,N=%d).\n',time(4:6),m,d,k,n); 
        % call corresponding causal discovery routine 
        % [P,MC,S,dP,dMC,dS] = BMRun_PC_MVG(D,mode);
        tic;        % track elapsed time
        switch BR.Alg
          case 'GPS_MML' 	% GPS with MAG score
            [P,Info] = Run_GPS(BMmdk.Cov,n,BR.Mode,[]);
          case {'GPS_SHD','GPSX_SHD','GPSGPS2_SHD','GPSGPS2X_SHD'}    % GPS with SHD score
            [P,Info] = Run_GPS([],[],BR.Mode,BMmdk.P);  % true PAG
          case 'GMS_MML'    % GMS with MAG score
            [P,Info] = Run_GMS(BMmdk.Cov,n,BR.Mode,[]);
          case 'GMS_SHD'    % GMS with SHD score
            [P,Info] = Run_GMS([],[],BR.Mode,BMmdk.P);  % true PAG
          case 'GPSGMS_SHD'    % GPS+GMS with SHD score
            [P,Info] = Run_GPSGMS([],[],BR.Mode,BMmdk.P);  % true PAG
%           case 'GPSGPS2_SHD'    % GPS+GMS with SHD score
% %            [P,Info] = Run_GPSGPS2([],[],BR.Mode,BMmdk.P);  % true PAG
%             [P,Info] = Run_GPS([],[],BR.Mode,BMmdk.P);  % true PAG
          case 'MMHC'       % MMHC with MAG score
            [P,Info] = Run_MMHC(BMmdk.Cov,n,BR.Mmode);
        end;
        Info.elapsedTime = toc;

        % collect single run output
        resBR{idx_m,idx_d,idx_k,idx_n}.P    = P;        % PAG
        resBR{idx_m,idx_d,idx_k,idx_n}.Info = Info;  	% counters/details of alg.run

        % evaluate and store next entry
        idx = idx + 1;
        % return accuracy/SHD/conf.mat of PAG found
        result = BM_eval_AccMcPAG(BMmdk,resBR{idx_m,idx_d,idx_k,idx_n});    % true vs. found
        resEvalBR(idx,:) = [idx,m,d,k,n,result];
        %     1   2 3 4 5   6    7    8   9-24
        % = [(idx,m,d,k,n,) accP,accS,SHD,confP(:)'];
        tmpCounts = Info.TotCounts.nMECs';  % layout transposed
        resCounts(idx,:) = [idx,m,d,k,n,tmpCounts(:)'];

      end; % for n (save after each model ?)

      % save intermediate result every 100 models?
      if (mod(idx,100) == 1), save(BRoutfile,'BR','resBR','resEvalBR'); end;

    end; % for m  
  end; % for k
  % make sure results are saved after each size model completed
  save(BRoutfile,'BR','resBR','resEvalBR');
end;  % for d

% stage 3: compute other global performance criteria (e.g. AUC)? not now
% not needed here: save(BRoutfile,'BRmode','resBR','resEvalBR');

end  % function


% ===========
% subfunction for algorithm defaults
% ===========
function mode = getAlgorithmMode(Alg)
  % configure default mode for type of algorithm
  mode.Alg = Alg;   % bit superfluous but hey
  switch(Alg)
    case 'GPS_MML'    % GPS with MAG Max.L score
      mode.Stages       = [1,1,1,1];    % one stage: [edge add, edge del, col -> noncol, noncol -> col]
                                  % (possibly as matrix for succesive stages)
      mode.MAGFlags     = [1,0];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'MML';    
      mode.tol          = 1e-4;         % default tolerance for RICF convergence
    case 'GPS_SHD'    % GPS with Struct.Hamming Distance score
      mode.Stages       = [1,1,1,1];    % [edge add, edge del, col -> noncol, noncol -> col]
      mode.MAGFlags     = [1,1];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'SHD';    
    case 'GPSX_SHD'    % GPS with Struct.Hamming Distance score
      mode.Stages       = [2,2,2,2];    % [edge add, edge del, col -> noncol, noncol -> col]
      mode.MAGFlags     = [1,1];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'SHD';    
    case 'GMS_MML'    % Greedy MAG Search with MAG ML score
      mode.Stages       = [1,1,1,1,1];  % flag [-X-, <--, -->, <->, ---]
      mode.MAGFlags     = [1,0];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'MML';    
      mode.tol          = 1e-4;         % default tolerance for RICF convergence
    case 'GMS_SHD'    % Greedy MAG Search with SHD score
      mode.Stages       = [1,1,1,1,1];  % flag [-X-, <--, -->, <->, ---]
      mode.MAGFlags     = [1,1];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'SHD';    
    case 'GPSGMS_SHD'    % GPS with Struct.Hamming Distance score
      mode.Stages       = [1,1,1,1,1];    % [edge add, edge del, col -> noncol, noncol -> col]
      mode.MAGFlags     = [1,1];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'SHD';    
    case 'GPSGPS2_SHD'    % GPS-GPSext with Struct.Hamming Distance score
      mode.Stages       = [1,1,1,1;     % [edge add, edge del, col -> noncol, noncol -> col]
                           2,2,2,2];    % 1=basic, then 2=extended
      mode.MAGFlags     = [1,1];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'SHD';    
    case 'GPSGPS2X_SHD'    % GPS-GPSext-GPS  with Struct.Hamming Distance score
      mode.Stages       = [1,1,1,1,1;     % [edge add, edge del, col -> noncol, noncol -> col]
                           2,2,2,2,1];    % 1=basic, then 2=extended
      mode.MAGFlags     = [1,1];        % [allow bidirected, allow undirected edges
      mode.ScoreType    = 'SHD';    
    case 'GMS_MMHC'    % Greedy MAG Search Sofia version
      mode.Stages       = [1,1,1,1,1];  % flag [-X-, <--, -->, <->, ---]
      mode.MAGFlags    = [1,0];        % [allow bidirected, allow undirected edges
      mode.tol         = 1e-4;         % default tolerance for RICF convergence
    case {'PC','CPC','FCI','FCI+'}  % PC/Conservative PC/FCI/FCI+ (MVG test)
      mode.maxTests = 1000000;    % max. nr of tests per stage per model
      mode.alpha    = 0.05;       % threshold p-value (adapt per data size?)
    case {'BPC','BCPC','BFCI','BFCI+'} % idem bootstrapped versions
      mode.maxTests   = 1000000;   % max. nr of tests per stage per model
      mode.alpha      = 0.05;     % threshold p-value (adapt per data size?)
      mode.nBootstrap = 100;      % nr. of bootstrap samples to use
      mode.Bmode      = 2;        % 1=bootstrap CI-test, 2=bootstrap models
    case {'PC_ORA','FCI_ORA','FCI+_ORA'} % idem Oracle versions
      mode.maxTests = 1000000;    % max. nr of tests per stage per model
    case 'BCCD' 
      mode.maxTests = 1000000;    % max. nr of tests per stage per model?
      mode.NO_SEL_BIAS = ~true;   % BCCD will/not assume 'no selection bias'
      mode.ensure_maxK = 1;       % ensure all relevant L statements (p(L)>minp_loci/2) tested on K variables
    otherwise
      fprintf('ERROR - BR_getInfo: Unknown algorithm %s',Alg); 
  end;  % switch algorithm
end  % sub-function getalgorithmmode
