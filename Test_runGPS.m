% Example script for running GPS in different modes
% Based on UAI-22 article 'Greedy Equivalence Search in the Presence of Latent
% Confounders' (Claassen & Bucur, 2022).
% Script overview:
% 1. set high level test case: model, score, and GPS run type
% 2. add detailed information (ignore, or tweak if needed)
% 3. run GPS
% 4. display top-scoring PAG

% =========================================================================
% 1: set high level test case: model, score, and GPS run type
% =========================================================================
global DEBUG
DEBUG = 1;          % 0=off (no logging), 1=basic, 2=detailed
rng(1234);          % set random seed for reproducibility (if needed)
% 1a: select GPS run type (baseline/extended/hybrid)
GPS_RunType = 0;    % 0=baseline, 1=extended, 2=hybrid
GPS_P0 = 'EMPTY';   % default start from empty graph P0
                    % optional: set starting graph P0 via external routine (e.g. FCI) in step 2c, below 

% 1a: select test graph (id 1=random, 2-5 predefined testgraphs)
graph_id = 4;       % id 4 = example Fig.1 GPS article
if (graph_id == 1),
  % configure random (ancestral) graph properties
  N     = 20;       % nr. of variables
  maxK  = 6;        % max. node degree
  avgK  = 3.0;      % graph density (avg. node degree)
  p_bi  = 0.20;     % chance edge is bidirected
  p_sel = 0.00;     % chance on (direct) selection per node
end;

% 1b: select GPS scoring metric (MML/SHD: pick one)
%Mode.ScoreType = 'MML';     % Gaussian MAG score (max.likelihood+BIC penalty)
Mode.ScoreType = 'SHD';     % Structural Hamming Distance vs. true PAG
% set score specific parameters
switch(Mode.ScoreType)
  case 'MML'
    Mode.RICF_tol = 1e-6;   % convergence threshold for RICF step 
    Mode.MAGFlags = [1,0];  % allow bidirected, not undirected
    nSamples = 10000;        % equivalent nr. of data points (1000/10000 gives reasonable performance)
    Sigma = [];             % optional: include predefined covariance matrix Sigma 
                            % (if empty then random MVG model generated for chosen ancestral graph in step 2b, below)
  case 'SHD'
    Mode.MAGFlags = [1,1];   % allow bidirected + undirected
end;


% =======================================================================
% 2: add detailed information 
% =======================================================================
% 2a: get ground truth ancestral graph (MAG/PAG, predefined or random)
switch (graph_id)
  case 1  % random generated ancestral graph (used in GPS experiments)
    % properties (adjust above as needed)
    % N     = 20;       % nr. of variables
    % maxK  = 4;        % max. node degree
    % avgK  = 2.5;      % graph density (avg. node degree)
    % p_bi  = 0.15;     % chance edge is bidirected
    % p_sel = 0.0;      % chance on (direct) selection per node
    [AG,Reach,Vtime,Sel] = mk_ag(N,maxK,avgK,p_bi,p_sel);
  case 2   % x1-->x2<->x3 (v-structure)
    AG = [0,1,0;
          2,0,2;
          0,2,0];
  case 3   % x1-->x2<->x3 + x2-->x4 (y-structure, invariant causal x2-->x4)
    AG = [0,1,0,0;
          2,0,2,1;
          0,2,0,0;
          0,0,0,2];
  case 4   % Example Fig1 GPS article: A-->B<->C<->D<->E + B/C-->E
    % SHD: both baseline + extended easily find true graph
    % MML: baseline finds correct graph when starting from rng(1234); 10000 samples 
    % but often fails to converge to true graph if edges are too weak
    AG = [0,1,0,0,0;  %A
          2,0,2,0,1;  %B
          0,2,0,2,1;  %C
          0,0,2,0,2;  %D
          0,2,2,2,0]; %E
  case 5   % larger graph over 15 nodes
    % baseline reaches local minimum (SHD=8) after 53 steps, extended
    % reaches true PAG (SHD=0) after 43 steps
    AG = [0,0,0,0,0,0,0,0,2,0,0,2,0,0,0;
          0,0,0,0,1,1,0,0,1,0,1,0,0,0,1;
          0,0,0,2,0,2,0,0,0,0,2,0,0,0,2;
          0,0,1,0,0,0,1,2,0,0,2,2,2,1,0;
          0,1,0,0,0,1,0,0,0,1,1,0,0,0,1;
          0,2,2,0,2,0,0,0,2,0,1,2,0,0,2;
          0,0,0,2,0,0,0,2,2,0,0,2,2,0,0;
          0,0,0,2,0,0,1,0,0,0,1,0,0,0,2;
          1,1,0,0,0,1,1,0,0,0,1,0,0,0,1;
          0,0,0,0,2,0,0,0,0,0,0,2,2,0,0;
          0,2,2,2,2,2,0,2,2,0,0,2,0,0,2;
          1,0,0,2,0,1,1,0,0,2,1,0,0,0,1;
          0,0,0,2,0,0,1,0,0,2,0,0,0,2,0;
          0,0,0,2,0,0,0,0,0,0,0,0,2,0,0;
          0,2,2,0,2,1,0,1,2,0,1,2,0,0,0];  
end;  % switch(graph)
% get corresponding (true) MAG+MEC+cPAG
G = ag_to_mag(AG);
M = mag_to_mec(G);
P = mag_to_cpag(G); % or equivalently: P = mec_to_cpag(M);
N = size(P,1);      % nr. of variables
% display 
figure(1); clf; draw_cpmag(G); title('Ground truth MAG G');
figure(2); clf; draw_cpmag(P); title('Ground truth completed PAG P');

% =======================================================================
% 2b: details per score
switch(Mode.ScoreType)
  case 'MML'
    % get covariance matrix Sigma (if not specified above)
    if isempty(Sigma)
      % create random MVG model for ancestral graph AG 
      MVG.b_range = [-1,+1];        % range of edge strengths
      MVG.s_range = [0.01,0.2];     % range of std.dev. of noise terms
      MVG.minabs_b = 0.1;           % avoid/include very weak edges 
      % generate Sigma for MVG with random edge strengths B + intrinsic noise with variance V
      [Sigma,Mu,B,V] = ag_to_random_MVG_PDF(AG,MVG);   % Mu = zeros(N,1);
    end;
    % Note: it is also possible to sample a cov.matrix e.g. via data or a Wishart,
    % but here we focus on comparing the perfomance of various discovery algorithms,
    % and so opt to remove that contribution to the variability by starting
    % from the (ground truth) covariance matrix Sigma, and assume it has been
    % obtained from nSamples observations.
  case 'SHD'
    % no extra action needed
end; 

% =======================================================================
% 2c: set Stages for GPS runtype (basic/extended/hybrid) + starting
switch(GPS_RunType)
  case 0 % 'baseline'
    Mode.Stages = [1,1,1,1,1];
  case 1 % 'extended'
    Mode.Stages = [2,2,2,2,1];
  case 2 % 'hybrid', baseline, try extended if stuck (then back to baseline)
    Mode.Stages = [1,1,1,1,1;
                   2,2,2,2,1]; % note: final '1' reverts back to first stage (baseline, at row 1) 
end;
% set starting graph (PAG)
switch(GPS_P0)
  case 'EMPTY'
    P0 = zeros(size(P));  % default empty starting graph
%  optional: obtain starting PAG from favourite FCI implementation 
%  on Gaussian data after matching cov.matrix Sigma generated in step 2b
%   case 'FCI',
%     P0 = Run_FCI(Sigma,nSamples,alpha);  % using part.corr. test with decision threshold alpha
end;

% =======================================================================
% 3: run GPS
time = clock;
fprintf('[%d:%d:%2.1f] => Running model %d (N=%d,Mode=%d,Score=%s).\n',time(4:6),graph_id,N,GPS_RunType,Mode.ScoreType); 
% call GPS with chosen score 
% note: if 'debug >= 1', then 'best PAG' progression in Run_GPS shown in 'figure(4)'
tic;        % track elapsed time
switch Mode.ScoreType    % (switch not needed: only to make difference explicit)
  case 'MML'        % GPS with Gaussian MAG score
    [P1,Info] = Run_GPS(Sigma,nSamples,Mode,P0,[]);
  case 'SHD'        % GPS with SHD score
    [P1,Info] = Run_GPS([],[],Mode,P0,P);  % P = true PAG, P1 = PAG found
end;
Info.elapsedTime = toc;

% show result: final score + output graph
time = clock;
fprintf('[%d:%d:%2.1f] GPS finished with score %d\n',time(4:6),Info.topScore);
figure(3); clf; draw_cpmag(P1);  title(sprintf('GPS - Final PAG'));
 
 
% =======================================================================
% END OF SCRIPT TEST_RUNGPS
% =======================================================================
