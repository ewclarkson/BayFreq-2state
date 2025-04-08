%% Analyse simulated 2D data with ns-HMM

% Neglect localisation error and motion blur.


% ------------------ initialise variables -----------------------

nTraj = 10; % number of trajectories
nSteps = 3E2; % number of steps
tau = 15E-3; % sampling time
D1 = 1.5; % diffusion constant in state 1
D2 = 0.2; % diffusion constant in state 2
p12 = 0.15; % transition probability 1-->2
p21 = 0.05; % transition probability 2-->1
sigmaB1 = sqrt(2*D1*tau); % standard deviation of Gaussian displacements
sigmaB2 = sqrt(2*D2*tau); % standard deviation of Gaussian displacements

% generate trajectories
for idxTrack = 1:nTraj  
    stateVec = utilF.twoState_Markov(p12, p21, nSteps); % random state sequence
    dispCell{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
end


%% run nested sampling based on ideal tracks

% priors
D1_max = 3.0; % known upper limit of D1
D1_min = 0.30; % known lower limit of D1
D2_max = 0.29; % known upper limit of D2
D2_min = 0.01; % known lower limit of D2
p12_max = 0.5; % known upper limit of p12
p12_min = 0.10; % known lower limit of p12
p21_max = 0.2; % known upper limit of p21
p21_min = 0.01; % known lower limit of p21

% algorithm parameters
nLive = 100; % number of live points
StopRatio = 1E-4; % stop criterion for evidence
priorLimits = [D1_min,D1_max;D2_min,D2_max;p12_min,p12_max;p21_min,p21_max];

% run nested sampling
[finalSeq, thetaMLE, logZ] = utilB.nestedsampling2D(nLive, StopRatio, priorLimits, tau, dispCell);
seqLen = length(finalSeq);

[logZ_error, thetaBayes, errorBayes] = utilB.bayesianestimate(finalSeq,logZ,nLive,false);

thetaBayes % parameter estimate









