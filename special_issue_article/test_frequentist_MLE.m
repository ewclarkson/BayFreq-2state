%% Frequentist MLE as an alternative to the method of moments

% Simulation parameters
N = 2E3;              % Number of sampling time intervals for "actual" data
sigmaB1 = 0.1;      % Std dev for dx and dy when in state 1
sigmaB2 = 0.5;      % Std dev for dx and dy when in state 2 % NOTE: switch so that D1 > D2, as we do for the Bayesian analysis
p12 = 0.10;         % transition probability from state 1 to state 2
p21 = 0.05;         % transition probability from state 2 to state 1  

% Analysis parameters
NDensityEst = 1E6;    % Number of sampling time intervals for estimation of 
                      % density of level-crossings for one-state diffusion
                      % (NDensityEst should be large)
XStdDevs = 2;         % we determine the length threshold for the spin flip
                      % procedure  using a mean + X std dev rule

nIter = 1E1;        % number of trajectories to analyze for each p12
nP = 1E3;       % number of p12-values


%% Given trajectory

fTrue = p21/(p12+p21);   % fraction of time when we are in state 1
stateVec0 = utilF.twoState_Markov(p12, p21, N); % Erik's code
% [~,~,n21,~] = utilF.calculate_nij(stateVec);
% p12Est = n12/(n11+n12);
data0 = utilF.brownian_displacements_2d(stateVec0,sigmaB1,sigmaB2);


%% Analysis of given trajectory, using only thresholding

n12Vec = zeros(1,nP);
sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % optimal threshold
p12Vec = linspace(0.01,0.99,nP);

for idx = 1:nP % do for p12-values
    n12Temp = zeros(1,nIter); % values from different realisations
    for idy = 1:nIter % re-do analysis
        % generate a trajectory
        stateVec = utilF.twoState_Markov(p12Vec(idx), fTrue/(1-fTrue)*p12Vec(idx), N);
        data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
        
        % threshold and categorize into states
        stateVecEst = zeros(1,length(data));
        stateVecEst(data <= sThresh) = 1;
        stateVecEst(data > sThresh) = 2;

        % Estimate n12, n21, n11 and n22 for the state sequence
        [~,n12Est,~,~] = utilF.calculate_nij(stateVecEst);
        n12Temp(idy) = n12Est;  
    end
    n12Vec(idx) = mean(n12Temp);
end

histogram(n12Vec) % appears almost uniform
histogram(n12Temp) % gaussian

%% Analysis of given trajectory, using state flipping

n12Vec = zeros(1,nP);
sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % optimal threshold
% Find mean and variance of miscategorized segments
[ ~ , ~ , mean2X , std2X ] = ...
    utilF.determine_segment_length_mean_and_std( sigmaB1 , sThresh , N   );
lThresh2X = ceil(mean2X + XStdDevs*std2X);      
                    % length threshold for regions 2X = '2','22','222', etc

[mean1X , std1X, ~ , ~ ] = ...
    utilF.determine_segment_length_mean_and_std( sigmaB2 , sThresh , N   );
lThresh1X = ceil(mean1X + XStdDevs*std1X);     
                    % length threshold for regions 1X = '1','11','111', etc
p12Vec = linspace(0.01,0.99,nP);

for idx = 1:nP % do for p12-values
    n12Temp = zeros(1,nIter); % values from different realisations
    for idy = 1:nIter % re-do analysis
        % generate a trajectory
        stateVec = utilF.twoState_Markov(p12Vec(idx), fTrue/(1-fTrue)*p12Vec(idx), N);
        data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
        
        % threshold and categorize into states
        stateVecEst = zeros(1,length(data));
        stateVecEst(data <= sThresh) = 1;
        stateVecEst(data > sThresh) = 2;
        
        % Apply spin flip procedure to reduce miscategorized segments
        flippedStateVec  = utilF.apply_spin_flip(stateVecEst, lThresh1X, lThresh2X);
        
        % Estimate n12, n21, n11 and n22 for the flipped state sequence
        [~,n12Est,~,~] = utilF.calculate_nij(stateVecEst);
        n12Temp(idy) = n12Est;  
    end
    n12Vec(idx) = mean(n12Temp);
end

histogram(n12Vec) % is it gaussian?
















