%% Method of moments using a consistency equation

% INPUT:
% Input parameters
% data - given trajectory of displacements
% fTrue - given value of p21/(p12+p21)
% sigmaB1 - Std dev for dx and dy when in state 1
% sigmaB2 - Std dev for dx and dy when in state 2 % NOTE: switch so that D1 > D2, as we do for the Bayesian analysis
% lThresh - 2 % NOTE: this should be replaced by a p-value
%

% --------- old parameters -----------------------
% Nsim -              number of simulations to determine segment length 
                      % statistics (zero crossing intervals)
% NDensityEst - Number of sampling time intervals for estimation of 
                      % density of level-crossings for one-state diffusion
                      % (NDensityEst should be large)
% XStdDevs - we determine the length threshold for the spin flip
                      % procedure  using a mean + X std dev rule

% OUTPUT:
% p12DoublyCorrected - estimate of transition probability p12
% p21DoublyCorrected - estimate of transition probability p21

% Dependencies: a lot...

function [p12DoublyCorrected, p21DoublyCorrected] = method_of_moments(data, fTrue, sigmaB1, sigmaB2, lThresh)

    N = length(data); % time series length

    % Determine optimal threshold 
    sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2);
                       % threshold for the displacement (everything below 
                       % the threshold is initially deemed to be in state 1, and 
                       % everything above is deemed to be in state 2)
    stateVecEst = zeros(1,length(data));
    stateVecEst(data <= sThresh) = 1;
    stateVecEst(data > sThresh) = 2;
    

    % Determine optimal length threshold
%     lThresh = 2; % TODO: set this based on p-value

%     % Find mean and variance of miscategorized segments
%     [ ~ , ~ , mean2X , std2X ] = ...
%         utilF.determine_segment_length_mean_and_std( sigmaB1 , sThresh , Nsim   );
%     lThresh2X = ceil(mean2X + XStdDevs*std2X);      
%                         % length threshold for regions 2X = '2','22','222', etc
%     
%     [mean1X , std1X, ~ , ~ ] = ...
%         utilF.determine_segment_length_mean_and_std( sigmaB2 , sThresh , Nsim   );
%     lThresh1X = ceil(mean1X + XStdDevs*std1X);     
                        % length threshold for regions 1X = '1','11','111', etc
    
    % Apply spin flip procedure to reduce miscategorized segments
    flippedStateVec = stateVecEst;
    flippedStateVec  = utilF.apply_spin_flip(flippedStateVec, lThresh, lThresh);
    
    % Estimate n12, n21, n11 and n22 for the flipped version of the state sequence
    [~,~,n21Est,~] = utilF.calculate_nij(flippedStateVec);
    %n21Est + n12Est + n11Est + n22Est  % should = N-1, i.e., all nearest neighbour 
                           % pairs: (1,2), (2,3), (3,4), ..., (N-1,N). 

%     % Estimate level crossings for one-state systems
%     [density21State2] = utilF.density_level_crossings(sigmaB2 , sThresh , NDensityEst , ...
%                                               lThresh1X, lThresh2X);
%     [density21State1] = utilF.density_level_crossings(sigmaB1,sThresh,NDensityEst, ...
%                                               lThresh1X, lThresh2X);
    
    % Estimate n12 and n21
%     n21FalsePositivesEst = density21State1*fTrue*N + density21State2*(1-fTrue)*N;
    q = raylcdf(sThresh,sigmaB1);
    % q2 = raylcdf(sThresh,sigmaB2);
    n12FP = q*(1-q)*(N-1)*(1-q)^lThresh;
    n21FP = n12FP;
    
    % Solve non-linear equation for p21. in case there are several solutions 
    % (empirically there seems to be two solutions in the range [0,1])
    % we choose the smallest solution (NOT CLEAR WHY...). Here, we use a
    % very simplistic approach for finding the root.
    D12 = (n21Est - n21FP)/(fTrue*N); 
    func = @(p) utilF.func_self_consistency_eq(p, fTrue, D12, lThresh, lThresh);
    p12Vec = 0:0.0001:1;
    fSelfCon = func(p12Vec);
    % idx = find(diff(sign(fSelfCon)) ~= 0);
    % p12DoublyCorrected = (p12Vec(min(idx)) + p12Vec(min(idx)+1))/2;
    % figure
    % plot(p12Vec,fSelfCon)
    [~,idx] = min(fSelfCon); % index of minimum
    p0 = [1E-4,p12Vec(idx)]; % search interval for solution
    p12DoublyCorrected = fzero(func,p0); % root-finding
    p21DoublyCorrected = (fTrue/(1-fTrue))*p12DoublyCorrected;
end
    
    
    




