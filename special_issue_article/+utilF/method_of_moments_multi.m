%% Method of moments using a consistency equation with support for several trajectories

% INPUT:
% data - given trajectoris of displacements as a 1D cell array
% fTrue - given value of p21/(p12+p21)
% sigmaB1 - std dev for dx and dy when in state 1
% sigmaB2 - std dev for dx and dy when in state 2 % NOTE: switch so that D1 > D2, as we do for the Bayesian analysis

% OUTPUT:
% estimated transition probabilities


% Dependencies: a lot...

function [p12DoublyCorrected, p21DoublyCorrected] = method_of_moments_multi(dataMul, fTrue, sigmaB1, sigmaB2)

    Ntraj = length(dataMul); % number of trajectories
    n21Vec = zeros(1,Ntraj); % hold all FP-corrected estimates of n21
    lengthVec = zeros(1,Ntraj); % hold the lengths of all trajectories

    % Determine optimal threshold 
    sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % displacement threshold
%     q = raylcdf(sThresh,sigmaB1); % probability of a single displacement being correctly labelled in the thresholding
    q = raylcdf(sThresh,sigmaB2); % use this if D1>D2, otherwise switch to sigmaB1 in the expression
    lThresh = round(log(0.03)/log(1-q)); % set threshold based on a p-value of 0.03  

    for idx = 1:Ntraj % do for every trajectory

        data = dataMul{idx}; % current trajectory
        N = size(data,2); % time series length
            
        stateVecEst = zeros(1,length(data));
        stateVecEst(data <= sThresh) = 1;
        stateVecEst(data > sThresh) = 2;
              
        % Apply spin flip procedure to reduce miscategorized segments
        flippedStateVec = stateVecEst;
        flippedStateVec  = utilF.apply_spin_flip(flippedStateVec, lThresh, lThresh);
        
        % Estimate n12, n21, n11 and n22 for the flipped version of the state sequence
        [~,~,n21Est,~] = utilF.calculate_nij(flippedStateVec); % observed number of downcrossings
        
        % n12FP = q*(1-q)*(N-1)*(1-q)^lThresh/(1-2*(1-q));
        n12FP = q*(1-q)*(N-1)*(1-q)^lThresh; % assuming long trajectories
        % n12FP = q*(1-q)*(N-1)*(1-(1-(1-q)^lThresh)/(1-(1-q)^N)); % with proper length normalisation
        n21FP = n12FP;
        n21Vec(idx) = n21Est-n21FP; % store this value
        lengthVec(idx) = N;
    end

    % Solve non-linear equation for p21
    D12 = sum(n21Vec)/(fTrue*sum(lengthVec)); 
    func = @(p) utilF.func_self_consistency_eq(p, fTrue, D12, lThresh, lThresh, sum(lengthVec));
    p12Vec = 0:0.0001:1;
    fSelfCon = func(p12Vec);
%     funcNew = @(p) utilF.func_self_consistency_eq_new(p, fTrue, D12, lThresh, lThresh, sum(lengthVec));
%     fSelfConNew = funcNew(p12Vec);
%     plot(p12Vec,fSelfCon)
%     hold on
%     plot(p12Vec,fSelfConNew)

    [~,idx] = min(fSelfCon); % index of minimum
    p0 = [1E-4,p12Vec(idx)]; % search interval for solution
    p12DoublyCorrected = fzero(func,p0); % root-finding
    p21DoublyCorrected = (fTrue/(1-fTrue))*p12DoublyCorrected;


end
    
    
    




