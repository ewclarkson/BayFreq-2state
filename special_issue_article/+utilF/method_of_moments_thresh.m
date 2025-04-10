%% Method of moments using a consistency equation with support for several trajectories

% INPUT:
% sThresh - displacement threshold for binarisation
% q - probability of correctly labelling a displacement
% data - given trajectoris of displacements as a 1D cell array
% fTrue - given value of p21/(p12+p21)
% sigmaB1 - std dev for dx and dy when in state 1
% sigmaB2 - std dev for dx and dy when in state 2 

% OUTPUT:
% estimated transition probabilities using various corrections

% Dependencies:
% determine_optimal_thresh_displ
% apply_spin_flip
% calculate_nij
% func_self_consistency_eq

function [p12DoublyCorrected, p21DoublyCorrected, p12NoCorrection, p21NoCorrection,...
    p12FP, p21FP, p12FN, p21FN] = ...
    method_of_moments_thresh(sThresh,q,lThresh, dataMul, fTrue)

    Ntraj = length(dataMul); % number of trajectories
    n21Vec = zeros(1,Ntraj); % hold all FP-corrected estimates of n21
    lengthVec = zeros(1,Ntraj); % hold the lengths of all trajectories

    % Determine optimal threshold 
    % sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % displacement threshold
    % q = raylcdf(sThresh,sigmaB2); % NOTE: use this if D1>D2, otherwise switch to sigmaB1 in the expression

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
       
        n21Vec(idx) = n21Est;
        lengthVec(idx) = N;
    end

    % Solve non-linear equation for p21
    p12NoCorrection = sum(n21Vec)/(fTrue*sum(lengthVec));
    p21NoCorrection = (fTrue/(1-fTrue))*p12NoCorrection;

    p12FP = (sum(n21Vec)-q*(sum(lengthVec)-2)*(1-q)^(lThresh+1))/(fTrue*sum(lengthVec));
    p21FP = (fTrue/(1-fTrue))*p12FP;

    func = @(p) utilF.func_self_consistency_eq(p, fTrue, p12FP, lThresh, lThresh);
    p12Vec = 0:0.0001:1;
    fSelfCon = func(p12Vec);
    [~,idx] = min(fSelfCon); % index of minimum
    p0 = [1E-4,p12Vec(idx)]; % search interval for solution
    p12DoublyCorrected = fzero(func,p0); % root-finding
    p21DoublyCorrected = (fTrue/(1-fTrue))*p12DoublyCorrected;

    funcFN = @(p) utilF.func_self_consistency_eq(p, fTrue, p12NoCorrection, lThresh, lThresh);
    fSelfConFN = funcFN(p12Vec);
    [~,idx] = min(fSelfConFN); % index of minimum
    p0 = [1E-4,p12Vec(idx)]; % search interval for solution
    p12FN = fzero(funcFN,p0); % root-finding
    p21FN = (fTrue/(1-fTrue))*p12FN;

end
    
    
    




