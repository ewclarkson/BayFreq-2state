function [rho] = density_level_crossings(sigmaB , sThresh , N , ...
                                  lThresh1X, lThresh2X )
    %
    %vDensity of crossings of a level set by threshQ for Brownian
    % displacements with standard deviation = sigmaB
    %
    % Input:
    % sigmaB = std of the Brownian steps
    % sThresh = threshold as applied to the Brownian displacement, s
    % N = length of the simulation
    % lThresh1X = all segments of type '1', '11', '111' which are 
    %             of length <= lThres1X are "flipped" to become segments
    %             '2', '22', '222' etc. Set lThresh1X = 0 to avoid spin
    %             flip.
    % lThresh2X = all segments of type '2', '22', '222' which are 
    %             of length <= lThres2X are "flipped" to become segments
    %             '1', '11', '111' etc. Set lThresh1X = 0 to avoid spin
    %             flip.
    % nRecurSpinFlip = number of recursions in the spin flip procedure 
    %                  (set = 0 if you do not want to apply the spin flip
    %                  procedure)
    %
    % Output:
    % rho = level crossing density
    %
    import('utilF.brownian_displacements_2d')
    import('utilF.apply_spin_flip')

    nSim = 1E2;
    n12FPiVec = zeros(1,nSim); % before spin-flip
    n12FPfVec = zeros(1,nSim); % after spin-flip

    stateVec = ones(1,N);
    for i = 1:nSim % number of repetitions of experiment
        % Simulate a one-state diffusion process
        data = brownian_displacements_2d(stateVec,sigmaB,sigmaB);
           
        % Apply a threshold
        stateVecEst = zeros(1,length(data));
        stateVecEst(data <= sThresh) = 1;
        stateVecEst(data > sThresh) = 2;
     
        diffX = diff(stateVecEst);
        n12Est = sum(diffX == 1);
        n12FPiVec(i) = n12Est; % store value

        % Apply spin flip to remove short segments
        flippedStateVec  = apply_spin_flip( stateVecEst , lThresh1X, lThresh2X );
   
        % Count the number of 1-2 crossings
        diffX = diff(flippedStateVec);
        n12Est = sum(diffX == 1); %length(find(diffX == 1));
        n12FPfVec(i) = n12Est; % store value

    end
%     q = raylcdf(sThresh,sigmaB);
%     n12FPiTheory = q*(1-q)*(N-1)
%     n12FPfTheory = q*(1-q)*(N-1)*(1-q)^lThresh1X
%     n12FPi = mean(n12FPiVec);
    n12FPf = mean(n12FPfVec);
    
%     n12FPf/n12FPi
%     (1-q)^lThresh1X

    rho = n12FPf/(N-1); % density of up-crossings

end

