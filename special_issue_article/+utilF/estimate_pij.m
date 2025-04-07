function [ p12Corrected , p21Corrected ] = ...
    estimate_pij(sigmaB1,sigmaB2,threshQ, flippedStateSeq , NDensityEst) 
    %
    % Calculates a corrected expression for the transition probabilities
    % p12, and p21.
    %
    % Input parameters:
    % NDensityEst = Number of sampling time intervals for estimation of 
                        % density of level-crossings for one-state diffusion
                        % (NDensityEst should be large)
    % sigmaB1 = std dev for dx and dy when in state 1
    % sigmaB2 = std dev for dx and dy when in state 2
    % threshQ = threshold for the displacement (everything below 
    %             the threshold is initiall deemed to be in state 1, and 
    %              everything above is deemed to be in state 2)
    % 
    % Output parameters:
    % p12Corrected = estimate for p12, i.e., the transition probability 
    %                between state 1 and state 2
    % p21Corrected = estimate for p21, , i.e., the transition probability 
    %                between state 1 and state 2
    %
    % Dependencies: None.



    n21FalsePositivesEst = density21State1*fTrue*N + density21State2*(1-fTrue)*N;
    n21Corrected = n21Est - n21FalsePositivesEst;

    n11Corrected = n11Est + fTrue*N*(1 - density11State1)...
             - (1-fTrue)*N*(1 - density22State2 - 2*density21State2); 

    n22Corrected = n22Est - fTrue*N*(1 - density11State1 - 2*density21State1 )...
        + (1-fTrue)*N*(1 - density22State2); 



    disp(['n_{11}, true = ',num2str(n11),', estimated = ',num2str(n11Est),...
        ', corrected = ',num2str(round(n11Corrected))])
    disp(['n_{21}, true = ',num2str(n21),', estimated = ',num2str(n21Est),...
        ', corrected = ',num2str(round(n21Corrected))])
    disp(['n_{22}, true = ',num2str(n22),', estimated = ',num2str(n22Est),...
        ', corrected = ',num2str(round(n22Corrected))])

    % Double-check: n11 + n22 + n21 + n12 = N - 1
    n11Corrected + n22Corrected + 2*n21Corrected

    fCorrected = (n11Corrected + n21Corrected)/(N-1)



    % Bring the estimated values for n_{ij} within the allowed range
    % [0 , N-1]
    n11Corrected = min(max(0,n11Corrected),N-1);
    n21Corrected = min(max(0,n21Corrected),N-1);
    n22Corrected = min(max(0,n22Corrected),N-1);


    % 
     p12Corrected = n21Corrected/(n11Corrected+n21Corrected);
     p21Corrected = n21Corrected/(n21Corrected+n22Corrected);
     % p12Corrected = n21Corrected/(fTrue*N)
    % p21Corrected = n21Corrected/((1-fTrue)*N)

     disp(['p_{12}, true = ',num2str(p12),', estimated and corrected = ',num2str(p12Corrected)])
    disp(['p_{21}, true = ',num2str(p21),', estimated and corrected = ',num2str(p21Corrected)])




end

