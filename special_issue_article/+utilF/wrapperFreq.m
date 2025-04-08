
function [p12DoublyCorrected, p21DoublyCorrected, p12NoCorrection, p21NoCorrection,...
    p12FP, p21FP, p12FN, p21FN] = wrapperFreq(data, fTrue, sigmaB1, sigmaB2) 

    % Estimates transition probabilities from two-state "emission" data
    % 
    % Input:
    % data = trajectories of displacements as a 1D cell array
    % fTrue = p21/(p12+p21)
    % sigmaB1 - std dev for dx and dy when in state 1
    % sigmaB2 - std dev for dx and dy when in state 2 
    % 
    % Output:
    % p12DoublyCorrected = estimated p12 transition probability
    % p21DoublyCorrected = estimated p21 transition probability
    % 
    % Dependencies: 
    % method_of_moments_thresh
    % 

    sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % displacement threshold
    q = raylcdf(sThresh,sigmaB2); % NOTE: use this if D1>D2, otherwise switch to sigmaB1 in the expression
    lThresh = round(log(0.03)/log(1-q)); % set threshold based on a p-value of 0.03
    
    flag = true; % true until a solution has been found 
    while flag && lThresh >=0
        flag = false;
        try
            [p12DoublyCorrected, p21DoublyCorrected, p12NoCorrection, p21NoCorrection,...
                p12FP, p21FP, p12FN, p21FN] = utilF.method_of_moments_thresh(...
                sThresh,q,lThresh,data,fTrue);
        catch
            lThresh = lThresh-1; % lower segment length threshold if necessary
            flag = true;
        end
    end

    if lThresh == -1
        error('frequentist method did not converge')
    end

end