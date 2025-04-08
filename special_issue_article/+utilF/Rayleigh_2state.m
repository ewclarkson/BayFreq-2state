

function [est, errorEst] = Rayleigh_2state(D1_guess, D2_guess, pi_guess, tau, trackCell)

    % Estimate diffusion constants D1&D2 and the fraction of time, pi, spent
    % in state 1, with maximum likelihood estimation.
    % 
    % Input:
    % D1_guess = guessed value of D1 
    % D2_guess = guessed value of D2
    % pi_guess = guessed value of pi
    % tau = sampling time
    % trackCell = cell array of trajectory displacements
    % 
    % Output:
    % est = array of estimated parameters [pi,D1,D2]
    % 
    % Dependencies: w_rayleigh_pdf
    
    nTracks = length(trackCell);
    
    % Extract step sizes
    StepSizes = []; % initialise holder for all estimated diffusion constants

    for k = 1:nTracks % do for every track
        for j = 1:length(trackCell{k}) % do for every displacement

            dist = trackCell{k}(j); % displacement
            StepSizes(end+1) = dist; % store displacement
        end
    end

    % Compute a maximum-likelihood estimates
    in_guess = [pi_guess, sqrt(2*tau*D1_guess), sqrt(2*tau*D2_guess)]; % f, th1, th2
    mypdf = @(data,f,th1,th2) utilF.w_rayleigh_pdf(data, f, th1, th2); % define function

    est = mle(StepSizes, 'pdf', mypdf, 'start', in_guess);  
%     acov = mlecov(est,StepSizes,'pdf',mypdf);
%     errorEst = sqrt(diag(acov));  % estimated errors for each of the parameters
    errorEst = NaN;

    % plotting
%     xVec = linspace(0, max(StepSizes), 1000); % create x-vector for plotting
%     histogram(StepSizes, 'Normalization', 'pdf');
%     xlabel('Step size (\mum^2)')
%     ylabel('Normalised frequency')
%     hold on
%     plot(xVec, utilF.w_rayleigh_pdf(xVec, est(1), est(2), est(3)), 'linewidth', 2.0')

    % output
    est(2) = est(2)^2/(2*tau); est(3) = est(3)^2/(2*tau); % convert to diffusion constants

end








