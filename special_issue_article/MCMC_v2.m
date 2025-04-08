function [averaged, errorUpper, errorLower, s1, s2, s3, s4] = MCMC_v2(num_MC_steps, burnInTime, scaleArr, guessArr, tau, data)
    % 
    % Optimises a function by using the Metropolis-Hastings algorithm.
    % 
    % Input:
    % num_MC_steps = total number of Mone Carlo steps to be taken
    % burnInTime = number of thermalisation steps until deemed equilibrium
    % scaleArr = std of proposal distribution, for steps in each of the
    %            parameter directions: [D1, D2, p12, p21]
    % guessArr = initial guess of parameter values, on the form 
    %            [D1, D2, p12, p21]
    % tau = sampling time (time between frames)
    % data = cell array containing all trajectory displacements
    % 
    % Output:
    % averaged = parameter estimations as averages after thermalisation
    % errorUpper - upper confidence levels on estimated parameters
    % errorLower - lower confidence levels on estimated parameters
    % s1, s2, s3, s4 = acceptance rates for proposals along parameteres
    %                  D1, D2, p12, p21, respectively
    % 
    % Dependencies:
    % Ens_likelihood
    % 

    import('utilB.EnsLikelihood2D')

    D1_min = guessArr(1,1); D1_max = guessArr(1,2);
    D2_min = guessArr(2,1); D2_max = guessArr(2,2);
    p12_min = guessArr(3,1); p12_max = guessArr(3,2);
    p21_min = guessArr(4,1); p21_max = guessArr(4,2);

    init_D1 = D1_min + (D1_max-D1_min)*rand; % initial D1-coordinates for every live point
    init_D2 = D2_min + (D2_max-D2_min)*rand; % initial D2-coordinates for every live point
    init_p12 = p12_min + (p12_max-p12_min)*rand; % initial p12-coordinates for every live point
    init_p21 = p21_min + (p21_max-p21_min)*rand; % initial p21-coordinates for every live point

    paramArr = [init_D1, init_D2, init_p12, init_p21]; % initial random position
    
    l_logLike = EnsLikelihood2D([paramArr(1), paramArr(2), paramArr(3), paramArr(4)], tau, data);

    n = num_MC_steps;
    l_theta = paramArr;
    
    all_D1s = zeros(1,num_MC_steps); % store all D1s
    all_D2s = zeros(1,num_MC_steps); % store all D2s
    all_p12s = zeros(1,num_MC_steps); % store all p12s
    all_p21s = zeros(1,num_MC_steps); % store all p21s
    
    all_D1s(1) = paramArr(1); all_D2s(1) = paramArr(2); all_p12s(1) = paramArr(3); all_p21s(1) = paramArr(4); % store initial parameter values
    
    acMoves = [0,0,0,0]; % store the number of accepted moves for each parameter (D1, D2, p12, p21)
    prMoves = [0,0,0,0]; % store the total number of proposed moves for each parameter (D1, D2, p12, p21)
    
    
    for i = 1:round(n/4)
       for k = 1:4 
           
           l = 4*(i-1)+k; % keeps count of the current iteration
           
           deltaTheta = normrnd(0,sqrt(scaleArr(k))); % propose a displacement for the kth parameter
           propTheta = l_theta; propTheta(k) = propTheta(k) + deltaTheta; % propose a new state
           
           prMoves(k) = prMoves(k) + 1; % add one count to the number of proposed moves along this axis 
           
           prop_logLike = EnsLikelihood2D([propTheta(1), propTheta(2), propTheta(3), propTheta(4)], tau, data); % log-likelihood for proposed parameters
           
           if prop_logLike >= l_logLike && propTheta(3)< 1 && propTheta(4) < 1 % then accept the proposed state
               
               l_theta = propTheta; % assign the current parameter set to be the proposed one
               l_logLike = prop_logLike; % assign the current log-likelihood to be the proposed one
               
               acMoves(k) = acMoves(k) + 1; % add one count to the number of accepted moves along parameter axis k
               
               all_D1s(l+1) = propTheta(1); all_D2s(l+1) = propTheta(2); all_p12s(l+1) = propTheta(3); all_p21s(l+1) = propTheta(4); % store new parameter values
          
           else % make a new trial for whether or not to accept the proposed state
               
               u = rand();
               if log(u) <= prop_logLike - l_logLike && propTheta(3)< 1 && propTheta(4) < 1 % then accept the proposed state
                   
                   l_theta = propTheta; % assign the current parameter set to be the proposed one
                   l_logLike = prop_logLike; % assign the current log-likelihood to be the proposed one
                   
                   acMoves(k) = acMoves(k) + 1; % add one count to the number of accepted moves along parameter axis k
                   
                   all_D1s(l+1) = propTheta(1); all_D2s(l+1) = propTheta(2); all_p12s(l+1) = propTheta(3); all_p21s(l+1) = propTheta(4); % store new values
                   
               else
                    all_D1s(l+1) = l_theta(1); all_D2s(l+1) = l_theta(2); all_p12s(l+1) = l_theta(3); all_p21s(l+1) = l_theta(4); % store "new" parameter values
               end
           end
       end
    end
  
    s1 = acMoves(1)/prMoves(1); s2 = acMoves(2)/prMoves(2); s3 = acMoves(3)/prMoves(3); s4 = acMoves(4)/prMoves(4); % store acceptance rates for every parameter

    D1s = all_D1s(burnInTime:end); D2s = all_D2s(burnInTime:end); p12s = all_p12s(burnInTime:end); p21s = all_p21s(burnInTime:end); % relevant values
    est_D1 = mean(D1s); est_D2 = mean(D2s); est_p12 = mean(p12s); est_p21 = mean(p21s); % calculate averages after thermalisation 
  
    averaged = [est_D1, est_D2, est_p12, est_p21]; % output the averaged values after thermalisation

    % Compute error bars from coverage intervals
    % sort parameter values
    D1sorted = sort(D1s); D2sorted = sort(D2s); p12sorted = sort(p12s); p21sorted = sort(p21s);
    % number of elements to be removed from each side of the interval
    Nrem = round(numel(D1sorted)*0.025);

    errorLower = [D1sorted(1+Nrem) D2sorted(1+Nrem) p12sorted(1+Nrem) p21sorted(1+Nrem)];
    errorUpper = [D1sorted(end-Nrem) D2sorted(end-Nrem) p12sorted(end-Nrem) p21sorted(end-Nrem)]; 

    % ------------------ plot MCMC walk --------------------
%     
%     MCMCstepVec = 0:1:num_MC_steps-1;
% 
%     figure
%     plot(MCMCstepVec, D1s, 'LineWidth', 1)
%     hold on
%     plot(MCMCstepVec, D2s, 'LineWidth', 1)
%     xlabel('MCMC steps')
%     ylabel('Diffusion constant (\mum^2/s)')
%     legend('D_{1}', 'D_{2}') 
%     hold off
% 
%     figure
%     plot(MCMCstepVec, p12s, 'LineWidth', 1)
%     hold on
%     plot(MCMCstepVec, p21s, 'LineWidth', 1)
%     xlabel('MCMC steps')
%     ylabel('Transition probabilities')
%     legend('p_{12}', 'p_{21}')
%     hold off
end