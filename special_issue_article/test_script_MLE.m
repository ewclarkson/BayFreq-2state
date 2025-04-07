%% Simulate tracks and analyse them with a HMM 2-state model

% ----------------------- initialise constants ---------------------------

nTraj = 10;
N = 3E2;
D1  = 10;
D2 = 3;
p12 = 0.05;
p21 = 0.05;
tau = 5E-3;
sigmaB1 = sqrt(2*D1*tau); 
sigmaB2 = sqrt(2*D2*tau);

% ---------------- generate state sequences and tracks --------------------

% generate trajectories
for idxTrack = 1:nTraj  
    stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
    data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
end
% ------------------ analyse a track with HMM ---------------------------


% loglikeli_acc = 0; % to hold accumulated log likelihoods for all tracks together
% for n = 1:nTracks
%     
%     trackDisps = zeros(1,nSteps); % initialise array of all nSteps displacements
%     for j = 1:nSteps
%         trackDisps(j) = Distance(trackCell{n}(j,1), trackCell{n}(j,2), trackCell{n}(j+1,1), trackCell{n}(j+1,2));
%     end
%     
%     loglikeli_M = Markov_likelihood_v2(D1, D2, p12, p21, deltaT, trackDisps); % the log-likelihood for track n 
%     loglikeli_acc = loglikeli_acc + loglikeli_M; % add this log-likelihood to the sum of all such
% end

% loglikeli_acc2 = Ens_likelihood([D1, D2, p12, p21], deltaT, dispCell);

% --------------- optimise with an MCMC metropolis scheme ----------------
close all

num_MC_steps = 1E4+1;
burnInTime = 1E3;
scaleArr = [1.1, 0.14, 0.0008, 0.0007]; % set this to achieve acceptance rates of 20-40%
D1_a = 7; D1_b = 14; D2_a = 0.3; D2_b = 3; p12_a = 0.01; p12_b = 0.3; p21_a = 0.01; p22_b = 0.3;
guessArr = [D1_a, D1_b; D2_a, D2_b; p12_a, p12_b; p21_a, p22_b]; % hold limits on random guess for initial parameter values  

% [latest_est, mean_est, D1_est, D2_est, p12_est, p21_est, s1, s2, s3, s4] = MCMC(num_MC_steps, burnInTime, scaleArr, guessArr, tau, data);
[averaged, errorEst, s1, s2, s3, s4] = MCMC_v2(num_MC_steps, burnInTime, scaleArr, guessArr, tau, data)

% MCMCstepVec = 0:1:num_MC_steps-1;
% plot(MCMCstepVec, D1_est)
% hold on
% plot(MCMCstepVec, D2_est)
% %legend('D1', 'D2')
% hold on
% plot(MCMCstepVec, p12_est)
% hold on
% plot(MCMCstepVec, p21_est)
% %legend('p12', 'p21')
% legend('D1', 'D2', 'p12', 'p21')
% hold off

