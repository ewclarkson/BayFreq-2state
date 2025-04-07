%
% 
% Estimate transition probabilities, p12 and p21,
% from a noisy Markov chain process using a frquentist-approach.
%
% import('twoState_Markov')
% addpath(genpath('Programs/special_issue_article'))

close all
clear all

rng(4) % for reproducable comparisons

%%
% Input parameters

% Simulation parameters
N=1E3;              % Number of sampling time intervals for "actual" data
sigmaB1 = 0.2;      % Std dev for dx and dy when in state 1
sigmaB2 = 0.8;      % Std dev for dx and dy when in state 2 % NOTE: switch so that D1 > D2, as we do for the Bayesian analysis
p12 = 0.05;         % transition probability from state 1 to state 2
p21 = 0.05;         % transition probability from state 2 to state 1  
p11 = 1-p12;
p22 = 1-p21;

%%
% Generate a Markov chain time series 
%
stateVec = utilF.twoState_Markov(p12, p21, N); % Erik's code

fTrue = p21/(p12+p21);   % fraction of time when we are in state 1

% Maximum likelihood estimates for p12 and p21 time series
% [n11,n12,n21,n22] = utilF.calculate_nij(stateVec);
% p12Est = n12/(n11+n12)



%%

% Generate actual data for a 2d Brownian walk
data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
figure
plot(data,'k-'); hold on;
plot(stateVec-2.1,'b-','linewidth',3)


%%
% Now turn the time-series into a noisy version 
% of the "two-state" process by applying the threshold
% above to the displacements, s.

% Determine optimal threshold 
sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2);
                   % threshold for the displacement (everything below 
                   % the threshold is initially deemed to be in state 1, and 
                   % everything above is deemed to be in state 2)
stateVecEst = zeros(1,length(data));
stateVecEst(data <= sThresh) = 1;
stateVecEst(data > sThresh) = 2;

plot(stateVecEst-4.1,'r-','linewidth',3)

%% Determine optimal length threshold

lThresh = 2; % TODO: set this so as to minimise corrections

% minimise the sum of corrections:
% q*(N-1)*(1-q)^(lThresh+1)+()

% lThresh2X = 2;
% lThresh1X = 2;

%% Apply spin flip procedure

flippedStateVec = stateVecEst;
flippedStateVec  = utilF.apply_spin_flip(flippedStateVec, lThresh, lThresh);
% if nRecurSpinFlip > 0
%     for idxSpinFlip = 1:nRecurSpinFlip
%         flippedStateVec  = utilF.apply_spin_flip( flippedStateVec, lThresh1X, lThresh2X );
%     end
% end
plot(flippedStateVec-6.1,'m-','linewidth',3)

%%
% Estimate n12, n21, n11 and n22 for the flipped version of the state sequence
[n11Est,n12Est,n21Est,n22Est] = utilF.calculate_nij(flippedStateVec);
%n21Est + n12Est + n11Est + n22Est  % should = N-1, i.e., all nearest neighbour 
                       % pairs: (1,2), (2,3), (3,4), ..., (N-1,N). 



%%
% Estimate level crossings for one-state systems
% (used to correct for false positives below)
[density12State] = utilF.density_level_crossings(sigmaB1, sThresh, 1E5, ...
                                          lThresh, lThresh);
[density21State] = utilF.density_level_crossings(sigmaB2, sThresh, 1E5, ...
                                          lThresh, lThresh);


%%

% Estimate p12 and p21 by correcting for 
% the expected number of false positive 
% level-crossing events left after state flipping, 
% and for accidently flipped crossings.

n21FalsePositivesEst = density12State*fTrue*N + density21State*(1-fTrue)*N
% n21Corrected = n21Est - n21FalsePositivesEst; % not used

% NEW METHOD:
q = raylcdf(sThresh,sigmaB1);
% q2 = raylcdf(sThresh,sigmaB2);
n12False = q*(1-q)*(N-1)*(1-q)^lThresh
% n21False = q2*(1-q2)*(N-1)*q2^lThresh*(1-fTrue);
n21False = n12False;
% n21FalsePositivesEstNew = q*(1-q)*(N-1)*(1-q)^lThresh;
% abs(n12False-n21FalsePositivesEst)/n12False % relative error

% Solve non-linear equation for p21. in case there are several solutions 
% (empirically there seems to be two solutions in the range [0,1])
% we choose the smallest solution (NOT CLEAR WHY...). Here, we use a
% very simplistic approach for finding the root.

D12 = (n21Est - n21FalsePositivesEst)/(fTrue*N); % simple estimate of p12
% q = lThresh + fTrue/(1-fTrue)*lThresh; % for analysis of solution

% [n11obs,n12obs,n21obs,n22obs] = utilF.calculate_nij(stateVecEst);
% D12 = (n21obs - n21FalsePositivesEst)/(fTrue*N);

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

disp(['p_{12}, true = ',num2str(p12),', estimated and doubly corrected = ',num2str(p12DoublyCorrected)])
disp(['p_{21}, true = ',num2str(p21),', estimated and doubly corrected = ',num2str(p21DoublyCorrected)])

legend('s(t)','true state sequence','noisy state sequance','spin-flipped state sequence')

%% 
% OLD
% 
% % estimated number of true 2-1- events which are missed due to the spin
% % flip procedure, see test_script_pure_markov_chain, then correct for this
% q =  lThresh1X/(fTrue*N) + lThresh2X/((1-fTrue)*N);
%         
% n21DoublyCorrectedOld = n21Corrected + q*n21Corrected^2;
% 
% n21DoublyCorrected = (1/(2*q))*(1 - sqrt(1 - 4*q*n21Corrected ));
% 
% 
% % n11Corrected = n11Est + fTrue*N*(1 - density11State1)...
% %          - (1-fTrue)*N*(1 - density22State2 - 2*density21State2); 
% % 
% % n22Corrected = n22Est - fTrue*N*(1 - density11State1 - 2*density21State1 )...
% %     + (1-fTrue)*N*(1 - density22State2); 
% % 
% 
% 
% % disp(['n_{11}, true = ',num2str(n11),', estimated = ',num2str(n11Est),...
% %     ', corrected = ',num2str(round(n11Corrected))])
% disp(['n_{21}, true = ',num2str(n21),', estimated = ',num2str(n21Est),...
%     ', corrected = ',num2str(round(n21Corrected)), ...
%     ' doubly corrected = ',num2str(round(n21DoublyCorrected))])
% % disp(['n_{22}, true = ',num2str(n22),', estimated = ',num2str(n22Est),...
% %     ', corrected = ',num2str(round(n22Corrected))])
% 
% % % Double-check: n11 + n22 + n21 + n12 = N - 1
% % n11Corrected + n22Corrected + 2*n21Corrected
% % 
% % fCorrected = (n11Corrected + n21Corrected)/(N-1)
% 
% 
% 
% % Bring the estimated values for n_{ij} within the allowed range
% % [0 , N-1]
% % n11Corrected = min(max(0,n11Corrected),N-1);
% %n21Corrected = min(max(0,n21Corrected),N-1);
% % n22Corrected = min(max(0,n22Corrected),N-1);
% 
% 
% % 
% %  p12Corrected = n21Corrected/(n11Corrected+n21Corrected);
% %  p21Corrected = n21Corrected/(n21Corrected+n22Corrected);
% p12DoublyCorrected = n21DoublyCorrected/(fTrue*N);
% p21DoublyCorrected = n21DoublyCorrected/((1-fTrue)*N);
% 
% disp(['p_{12}, true = ',num2str(p12),', estimated and doubly corrected = ',num2str(p12DoublyCorrected)])
% disp(['p_{21}, true = ',num2str(p21),', estimated and doubly corrected = ',num2str(p21DoublyCorrected)])
% 
% 
% 
%     
