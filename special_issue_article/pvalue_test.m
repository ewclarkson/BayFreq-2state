%% plot the decision boundary


% Simulation parameters
nTraj = 10;
N = 300;              % Number of sampling time intervals for "actual" data
sigmaB1 = .8;      % Std dev for dx and dy when in state 1
sigmaB2 = .4;      % 

p12 = 0.1;         % transition probability from state 1 to state 2
p21 = 0.1;         % transition probability from state 2 to state 1  
p11 = 1-p12;
p22 = 1-p21;

stateVec = utilF.twoState_Markov(p12, p21, N); % Erik's code
fTrue = p21/(p12+p21);   % fraction of time when we are in state 1
% data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
data = cell(1,nTraj);
for idxTrack = 1:nTraj  
    stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
    data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
end


%%
sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % displacement threshold
% q = raylcdf(sThresh,sigmaB1); % probability of a single displacement being correctly labelled in the thresholding
q = raylcdf(sThresh,sigmaB2); % use this if D1>D2, otherwise switch to sigmaB1 in the expression
pThresh = 0.03;
% log(pThresh)/log(1-q)
lThresh = round(log(pThresh)/log(1-q)); % set threshold based on a p-value of 0.03  

[p12DoublyCorrected, p21DoublyCorrected] = utilF.method_of_moments_thresh(lThresh, data, fTrue, sigmaB1, sigmaB2)

%%
pThresh = 0.01;
% log(pThresh)/log(1-q)
lThresh = round(log(pThresh)/log(1-q)); % set threshold based on a p-value of 0.03  

[p12DoublyCorrected, p21DoublyCorrected] = utilF.method_of_moments_thresh(lThresh, data, fTrue, sigmaB1, sigmaB2)
%%

% [p12DoublyCorrected, p21DoublyCorrected] = utilF.wrapperFreq({data}, fTrue, sigmaB1, sigmaB2)

%%
% % Generate actual data for a 2d Brownian walk
% data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
% sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2);
% 
% histogram(data,'Normalization','pdf','LineStyle','none')
% hold on
% xline(sThresh)
% 

