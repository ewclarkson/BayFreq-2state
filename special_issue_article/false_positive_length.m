%% Find the PDF of false segment lengths

% Find out the PDF of the number false positives (1->2 events, number of segments)
close all
% Simulation parameters
N=1E4;              % Trajectory length
sigmaB1 = 0.5;      % Std dev for dx and dy when in state 1
sigmaB2 = 0.8;      % Std dev for dx and dy when in state 2 % NOTE: switch so that D1 > D2, as we do for the Bayesian analysis
% p12 = 0.1;
% p21 = 0.3;
% p22 = 1-p21; % follows from normalisation

nSim = 1E2; % number of different "experiments"

%% Find the PDF of segment lengths numerically

lengthVec = []; % vector of segment lengths of 2's
stateVec = ones(1,N);
sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % find threshold
q = raylcdf(sThresh,sigmaB1); % probability of being correctly labelled
for i = 1:nSim
%     stateVec = utilF.twoState_Markov(p12, p21, N); % state sequence
    data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2); % generate displacements
    % Threshold
    stateVecEst = zeros(1,length(data));
    stateVecEst(data <= sThresh) = 1;
    stateVecEst(data > sThresh) = 2;
    % Count the number of 1->2 events
    ind = find(stateVecEst==2); % indices of all 2's
    diffInd = diff(ind); %
    tempVec = [];
    tempLen = 0;
    for k = 1:length(ind)-1
        if ind(k+1) == ind(k)+1
            tempLen = tempLen+1; % add to this segment's length
        else
            if tempLen == 0 % a segment of length 1
%                 tempLen = tempLen+1; % add the final 2
                tempVec(end+1) = 1;
            else % end of a segment longer than 1
                tempLen = tempLen+1;
                tempVec(end+1) = tempLen;
                tempLen = 0; % reset
            end
        end
    end
    % include also the final segment
%     if stateVec(end) == 2
%         
%     end
    lengthVec = [lengthVec,tempVec]; % add segment lenghts
end

%% Plot the distribution and theory

% histogram(lengthVec,'Normalization','pdf');
histogram(tempVec,'Normalization','pdf');
hold on
lvec = 1:max(tempVec);
lenTheory = (1-q).^(lvec-1)*q;
plot(lvec,lenTheory,'LineWidth',1)


% 

















