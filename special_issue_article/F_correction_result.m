%% Compare various estimations using various corrections within tpe-RMM


% ----------- Set data parameters -------------------------
numExp = 4; % number of data sets to analyse
nTraj = 100;         % number of trajectories
N = 3E2;            % number of displacements
tau = 5E-3;
D1vec = [10,10,10,10]; %linspace(10,10,numExp);
D2vec = [1,3,1,3]; %linspace(1,4,numExp);
p12vec = [0.01 0.01 0.05 0.05]; %linspace(0.01,0.01,numExp);          
p21vec = [0.01 0.01 0.05 0.05]; %linspace(0.01,0.01,numExp);
% fTrueVec = p21vec./(p12vec+p21vec);
data = cell(1,nTraj);
% ----------------------------------------------------------

freqEst = cell(1,numExp); % using all corrections
freqEstNo = cell(1,numExp); % using no corrections
freqEstFP = cell(1,numExp); % using only the false positive correction
freqEstFN = cell(1,numExp); % using only the false negative correction

for idx = 1:numExp
   
    D1 = D1vec(idx); D2 = D2vec(idx);
    p12 = p12vec(idx); p21 = p21vec(idx); fTrue = p21/(p12+p21);
    sigmaB1 = sqrt(2*D1*tau); sigmaB2 = sqrt(2*D2*tau);

    % generate a data set
    for idxTrack = 1:nTraj  
        stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
        data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
    end

    D1_max = 20; %D1+0.5*D1; % known upper limit of D1
    D1_min = 3; %max(D1-0.5*D1,0.01); % known lower limit of D1
    D2_max = 8; %D2+0.5*D2; % known upper limit of D2
    D2_min = 0.3; %max(D2-0.5*D2,0.01); % known lower limit of D2
    p12_max = 0.2; %min(5*p12,1); % known upper limit of p12
    p12_min = 0.005; %min(1/5*p12,0.01); % known lower limit of p12
    p21_max = 0.2; %min(5*p21,1); % known upper limit of p21
    p21_min = 0.005; %min(1/5*p21,0.01); % known lower limit of p21 

    % ----------- do frequentist analysis ---------------
    tic

    [estInitial,~] = utilF.Rayleigh_2state((D1_min+D1_max)/2, (D2_min+D2_max)/2, (p21_max/(p21_max+p12_min)+p21_min/(p21_min+p12_max))/2, tau, data); 
    freq_fEst = estInitial(1); freq_D1Est = estInitial(2); freq_D2Est = estInitial(3); % estimates of f,D1,D2
    [freq_p12Est, freq_p21Est, p12No, p21No, p12FP, p21FP, p12FN, p21FN] =...
       utilF.wrapperFreq(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau));

    freqEst{idx} = [freq_D1Est,freq_D2Est,freq_p12Est,freq_p21Est]; % store parameter estimates
    freqEstNo{idx} = [freq_D1Est,freq_D2Est,p12No,p21No]; % store parameter estimates
    freqEstFP{idx} = [freq_D1Est,freq_D2Est,p12FP,p21FP];
    freqEstFN{idx} = [freq_D1Est,freq_D2Est,p12FN,p21FN];

    toc
end


%%

f = figure('Position',[500 100 580 640]);
tiledlayout(numExp,2,'TileSpacing','compact','Padding','compact');
panelNames = ["A." "B." "C." "D."];

for idx = 1:numExp
    for idy = 3:4 % plotting the result of the first experiment
    
        nexttile
        if idy == 3
            xlabel('p_{12}')
            yticks([])
            hold on
            xline(p12vec(idx),'--','LineWidth',1.5,'Color',[0,0,0])
            hold on
            xline(freqEst{idx}(3),'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
            hold on
            xline(freqEstNo{idx}(3),'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
            hold on
            xline(freqEstFP{idx}(3),'-','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330])
            hold on
            xline(freqEstFN{idx}(3),'-','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880])
            if idx == 1
                legend({'Ground truth','Full correction', 'No correction', 'FP correction',...
                    'FN correction'}, 'Location', 'Northwest')
            end
            text(-25,102,{panelNames(idx)},'FontSize',12,'unit','pixel'); % text box for indicating panel A,B,C

        else
            xlabel('p_{21}')
            yticks([])
            hold on
            xline(p21vec(idx),'--','LineWidth',1.5,'Color',[0,0,0])
            hold on
            xline(freqEst{idx}(4),'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
            hold on
            xline(freqEstNo{idx}(4),'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
            hold on
            xline(freqEstFP{idx}(4),'-','LineWidth',1.5,'Color',[0.3010 0.7450 0.9330])
            hold on
            xline(freqEstFN{idx}(4),'-','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880])
        end
    end
end


