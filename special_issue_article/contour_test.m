%% script for testing to plot an iso-contour


numExp = 1;
nTraj = 10;         % number of trajectories
N = 3E2;            % number of displacements
tau = 5E-3;
D1 = 10;
D2 = 1;
% sigmaB1vec = sqrt(2*D1vec*tau);
% sigmaB2vec = sqrt(2*D2vec*tau);
p12 = 0.05;          
p21 = 0.05;     
data = cell(1,nTraj);

% D1BayesEst = [0]; % hold bayesian parameter estimates
% BayesError = cell(1,numExp); % hold error estimates
% BayesError2std = cell(1,numExp); % hold 2 standard deviation error estimates
% BayesError3std = cell(1,numExp); % hold 2 standard deviation error estimates

for idx = 1:numExp

    sigmaB1 = sqrt(2*D1*tau); sigmaB2 = sqrt(2*D2*tau);

    % generate a data set
    for idxTrack = 1:nTraj  
        stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
        data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
    end
    
    % ----------------- do bayesian analysis ------------------------
    D1_max = D1+0.5*D1; % known upper limit of D1
    D1_min = max(D1-0.5*D1,0.01); % known lower limit of D1
    D2_max = D2+0.5*D2; % known upper limit of D2
    D2_min = max(D2-0.5*D2,0.01); % known lower limit of D2
    p12_max = min(5*p12,1); % known upper limit of p12
    p12_min = min(1/5*p12,0.01); % known lower limit of p12
    p21_max = min(5*p21,1); % known upper limit of p21
    p21_min = min(1/5*p21,0.01); % known lower limit of p21 
    % algorithm parameters
    nLive = 300; % number of live points
    StopRatio = 1E-4; % stop criterion for evidence
    priorLimits = [D1_min,D1_max,D2_min,D2_max,p12_min,p12_max,p21_min,p21_max];

    [finalSeq, thetaMLE, logZ] = utilB.nestedsampling2D(nLive, StopRatio, priorLimits, tau, data);
    [logZ_error, thetaBayes, errorBayes] = utilB.bayesianestimate(finalSeq,logZ,nLive,false);
    BayesEst = thetaBayes; % store parameter estimates
    BayesError = errorBayes; % store error estimates

end

%%

seqLen = length(finalSeq);
nBins = round(sqrt(seqLen));
% extract data1 and data2, i.e. all values of D1 and D2
dataT = zeros(seqLen,4); % all coordinate values
w = zeros(1,seqLen); % corresponding weights
logLVec = zeros(seqLen,1); % likelihoods
for i = 1:seqLen
    
    dataT(i,:) = finalSeq(i).pos; % parameter values
    w(i) = finalSeq(i).postWt; % posterior weights
    logLVec(i) = finalSeq(i).logL; % log-likelihood values
end
data1 = dataT(:,1); % 1st coordinate
data2 = dataT(:,2); % 2nd coordinate
data3 = dataT(:,3); % 3rd coordinate
data4 = dataT(:,4); % 4th coordinate

% nBins = 100; % for testing, otherwise pick as above
[xbar,binCounts] = utilB.bincounts(nBins, BayesEst(1), BayesError(1), data1, w);
bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram

figure;
[ybar,binCounts] = utilB.bincounts(nBins, BayesEst(2), BayesError(2), data2, w);
bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram


% [xbar,binCounts] = utilB.bincounts(nBins, BayesEst(3), BayesError(3), data3, w);
% bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
% [xbar,binCounts] = utilB.bincounts(nBins, BayesEst(4), BayesError(4), data4, w);
% bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram


%% posterior contour from histogram

plot1 = errorbar(BayesEst(2),BayesEst(1),BayesError(1),BayesError(1),BayesError(2),BayesError(2),...
    'o','Color','#0072BD','MarkerSize',3,'MarkerFaceColor','#0072BD','LineWidth',1.0);
hold on

plot0 = plot(D2,D1,'o','MarkerSize',4,'MarkerFaceColor','#D95319','MarkerEdgeColor','none');
hold on

[xbar1,xbar2,binCounts] = utilB.isocontour(nBins, [BayesEst(1),BayesEst(2)], [BayesError(1),BayesError(2)], data1, data2, w);
% bar3(binCounts) % 2-parameter histogram
bc = binCounts;
% bc = binCounts(:,1:40);
% bar3(bc)
Pmax = max(binCounts,[],'all'); % maximum of posterior
lv = Pmax*exp(-1/2);
% lv = Pmax*0.4;
% contour(xbar2, xbar1, bc, [lv,lv])
% contour(xbar2, xbar1, bc)

% interpPoints = 50;
% [xq,yq] = meshgrid(linspace(min(xbar2),max(xbar2),interpPoints), linspace(min(xbar1),max(xbar1),interpPoints));
% Pq = griddata(xbar2,xbar1,binCounts,xq,yq, 'cubic');
% mesh(xq,yq,Pq)
% contour(xq, yq, Pq, [lv,lv])

bcGF = imgaussfilt(bc,1.5);
PmaxGF = max(bcGF,[],'all'); % define maximum after smoothening
% contour(xbar2,xbar1,bcGF)
lv1 = PmaxGF*exp(-1/2);
lv2 = PmaxGF*exp(-2/2);
lv3 = PmaxGF*exp(-3/2);
contour(xbar2,xbar1,bcGF, [lv1,lv1],'LineWidth',1.0,'LineColor',[0 0 1]);
hold on
contour(xbar2,xbar1,bcGF, [lv2,lv2],'LineWidth',1.0,'LineColor',[0 .5 1])
hold on
contour(xbar2,xbar1,bcGF, [lv3,lv3],'LineWidth',1.0,'LineColor',[0 .7 1])


% surf(bcGF)

%% likelihood contour from formula - ignoring marginalization
% 
% plot1 = errorbar(BayesEst(2),BayesEst(1),BayesError(1),BayesError(1),BayesError(2),BayesError(2),...
%     'o','Color','#0072BD','MarkerSize',3,'MarkerFaceColor','#0072BD','LineWidth',1.0);
% hold on
% 
% nPoints = 100; % increase this later
% x = linspace(D2_min,D2_max,nPoints);
% y = linspace(D1_min,D1_max,nPoints);
% Z = zeros(nPoints,nPoints);
% 
% for idx = 1:nPoints
%     for idy = 1:nPoints
% 
%         Z(idx,idy) = utilB.EnsLikelihood2D([x(idx) y(idy) BayesEst(3) BayesEst(4)], tau, data);
%     end
% end
% % contour(x,y,Z)
% 
% Lmax = max(Z,[],'all'); % maximum of likelihood % NOTE: this is log-likelihood!
% % lv = Lmax*exp(-1/2); % 1 std level
% lv = Lmax-1/2; % 1 std in log-space
% contour(x,y,Z,[lv,lv])
% hold on
% lv2 = Lmax-2/2; % 2 std in log-space
% contour(x,y,Z,[lv2,lv2])
% hold on
% lv3 = Lmax-3/2; % 3 std in log-space
% contour(x,y,Z,[lv3,lv3])
% legend({'1\sigma','2\sigma','3\sigma'})


%% posterior contour from formula, maximum a posteriori, ignoring marginalization

% plot1 = errorbar(BayesEst(2),BayesEst(1),BayesError(1),BayesError(1),BayesError(2),BayesError(2),...
%     'o','Color','#0072BD','MarkerSize',3,'MarkerFaceColor','#0072BD','LineWidth',1.0);
% hold on
% 
% nPoints = 100; % increase this later
% x = linspace(D2_min,D2_max,nPoints);
% y = linspace(D1_min,D1_max,nPoints);
% Z = zeros(nPoints,nPoints);
% 
% for idx = 1:nPoints
%     for idy = 1:nPoints
% 
%         Z(idx,idy) = utilB.EnsLikelihood2D([x(idx) y(idy) BayesEst(3) BayesEst(4)], tau, data) +...
%             log(1/(D1_max-D1_min)*1/(D2_max-D2_min)*1/(p12_max-p12_min)*1/(p21_max-p21_min))-logZ;
%     end
% end
% % contour(x,y,Z)
% 
% Pmax = max(Z,[],'all'); % maximum of likelihood % NOTE: this is log-likelihood!
% % lv = Lmax*exp(-1/2); % 1 std level
% lv = Pmax-1/2; % 1 std in log-space
% contour(x,y,Z,[lv,lv])
% hold on
% lv2 = Pmax-2/2; % 2 std in log-space
% contour(x,y,Z,[lv2,lv2])
% hold on
% lv3 = Pmax-3/2; % 3 std in log-space
% contour(x,y,Z,[lv3,lv3])
% legend({'1\sigma','2\sigma','3\sigma'})
% 




%%
% % iso-curve of posterior
% hold on
% % figure;
% Pmax = max(binCounts,[],'all'); % maximum of posterior
% bcNew = bc;
% bcNew(bc<Pmax*exp(-1/2)) = 0;
% % bcNew(abs(bcNew-Pmax*exp(-1/2)) > 0.0002) = 0; % include a numerical error tolerance % NOTE: tolerance needs to be chosen carefully
% 
% % bar3(binCountsNew) % 2-parameter histogram
% [row,col] = find(bcNew~=0);
% % scatter(xbar1(row),xbar2(col))
% scatter(xbar2(col),xbar1(row))
% 
% % plot(xbar2(col),xbar1(row))






