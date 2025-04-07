


%
% Segment-length histograms for a pure Markov chain 
%

close all

%%
% Input parameters

                                            
N = 2E5;                % Number of sampling time intervals for "actual" data
p12 = 0.05;
p21 = 0.02;
lThresh1X = 5;        % length threshold for regions 1X = '1','11','111', etc
lThresh2X = 5;        % length threshold for regions 2X = '2','22','222', etc
nRecurSpinFlip = 10;  % number of times we run the spin flip procedure
                      % (alternatively, we could simply run this procedure 
                      % until "convergence")


% Derived parameters
p11 = 1 - p12;
p22 = 1 - p21;

%%
% Generate a Markov chain time series 
%

stateVec = utilF.twoState_Markov(p12, p21, N); % Erik's code

% % n_i = number of time step when we are in state i (i=1,2)
n1 = length(find(stateVec==1));
n2 = length(find(stateVec==2));
%n1+n2 % should equal N
fEst = n1/(n1+n2) % fraction of time when we are in state 1
fTrue = p21/(p12+p21)

%%
% Find all segments of type 1X = '1','11, ..., '11111' etc and 
% store their start positions and  lengths 
%
[segment1XStartPos, segment1XLengths] = utilF.find_1X_segments(stateVec);

% Make histograms for segment lengths of type 1X and compare to theory.
binEdges1X = 0.5:1:max(segment1XLengths)+0.5;
H1X  = histcounts(segment1XLengths,binEdges1X);
binCenters1X = binEdges1X(1:end-1) + 0.5;
deltaBin= binEdges1X(2)-binEdges1X(1);

figure
bar(binCenters1X,H1X);
hold on;

nStop = 100;
n=1:1:nStop;

% Theory
histTheory1X = deltaBin*sum(H1X)*(1-p11)*p11.^(n-1);

plot(n,histTheory1X,'linewidth',2);
title('1X segments')
xlabel('segment lengths')
ylabel('counts')


%%
% Find all segments of type 2X = '2','22, ..., '22222' etc and 
% store their start positions and  lengths 
%
[segment2XStartPos, segment2XLengths] = utilF.find_2X_segments(stateVec);

% Make histograms for segment lengths of type 2X and compare to theory.
binEdges2X = 0.5:1:max(segment2XLengths)+0.5;
H2X = histcounts(segment2XLengths,binEdges2X);
binCenters2X = binEdges2X(1:end-1) + 0.5 ;
deltaBin= binEdges2X(2)-binEdges2X(1);

figure
bar(binCenters2X,H2X);
hold on;

nStop = 100;
n=1:1:nStop;

% Theory
histTheory2X = deltaBin*sum(H2X)*(1-p22)*p22.^(n-1);

plot(n,histTheory2X,'linewidth',2);
title('2X segments')
xlabel('segment lengths')
ylabel('counts')


% 

%% 
% Apply spin flip procedure. See how the histograms change

flippedStateVec = stateVec;
if nRecurSpinFlip > 0
    for idxSpinFlip = 1:nRecurSpinFlip
        flippedStateVec  = utilF.apply_spin_flip( flippedStateVec , lThresh1X, lThresh2X );
    end
end

%
% Find all segments of type 2X = '2','22, ..., '22222' etc and 
% store their start positions and  lengths 
%
[segment1XStartPosFlip, segment1XLengthsFlip] = utilF.find_1X_segments(flippedStateVec);

% Make histograms for segment lengths of type 1X and compare to theory.
binEdges1XFlip = 0.5:1:max(segment1XLengthsFlip)+0.5;
H1XFlip  = histcounts(segment1XLengthsFlip,binEdges1XFlip);
binCenters1XFlip = binEdges1XFlip(1:end-1) + 0.5;

figure
bar(binCenters1XFlip,H1XFlip);
hold on;

plot(n,histTheory1X,'linewidth',2);
title('1X segments after spin flip')
xlabel('segment lengths')
ylabel('counts')


%
% Find all segments of type 2X = '2','22, ..., '22222' etc and 
% store their start positions and  lengths 
%
[segment2XStartPosFlip, segment2XLengthsFlip] = utilF.find_2X_segments(flippedStateVec);

% Make histograms for segment lengths of type 2X and compare to theory.
binEdges2XFlip = 0.5:1:max(segment2XLengthsFlip)+0.5;
H2XFlip = histcounts(segment2XLengthsFlip,binEdges2XFlip);
binCenters2XFlip = binEdges2XFlip(1:end-1) + 0.5 ;

figure
bar(binCenters2XFlip,H2XFlip);
hold on;

plot(n,histTheory2X,'linewidth',2);
title('2X segments after spin flip')
xlabel('segment lengths')
ylabel('counts')


%% Plot stateVec, before and after spin flip

figure
plot(stateVec(1:1000),'k-','linewidth',3); hold on;
plot(flippedStateVec(1:1000)-2.1,'b-','linewidth',3)

%% Check n12 before and after spin flip

[n11Est,n12Est,n21Est,n22Est] = utilF.calculate_nij(stateVec);
disp(['n21 before spin flip = ',num2str(n21Est)])

[n11Flip,n12Flip,n21Flip,n22Flip] = utilF.calculate_nij(flippedStateVec);
disp(['n21 after spin flip = ',num2str(n21Flip)])

n1X = length(find(segment1XLengths <= lThresh1X));
disp(['number of short 1X segments before spin flip = ', num2str(n1X)])

n2X = length(find(segment2XLengths <= lThresh2X));
disp(['number of short 2X segments before spin flip = ', num2str(n2X)])

n1XFlip = length(find(segment1XLengthsFlip <= lThresh1X));
disp(['number of short 1X segments after spin flip = ', num2str(n1XFlip)])

n2XFlip = length(find(segment2XLengthsFlip <= lThresh2X));
disp(['number of short 2X segments after spin flip = ', num2str(n2XFlip)])

% Theoretical estimate for the number of short segments:
% sum_{n=1}^lThresh = deltaBin*sum(H1X)*(1-p11)*p11.^(n-1) 
% = deltaBin*sum(H1X)*(1-p11^lThresh) = deltaBin*sum(H1X)*(1-(1-p12)^lThresh)
% = / assume p12 << 1 /
% = deltaBin*sum(H1X)*(1-(1- lThresh* p12)) = deltaBin*sum(H1X)*lThresh* p12
% deltaBin = 1. Now, sum(H1X) and sum(H2X) are simply equal to the 
% number of state transition in the time series, i.e., = n21.

[n11,n12,n21,n22 ] = utilF.calculate_nij(stateVec);

nShortTheory1X = n21*(1-p11^lThresh1X);
nShortTheory1XApprox = n21*lThresh1X*p12;
disp(['theoretical estimate for number of short 1X segments = ', num2str(nShortTheory1X)])
disp(['approximate theoretical estimate for number of short 1X segments = ', num2str(nShortTheory1XApprox)])


nShortTheory2X = n21*(1-p22^lThresh2X);
nShortTheory2XApprox = n21*lThresh2X*p21;
disp(['theoretical estimate for number of short 2X segments = ', num2str(nShortTheory2X)])
disp(['approximate theoretical estimate for number of short 2X segments = ', num2str(nShortTheory2XApprox)])

n21Diff = n21Est - n21Flip
% We here assume that for every missing 1X or 2X segment there is a missing
% '2-1' event in the flipped time series (we also assume that p21, p12 <<1)
n21DiffTheory = nShortTheory1XApprox + nShortTheory2XApprox

