function [mean1X , std1X, mean2X, std2X] = ...
    determine_segment_length_mean_and_std( sigmaB , sThresh , N   )

    %
    % Mean and variance for length of thresholded displacements
    % for 2D Brownian motion 
    %
    % Input:
    % sigmaB = variance of the Brownian steps
    % threshQ = threshold as applied to the Brownian displacement
    % N = length of the simulation
    %
    % Output:
    % mean1X = mean length of the 1X segments
    %         (1X = segments of type '1',
    %           '11', '111', etc)
    % std1X = std in length of the 1X segments
    % mean2X = mean length of the 2X segments
    %         (2X = segments of type '2',
    %           '22', '222', etc)
    % std2X = std in length of the 2X segments
    %
    % Dependencies: brownian_displacements_2d, find_1X_segments, find_2X_segments
    %
    import('utilF.brownian_displacements_2d')
    import('utilF.find_1X_segments')
    import('utilF.find_2X_segments')
    
    % Simulate a one-state diffusion process
    stateVec = ones(1,N);
    data = brownian_displacements_2d(stateVec,sigmaB,sigmaB);
    
       
    % Apply a threshold
    stateVecEst = zeros(1,length(data));
    stateVecEst(data <= sThresh) = 1;
    stateVecEst(data > sThresh) = 2;
 
  
   %%
    % Find all segments of type 1X = '1','11, ..., '11111' etc 
    % calculate the mean and variance for the lengths
    %
    [segment1XStartPos, segment1XLengths] = find_1X_segments(stateVecEst);
    mean1X = mean(segment1XLengths);
    std1X = std(segment1XLengths);


    %%
    % Find all segments of type 2X = '2','22, ..., '22222' etc and 
    %  calculate the mean and variance for the lengths
    %
    [segment2XStartPos, segment2XLengths] = find_2X_segments(stateVecEst);
    mean2X = mean(segment2XLengths);
    std2X = std(segment2XLengths);


    %% Plot histograms 
% 
%     % 1X
% 
%     % Make histograms for segment lengths of type 1X and compare to theory.
%     binEdges1X = 0.5:1:max(segment1XLengths)+0.5;
%     H1X  = histcounts(segment1XLengths,binEdges1X);
%     binCenters1X = binEdges1X(1:end-1) + 0.5;
%     deltaBin= binEdges1X(2)-binEdges1X(1);
% 
%     figure
%     bar(binCenters1X,H1X);
%     title('1X segments')
%     xlabel('segment lengths')
%     ylabel('counts')
% 
% 
%     % 2X segments
% 
%     % Make histograms for segment lengths of type 2X and compare to theory.
%     binEdges2X = 0.5:1:max(segment2XLengths)+0.5;
%     H2X = histcounts(segment2XLengths,binEdges2X);
%     binCenters2X = binEdges2X(1:end-1) + 0.5 ;
%     deltaBin= binEdges2X(2)-binEdges2X(1);
%     figure
%     bar(binCenters2X,H2X);
%     title('2X segments')
%     xlabel('segment lengths')
%     ylabel('counts')


  
end

