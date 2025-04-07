function [n11,n12,n21,n22 ] = calculate_nij(binaryStateSeq)
    %
    % Estimates the number of transitions, n_{ij} (with i,j=1,2),
    % for a given state sequence
    % 
    % Input parameters:
    % binaryStateSeq = binary state sequence of type [1 2 1 2 1 1 2 2 ... ]
    %
    % Output parameters:
    % n11 = number of "dimers" '11' in the state sequence
    % n12 = number of "dimers" '12' in the state sequence
    % n21 = number of "dimers" '21' in the state sequence
    % n22 = number of "dimers" '22' in the state sequence
    %
    % Comment: using n_{ij} we can estimate the transition probabilities
    % p12, and p21: the maximum likelihood estimates for p12 and p21 
    % are obtained using Eq. (10) in 
    % Singer, Philipp, et al. "Detecting memory and structure in human 
    % navigation patterns using Markov chain models of varying order." 
    % PloS one 9.7 (2014): e102070.
    % 

    diffX = diff(binaryStateSeq);
    n12 = sum(diffX == -1); % number of 1-2 transition in time series
    n21 = sum(diffX == 1); % number of 2-1 transition in time series

    sumX = binaryStateSeq(1:end-1)+binaryStateSeq(2:end);
    n11 = sum(sumX == 2); % number of 1-1 transition in time series
    n22 = sum(sumX == 4); % number of 2-2 transition in time series



end

