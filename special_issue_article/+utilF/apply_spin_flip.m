function [ flippedStateSeq ] = apply_spin_flip( binaryStateSeq, lThresh1X, lThresh2X )

    %
    %  Apply a spin flip procedure to a state sequence 
    %  of type 12112222111122 until convergence, i.e. all
    %  small segments are eliminated.
    % 
    % Input parameters:
    % binaryStateSeq = binary state sequence of type [1 2 1 2 1 1 2 2 ... ]
    % lThresh1X = all segments of type '1', '11', '111' which are 
    %             of length <= lThres1X are "flipped" to become segments
    %             '2', '22', '222' etc. Set lThresh1X = 0 to avoid spin
    %             flip.
    % lThresh2X = all segments of type '2', '22', '222' which are 
    %             of length <= lThres2X are "flipped" to become segments
    %             '1', '11', '111' etc. Set lThresh1X = 0 to avoid spin
    %             flip.
    %
    % Output:
    % flippedStateSeq = state sequence where all small segments 
    %                   have been "flipped"
    %             
    % Dependencies: find_1X_segments, find_2X_segments 
    %
    import('utilF.find_1X_segments')
    import('utilF.find_2X_segments')
     
    flippedStateSeq = binaryStateSeq;
    stateDiff = true; % true if there are small segments left, false otherwise

    %  Find all segments of type 1X = '1','11, ..., '11111'
    [segment1XStartPos, segment1XLengths] = find_1X_segments(binaryStateSeq);
        
    % Find all segments of type 2X = '2', '22', '222'  
    [segment2XStartPos, segment2XLengths] = find_2X_segments(binaryStateSeq);

    if isempty(segment1XLengths) || isempty(segment2XLengths) % then no flipping is required
        return % return thresholded state sequence
    end

    while stateDiff % do flipping until there are no small segments left

        % Apply a "spin flip" to all 1X segments of length <= lThresh1X
        %
        if segment1XLengths(1) <= lThresh1X
            if segment1XStartPos(1) == 0 % should be == 1 ? Why do the first step separately?
                flippedStateSeq(1:segment1XLengths(1)) = 2;
            else
                flippedStateSeq(segment1XStartPos(1):segment1XStartPos(1)+segment1XLengths(1)-1) = 2;
            end
        end
        for k = 2:length(segment1XStartPos)
            if segment1XLengths(k) <= lThresh1X
               flippedStateSeq(segment1XStartPos(k):segment1XStartPos(k)+segment1XLengths(k)-1) = 2;
            end
        end
    
        if all(flippedStateSeq==2) && all(segment1XLengths<=lThresh1X) && all(segment2XLengths<=lThresh2X) % all states are 2's
            return % all segment lengths are shorter than the threshold and so flipping will not converge. Return sequence of only 2's. 
        end
    %     plot(binaryStateSeq,'k-'); hold on;
    %     plot(flippedStateSeq-2.1,'r-','linewidth',3)
    
        %               
        % Apply a "spin flip" to all 2X segments of length <= lThresh2X
        %
        if segment2XLengths(1) <= lThresh2X
            if segment2XStartPos(1) == 0 % should be == 1 ?
                flippedStateSeq(1:segment2XLengths(1)) = 1;
            else
                flippedStateSeq(segment2XStartPos(1):segment2XStartPos(1)+segment2XLengths(1)-1) = 1;
            end
        end
        for k = 2:length(segment2XStartPos)
            if segment2XLengths(k) <= lThresh2X
               flippedStateSeq(segment2XStartPos(k):segment2XStartPos(k)+segment2XLengths(k)-1) = 1;
            end
        end
        
        % Check convergence of flipping
        [segment1XStartPos, segment1XLengths] = find_1X_segments(flippedStateSeq);
        [segment2XStartPos, segment2XLengths] = find_2X_segments(flippedStateSeq);

        if isempty(segment1XLengths(segment1XLengths<=lThresh1X)) && isempty(segment2XLengths(segment2XLengths<=lThresh2X)) 
            stateDiff = false; % convergence reached: no small segments left
        end
    end

%     plot(flippedStateSeq-4.1,'b-','linewidth',3)
%     hold off 
%     

end
