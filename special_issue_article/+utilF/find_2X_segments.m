

function [segment2XStartPos, segment2XLengths] = find_2X_segments(binaryStateSeq)

    % Find all segments of type 2X = '2','22', ..., '22222' etc 
    % in a binary state sequense of type 121122111222. Store 
    % their start positions and lengths of the segments. 
    %
    % Input parameters:
    % binaryStateSeq = binary state sequence of type [1 2 1 2 1 1 2 2 ... ]
    %
    % Output parameters:
    % segment2XStartPos = start positions of the segments
    % segmentXLengths = lengths of the segments
    %

    diffX = diff(binaryStateSeq);
    idxUp = find(diffX == 1);
    idxDown = find(diffX == -1);

    % if isempty(idxUp) || isempty(idxDown) % in case of no state switches
    %     segment2XStartPos = [];
    %     segment2XLengths = [];
    %     return
    % end

    if isempty(idxUp) && isempty(idxDown) % in case of no state switches
        segment2XStartPos = [];
        segment2XLengths = [];
        return
    end

    if isempty(idxUp) 
        if ~isempty(idxDown) % in case of no state switches
            segment2XStartPos = 1;
            segment2XLengths = idxDown;
            return
        end
    elseif isempty(idxDown)
        if ~isempty(idxUp) % in case of no state switches
            segment2XStartPos = idxUp+1;
            segment2XLengths = length(binaryStateSeq)-idxUp;
            return
        end
    end

    if idxUp(1) < idxDown(1) % first level-crossing is a down-crossing
                             % (level crossing  = crossing of the value 3/2, 
                             % i.e. a transition from 1 to 2 = down-crossing
                             % or from 2 to 1 = up-crossing)

        segment2XStartPos = idxUp + 1;

        segment2XLengths = zeros(1,length(segment2XStartPos));
        if idxUp(end) < idxDown(end) % last level-crossing is an down-crossing       
            segment2XLengths = idxDown-idxUp;        
        else % last level-crossing is an up-crossing       
            segment2XLengths(1:end-1) = idxDown-idxUp(1:end-1);
            segment2XLengths(end) = length(binaryStateSeq) - idxUp(end);
        end      

    else  % idxUp(1) > idxDown(1), i.e., first level-crossing is a down-crossing

        segment2XStartPos = [1 , idxUp + 1];

        segment2XLengths = zeros(1,length(segment2XStartPos));
        if idxUp(end) < idxDown(end) % last level-crossing is an down-crossing       
           segment2XLengths(1) = idxDown(1);
           segment2XLengths(2:end) = idxDown(2:end) - idxUp;    
        else % last level-crossing is an up-crossing      
            segment2XLengths(1) = idxDown(1);
            segment2XLengths(2:end-1) = idxDown(2:end)-idxUp(1:end-1);
            segment2XLengths(end) = length(binaryStateSeq) - idxUp(end);          
        end

    end

end

