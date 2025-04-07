
function [segment1XStartPos, segment1XLengths] = find_1X_segments(binaryStateSeq)
   
    %  Find all segments of type 1X = '1','11, ..., '11111' etc 
    % in a binary state sequense of type 121122111222. Store 
    % their start positions and lengths of the segments. 
    %
    % Input parameters:
    % binaryStateSeq = binary state sequence of type [1 2 1 2 1 1 2 2 ... ]
    %
    % Output parameters:
    % segment1XStartPos = start positions of the segments
    % segment1XLengths = lengths of the segments
    %

    diffX = diff(binaryStateSeq);
    idxUp = find(diffX == 1);
    idxDown = find(diffX == -1);

    if isempty(idxUp) && isempty(idxDown) % in case of no state switches
        segment1XStartPos = [];
        segment1XLengths = [];
        return
    end

    if isempty(idxUp) 
        if ~isempty(idxDown) % in case of no state switches
            segment1XStartPos = idxDown+1;
            segment1XLengths = length(binaryStateSeq)-idxDown;
            return
        end
    elseif isempty(idxDown)
        if ~isempty(idxUp) % in case of no state switches
            segment1XStartPos = 1;
            segment1XLengths = idxUp;
            return
        end
    end

    if idxDown(1) < idxUp(1) % first level-crossing is a down-crossing
                             % (level crossing  = crossing of the value 3/2, 
                             % i.e. a transition from 1 to 2 = down-crossing
                             % or from 2 to 1 = up-crossing)

        segment1XStartPos = idxDown + 1;

        segment1XLengths = zeros(1,length(segment1XStartPos));
        if idxDown(end) < idxUp(end) % last level-crossing is an up-crossing       
            segment1XLengths = idxUp-idxDown;        
        else % last level-crossing is a down-crossing       
            segment1XLengths(1:end-1) = idxUp-idxDown(1:end-1);
            segment1XLengths(end) = length(binaryStateSeq) - idxDown(end);
        end      

    else  % idxDown(1) > idxUp(1), i.e., first level-crossing is an up-crossing

        segment1XStartPos = [1 , idxDown + 1];

        segment1XLengths = zeros(1,length(segment1XStartPos));
        if idxDown(end) < idxUp(end) % last level-crossing is an up-crossing       
           segment1XLengths(1) = idxUp(1);
           segment1XLengths(2:end) = idxUp(2:end) - idxDown;    
        else % last level-crossing is a down-crossing      
            segment1XLengths(1) = idxUp(1);
            segment1XLengths(2:end-1) = idxUp(2:end)-idxDown(1:end-1);
            segment1XLengths(end) = length(binaryStateSeq) - idxDown(end);          
        end

    end

    % segment1XStartPos
    % segment1XLengths
    % 

    % figure
    % plot(stateVec,'b-','linewidth',3)


end

