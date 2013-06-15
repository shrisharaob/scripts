function [decodedPos, decodedBin] = decodedPosMAP(gt,posterior)
    
    nDims = ndims(posterior);
    if nDims == 3
        nBins = size(posterior,3);
        
        for kBin = 1 : nBins 
            curPos = posterior(:,:,kBin);
            [~, linIdx] = max(curPos(:));
            [i,j]=ind2sub(size(curPos), linIdx);
            decodedBin(kBin, :)  =  [i,j];
            maxRateXLoc(kBin) = gt.pfObject.xBin(i);
            maxRateYLoc(kBin) = gt.pfObject.yBin(j);
        end
        decodedPos = [maxRateXLoc', maxRateYLoc'];
    elseif nDims == 2
        nBins = size(posterior, 2);
        for kBin = 1 : nBins 
            curPos = posterior(:, kBin);
            [~, decodedBin(kBin)] = max(curPos(:));
            %            [~, linpos] = cart2pol(decodedBin);
            %            defPosRange = [min(linPos), max(linPos)];
            %decodedBin(kBin, :)  =  [i,j];
            %maxRateXLoc(kBin) = gt.pfObject.xBin(i);
            %maxRateYLoc(kBin) = gt.pfObject.yBin(j);
            decodedPos = [];
        end
    end
end