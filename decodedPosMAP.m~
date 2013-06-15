function decodedPos = decodedPosMAP(gt,posterior)
nBins = size(posterior,3);
for kBin = 1 : nBins 
	curPos = posterior(:,:,kBin);
	 [~, linIdx] = max(curPos(:));
        [i,j]=ind2sub(size(curPos), linIdx);
        maxRateXLoc(kBin) = gt.pfObject.xBin(i);
        maxRateYLoc(kBin) = gt.pfObject.yBin(j);
end
decodedPos = [maxRateXLoc', maxRateYLoc'];
end