function out = CCG(gt, varargin)
    % out = CCG(gt, varargin)
    % compute ccg for specified cell pairs in  the trial
     
    [cluId, roi, arena, binSize, maxTimeLag, ccgSmthFactor] = DefaultArgs(varargin, {[], 'CA3', 'bigSquare', 10e-3, 200e-3, 0.03});

    if isempty(cluId), load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']); end
    if isempty(gt.clu), gt = gt.LoadCR; end
    [res, clu] = gt.LoadStateRes('RUN', 1);
     if length(commonClus) > 1
         cellPairs = nchoosek(commonClus, 2);
     else 
         out = [];
         return;
     end
     nPairs = size(cellPairs, 1);
     halfBins = round(maxTimeLag / binSize); % number of bins on each side of zero lag
     binSize = round(binSize * gt.sampleRate);
     for lPair = 1 : nPairs
         pRes = [];
         pClu = [];
         for mClu = 1 : 2
             pRes = [pRes; res(clu == cellPairs(lPair, mClu))];
             pClu = [pClu; cellPairs(lPair, mClu) * ones(length(res(clu == cellPairs(lPair, mClu))) , 1)];
         end

         [ccgOut, ccgTimeAx, pp] = myCCG(pRes, pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count');
         ccg.TimeAx = ccgTimeAx;
         %  gaussian with std ccgSmthFactor
         tt = ccgTimeAx(1):.1:ccgTimeAx(end);
         xBins = [-halfBins : halfBins] / halfBins  / 2;
         gw = exp(-power(xBins, 2) ./ (2 * ccgSmthFactor ^ 2));
         gw = gw ./ sum(gw);
         yy = conv(ccgOut(:, 1, 2), gw, 'same');
         ccgSmooth(:, lPair) = spline(ccgTimeAx',yy,tt);
         [T(lPair), offset(lPair), firstPeak(lPair)] = FindCCGPars(ccgSmooth(:,lPair), tt);
     end
     out.period = T;
     out.offset = offset;
     out.firstPeak = firstPeak;
     close all;
end
        
