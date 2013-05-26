function out = CCGPars(gt, varargin)
% out = CCG(gt, varargin)
% compute ccg for specified cell pairs in  the trial

    [cluId, roi, arena, binSize, maxTimeLag, ccgSmthFactor, jitterWinSiz, nResamples] = ...
        DefaultArgs(varargin, {[], 'CA3', 'bigSquare', 10e-3, 200e-3, 0.03, 10, 1e3});

    if isempty(cluId), load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
    else commonClus = cluId; end
    if isempty(gt.clu), gt = gt.LoadCR; end
%     [res, resInd] = SelectPeriods(gt.res, ConvertFs(gt.goodPosPeriods, gt.trackingSampleRate, ), gt.sampleRate, 'd');
%     clu = gt.clu(resInd);
%     nClu = length(unique(clu));
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

         [ccgOut, ccgTimeAx, pp] = myCCG(pRes, pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count', [], 1);
         %         ccg.Out(:,:,:,lPair) = ccgOut;
         ccg.TimeAx = ccgTimeAx;
         %  gaussian with std ccgSmthFactor
         tt = ccgTimeAx(1):.1:ccgTimeAx(end);
         xBins = [-halfBins : halfBins] / halfBins  / 2;
         gw = exp(-power(xBins, 2) ./ (2 * ccgSmthFactor ^ 2));
         gw = gw ./ sum(gw);
         yy = conv(ccgOut(:, 1, 2), gw, 'same');
         ccgSmooth(:, lPair) = spline(ccgTimeAx',yy,tt);
         subplot(2, 2, 1);
         hold on;
         plot(tt, ccgSmooth(:, lPair), 'r');
         [T(lPair), offset(lPair), firstPeak(lPair)] = FindCCGPars(ccgSmooth(:,lPair), tt);

         % emperical sampling distribution
         pRes1 = pRes(pClu == cellPairs(lPair, 1));
         pRes2 = pRes(pClu == cellPairs(lPair, 2));
         jitter = @(jitterWinSiz, x) round((- jitterWinSiz / 2) + round(jitterWinSiz) .* rand(size(x))); 
        for kResample  = 1 : nResamples
            jitteredRes1 = pRes1 + jitter(jitterWinSiz, pRes1);
            jitteredRes2 = pRes2 + jitter(jitterWinSiz, pRes2);
            [shuffledCCGOut, ccgTimeAx, ~] = myCCG([jitteredRes1; jitteredRes2], pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count');
            %  gaussian with std ccgSmthFactor
            yy = conv(shuffledCCGOut(:, 1, 2), gw, 'same');
            shuffledCCGSmooth(:, lPair) = spline(ccgTimeAx',yy,tt);
%             subplot(2, 2, 1);
%             hold on;
%             plot(tt, ccgSmooth(:, lPair), 'r');
%             keyboard;
           [~, sOffset(kResample, lPair), ~] = FindCCGPars(shuffledCCGSmooth(:,lPair), tt);
           kResample
        end            
        keyboard;
     end

     out.period = T;
     out.offset = offset;
     out.firstPeak = firstPeak;
     close all;
end
        
