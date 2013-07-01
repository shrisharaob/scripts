function out = CCGPars(gt, varargin)
% out = CCG(gt, varargin)
% compute ccg for specified cell pairs in  the trial

    [cluId, roi, arena, binSize, maxTimeLag, ccgSmthFactor, jitterWinSiz, nResamples, alpha] = ...
        DefaultArgs(varargin, {[], 'CA3', 'bigSquare', 120e-3, 3000e-3, 0.03, 10, 1e2, 5e-2});

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
     if isempty(gt.pfObject), gt = gt.LoadPF; end
     for lPair = 1 : nPairs
         pRes = [];
         pClu = [];
         for mClu = 1 : 2
             pRes = [pRes; res(clu == cellPairs(lPair, mClu))];
             pClu = [pClu; cellPairs(lPair, mClu) * ones(length(res(clu == cellPairs(lPair, mClu))) , 1)];
         end
         [ccgOut, ccgTimeAx, pp] = myCCG(pRes, pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count', [], 1);
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
         subplot(2,2,4);
         gt.pfObject.PlotRateMaps(1, 0, 1, [],[],[], cellPairs(lPair, :));
         [T(lPair), offset(lPair), firstPeak(lPair)] = FindCCGPars(ccgSmooth(:,lPair), tt);
         % emperical sampling distribution
         pRes1 = pRes(pClu == cellPairs(lPair, 1));
         pRes2 = pRes(pClu == cellPairs(lPair, 2));
         jitter = @(jitterWinSiz, x) round((- jitterWinSiz / 2) + round(jitterWinSiz) .* rand(size(x))); 
         lPair
         drawnow;
keyboard;
%         for kResample  = 1 : nResamples % use 1s scram
%             jitteredRes1 = pRes1 + jitter(jitterWinSiz, pRes1);
%             jitteredRes2 = pRes2 + jitter(jitterWinSiz, pRes2);
%             [shuffledCCGmat, ~, ~] = myCCG([jitteredRes1; jitteredRes2], pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count');
%             shuffledCCG(:, kResample) = sq(shuffledCCGmat(:, 1, 2));
%             %  gaussian with std ccgSmthFactor
%             yy = conv(shuffledCCG(:, kResample), gw, 'same');
%             shuffledCCGSmooth = spline(ccgTimeAx',yy,tt);
%            [~, sOffset(kResample, lPair), ~] = FindCCGPars(shuffledCCGSmooth, tt);
%            kResample
%         end            
%         out.pVal(lPair) = sum(sum(shuffledCCG > repmat(ccgOut(:, 1, 2), 1, kResample))) ./ kResample ;
     clf;
     end
     out.IS_SIGNF = out.pVal < alpha ./ (2 * halfBins + 1);
     out.period = T;
     out.offset = offset;
     out.firstPeak = firstPeak;
     save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
 end
        
