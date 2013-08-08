function out = CCGPars(gt, varargin)
% out = CCG(gt, varargin)
% compute ccg for specified cell pairs in  the trial

    [type, cluId, roi, arena, binSize, maxTimeLag, ccgSmthFactor, jitterWinSiz, nResamples, alpha] = ...
        DefaultArgs(varargin, {'display', [], 'CA3', 'bigSquare', 10e-3, 2000e-3, 0.03, 10, 0, 5e-2});

    if ~FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']), type = 'compute'; end
    switch type
      case 'compute'
        switch gt.datasetType
            case 'kenji'
              if isempty(cluId), load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
              else commonClus = cluId; end
              case 'MTA'
                commonClus = cluId;
        end
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
        if isempty(gt.pfObject), gt = gt.LoadPF; end
        for lPair = 1 : nPairs
            pRes = [];
            pClu = [];
            for mClu = 1 : 2
                pRes = [pRes; res(clu == cellPairs(lPair, mClu))];
                pClu = [pClu; cellPairs(lPair, mClu) * ones(length(res(clu == cellPairs(lPair, mClu))) , 1)];
            end
            [ccgOut, ccgTimeAx, pp] = myCCG(pRes, pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count', [], 0);
            out.ccg{lPair} = ccgOut;
            ccg.TimeAx = ccgTimeAx;
            %  gaussian with std ccgSmthFactor
            tt = ccgTimeAx(1) : .1 : ccgTimeAx(end);
            xBins = [-halfBins : halfBins] / halfBins  / 2;
            gw = exp(-power(xBins, 2) ./ (2 * ccgSmthFactor ^ 2));
            gw = gw ./ sum(gw);
            yy = conv(ccgOut(:, 1, 2), gw, 'same');
            ccgSmooth(:, lPair) = spline(ccgTimeAx',yy,tt);
            [T(lPair), offset(lPair), firstPeak(lPair)] = FindCCGPars(ccgSmooth(:,lPair), tt);
            % emperical sampling distribution
            pRes1 = pRes(pClu == cellPairs(lPair, 1));
            pRes2 = pRes(pClu == cellPairs(lPair, 2));
            jitter = @(jitterWinSiz, x) round((- jitterWinSiz / 2) + round(jitterWinSiz) .* rand(size(x))); 
            if nResamples
                for kResample  = 1 : nResamples % use 1s scram
                    jitteredRes1 = pRes1 + jitter(jitterWinSiz, pRes1);
                    jitteredRes2 = pRes2 + jitter(jitterWinSiz, pRes2);
                    [shuffledCCGmat, ~, ~] = myCCG([jitteredRes1; jitteredRes2], pClu, binSize, halfBins, gt.sampleRate, cellPairs(lPair, :), 'count');
                    shuffledCCG(:, kResample) = sq(shuffledCCGmat(:, 1, 2));
                    % gaussian with std ccgSmthFactor
                    yy = conv(shuffledCCG(:, kResample), gw, 'same');
                    shuffledCCGSmooth = spline(ccgTimeAx',yy,tt);
                    [~, sOffset(kResample, lPair), ~] = FindCCGPars(shuffledCCGSmooth, tt);
                end            
                out.pVal(lPair) = sum(sum(shuffledCCG > repmat(ccgOut(:, 1, 2), 1, kResample))) ./ kResample ;
            end
        end
        % out.IS_SIGNF = out.pVal < alpha ./ (2 * halfBins + 1);
        out.ccgTimeAx = ccgTimeAx;
        out.cellPairs = cellPairs;
        out.smthCCG = ccgSmooth;
        out.smthTAx  = tt;
        out.period = T;
        out.offset = offset;
        out.firstPeak = firstPeak;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'out');
      
      case 'display'

   load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'out');
keyboard;     

        figure;
        cellPairs = nchoosek(gt.pfObject.acceptedUnits(gt.pfObject.sparsity < .1), 2);
        idx = find(ismember(out.cellPairs, cellPairs, 'rows'));
        for kCellPair = 1 : length(idx);
            subplot(2, 2, 1);
            bar(out.ccgTimeAx, out.ccg{idx(kCellPair)}(:, 1, 2), 'FaceColor', 'k'); axis tight;
            hold on;
            plot(out.smthTAx, out.smthCCG(:, idx(kCellPair)), 'r'); axis tight;
            subplot(2, 2, 2);
            bar(out.ccgTimeAx, out.ccg{idx(kCellPair)}(:, 1, 1), 'FaceColor', 'k'); axis tight;
            subplot(2, 2, 3);
            bar(out.ccgTimeAx, out.ccg{idx(kCellPair)}(:, 2, 2), 'FaceColor', 'k'); axis tight;
            subplot(2,2,4);
            gt.pfObject.PlotRateMaps(1, 0, 1, [],[],[], out.cellPairs(idx(kCellPair), :));
            waitforbuttonpress;
            clf;
        end
    end
end

