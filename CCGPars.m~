function out = CCGPars(res, clu, sampleRate, varargin)

    if length(unique(clu)) > 1, defPairs = nchoosek(unique(clu), 2); end
    defOptions = struct('type', 'jitter', 'winSize', 20e-3, 'nResamples', 1e2);
    [pairs2Test, options, binSize, maxTimeLag, ccgSmthFactor, IF_RESAMPLE, IF_PLOT]  =  ...
        DefaultArgs(varargin, {defPairs, defOptions, 20e-3, 1000e-3, 0.03, 0, 1});
    
    halfBins = round(maxTimeLag / binSize); % number of bins on each side of zero lag
    binSize = round(binSize * sampleRate);
    jitter = @(jitterWinSiz, x) round((- jitterWinSiz / 2) + round(jitterWinSiz) .* rand(size(x))); 
    circShuffle = @(shuffle, x) circshift(x(:), [round(shuffle * rand), 0]);
    for lPair = 1 : size(pairs2Test, 1)
        res1 = res(clu == pairs2Test(lPair, 1));
        res2 = res(clu == pairs2Test(lPair, 2));
        lClu = [ones(length(res1), 1) * pairs2Test(lPair, 1); ones(length(res2), 1) * pairs2Test(lPair, 2)];
        [ccgOut, ccgTimeAx, pp] = ...
            myCCG([res1; res2], lClu, binSize, halfBins, sampleRate, pairs2Test(lPair, :), 'count', [], IF_PLOT);
        tt = ccgTimeAx(1):.1:ccgTimeAx(end);
        xBins = [-halfBins : halfBins] / halfBins  / 2;
        gw = exp(-power(xBins, 2) ./ (2 * ccgSmthFactor ^ 2));
        gw = gw ./ sum(gw);
        yy = conv(ccgOut(:, 1, 2), gw, 'same');
        ccgSmooth(:, lPair) = spline(ccgTimeAx',yy,tt);
        [T, offset, firstPeak] = FindCCGPars(ccgSmooth(:,lPair), tt);
        ccg{lPair}.ccg = ccgOut;
        ccg{lPair}.T  = T;
        ccg{lPair}.offset = offset;
        ccg{lPair}.firstPeak = firstPeak;
        ccg{lPair}.ccgSmooth = ccgSmooth(:, lPair);
        ccg{lPair}.smoothTAx = tt;
        %% Resample
        if IF_RESAMPLE
            for kResample  = 1 : options.nResamples % use 1s scram
                switch options.type
                  case 'jitter'
                    kRes1 = res1 + jitter(options.winSize * sampleRate, res1);
                    kRes2 = res2 + jitter(options.winSize * sampleRate, res2);
                  case 'circShift'
                    kRes1 = circShuffle(round(options.winSize * sampleRate), res1);
                    kRes2 = circShuffle(round(options.winSize * sampleRate), res2);
                end
                [shuffledCCGmat, ccgTimeAx, ~] = myCCG([kRes1; kRes2], lClu, binSize, halfBins, sampleRate, pairs2Test(lPair, :), 'count');
                shuffledCCG(:, kResample) = sq(shuffledCCGmat(:, 1, 2));
                yy = conv(shuffledCCG(:, kResample), gw, 'same'); %  gaussian with std ccgSmthFactor
                shuffledCCGSmooth = spline(ccgTimeAx',yy,tt);
                [sT(kResample, lPair), sOffset(kResample, lPair), sfirstPk(kResample, lPair)] = FindCCGPars(shuffledCCGSmooth, tt);
                kResample
                ss(:, kResample) = shuffledCCGmat(:, 1, 2);
                out.pVal(lPair) = sum(sum(shuffledCCG > repmat(ccgOut(:, 1, 2), 1, kResample))) ./ kResample ;  
            end
        end
    end
    out = ccg;
end