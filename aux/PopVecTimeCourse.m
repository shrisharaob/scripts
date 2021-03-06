 function [popVec, avgVector, dotProd] = PopVecTimeCourse(gt, varargin)
% [popVec, avgVector, dotProd] = PopVecTimeCourse(gt, varargin)
% [IF_COMPUTE, IF_CHUNKS, nChunks, ThPh, commonClus, roi, arena, IF_REPORTFIG, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins, nThCycles, rmSmoothFactor]

% [], [], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50], 1
% popVec nDims-by-nThetaCycles
% population vector time course for the entire filebase, considers only the common clus
    
    [IF_COMPUTE, IF_CHUNKS, nChunks, ThPh, commonClus, roi, arena, IF_REPORTFIG, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins, nThCycles, rmSmoothFactor] = ...
        DefaultArgs(varargin, { 0, 0, 3, [], [], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50], 1, 0.02});

avgVector = [];    

    switch gt.datasetType
      case 'kenji'
        
      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
    end
    if isempty(gt.pfObject), gt.LoadPF; end
    dp = @(a, b) a' * b ./ (norm(a) * vnorm(b)); % normalized dot product
    if isempty(commonClus), load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']); end
    if ~IF_COMPUTE & ~IF_CHUNKS
        fprintf('loading  ...');
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
        avgVector = mean(popVec, 4);
        sRateMaps = sq(gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus)));
        dotProd = atan(dp(sRateMaps(:), popVec));
                     %popVec = pv;
        return;
    end
    %% COMPUTE AVG RATE MAPS FOR CHUNKS
    if IF_CHUNKS
        if ~IF_COMPUTE
            load([gt.paths.analysis, gt.filebase, '.' gt.trialName, GenFiletag(roi, arena), 'CHUNKS.', num2str(nChunks), '.',mfilename, '.mat']);
            if IF_REPORTFIG
                if isempty(gt.pfObject), gt = gt.LoadPF;end
                for kk = 1 : length(out.clu)
                    p1 = reshape(full(out.popVec), [50, 50, length(out.clu), nChunks]);
                    subplot(1, nChunks + 1, nChunks + 1);
                    imagesc(gt.pfObject.smoothRateMap(:,:,ismember(gt.pfObject.acceptedUnits, out.clu(kk))));
                    axis square;
                    clims = caxis;
                    for i = 1: nChunks
                        subplot(1, nChunks + 1, i);
                        imagesc(sq(p1(:, :, (kk), i)))
                        axis square;
                        caxis(clims);
                    end
                    reportfig(gcf, ['smthmap.CHUNKS.', num2str(nChunks)], 0, [gt.filebase, '.' gt.trialName '      clu id :: ' num2str(out.clu(kk))]);
                    clf
                end
            end
            popVec = out.popVec;
            dotProd = out.dotProd;
            return;
        end

        if isempty(gt.pfObject), gt = gt.LoadPF; end
        if isempty(gt.res), gt = gt.LoadCR; end
        switch gt.datasetType
          case 'kenji'
            chunkBoundaries = gt.trialPeriods(1) : floor(diff(gt.trialPeriods) / nChunks) : gt.trialPeriods(2); 
        end
        chunkPeriods = ConvertFs([chunkBoundaries(1: end - 1)', chunkBoundaries(2 : end)'], gt.lfpSampleRate, gt.trackingSampleRate);
        [res, clu, pos] = gt.LoadStateRes('RUN', 1, gt.trackingSampleRate, chunkPeriods);
        nClus = length(commonClus);
        sRateMaps = sq(gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus)));
        for kChunk = 1 : nChunks
            str = sprintf('chunck %d of %d \n', kChunk, nChunks);
            fprintf(str);
            kChunkRes = res{kChunk};
            for kClu = 1 : nClus
                rm{kChunk, kClu} = GenericPF.ComputeRateMap(gt, kChunkRes(clu{kChunk} == commonClus(kClu)), pos{kChunk}, [], 0.03);
                srm(:, kClu) = Mat2Vec(SmoothSurface(rm{kChunk, kClu}, rmSmoothFactor));
            end
            popVec(:, kChunk) = srm(:);
        end
        out.popVec = sparse(popVec);
        out.clu = commonClus;
        out.rateMap = rm;
        out.dotProd = atan(dp(sRateMaps(:), out.popVec)) ; %  ./ nClus);  % normalize dp by the number of cells 
        save([gt.paths.analysis, gt.filebase, '.' gt.trialName, GenFiletag(roi, arena), 'CHUNKS.', num2str(nChunks), '.', mfilename, '.mat'], 'out', '-v7.3');
        %save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'CHUNKS.', num2str(nChunks), '.', mfilename, '.mat'], 'popVec', 'dotProd','-v7.3');
        popVec = out.popVec;
        avgVector = [];
        dotProd = out.dotProd;
        return;
    end
    %% COMPUTE PV FOR THETA CYCLES
    if isempty(gt.clu), gt = gt.LoadCR; end
    if isempty(gt.pfObject), gt = gt.LoadPF; end
    if isempty(ThPh)
        switch gt.datasetType
          case 'kenji'
            load([gt.paths.data, gt.filebase, '.thpar.mat']);
          case 'MTA'
            [ThPh, ~] = AnalyticTheta(gt, [], [], 0);
        end
    end

    nClus = length(commonClus);
    res = gt.res(ismember(gt.clu, commonClus)); % load res only for the units in roi
    clu = gt.clu(ismember(gt.clu, commonClus));
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd', 1, 1);
    clu = clu(resIdx);
    [res, resIdx] = SelectPeriods(res, ConvertFs(gt.goodPosPeriods, gt.trackingSampleRate, gt.lfpSampleRate), 'd');
    clu = clu(resIdx);
    [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
    nDims = nRows * nClmns;
    refVector = zeros(nClus, nRows, nClmns); % each clm of this matrix is a pop vector at one of the spatial bins
                                             % indices of thetapeak in lfp sample rate thetaBoundaries
    if strcmp(gt.datasetType, 'kenji')
        markerNo = 1;
    elseif strcmp(gt.datasetType, 'MTA')
        markerNo = 7;
    end
    oldStr = [];
    [binnedPos, coverageMask] = Coverage(gt, roi, arena, markerNo, 0, 0);
    inout = InOut(~isnan(binnedPos(:,1)), [0, 1]); % select res from only when rat is in a good pf
    inout = inout{1};
    thetaPeriods = gt.ThetaBoundaries(ThPh, commonClus, [], [], nThCycles, ConvertFs(inout, gt.trackingSampleRate, gt.lfpSampleRate));
    nCycles = size(thetaPeriods, 1);
    popVec = zeros(nClus, nRows, nClmns, nCycles);
    [res, resIdx] = SelectPeriods(res, ConvertFs(inout, gt.trackingSampleRate, gt.lfpSampleRate), 'd');
    clu = clu(resIdx);
    fprintf('\n computing popvector...');
    for kPopVec = 1 : nCycles
        str = sprintf(['#', num2str(kPopVec) ' of ' num2str(nCycles)]);
        fprintf([repmat('\b', 1, length(oldStr)), str]);
        oldStr = str;
        [curRes, curResId] = SelectPeriods(res(ismember(clu, commonClus)), thetaPeriods(kPopVec, :), 'd');
        curClu = clu(curResId);
        if ~isempty(curRes)% & length(curRes) > 5 
            tPos = SelectPeriods(binnedPos, ConvertFs(thetaPeriods(kPopVec, :), gt.lfpSampleRate, gt.trackingSampleRate), 'c');
            tPos = round(nanmean(tPos, 1));
            spkCnt = zeros(nClus, 1);
            if ~(all(isnan(tPos(:))))
                for mClu = 1 : nClus
                    spkCnt(mClu) = sum(curClu == commonClus(mClu));
                end
                popVec(:, tPos(1), tPos(2), kPopVec) = spkCnt; 
            end
        end
    end
    avgVector = mean(popVec, 4);
    popVec = sparse(reshape(popVec, nClus * nDims, nCycles));
    sRateMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
    fprintf('  done !!! \n');
    dotProd = atan(dp(sRateMaps(:), popVec) / nClus); % Fisher z - var equalizing for the corr estimator
    dotProd(isnan(dotProd)) = 0;
    if IF_COMPUTE
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'popVec', 'dotProd','-v7.3');
    end
end 
