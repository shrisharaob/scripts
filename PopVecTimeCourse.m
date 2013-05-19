function [popVec, avgVector, dotProd] = PopVecTimeCourse(gt, varargin)
% [popVec, avgVector, dotProd] = PopVecTimeCourse(gt, varargin)
% [ThPh, commonClus, roi, arena, IF_COMPUTE, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins, nThCycle]s
% [], [], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50], 1
% population vector time course for the entire filebase, considers only the common clus
  
  [ThPh, commonClus, roi, arena, IF_COMPUTE, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins, nThCycles] = ...
        DefaultArgs(varargin, {[], [], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50], 1});
    
    switch gt.datasetType
      case 'kenji'
          
      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
    end
    if ~IF_COMPUTE
        fprintf('loading  ...');
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
        avgVector = mean(popVec, 4);
        return;
    end
    if isempty(gt.clu), gt = gt.LoadCR; end
    if isempty(gt.pfObject), gt = gt.LoadPF; end
    if isempty(ThPh), 
        switch gt.datasetType
          case 'kenji'
            load([gt.paths.data, gt.filebase, '.thpar.mat']);
          case 'MTA'
            [ThPh, ~] = AnalyticTheta(gt, [], [], 0);
        end
    end
    if isempty(commonClus), load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']); end
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
        if ~isempty(curRes) & length(curRes) > 5
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
    sRateMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
    %    sRateMaps = reshape(permute(sRateMaps, [3 1 2]), [], nClus)
    fprintf('  done !!! \n');
    dp = @(a, b) a' * b ./ (norm(a) * vnorm(b)); % normalized dot product
    dotProd = dp(sRateMaps(:), reshape(popVec, [], size(popVec,4)));
    dotProd(isnan(dotProd)) = 0;
    if IF_COMPUTE
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'popVec', 'dotProd','-v7.3');
    end
end 
