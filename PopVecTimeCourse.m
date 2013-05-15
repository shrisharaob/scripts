function [popVec, avgVector, dotProd] = PopVecTimeCourse(gt, varargin)
% population vector time course for the entire filebase, considers only the common clus
    [ThPh, commonClus, roi, arena, IF_COMPUTE, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins, nThCycles] = ...
        DefaultArgs(varargin, {[], [], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50], 4});
    
    switch gt.datasetType
      case 'kenji'
          
      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
    end

    if ~IF_COMPUTE
        fprintf('loading  ...');
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
        return;
    end
    if isempty(gt.clu), gt = gt.LoadCR; end
    if isempty(gt.pfObject), gt = gt.LoadPF; end
    if isempty(ThPh), load([gt.paths.data, gt.filebase, '.thpar.mat']); end
    if isempty(commonClus), load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']); end
    nClus = length(commonClus);
    res = gt.res(ismember(gt.clu, commonClus)); % load res only for the units in roi
    clu = gt.clu(ismember(gt.clu, commonClus));
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd', 1, 1);
    clu = clu(resIdx);
    %    [res, resIdx] = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ./ gt.trackingSampleRate) + 1, 'd');
    % clu = clu(resIdx);
    [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
    nDims = nRows * nClmns;
    refVector = zeros(nClus, nRows, nClmns); % each clm of this matrix is a pop vector at one of the spatial bins
    % indices of thetapeak in lfp sample rate thetaBoundaries
    thetaPeriods = gt.ThetaBoundaries(ThPh, commonClus, [], [], nThCycles); 
    if strcmp(gt.datasetType, 'kenji')
        markerNo = 1;
    elseif strcmp(gt.datasetType, 'MTA')
        markerNo = 7;
    end
    nCycles = size(thetaPeriods, 1);
    popVec = zeros(nClus, nRows, nClmns, nCycles);
    oldStr = [];
    fprintf('\n computing popvector...');
    [binnedPos, coverageMask] = Coverage(gt, roi, arena, markerNo, 0, 0);
    inout = InOut(~isnan(binnedPos(:,1)), [0, 1]); % select res from only when rat is in a good pf
    inout = inout{1};
    [res, resIdx] = SelectPeriods(res, round(inout * gt.lfpSampleRate ./ gt.trackingSampleRate) + 1, 'd');
    clu = clu(resIdx);
    for kPopVec = 1 : nCycles
        str = sprintf(['#', num2str(kPopVec) ' of ' num2str(nCycles)]);
        fprintf([repmat('\b', 1, length(oldStr)), str]);
        oldStr = str;
        [curRes, curResId] = SelectPeriods(res(ismember(clu, commonClus)), thetaPeriods(kPopVec, :),'d');
        curClu = clu(curResId);
        if ~isempty(curRes) & length(curRes) > 5
            tPos = SelectPeriods(binnedPos, round(thetaPeriods(kPopVec, :) .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1, 'c');
            tPos = round(nanmean(tPos, 1));
            spkCnt = zeros(nClus, 1);
            if ~(all(isnan(tPos(:))))
                for mClu = 1 : nClus
                    spkCnt(mClu) = sum(curClu == commonClus(mClu));
                end
                %                xyBin = Sub2Ind(size(coverageMask), tPos);
                popVec(:, tPos(1), tPos(2), kPopVec) = spkCnt; 
            end
        end
    end
    avgVector = mean(popVec, 4);
    fprintf('  done !!! \n');
    dotProd = avgVector(:)' * reshape(popVec, [], size(popVec,4));
    if IF_OVERWRITE
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'popVec', 'dotProd','-v7.3');
    end
end 
