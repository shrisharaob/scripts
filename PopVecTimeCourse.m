function [popVec, refVector, dotProd] = PopVecTimeCourse(gt, ThPh, varargin)
% population vector time course for the entire filebase, considers only the common clus
    [commonClus, roi, arena, IF_COMPUTE, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins, nThCycles] = ...
        DefaultArgs(varargin, {[], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50], 3});
    
    if ~IF_COMPUTE
        fprintf('loading  ...');
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
        return;
    end
    nClus = length(commonClus);
    res = gt.res(ismember(gt.clu, commonClus)); % load res only for the units in roi
    clu = gt.clu(ismember(gt.clu, commonClus));
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd');
    clu = clu(resIdx);
    [res, resIdx] = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ./ gt.trackingSampleRate) + 1 + gt.trialPeriods(1,1),'d', 1, 1);
    clu = clu(resIdx);
    % in lfp sample rate thetaBoundaries
    % indices of thetapeak
    [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
    nDims = nRows * nClmns;
    refVector = zeros(nClus, nDims); % each clm of this matrix is a pop vector at one of the spatial bins
    for mClu = 1 : length(commonClus)
        refVector(mClu, :) = Mat2Vec(gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus(mClu))), 0);
    end
    thetaPeriods = gt.ThetaBoundaries(ThPh, commonClus, [], [], nThCycles); 
    if strcmp(gt.datasetType, 'kenji')
        markerNo = 1;
    elseif strcmp(gt.datasetType, 'MTA')
        markerNo = 7;
    end
    nCycles = size(thetaPeriods, 1);
    popVec = zeros(nClus, nDims, nCycles);
    oldStr = [];
    fprintf('\n computing popvector...');
    [pos, ~] = Coverage(gt, roi, arena, markerNo, 0, 0);
    keyboard;
    for kPopVec = 1 : nCycles
        str = sprintf(['#', num2str(kPopVec) ' of ' num2str(nCycles)]);
        fprintf([repmat('\b', 1, length(oldStr)), str]);
        oldStr = str;
        [curRes, curResId] = SelectPeriods(res(ismember(clu, commonClus)), thetaPeriods(kPopVec, :),'d',1,1);
        curClu = clu(curResId);
        if ~isempty(curRes)
            pos = SelectPeriods(pos, round(thetaPeriods(kPopVec, :) .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1);
            pos = nanmean(pos, 1);
            spkCnt = zeros(nClus, 1);
            if ~(all(isnan(pos(:))))
                for mClu = 1 : nClus
                    spkCnt(mClu) = sum(curRes == commonClus(mClu));
                end
                %xyBin = sub2ind(spatialBins, find(IsEqual(pos(1), gt.pfObject.xBin, 2), 1), find(IsEqual(pos(2), gt.pfObject.yBin, 2), 1));
                xbin = find(histc(pos(1), gt.pfObject.xBin), 1);
                if isempty(xbin)
                    xbin = 1;
                end
                ybin = find(histc(pos(2), gt.pfObject.yBin), 1);
                if isempty(ybin)
                    ybin = 1;
                end
                xyBin = sub2ind(spatialBins, xbin, ybin);
                popVec(:, xyBin, kPopVec) = spkCnt; 
            end
        end
    end
    fprintf('  done !!! \n');       
    dotProd = refVector(:)' * reshape(popVec, [], size(popVec,3));
    save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat'], 'refVector', 'popVec', 'dotProd','-v7.3');
end 
