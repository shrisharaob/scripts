function [popVec, refVector] = RemappingTimeCourse(gt, varargin)
    % RemappingTimeCourse(gt, varargin)

    [binSize, roi, tolerence, IF_OVERWRITE] = DefaultArgs(varargin, { 10,  {'CA3'}, 1e-1,1});
    %binSize in degrees
    % .1 ~= 5 degrees
    sK.roi = roi;
    if ~SearchKenji(sK, gt.filebase, gt.trialName), return; end;
    if isempty(gt.res)
        gt =gt.LoadCR;
    end
    roiClus = cell2mat(gt.GetRegionClu(roi)); % shd also add
    roiClus= roiClus(ismember(roiClus, gt.pyrCluIdx(gt.pfObject.acceptedUnits))); % clus in roi & with PFs
    if isempty(roiClus), return; end; % no clusters in the roi
    if isempty(gt.pfObject)
        gt = gt.LoadPF;
    end
    roiPFClu = unique(gt.clu(CompareVectors(gt.clu, roiClus)));
    validClu = roiPFClu(~cellfun(@isempty,gt.pfObject.rateMap(roiPFClu))); % remove clus for with no rate maps
    validClu = validClu(ismember(validClu, gt.pyrCluIdx(gt.pfObject.idealPFUnits)));
    res = gt.res(ismember(gt.clu, validClu)); % load res only for the units in roi
    clu = gt.clu(ismember(gt.clu, validClu));
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd');
    clu = clu(resIdx);
    [res, resIdx] = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ./ gt.trackingSampleRate) + 1 + gt.trialPeriods(1,1),'d', 1, 1);
    clu = clu(resIdx);
    fprintf('\n loading thpar...');
    load ([gt.paths.data, gt.filebase '.thpar.mat'] );
    fprintf('done !! \n'); 
    [trialThPh, thIdx] = SelectPeriods(ThPh, gt.trialPeriods, 'c', 1);
    nBins = 360 / binSize;
    binEdges = linspace(-pi, pi, nBins);
    [count, binIdx] = histc(trialThPh(res), binEdges);
    xx = [binEdges, binEdges(2 : end) + 2*pi];
    [~, minCIdx] = min(count);
    minPh = xx(minCIdx);
    thetaBoundaries = find(IsEqual(trialThPh, minPh, tolerence, 0)); % in lfp sample rate 
    [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
    nDims = nRows * nClmns;
    popVec = zeros(nDims, length(thetaBoundaries) - 1);
    refVector = zeros(nDims, 1);
    % refVector = sum(reshape(Nan2Zero(cell2mat(gt.pfObject.rateMap(roiClus))), nRows, nClmns, length(roiClus)), 3);
    refVector = nansum(gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, validClu)), 3);
    refVector = refVector(:);
    refVector = refVector ./ norm(refVector);
    fprintf('\n computing popvector...');
    thetaPeriods = [thetaBoundaries(1 : end - 1)- 1, thetaBoundaries(2 : end)]; 
    if strcmp(gt.datasetType, 'kenji')
        markerNo = 1;
    elseif strcmp(gt.datasetType, 'MTA')
        markerNo = 7;
    end
    for kPopVec = 1 : length(thetaBoundaries) - 1
        % population vector for each theta cycle, nDims -by- nClycles 
        for kClu = 1 : length(validClu)
            curRes = SelectPeriods(res(clu == validClu(kClu)), thetaPeriods(kPopVec, :),'d',1,1);
            if ~isempty(curRes)
                pos = SelectPeriods(gt.position(:, markerNo, :), round(thetaPeriods(kPopVec, :) .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1);
                if ~(all(isnan(pos(:))))
                    fprintf('*');
                    popVec(:, kPopVec) = popVec(:, kPopVec) + SpkCntAtPos(gt, curRes, pos); 
                end
            end
        end
    end
    fprintf('done !! \n');
    dotProd = popVec' * refVector;
    save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'refVector', 'popVec', 'dotProd');
    %% figures
    figure;
    %    bar(xx * 180 /pi, [count', fliplr(count(2:end))'],'FaceColor', 'k');
    %    hold on;
    % line([minPh, minPh], ylim, 'Color', 'g');
    % line([minPh, minPh]+2*pi, ylim, 'Color', 'g');
    figure;
    plot(linspace(0,sum(diff(gt.goodPosPeriods, 1, 2)) / gt.trackingSampleRate, length(dotProd)), dotProd);   
end    


