function [popVec, refVector] = RemappingTimeCourse(gt, varargin)
    % RemappingTimeCourse(gt, varargin)

    [binSize, roi, tolerence, IF_OVERWRITE] = DefaultArgs(varargin, { 10,  {'CA3'}, 1e-1,1});
    %binSize in degrees
    % .1 ~= 5 degrees
    fprintf('\n loading thpar...');
    load ([gt.paths.data, gt.filebase '.thpar.mat'] );
    fprintf('done !! \n');
    if isempty(gt.res)
        gt =gt.LoadCR;
    end
    if isempty(gt.pfObject.rateMap)
        gt = gt.LoadPF;
    end
    roiClus = cell2mat(gt.GetRegionClu(roi)); % shd also add
    roiClus= roiClus(ismember(roiClus, gt.pyrCluIdx(gt.pfObject.acceptedUnits)));
    res = gt.res(CompareVectors(gt.clu, roiClus)); % load res only for the units in roi
    res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
    [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd');
    clu = gt.clu(resIdx);
    [res, resIdx] = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ./ gt.trackingSampleRate) + 1 + gt.trialPeriods(1,1),'d', 1, 1);
    clu = clu(resIdx);
    [trialThPh, thIdx] = SelectPeriods(ThPh, gt.trialPeriods, 'c', 1);
    nBins = 360 / binSize;
    bins = linspace(-pi, pi, nBins);
    [count, binIdx] = hist(trialThPh, bins);
    xx = [bins, bins(2:end)+2*pi, bins(end) + 2*pi + 2 * pi / nBins];
    [~, minCIdx] = min(count);
    minPh = xx(minCIdx);
    thetaBoundaries = find(IsEqual(trialThPh, minPh, tolerence, 0)); 
    % in lfp sample rate thetaBoundaries = find(IsEqual(trialThPh,pi,
    % tolerence, 0));% indeices of thetapeak
    [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
    nDims = nRows * nClmns;
    popVec = zeros(nDims, length(thetaBoundaries) - 1);
    refVector = zeros(nDims, 1);
    validClu = roiClus(~cellfun(@isempty,gt.pfObject.rateMap(roiClus))); % remove clus for with no rate maps 
    refVector = sum(reshape(cell2mat(gt.pfObject.rateMap(roiClus)), nRows, nClmns, length(roiClus)), 3);
    refVector = refVector(:);
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
            curRes = SelectPeriods(res(clu == validClu(kClu)), thetaPeriods(kPopVec, :),'d');
            if ~isempty(curRes)
                kClu
                tic
                position = SelectPeriods(gt.position(:, markerNo, :), round(thetaPeriods(kPopVec, :) .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1);
                toc
                if ~(all(isnan(pos(:))))
                    fprintf('*');
                    popVec(:, kPopVec) = popVec(:, kPopVec) + SpkCntAtPos(gt, curRes, pos); 
                end
            end
        end
    end
    fprintf('done !! \n');
    dotProd = popVec' * refVector;
    save([gt.paths.analysis, gt.filebase, gt.trialName, mfilename, '.mat'], 'refVector', 'popVec', 'dotProd');
    %% figures
    figure;
    bar(xx, [count, count],'FaceColor', 'k');
    hold on;
    line([minPh, minPh], ylim, 'Color', 'g');
    line([minPh, minPh]+2*pi, ylim, 'Color', 'g');
    figure;
    plot(linspace(0,sum(diff(gt.goodPosPeriods, 1, 2)) / gt.trackingSampleRate, length(dotProd)), dotProd);   
    


