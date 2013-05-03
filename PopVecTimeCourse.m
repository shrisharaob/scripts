function [popVec, refVector, dotProd] = PopVecTimeCourse(filebase, varargin)
    % population vector time course for the entire filebase, considers only the common clus
    [commonClus, roi, arena, IF_COMPUTE, trialName, binSize, tolerence, IF_OVERWRITE,  spatialBins] = ...
        DefaultArgs(varargin, {[], {'CA3'},  {'bigSquare'}, 0, [], 10, 1e-1, 1, [50, 50]});
  
   if ~IF_COMPUTE
       gt = GenericTrial(filebase, trialName);
       fprintf('loading  ...');
       load([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), mfilename, '.mat']);
       return;
   end
   
   filetag = GenFiletag(roi, arena);
   if isempty(commonClus)
        load(['~/data/analysis/kenji/', filebase, '/', filebase, filetag 'commonClus.mat']);
    end
    if  isempty(commonClus), return; end
    nClus = length(commonClus);
    sK.roi = roi;
    sK.arena = arena;
    filebases = SearchKenji(sK);
    workingFb = filebase;
    trialNames = filebases(strcmp(filebases(:,1), workingFb), 2);
    if isempty(trialNames), return; end;
    for kTr = 1 : length(trialNames)
        gt = GenericTrial(workingFb, trialNames{kTr});
        if isempty(gt.res)
            if kTr == 1
                gt =gt.LoadCR;
                allRes = gt.res;
                allClu = gt.clu;
            else
                gt.clu = allClu;
                gt.res = allRes;
            end
        end
        if isempty(gt.pfObject)
            gt = gt.LoadPF;
        end
        res = gt.res(ismember(gt.clu, commonClus)); % load res only for the units in roi
        clu = gt.clu(ismember(gt.clu, commonClus));
        res = round(res .* gt.lfpSampleRate ./ gt.sampleRate) + 1; % convert res to lfp sample rate
        [res, resIdx] = SelectPeriods(res, gt.trialPeriods, 'd');
        clu = clu(resIdx);
        [res, resIdx] = SelectPeriods(res, round(gt.goodPosPeriods .* gt.lfpSampleRate ./ gt.trackingSampleRate) + 1 + gt.trialPeriods(1,1),'d', 1, 1);
        clu = clu(resIdx);
        if kTr == 1
            fprintf('\n loading thpar...');
            load ([gt.paths.data, gt.filebase '.thpar.mat'] );
            fprintf('done !! \n');
        end
        [trialThPh, thIdx] = SelectPeriods(ThPh, gt.trialPeriods, 'c', 1);
        nBins = 360 / binSize;
        binEdges = linspace(-pi, pi, nBins);
        [count, binIdx] = histc(trialThPh(res), binEdges);
        xx = [binEdges, binEdges(2 : end) + 2*pi];
        [~, minCIdx] = min(count);
        minPh = xx(minCIdx);
        thetaBoundaries = find(IsEqual(trialThPh, minPh, tolerence, 0)); 
        % in lfp sample rate thetaBoundaries
        % indices of thetapeak
        [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
        nDims = nRows * nClmns;
        refVector = zeros(nClus, nDims); % each clm of this matrix is a pop vector at one of the spatial bins
        for mClu = 1 : length(commonClus)
            refVector(mClu, :) = Mat2Vec(gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus(mClu))), 0);
        end
        % refVector = refVector(:);
        % refVector = refVector ./ norm(refVector);
        fprintf('\n computing popvector...');
        thetaPeriods = [thetaBoundaries(1 : end - 1)- 1, thetaBoundaries(2 : end)]; 
        if strcmp(gt.datasetType, 'kenji')
            markerNo = 1;
        elseif strcmp(gt.datasetType, 'MTA')
            markerNo = 7;
        end
        nCycles = size(thetaPeriods, 1);
        popVec = zeros(nClus, nDims, nCycles);
        oldStr = [];
        for kPopVec = 1 : nCycles
            str = sprintf(['#', num2str(kPopVec) ' of ' num2str(nCycles)]);
            fprintf([repmat('\b', 1, length(oldStr)), str]);
            oldStr = str;
            [curRes, curResId] = SelectPeriods(res(ismember(clu, commonClus)), thetaPeriods(kPopVec, :),'d',1,1);
            curClu = clu(curResId);
            if ~isempty(curRes)
                pos = SelectPeriods(gt.position(:, markerNo, :), round(thetaPeriods(kPopVec, :) .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1);
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
    %% figures
    % figure;
    % bar(xx * 180 /pi, [count', fliplr(count(2:end))'],'FaceColor', 'k');
    % hold on;
    % line([minPh, minPh], ylim, 'Color', 'g');
    % line([minPh, minPh]+2*pi, ylim, 'Color', 'g');
    % figure;
    % plot(linspace(0,sum(diff(gt.goodPosPeriods, 1, 2)) / gt.trackingSampleRate, length(dotProd)), dotProd);   
    


