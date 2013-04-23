function [popVec, refVector] = PopVecTimeCourse(filebase, varargin)
    % population vector time course for the entire filebase, considers only the common clus

    [commonClus, binSize, roi, tolerence, IF_OVERWRITE, arena] = DefaultArgs(varargin, {[],  10,  {'CA3'}, 1e-1,1, {'bigSquare'}});
    %binSize in degrees
    % .1 ~= 5 degrees
    %    if FileExists(
    filetag = [];
    for mRoi = 1 : length(roi)
        if mRoi == 1
            filetag = ['.', char(roi{1})];
        end
        if mRoi > 1
            filetag = [filetag, '.', char(roi{mRoi})];
        end
    end
    for lArena = 1 : length(arena)
        filetag = [filetag, '.', char(arena{lArena})];
    end
    if isempty(commonClus)
        load(['~/data/analysis/kenji/', filebase, '/', filebase, filetag '.commonClus.mat']);
    end
    if  isempty(commonClus), return; end
    sK.roi = roi;
    sK.arena = arena;
    %if ~SearchKenji(sK, gt.filebase, gt.trialName), return; end;
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
        %roiClus = cell2mat(gt.GetRegionClu(roi)); % shd also add
        %roiClus= roiClus(ismember(roiClus, gt.pyrCluIdx(gt.pfObject.acceptedUnits))); % clus in roi & with PFs
        %if isempty(roiClus), return; end; % no clusters in the roi
        if isempty(gt.pfObject)
            gt = gt.LoadPF;
        end
        % roiPFClu = unique(gt.clu(CompareVectors(gt.clu, roiClus)));
        % validClu = roiPFClu(~cellfun(@isempty,gt.pfObject.rateMap(roiPFClu))); % remove clus for with no rate maps
        % validClu = validClu(ismember(validClu, gt.pyrCluIdx(gt.pfObject.idealPFUnits)));
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
        % in lfp sample rate thetaBoundaries = find(IsEqual(trialThPh,pi,
        % tolerence, 0));% indices of thetapeak
        [nRows, nClmns] = size(gt.pfObject.rateMap{find(~cellfun(@isempty, gt.pfObject.rateMap), 1)});
        nDims = nRows * nClmns;
        
        refVector = zeros(nDims, 1);
        refVector = nansum(reshape(Nan2Zero(cell2mat(gt.pfObject.rateMap(commonClus))), nRows, nClmns, length(commonClus)), 3);
        % refVector = nansum(gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus)), 3);
        refVector = refVector(:);
        % refVector = refVector ./ norm(refVector);
        fprintf('\n computing popvector...');
        thetaPeriods = [thetaBoundaries(1 : end - 1)- 1, thetaBoundaries(2 : end)]; 
        if strcmp(gt.datasetType, 'kenji')
            markerNo = 1;
        elseif strcmp(gt.datasetType, 'MTA')
            markerNo = 7;
        end
        nCycles = size(thetaPeriods, 1);
        popVec = zeros(nDims, nCycles);
        oldStr = [];
        for kPopVec = 1 : nCycles
            str = sprintf(['#', num2str(kPopVec) ' of ' num2str(nCycles)]);
            fprintf([repmat('\b', 1, length(oldStr)), str]);
            oldStr = str;
            % population vector for each theta cycle, nDims -by- nClycles
            for kClu = 1 : length(commonClus)
                curRes = SelectPeriods(res(clu == commonClus(kClu)), thetaPeriods(kPopVec, :),'d',1,1);
                if ~isempty(curRes)
                    pos = SelectPeriods(gt.position(:, markerNo, :), round(thetaPeriods(kPopVec, :) .* gt.trackingSampleRate ./ gt.lfpSampleRate) + 1);
                    if ~(all(isnan(pos(:))))
                        
                        popVec(:, kPopVec) = popVec(:, kPopVec) + SpkCntAtPos(gt, curRes, pos);
                    end
                end
            end
        end
        fprintf('  done !!! \n');
%         dotProd = popVec' * refVector;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, filetag, mfilename, '.mat'], 'refVector', 'popVec');
    end 
    %% figures
    % figure;
    %    bar(xx * 180 /pi, [count', fliplr(count(2:end))'],'FaceColor', 'k');
    %    hold on;
    % line([minPh, minPh], ylim, 'Color', 'g');
    % line([minPh, minPh]+2*pi, ylim, 'Color', 'g');
    %    figure;
    % plot(linspace(0,sum(diff(gt.goodPosPeriods, 1, 2)) / gt.trackingSampleRate, length(dotProd)), dotProd);   
    


