function out = GetGoodPFUnits(filebase, varargin)

    [roi, arena, IF_SmthRM, IF_REPORTFIG] = DefaultArgs(varargin, { {'CA1'}, {'bigSquare'}, 1, 0 });
    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    nTrials = size(matches, 1);
    refTr = 1;
    nPairs = 0;
    out.commonPairs = [];
    out.goodUnits = [];
    if nTrials == 0, return; end
    for kTr = 1 : nTrials
        gt = GenericTrial(filebase, matches{kTr, 2});
        if kTr == 1
            if FileExists([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
                load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
            else
               FindCommonKenjiClu(roi, arena);
               load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
            end
        end

        if isempty(commonClus), return; end
        if length(commonClus) > 2
            cellPairs = nchoosek(commonClus, 2);
        else
            cellPairs = [];
        end
        fprintf(['\n trial :' gt.trialName ]);
        gta{kTr} = gt.LoadPF;
        commonClus = commonClus(ismember(commonClus, gt.pyrCluIdx(gta{kTr}.pfObject.idealPFUnits))); 
        sparsity(:, kTr) = gta{kTr}.pfObject.sparsity(ismember(gta{kTr}.pfObject.acceptedUnits, commonClus));
        roiPFPairs{kTr} = nchoosek(commonClus, 2);
    end
    for lClu = 1 : length(commonClus)
        if IF_SmthRM
            for mTr = 1 : nTrials
                srm = gta{mTr}.pfObject.smoothRateMap(:,:,ismember(gta{mTr}.pfObject.acceptedUnits, commonClus(lClu)));
                mapEntropy(lClu, mTr) = Entropy(srm);
                
            end
        end
    end
    goodUnitsIdx = mapEntropy < 9 | sparsity < .07;
    gu = true(size(goodUnitsIdx, 1), 1);
    for kk = 1 : size(goodUnitsIdx, 2)
       gu = gu & goodUnitsIdx(:, kk);
    end
    out.goodUnits = commonClus(gu);
    if nTrials == 1 | length(out.goodUnits) == 1 | isempty(out.goodUnits)
        nPairs = 0; 
        nCommonPairs = 0;
    else 
        commonPFPairs = roiPFPairs{refTr};
        commonPFPairs = commonPFPairs(ismember(commonPFPairs, nchoosek(out.goodUnits, 2), 'rows'), :);
        nPairs = size(commonPFPairs, 1);
        nCommonPairs = nchoosek(length(commonClus), 2);
        out.commonPairs = commonPFPairs;
    end
    nPairs
    nTrials
    fp = fopen(['~/data/analysis/kenji/', mfilename, GenFiletag(arena, roi)], 'a');
    fprintf(fp, ['\n', repmat('::', 1, 10), filebase, repmat('::', 1, 10)]);
    fprintf(fp, ['\n #Trials :' num2str(nTrials)]);
    fprintf(fp, ['\n # good PF units :' num2str(length(commonClus))]);
    fprintf(fp, ['\n # common pairs :' num2str(nCommonPairs)]);
    fprintf(fp, ['\n # overalpping pairs  :' num2str(nPairs)]);
    fclose(fp);
end


