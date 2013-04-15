function FindCommonKenjiClu(varargin)

    [roi, arena] = DefaultArgs(varargin, {{'CA3'}, {'bigSquare'}});
    searchStruct.roi = roi;
    searchStruct.arena = arena;
    list = SearchKenji(searchStruct);
    for i = 1 : length(list)
        trialNames = list(strcmp(list(:, 1), list(i,1)), 2);
        if ~isempty(trialNames)
            fprintf('\n ********* filebase: %s ************** \n', list{i});
            cnt = 1;
            for kTr = 1 : length(trialNames)
                fprintf('subtrial %d of %d \n', kTr, length(trialNames));
                try
                    gt = GenericTrial(list{i}, trialNames{kTr});
                    gt = gt.LoadPF;
                    kClu{kTr} = gt.pyrCluIdx(gt.pfObject.idealPFUnits);
                    if cnt == 1, commonClus = kClu{1}; end
                    commonClus = intersect(commonClus, kClu{cnt});
                    cnt = cnt + 1;
                catch err
                    fprintf('error  !!!!! \n');
                end
            end
            %          for mTr = 1 : cnt - 1
            %   sIdx{mTr} = ismember(kClu{mTr}, commonClus);
            % end
            save(['~/data/analysis/kenji/', list{i, 1},'.', char(roi), '.commonClus.mat'], 'commonClus');
        end
    end
end

