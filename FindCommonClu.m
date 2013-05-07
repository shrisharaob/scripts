function FindCommonClu(varargin)
% this script returns clusters which are active acros trials in all filebases
% [filebases, datasetType,roi, arena, sparsityThresh]
% filebases - cell array containing list of fbs    
    
    [filebases, datasetType,roi, arena, sparsityThresh] = ...
        DefaultArgs(varargin, {{}, 'kenji', {'CA3'}, {'bigSquare'}, 0.35});

    switch datasetType
      case 'kenji'
        searchStruct.roi = roi;
        searchStruct.arena = arena;
        list = SearchKenji(searchStruct);
        filebases = unique(list(:, 1));
      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
        if ~iscell(filebases), filebases = cellstr(filebases); end
    end
    
    for i = 1 : length(filebases)
        switch datasetType 
          case 'kenji'
            trialNames = list(strcmp(list(:, 1), filebases(i)), 2);
          case 'MTA'
            trialNames = GenericTrial.MTATrialNames(filebases{i});
        end
        if ~isempty(trialNames)
            fprintf('\n ********* filebase: %s ************** \n', filebases{i});
            cnt = 1;
            for kTr = 1 : length(trialNames)
                fprintf('subtrial %d of %d \n', kTr, length(trialNames));
                try
                    gt = GenericTrial(filebases{i}, trialNames{kTr});
                    gt = gt.LoadPF;
                    kClu{kTr} = gt.pfObject.acceptedUnits(gt.pfObject.sparsity < sparsityThresh);
                    if cnt == 1, commonClus = kClu{1}; end
                    commonClus = intersect(commonClus, kClu{cnt});
                    cnt = cnt + 1;
                catch err
                    fprintf('error  !!!!! \n');
                end
            end
            filetag = GenFiletag(roi, arena);
            save([gt.paths.analysis, gt.filebase, filetag, '.commonClus.mat'], 'commonClus');
        end
    end
end

