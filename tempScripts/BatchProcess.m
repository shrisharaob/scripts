function out = BatchProcess(funcHandle, varargin)
% function BatchKenji(funcHandle, varargin)
% [datasetType, roi, arena, IF_LOAD_GT, funcArgs, type, IF_SAVE, options]
% '', {'CA3'}, {'bigSquare'}, 0, {}, 'passFb', 0
% evaluate func for all filebases

    [datasetType, roi, arena, IF_LOAD_GT, funcArgs, type, IF_SAVE, options] = ...
        DefaultArgs(varargin, {'', {'CA3'}, {'bigSquare'}, 0, {}, 'passFb', 0, []});

    switch datasetType
      case 'kenji'
        searchStruct.roi = roi;
        searchStruct.arena = arena;
        matches = SearchKenji(searchStruct);
        filebases = unique(matches(:, 1));
        %   matches = matches(strcmp(matches(:, 1), filebase), :);
        % trialNames = matches(:, 2);

      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
        %if ~iscell(filebases), filebases = cellstr(filebases); end
        load('~/data/analysis/MTA_fb_list.mat');
      otherwise
        fprintf('\n data set type not specified \n');
        return;
    end
    
    if ~isempty(options), 
        poolVar = options.poolVar; 
        poolArray = [];
        poolArrayId = {};
        poolCounter = 0;
    end
    
    for i = 1 : length(filebases);
        switch datasetType
          case 'kenji'
            trialNames = matches(strcmp(matches(:,1), filebases(i)), 2);
          case 'MTA'
            trialNames = MTATrialNames(filebases{i});
        end

        switch type
          case 'allTrials'
            fprintf(['\n ********* filebase: %s ************** \n'], filebases{i});
            for lTr = 1 : length(trialNames)
                if ~isempty(trialNames{lTr})
                    fprintf('\n subtrial %d of %d \n', lTr, length(trialNames));
                    if strcmp(datasetType, 'kenji')
                        if FileExists(['~/data/kenji/' filebases{i} '/' trialNames{lTr} '.whl']) && FileExists(['~/data/kenji/whl/' filebases{i} '.eegTime'])
                        else
                            fprintf('\n error !!!')
                        end
                    end
                    try
                        if IF_LOAD_GT
                            gt = GenericTrial(filebases{i}, trialNames{lTr});
                            fout = feval(funcHandle, gt, funcArgs{:});
                        else
                            feval(funcHandle, filebases{i}, funcArgs{:});
                        end
                        %if lTr == 1
                        %    fp = fopen(['~/data/analysis/kenji/', func2str(funcHandle)], 'a');
                        %   fprintf(fp, ['\n', repmat('*',1 ,60), '\n' gt.filebase ]);
                        %end
                        %fprintf(fp, ['\n ::: '  gt.trialName]);
                        %fclose(fp)
                    catch err
                        fprintf('error');
                        keyboard;
                    end
                end
                if IF_SAVE
                    save([gt.paths.analysis, gt.filebase, '.' gt.trialName, GenFiletag(roi, arena), func2str(funcHandle),  '.mat'], 'fout');
                end
            end
            
          case 'allBases'
            fprintf(['\n ********* filebase: %s ************** \n'], filebases{i});
            for lTr = 1 : length(trialNames)
                if ~isempty(trialNames{lTr})
                    if IF_LOAD_GT
                        if FileExists(['~/data/kenji/' filebases{i} '/' trialNames{lTr} '.whl']) && FileExists(['~/data/kenji/whl/' filebases{i} '.eegTime'])
                            fprintf('\n subtrial %d of %d \n', lTr, length(trialNames));
                            try
                                gt = GenericTrial(filebases{i}, trialNames{lTr});
                                feval(funcHandle, gt, funcArgs{:});
                                if lTr == 1
                                    fp = fopen(['~/data/analysis/kenji/', func2str(funcHandle)], 'a');
                                    fprintf(fp, ['\n', repmat('*',1 ,60), '\n' gt.filebase ]);
                                end
                            catch err
                                fprintf('error ........ !');
                            end
                        end
                    else
                        try 
                            fout = feval(funcHandle, filebases{i}, funcArgs{:});
                            fp = fopen(['~/data/analysis/kenji/', func2str(funcHandle)], 'a');
                            fprintf(fp, ['\n ::: '  filebases{i}]);
                            fclose(fp);
                        catch err
                        end
                    end
                end
            end
          case 'passFb' % pass the filebase as an argument and all the subtrials are taken care 
            fprintf(['\n ********* filebase: %s ************** \n'], filebases{i});
            try
                fout =feval(funcHandle, filebases{i}, funcArgs{:});
            catch err
                
                fprintf('error .......  !!!!!!');
            end
    
          case 'pool' % concatinate the returned variable
            fprintf(['\n ********* filebase: %s ************** \n'], filebases{i});
            for lTr = 1 : length(trialNames)
                if ~isempty(trialNames{lTr})
                    fprintf('\n subtrial %d of %d \n', lTr, length(trialNames));
                    trialNames{lTr}
                    if strcmp(datasetType, 'kenji')
                        if FileExists(['~/data/kenji/' filebases{i} '/' trialNames{lTr} '.whl']) && FileExists(['~/data/kenji/whl/' filebases{i} '.eegTime'])
                        else
                            fprintf('\n error !!!')
                        end
                    end
                    try
                        if IF_LOAD_GT
                            gt = GenericTrial(filebases{i}, trialNames{lTr});
                            fout = feval(funcHandle, gt, funcArgs{:});
                            poolCounter = poolCounter + 1;
                        else
                            fout = feval(funcHandle, filebases{i}, funcArgs{:});
                            poolCounter = poolCounter + 1;
                        end
                        if isempty(eval(['fout.' poolVar])), continue; end
                        nRowsOld = size(poolArray, 1);
                        eval(['poolArray = [poolArray; fout.' poolVar '];']);
                        eval(['nRows = size(fout.', poolVar, ', 1);']);
                        if strcmp(datasetType, 'kenji'), tArena = SearchKenji(gt.trialName); end
                        tArena = tArena{2};
                        poolArrayId = [poolArrayId; {filebases{i}, trialNames{lTr}, tArena ,[1 : nRows] + nRowsOld}];
                       %  if lTr == 1
%                             fp = fopen(['~/data/analysis/kenji/', func2str(funcHandle)], 'a');
%                             fprintf(fp, ['\n', repmat('*',1 ,60), '\n' gt.filebase ]);
%                         end
%                         fprintf(fp, ['\n ::: '  gt.trialName]);
%                         fclose(fp)
                    catch err
                        keyboard;
                        fprintf('error');
                    end
                end
            end
        end
        
    end
    if strcmp(type, 'pool')
        out.poolArrayId = poolArrayId;
        out.poolArray = poolArray;
    end
    if IF_SAVE
        save(['~/data/analysis/', datasetType, func2str(funcHandle), GenFiletag(arena, roi) 'mat'], 'fout');
    end
%     if IF_SAVE
%         save([gt.paths.analysis, func2str(funcHandle), GenFiletag(roi, arena),   'mat'], '');
%     end

    fclose('all');
end

