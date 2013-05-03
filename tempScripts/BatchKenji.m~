function BatchKenji(funcHandle, varargin)
    %% evaluate func for all kenji
    [roi, arena, IF_LOAD_GT, funcArgs, type, IF_SAVE] = ...
        DefaultArgs(varargin, {{'CA1'}, {'bigSquare'}, 0, {}, 'passFb', 0});
    %   list = importdata('~/data/kenji/list');
    %   load('~/data/kenji/Beh_time_ind.mat');
    sK.roi = roi;
    sK.arena = arena;
    matches = SearchKenji(sK);
    filebases = unique(matches(:,1));
    for i = 1 : length(filebases);
        trialNames = matches(strcmp(matches(:,1), filebases(i)), 2);
        switch type
            case 'allTrials'
                fprintf(['\n ********* filebase: %s ************** \n'], filebases{i});
                for lTr = 1 : length(trialNames)
                    if ~isempty(trialNames{lTr})
                        fprintf('\n subtrial %d of %d \n', lTr, length(trialNames));
                        if FileExists(['~/data/kenji/' filebases{i} '/' trialNames{lTr} '.whl']) && FileExists(['~/data/kenji/whl/' filebases{i} '.eegTime'])
                            try
                                if IF_LOAD_GT
                                    gt = GenericTrial(filebases{i}, trialNames{lTr});
                                    feval(funcHandle, gt, funcArgs{:});
                                else
                                    feval(funcHandle, filebases{i}, funcArgs{:});
                                end
                                if lTr == 1
                                    fp = fopen(['~/data/analysis/kenji/', func2str(funcHandle)], 'a');
                                    fprintf(fp, ['\n', repmat('*',1 ,60), '\n' gt.filebase ]);
                                end
                                fprintf(fp, ['\n ::: '  gt.trialName]);
                                fclose(fp)
                            catch err
                                
                            end
                        else
                            fprintf('\n error !!!')
                        end
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
                            feval(funcHandle, filebases{i}, funcArgs{:});
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
                %                goodUnits(i).filebase = filebases{i};
                try
                    fout =feval(funcHandle, filebases{i}, funcArgs{:});
                    %                    goodUnits(i).clu = fout;
                catch err
                    fprintf('error .......  !!!!!!');
                end
                
        end
    end
    if IF_SAVE
        save(['~/data/analysis/kenji/', func2str(funcHandle), GenFiletag(arena, roi) 'mat'], 'goodUnits');
    end
    fclose('all');
end

