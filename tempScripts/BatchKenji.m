function BatchKenji(funcHandle, varargin)
%% compute funch
        list = importdata('~/data/kenji/list');
        load('~/data/kenji/Beh_time_ind.mat');
        for i = 1 : length(list);
            nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) , 4);
            if ~isempty(nonSleepTrials{i})
                fprintf('\n ********* filebase: %s ************** \n', list{i});
                subTrialNames = nonSleepTrials{i};
                
                for kSubTr = 1 : length(subTrialNames)
                    fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
                    if FileExists(['~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl']) && FileExists(['~/data/kenji/whl/' list{i} '.eegTime'])
                        try
                            gt = GenericTrial(list{i}, subTrialNames{kSubTr});
                            gt = gt.Load({{'PF'}});
                            feval(funcHandle, gt, varargin);                            
                            if kSubTr == 1
                                fp = fopen(['~/data/analysis/kenji/', func2str(funcHandle)], 'a');
                                fprintf(fp, ['\n', repmat('*',1 ,60), '\n' gt.filebase ]);
                            end
                            fprintf(fp, ['\n ::: '  gt.trialName]); 
                            fclose(fp)
                        catch err
                        
                        end
                    else
                        fprintf('\n whl |& eegTime not found !!!') 
                    end
                end
            end
        end
        fclose('all');
    end
