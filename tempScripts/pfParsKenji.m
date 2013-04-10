        %% compute pf pars kenji
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
                            fprintf('\n computing pfPars\n');
                            gt.pfObject.FindPFPars(1:length(gt.pfObject.rateMap),[],[],1);
                            if kSubTr == 1
                                fp = fopen('~/data/analysis/kenji/pfParsLogFile','a');
                                fprintf(fp, ['\n', repmat('*',1 ,60), '\n' gt.filebase ]);
                                fclose(fp);
                            end
                            %                          PlotPlaceFields(gt.pfObject, [],[],1);
                            %set(gcf, 'Name', [list{i}, ' ' gt.maze.name]);
                            fprintf(fp, ['\n ::: '  gt.trialName]); 
                        catch
                            fprintf('\n whl |& eegTime not found !!!') 
                        end
                    end
                end
            end
            fclose('all');
end