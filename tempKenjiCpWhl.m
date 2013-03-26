% 
% 
% %% copy whl files 
% list = importdata('~/data/kenji/list');
% load('~/data/kenji/Beh_time_ind.mat');
% % for i = 1 : length(list)
% % nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) , 4);
% % if ~isempty(nonSleepTrials{i})
% %     subTrialNames = nonSleepTrials{i};
% %     for kSubTr = 1 : length(subTrialNames)
% %         system(['cp -v /data/homes/antsiro/data/blab/kenji/whl/', subTrialNames{kSubTr} '.whl.gz ~/data/kenji/' list{i} '/']);
% %         system(['gzip -d ~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl.gz']);
% % %         system(['rm -v ~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl.gz']);
% %     end
% % end
% % end
% 
% %%
% 
% % for i = 1 : length(list)
% %     nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) , 4);
% %     if ~isempty(nonSleepTrials{i})
% %         for kSubtrial = 1 : length(whos)
% %             subTrialNames = nonSleepTrials{i};
% %             duration = Beh(strcmp(Beh(:,2),list{i}), 7);
% %             trialPeriods =[0; cumsum(str2num(char(Beh(strcmp(Beh(:,2),list{i}), 7))))];
% %             trialPeriods = [trialPeriods(1:end-1), trialPeriods(2:end)];
% %             trialPeriods = trialPeriods .* genericTrial.lfpSampleRate; % trialPeriods @lfp fs
% %         end
% %     end
% % end
% 
% 
% 
% %%
% for i = 1 : length(list)
%    files = dir(['~/data/kenji/' list{i}]);
%     nWhlFiles = 0;
%     for kFile = 1 : length(files)
%         if ~files(kFile).isdir
%             [~,~,ext] = fileparts(files(kFile).name);
%             if strcmp(ext, '.whl')
%                 nWhlFiles = nWhlFiles + 1;
%                 whlFileNames{nWhlFiles} = files(kFile).name;
%             end
%         end
%     end
%     whlFiles{i} = whlFileNames;
%     
%     
% end
% 
% %% copute place fields for all kenji
% list = importdata('~/data/kenji/list');
% load('~/data/kenji/Beh_time_ind.mat');
%     for i = 1 : length(list)
%         nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) , 4);
%         if ~isempty(nonSleepTrials{i})
%             fprintf('\n ********* filebase: %s ************** \n', list{i});
%             subTrialNames = nonSleepTrials{i};
%             for kSubTr = 1 : length(subTrialNames)
%                 fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
%                 if FileExists(['~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl'])
%                     try 
%                         gt = GenericTrial(list{i}, subTrialNames{kSubTr});
%                         gt = gt.Load({{'PF', [], [],1}});
%                     catch err
%                         fp = fopen('~/data/analysis/kenji/logFile','a');
%                         str = sprintf(['\n' list{i} : subTrialNames{kSubTr} 'not done']);
%                         str = [str '\n' err.message];
%                         fwrite(fp, str);
%                         fclose(fp);
%                     end
%                 end
%                 fclose('all');
%             end
%         end
%     end
%     
%     
%     %% borwse kenji data
%     list = importdata('~/data/kenji/list');
%     load('~/data/kenji/Beh_time_ind.mat');
%     for i = 1 : length(list)
%          nonSleepTrials{i} = Beh(strcmp(Beh(:,5) ,'bigSquare') | strcmp(Beh(:,5) ,'midSquare') | strcmp(Beh(:,5) ,'linear') & strcmp(Beh(:,2),list{i}) , 4);
%         if ~isempty(nonSleepTrials{i})
%             fprintf('\n ********* filebase: %s ************** \n', list{i});
%             subTrialNames = nonSleepTrials{i};
%             for kSubTr = 1 : length(subTrialNames)
%                 fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
%                 if FileExists(['~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl'])
%                     try 
%                         gt = GenericTrial(list{i}, subTrialNames{kSubTr});
%                         gt = gt.Load({{'PF'}});
%                         PlotPlaceFields(gt.pfObject, [],[],1);
%                         set(gcf, 'Name', [list{i}, ' ' gt.maze.name]);
%                     catch
%                         fprintf(' ********* ')
%                     end
%                 end                   
%             end
%         end
%      end
%      
     %% kenji sum
        list = importdata('~/data/kenji/list');
    load('~/data/kenji/Beh_time_ind.mat');
    arenas = {'bigSquare', 'midSquare', 'linear'};
     fp = fopen('~/data/analysis/kenji/stats','a');
    fwrite(fp, ['********* stats for ' arenas{1}, ' ', arenas{2}, ' ', arenas{3} '****** \n']); 
    fclose(fp);
    for i = 1 : length(list)
        nonSleepTrials{i} = Beh(strcmp(Beh(:,5) ,'bigSquare') | strcmp(Beh(:,5) ,'midSquare') | strcmp(Beh(:,5) ,'linear') & strcmp(Beh(:,2),list{i}) , 4);
        if ~isempty(nonSleepTrials{i})
            fprintf('\n ********* filebase: %s ************** \n', list{i});
            subTrialNames = nonSleepTrials{i};
            for kSubTr = 1 : length(subTrialNames)
                
                fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
                gt = GenericTrial(list{i}, subTrialNames{kSubTr});
                files = dir(gt.paths.data);
                nCluFiles = 0;
                for kFile = 1 : length(files)
                    if ~files(kFile).isdir
          %             if ~isempty(regexp(files(kFile).name, [arg.filebase '\.clu\.']))
                        if FileExists([gt.paths.data, gt.filebase, '.clu.' num2str(nCluFiles + 1)])
                            nCluFiles = nCluFiles + 1;
                            [clu, ~] = LoadClu([gt.paths.data, gt.filebase, '.clu.', num2str(nCluFiles)]);
                            clu = SelectPeriods(clu, round(gt.trialPeriods .* gt.sampleRate ./ gt.lfpSampleRate)+1,'d');
                            if kSubTr == 1, uniqueClu = unique(clu);
                            else uniqueClu = intersect(uniqueClu, unique(clu)); end
                        end
                    end
                end
    
                nClus = NClusters(gt);
            end
                fp = fopen('~/data/analysis/kenji/stats','a');
                str = ['\n fielbase : ' list{i} '\n'];
                str = [str, 'nCluster = ' nClus '\n common units :' uniqueClu '\n'];
                fwrite(fp, str);
                fclose(fp);  
                
        end
    end
%                 
%      
%      
%     
%     