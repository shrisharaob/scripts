list = importdata('~/data/kenji/list');
load('~/data/kenji/Beh_time_ind.mat');
    for i = 1 : length(list)
        nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) , 4);
        if ~isempty(nonSleepTrials{i})
            fprintf('\n ********* filebase: %s ************** \n', list{i});
            subTrialNames = nonSleepTrials{i};
            for kSubTr = 1 : length(subTrialNames)
                fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
                if FileExists(['~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl'])
                    try 
                        gt = GenericTrial(list{i}, subTrialNames{kSubTr});
%                         gt = gt.Load({{'PF'}});
%                         pfPars = FindPFPars(gt.pfObject, 1: length(gt.pfObject.rateMap));
                        system(['rm ' gt.paths.analysis, gt.filebase, gt.trialName, '.FindPFPars.mat']);
                        save([gt.paths.analysis, gt.filebase, gt.trialName, '.FindPFPars.mat']);
                    catch err
                        fp = fopen('~/data/analysis/kenji/logFilePFPar','a');
                        str = sprintf(['\n' list{i} : subTrialNames{kSubTr} 'not done']);
                        str = [str '\n' err.message];
                        fwrite(fp, str);
                        fclose(fp);
                    end
                end
                fclose('all');
            end
        end
    end