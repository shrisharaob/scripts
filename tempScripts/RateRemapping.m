list = importdata('~/data/kenji/list');
load('~/data/kenji/Beh_time_ind.mat');

    for i = 1 : length(list)
        nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) & (strcmp(Beh(:,5) ,'bigSquare') | strcmp(Beh(:,5) ,'midSquare')) , 4);
        if ~isempty(nonSleepTrials{i})
            fprintf('\n ********* filebase: %s ************** \n', list{i});
            fp = fopen('~/data/analysis/kenji/logFilePF','a');
            fprintf(fp, '\n filebase: %s', list{i});
            subTrialNames = nonSleepTrials{i};
            for kSubTr = 1 : length(subTrialNames)
                fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
                if FileExists(['~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl'])
                    try 
                        gt = GenericTrial(list{i}, subTrialNames{kSubTr});
                        fprintf(fp, '\n trialname : %s', subTrialNames{kSubTr});fclose(fp);
                    catch err
                        fp = fopen('~/data/analysis/kenji/logFilePF','a');
                        str = sprintf(['\n' subTrialNames{kSubTr}: 'not done']);
                        str = [str '\n' err.message];
                        fwrite(fp, str);
                        fclose(fp);
                    end
                end
 
            end
        fp = fopen('~/data/analysis/kenji/logFilePF','a');
        fprintf(fp, ['\n' sprintf( repmat('-+', 1, 30))]);
        fclose(fp);
        fclose('all');
        end

    end