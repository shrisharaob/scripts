function genericTrial = ProcessKenji(genericTrial)


    % find all the tracking files 
    files = dir(genericTrial.paths.data);
    nWhlFiles = 0;
    for kFile = 1 : length(files)
        if ~files(kFile).isdir
            [~,~,ext] = fileparts(files(kFile).name);
            if strcmp(ext, '.whl')
                nWhlFiles = nWhlFiles + 1;
                whlFileNames{nWhlFiles} = files(kFile).name;
            end
        end
    end
    %
    load('~/data/kenji/Beh_time_ind.mat');
    trialNames = Beh(strcmp(Beh(:,2),genericTrial.filebase), 4);
    rowNo  = find(~cellfun(@isempty, regexp(genericTrial.trialName, trialNames)));
    times = load(['~/data/kenji/whl/', genericTrial.filebase, '.eegTime']);
    trialPeriods = [0; round(times(:,2))];
%     trialPeriods =[0; cumsum(str2num(char(Beh(strcmp(Beh(:,2),genericTrial.filebase), 7))))];
  trialPeriods = [trialPeriods(1:end-1), trialPeriods(2:end)-1];
%     trialPeriods = trialPeriods .* genericTrial.lfpSampleRate; % trialPeriods @lfp fs
%     
    genericTrial.trialPeriods = trialPeriods(rowNo,:);
    
end
        

