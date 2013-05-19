function genericTrial  = Load(genericTrial, loadProperty, varargin)
    % general load function for Generictrial class
    % returns GenericTrial object with specified fields loaded
    % loadProperty - cell array {{'<load property1 name>', '<arg1>', '<arg2>'}, {<load property2 name>', '<arg1>', '<arg2>'}}

    nPropsToLoad = length(loadProperty);
    for kProperty = 1 : nPropsToLoad
        curLoadProperty = loadProperty{kProperty};
        %keyboard;
        switch curLoadProperty{1}
            case 'PF'
                if strcmp(genericTrial.datasetType, 'MTA')
                    mtaPFObj = LoadMTAPFObject(genericTrial.filebase, genericTrial.trialName);
                    genericTrial.pfObject = GenericPF(mtaPFObj);
                else
% keyboard;
                if isempty(genericTrial.trialName)
                         error(['\n trial name not specified  !!! ' ...
                                '\n']);  
               end
                    if length(curLoadProperty) > 1
                        genericTrial.pfObject = GenericPF(genericTrial, curLoadProperty{2:end});
                    else
                        genericTrial.pfObject = GenericPF(genericTrial);
                    end
                end
            if FileExists([genericTrial.paths.analysis, genericTrial.filebase, 'SelectCells.mat'])
                load([genericTrial.paths.analysis, genericTrial.filebase, '.SelectCells.mat']);
            else
                genericTrial.pyrCluIdx = 1 : size(genericTrial.pfObject.rateMap, 2);
            end
            case 'CluRes'
                fprintf('\n loading CluRes ... ');
                if strcmp(genericTrial.datasetType, 'MTA')
                    if length(curLoadProperty) > 1
                        genericTrial = genericTrial.Convert2Generic(MTATrial(genericTrial.filebase, {{'CluRes', curLoadProperty{2:end}},genericTrial.trialName}));
                    else
                        genericTrial = genericTrial.Convert2Generic(MTATrial(genericTrial.filebase, {{'CluRes', genericTrial.sampleRate}},genericTrial.trialName));
                    end
                else
                   
                    [res, clu, map] = LoadCluRes([genericTrial.paths.data, genericTrial.filebase]);
               %      if strcmp(genericTrial.datasetType, 'kenji')
%                         [res, resIdx] = SelectPeriods(res, ConvertFs(genericTrial.trialPeriods, genericTrial.lfpSampleRate, genericTrial.sampleRate), 'd');
%                         clu = clu(resIdx);
%                     end
                    genericTrial.elClu = map(:,[2, 3]);
                    genericTrial.res = res;
                    genericTrial.clu = clu;
                    if length(curLoadProperty) > 1 % the second arg specifies if only clures for the trial is to be loaded
                        if curLoadProperty{2}
                            [genericTrial.res, resIdx] = SelectPeriods(res, genericTrial.trialPeriods .* genericTrial.sampleRate ./ genericTrial.lfpSampleRate, 'd', 1, 1);
                            genericTrial.clu = clu(resIdx);
                        end
                    end
                end
                fprintf(' done !! \n');       
            case 'ccgSegsA2B' % load ccg for A -> B & B -> A epochs
                if ~isempty(genericTrial.trialSubType)
                    filename = [genericTrial.paths.analysis, genericTrial.filebase, '.CCGofSegments.', genericTrial.trialName, '.', genericTrial.trialSubType '.mat'];
                else
                    filename = [genericTrial.paths.analysis, '/' genericTrial.filebase, '.CCGofSegments.', genericTrial.trialName, '.mat'];
                end
                if FileExists(filename)
                    fprintf('\n loading %s \n', filename);
                    ccg = load(filename);
                    genericTrial.ccg = ccg.ccgSegsA2B;
                else
                    error('\n ccg not computed \n')
                end
        end




    end
end
