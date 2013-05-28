classdef GenericPF
% general Place field class
% [datasetType, trialName, IF_OVERWRITE, state, markerNo, smoothSigma] 
% Defaults :  {'default', 'df', 0, 'RUN', 1, 0.03 px});  
%-------
    % History:
    %  Shrisha - Created
    
    properties (SetAccess = protected, GetAccess = public)
        % tag identifying the type of dataset
        datasetType;
       
        % dir name containing all the recording data
        filebase;
    end
    
    properties
        
        % file paths
        paths;
        
        trialName;
        
        trialSubType;
        
        % behavioural state
        state
        
        % trial begin & end times
        trialPeriods;
        
        % vaild position tracking periods 
        goodPosPeriods;
        
        % place field parameters
        

        rateMap;
        
        %
        occupancy;
        
        % spatial bins
        xBin;
        
        %
        yBin;
        
        %
        nBins;
        
        % arena specifications
        maze;

        %
        smoothRateMap;
        
        % center of mass
        com;
       
       % maximum peak of the smoothed rateMap
        ratePk; 
        
        % peak location 
        pkLoc;
        
        % distance of the peaks
        pkDist;
                
        % accepted units are the place cells which satisfy an absolute rate
        % threshold and overlap criterion
        acceptedUnits;
        
        % pairs of place cells which have overlapping fileds
        selectedPairs;
        
        % eucleadian distance of the com's of Place cells with overlapping
        % place fields
        comDist;
        
        % place cells stable/ideal in the sense that they satisfy sparsity
        % spatial coherence constraints
        idealPFUnits;
        
        % ideal place cells with overlapping fields
        idealPFPairs;
        
        % spatial coherence of firing rate
        spatialCoherence;
        
        % measure for quantifying the spatial localization of firing
        sparsity;

        % 
        jntEntropy;
    end
    
    
    methods
        
        function genericPF = GenericPF(arg, varargin)
            % class constructor
            
            [datasetType, trialName, IF_OVERWRITE, state, markerNo, smoothSigma] = DefaultArgs(varargin, {'default', 'df', 0, 'RUN', 1, 0.03});
            mazePars = [];
            if isempty(arg)
                genericPF.filebase = [];
            elseif isa(arg, 'GenericTrial')
                trialProps = properties(arg);
                pfProps = properties('GenericPF');
                for kProp = 1 : length(pfProps)
                    simPropIdx = strcmp(pfProps{kProp}, trialProps);
                   	if any(simPropIdx)
                       genericPF.(pfProps{kProp}) = arg.(trialProps{simPropIdx});
                    end
                end
                
            elseif isa(arg, 'GenericPF')
                propertyNames = properties(arg);
                for kProperty = 1 : length(propertyNames)
                    genericPF.(propertyNames{kProperty}) = arg.(propertyNames{kProperty});
                end
            elseif isa(arg, 'MTAPlaceField')
                genericPF = genericPF.Convert2GenericPF(arg);
            elseif isa(arg, 'MTATrial')
                genericPF = genericPF.Convert2GenericPF(arg);
                mazePars.name = arg.Maze.name;
                mazePars.boundaries = arg.Maze.boundaries;
            else
                % creat GenericPF object
                genericPF.datasetType = datasetType;
                genericPF.filebase = arg;
                genericPF.trialName = trialName;
                genericPF.paths.data = ['~/data/', genericPF.filebase, '/'];
                if ~DirExists('~/data/', genericPF.filebase)
                    % if filebase not found in ~/data/ search in ~/data/nlx/ folder
                    if DirExists('~/data/nlx/', genericPF.filebase)
                        genericPF.paths.data = ['~/data/nlx/', genericPF.filebase, '/'];
                        genericPF.datasetType = 'MTA';
                    end
                end
            end
            if isempty(genericPF.trialSubType)
                fileName = [genericPF.filebase, '.PF.', genericPF.trialName '.mat'];
            else
                fileName = [genericPF.filebase, '.PF.', genericPF.trialName '.' genericPF.trialSubType '.mat'];
            end
            
            if ~FileExists([genericPF.paths.analysis, fileName]) &&  ~strcmp(genericPF.datasetType, 'MTA')
                fprintf('\n rate maps not precomputed \n'); IF_OVERWRITE = 1; 
            end
            if IF_OVERWRITE % compute rateMaps and place field parameters
                if isa(arg, 'GenericTrial')
                    trial = arg;
                else
                    trial  = GenericTrial(genericPF.filebase,[], genericPF.trialName, mazePars);
                end
                trial = trial.Load({{'CluRes'}});

                switch trial.datasetType
                  case  'MTA'
                    mtaTrial = MTATrial(trial.filebase, [], trial.trialName);
                    if strcmp(state, 'RUN')
                        state = 'walk'; 
                    else
                        error('\n no %s state in filebase %s', state, trial.filebase); 
                    end
                    statePeriods = mtaTrial.statePeriods(state); %% state periods starts with 0
                    markerNo = 7;
                    posStatePeriods = round(statePeriods .* trial.trackingSampleRate ./ trial.lfpSampleRate) + 1;
                    statePeriods = round(statePeriods .* trial.sampleRate ./ trial.lfpSampleRate) + 1;
                    trialStartTime_pos = 0;
                  case 'default'
                    statePeriods = load([trial.paths.data, trial.filebase '.sts.', state]);
                    posStatePeriods = round(statePeriods .* trial.trackingSampleRate  ./ trial.lfpSampleRate) + 1;
                    statePeriods = round(statePeriods .* trial.sampleRate ./ trial.lfpSampleRate) + 1;
                    trialStartTime_pos = 0;
                    markerNo = 1;
                  case  'kenji'
                    statePeriods = load([trial.paths.data, trial.filebase '.sts.', state]); % @lfp fs
                    statePeriods = IntersectRanges(statePeriods, trial.trialPeriods);
                    posStatePeriods = round(statePeriods .* trial.trackingSampleRate  ./ trial.lfpSampleRate) + 1;
                    statePeriods = round(statePeriods .* trial.sampleRate ./ trial.lfpSampleRate) + 1;       
                    trialStartTime_pos =  round(trial.trialPeriods(1, 1) .*  trial.trackingSampleRate ./ trial.lfpSampleRate) + 1;
                    posStatePeriods = IntersectRanges(posStatePeriods, trial.goodPosPeriods + trialStartTime_pos); 
                    markerNo = 1;
                end


                genericPF.state = state;
                nClus= length(unique(trial.clu));
                genericPF.rateMap = cell(1, nClus);
                genericPF.occupancy = cell(1, nClus);


                pos = SelectPeriods(sq(trial.position(:,markerNo,:)), posStatePeriods - trialStartTime_pos, 'c');
                %convert spike times to to tracking sample rate
                res = round(trial.res .* trial.trackingSampleRate ./ trial.sampleRate) + 1;
                [kRes, resIdx] = SelectPeriods(res, posStatePeriods, 'd',1,1);
                for kClu = 1 : nClus
                    fprintf('\n computing rate maps for unit %d of %d units \n', kClu, nClus);
                    nSpikes = length(kRes(trial.clu(resIdx) == kClu));
                    if nSpikes > 10
                        cluKRes = trial.clu(resIdx);
                        %        pos(isnan(pos(:, 1)), :) = []; % !!!!!!!!!!!!!!!!!!! 
                        [genericPF.rateMap{kClu}, genericPF.occupancy{kClu}, xBin, yBin] = ...
                            GenericPF.ComputeRateMap(trial, kRes(cluKRes == kClu), pos, [], smoothSigma);
                        if ~isempty(xBin), genericPF.xBin = xBin; end 
                        if ~isempty(yBin), genericPF.yBin = yBin; end
                    end
                end
                rm = genericPF.rateMap;
                xBin = genericPF.xBin;
                yBin = genericPF.yBin;
                save([genericPF.paths.analysis, fileName], 'rm', 'xBin', 'yBin');
                clear rm;
                if isempty(genericPF.trialSubType)
                fileName = [genericPF.filebase, '.FindPFPars.', genericPF.trialName '.mat'];
                else
                fileName = [genericPF.filebase, '.FindPFPars.', genericPF.trialName '.' genericPF.trialSubType '.mat'];
                end
                fprintf('\n computing place field parameters\n');
                pfPars = GenericPF.FindPFPars(genericPF, 1 : nClus);
                save([genericPF.paths.analysis, fileName], 'pfPars')
            else
                if strcmp(genericPF.datasetType, 'MTA') && isempty(genericPF.rateMap) && isempty(genericPF.xBin)
                    fprintf('\n loading rate maps ... \n');
                    mtaPFObj = LoadMTAPFObject(genericPF.filebase, genericPF.trialName);
                    gPF = genericPF.Convert2GenericPF(mtaPFObj);
                    genericPF.rateMap = gPF.rateMap;
                    genericPF.xBin = gPF.xBin;
                    genericPF.yBin = gPF.yBin;
                    genericPF.nBins = gPF.nBins;
                    clear gPF;
                elseif strcmp(genericPF.datasetType,'kenji')
                    fprintf('\n loading rate maps ... \n');
                    load([genericPF.paths.analysis, fileName]);
                    genericPF.rateMap = rm;
                    genericPF.xBin = xBin;
                    genericPF.yBin = yBin;
                end
             if isempty(genericPF.trialSubType)
                 fileName = [genericPF.filebase, '.FindPFPars.', genericPF.trialName '.mat'];
             else
                 fileName = [genericPF.filebase, '.FindPFPars.', genericPF.trialName '.' genericPF.trialSubType '.mat'];
             end
             if FileExists([genericPF.paths.analysis, fileName])
                    load([genericPF.paths.analysis, fileName]);
             else                   
                    fprintf('\n computing place field parameters \n');
                    nClus = NClusters(genericPF);
                    pfPars = GenericPF.FindPFPars(genericPF, 1 : nClus);
                    save([genericPF.paths.analysis, fileName], 'pfPars')
             end
            end
            try
                genericPF.com = pfPars.com;
                genericPF.smoothRateMap = pfPars.smoothRateMap;
                genericPF.selectedPairs = pfPars.selectedPairs;
                genericPF.idealPFUnits = pfPars.idealPFUnits;
                genericPF.idealPFPairs = pfPars.idealPFPairs;
                genericPF.spatialCoherence = pfPars.spatialCoherence;
                genericPF.sparsity = pfPars.sparsity;
                genericPF.ratePk = pfPars.ratePk;
                genericPF.acceptedUnits = pfPars.acceptedUnits;
                genericPF.comDist = pfPars.comDist;
                genericPF.pkLoc = pfPars.pkLoc;
                genericPF.pkDist = pfPars.pkDist;
                if isfield(pfPars, 'entropy');
                    genericPF.jntEntropy = pfPars.entropy;
                end
               if isempty(genericPF.maze)
                   load([genericPF.paths.data, genericPF.filebase, '.miscPar.mat']);
                   genericPF.maze = miscPar.maze;
               end
            catch err
                
            end

        end % END of Class Constructor

        function genericPF = Convert2GenericPF(genericPF, anyPFObj)
                
                trClass = class(anyPFObj);
                switch trClass
                    case 'MTAPlaceField'
                        genericPF.datasetType = 'MTA';
                        % since MTAPlaceField has inconsistent terminology btw
                        % MTATrial MTAPlaceField, get the filebase
                        filebase = anyPFObj.filebase;
                        posOfDots = regexp(filebase,'\.');
                        filebase = filebase(1: posOfDots(1) -1);
                        genericPF.filebase = filebase;
                        genericPF.paths.data = ['~/data/', filebase, '/'];
                        if DirExists('~/data/nlx/', filebase)
                            genericPF.paths.data = ['~/data/nlx/', filebase, '/'];
                        else
                            error('\n %s filebase not found in either ~/data/ or ~/data/nlx directories \n', filebase);
                        end
                        genericPF.paths.analysis = ['~/data/analysis/', genericPF.filebase, '/'];
                        if ~DirExists('~/data/analysis', genericPF.filebase);
                            mkdir([genericPF.paths.analysis, filebase]);
                        end
                        % create a subfolder under analysis folder with  name
                        %  <filebase>, in case it does not exist
                        if ~DirExists('~/data/analysis', filebase);
                            mkdir([genericPF.paths.analysis, filebase]);
                        end
                        genericPF.trialName = anyPFObj.trialName;
                        genericPF.trialSubType = [];
                        genericPF.rateMap = anyPFObj.rateMap;
                        genericPF.xBin = anyPFObj.xbin;
                        genericPF.yBin = anyPFObj.ybin;
                        genericPF.nBins = anyPFObj.nbins;
                        if strcmp(anyPFObj.stateLabel, 'head.theta')
                            genericPF.state = 'walk'; 
                        end
                        genericPF.maze.name = anyPFObj.mazeName;

                    case 'MTATrial'
                        fprintf('\n loading rate maps ...');
                        mtaPFObj = LoadMTAPFObject(anyPFObj.name, anyPFObj.trialName);
                        genericPF = genericPF.Convert2GenericPF(mtaPFObj);
                end
        end
    end  
    
    methods(Static)
        %        [rateMap, occupancy, Bin1, Bin2, count, occCount] = ComputeRateMap(genericTrial, cluIdx, varargin);
        [rateMap, varargout] =  ComputeRateMap(genericTrial, cluIdx, varargin);
        pfPars = FindPFPars(arg, varargin);
    end
    
end

