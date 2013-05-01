classdef GenericTrial
    % class to structure general recording data
    %-------
    % History:
    %  Shrisha - Created
    
    properties (SetAccess = protected, GetAccess = public)
        % tag identifying the type of dataset
        datasetType;
    end
    
    
    properties
        
        % string identifying the filebase
        filebase;
        
        % data location on disk, all raw data files go here : ~/data/<filebase>
        % all analysis files go here : ~/data/analysis/<filebase>/
        % subfields paths.data & paths.analysis
        paths;
        
        % tag identifying the type of trial within a session
        trialName;
        
        % for kenji data set where several trials are concatinated into a
        % single filebase @lfp  fs
        trialPeriods;
        
        % tag to identify trial subtypes, if any
        trialSubType;
        
        % recording system sample rate
        sampleRate
        
        %
        lfpSampleRate;
        
        % lfp
        lfp;
        
        % statep
        states;
        
        % Maze parameters name, boundaries
        maze;
        
        % position tracking sample rate
        trackingSampleRate;
        
        % tracking data
        position;   % nSamples x nMarkers x  nSpatialDims
        
        % periods where the tracking is valid @tracking fs
        goodPosPeriods;
        
        % units classified as pyr
        pyrCluIdx;  %%%%%% result of SelectPyr cells stored here %%%
        
        % linear cluster indices for all electrodes
        clu;
        
        res;
        
        % each row of ElClu addressed by the linear Cluster index clu gives
        % the electrode# and Cluster# on that electrode
        elClu;
        
        %electrode position
        elPos;
        
        % cell array of placefields objects
        pfObject ;
        
        %
        ccg;
        
        
        
    end
    
    
    methods
        function genericTrial = GenericTrial(arg, varargin)
        % arg is the filebase name string
        % genericTrial = GenericTrial(arg.name);
            [trialName, datasetType, mazePar, trackingSampleRate, ...
             trialSubType, sampleRate, lfpSampleRate, state] = ...
                DefaultArgs(varargin, { [], 'default', [], [], [], [], [], {'RUN'}});
            
            if nargin < 1 
                return; % return empty object
            elseif isa(arg, 'GenericTrial')
                propertyNames = properties('GenericTrial');
                    for kProperty = 1 : length(propertyNames)
                        genericTrial.(propertyNames{kProperty}) = arg.(propertyNames{kProperty});
                    end
            elseif isa(arg, 'MTATrial')
                % convert MTATrial to GenericTrial class ;
                % add OR expressions here isa(arg, <new class name> to
                % extend to other class structures
                genericTrial = genericTrial.Convert2Generic(arg);
                return;
            else
                
                % creat the object
                genericTrial.datasetType = datasetType;
                fileBase = arg;
                genericTrial.filebase = fileBase;
                
                switch genericTrial.datasetType
                  case 'MTA'
                    genericTrial.paths.data = ['~/data/nlx/', fileBase, '/'];
                    genericTrial.paths.data = ['~/data/analysis/MTA/', fileBase, '/'];
                  case 'kenji'
                    genericTrial.paths.data = ['~/data/kenji/', fileBase, '/'];
                    genericTrial.paths.analysis = ['~/data/analysis/kenji/', fileBase, '/'];
                  otherwise
                    genericTrial.paths.data = ['~/data/', fileBase, '/'];
                    genericTrial.paths.analysis = ['~/data/analysis/', fileBase, '/'];
                    if ~DirExists('~/data/', fileBase)
                        % if filebase not found in ~/data/ search in ~/data/nlx/ folder
                        if DirExists('~/data/nlx/', fileBase)
                            genericTrial.paths.data = ['~/data/nlx/', fileBase, '/'];
                            genericTrial.paths.analysis = ['~/data/analysis/', fileBase, '/'];
                            genericTrial.datasetType = 'MTA';
                            if isempty(trialName)
                                mtaTrial = MTATrial(fileBase, [], 'all');
                                trialNames = mtaTrial.list_trialNames;
                                disp(trialNames);
                                trialName = input('trial name :' ,'s');
                                genericTrial = genericTrial.Convert2Generic(MTATrial(fileBase, [],trialName));
                            else
                                genericTrial = genericTrial.Convert2Generic(MTATrial(fileBase, [],trialName));
                            end
                        elseif DirExists('~/data/kenji/', fileBase)
                            genericTrial.datasetType = 'kenji';
                            genericTrial.paths.data = ['~/data/kenji/', fileBase, '/'];
                            genericTrial.paths.analysis = ['~/data/analysis/kenji/', fileBase, '/'];
                            genericTrial.trackingSampleRate = 39.0625;
                            % create a subfolder under analysis folder with  name
                            %  <filebase>, in case it does not exist
                            if ~DirExists('~/data/analysis/kenji/', fileBase);
                                mkdir([genericTrial.paths.analysis]);
                            end
                        elseif DirExists('~/data/', fileBase)
                            if ~DirExists('~/data/analysis/', fileBase);
                                mkdir(genericTrial.paths.analysis);
                            end
                        else
                            error('\n %s filebase not found in either ~/data/ or ~/data/nlx directories \n', fileBase);
                        end
                    end
                end
                
                par = LoadPar([genericTrial.paths.data, genericTrial.filebase]);
                genericTrial.trialName = trialName;
                genericTrial.trialSubType = trialSubType;
                genericTrial.sampleRate = par.SampleRate;
                genericTrial.lfpSampleRate = par.lfpSampleRate;
                genericTrial.lfp = {};
                if isempty(genericTrial.trackingSampleRate)
                    if FileExists([genericTrial.paths.data, genericTrial.filebase, '.miscPar.mat'])
                        load([genericTrial.paths.data, genericTrial.filebase, '.miscPar.mat']);
                        genericTrial.trackingSampleRate = miscPar.trackingSampleRate;
                    elseif ~isempty(trackingSampleRate)
                        genericTrial.trackingSampleRate = trackingSampleRate;
                    else
                        genericTrial.trackingSampleRate = input(' \n tracking sample rate : ');
                    end
                end
                miscPar.trackingSampleRate = genericTrial.trackingSampleRate;
                switch genericTrial.datasetType
                    % .whl contains the tracking data nSamples x (nClmns);
                    % nClmns = (x,y) x nMarkers = 2 x nMarkers
                    % position loaded as nSamples x nMarkers x  nSpatialDims
                  case 'default'
                    fprintf('\n loading position data \n');
                    xy = importdata([genericTrial.paths.data, genericTrial.filebase, '.whl']);
                  case 'kenji'
                    load('~/data/kenji/Beh_time_ind.mat'); % loads Beh

                    if isempty(genericTrial.trialName)
                        fprintf(['\n trialName not specified, choose from :\n \n']);
                        disp(Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),genericTrial.filebase) , [4, 5]));
                        trialName = input(['\n enter trial name  :'],'s');
                        trialNameStr = regexp(genericTrial.filebase, '\.', 'split');
                        genericTrial.trialName = [trialNameStr{1}, '.', trialName];
                        if ~FileExists(['~/data/kenji/whl/', genericTrial.filebase, '.eegTime'])
                            fprintf(['\n' repmat('*', 1, 50) '\n '])
                            fprintf(['\n ~/data/kenji/whl/',  genericTrial.filebase, '.eegTime not found !!! \n']);
                            fprintf(['\n' repmat('*', 1, 50) '\n '])
                            return;
                        end
                        %        return;
                    end
                    try
                        xy = importdata([genericTrial.paths.data, genericTrial.trialName '.whl']);
                        inValidIdx = xy(:,1) == -1;
                    catch err
                        fprintf(['\n' repmat('*', 1, 50) '\n '])
                        fprintf([genericTrial.paths.data, genericTrial.trialName '.whl not found']);
                        fprintf(['\n' repmat('*', 1, 50) '\n '])
                        return;
                    end
                end
            end
            genericTrial.clu = [];
            genericTrial.res = [];
            genericTrial.pfObject = {};
            genericTrial.ccg = {};
            switch genericTrial.datasetType
              case 'default'
                if FileExists([genericTrial.paths.data, genericTrial.filebase, '.miscPar.mat'])
                    load([genericTrial.paths.data, genericTrial.filebase, '.miscPar.mat']);
                    genericTrial.maze = miscPar.maze;
                elseif isempty(genericTrial.maze)
                    genericTrial.maze.name = input('maze name');
                    genericTrial.maze.dimsInCm = input('maze boundries: enter as [width, length]');
                    markerNo = 1;
                else
                    genericTrial.maze = mazePar;
                end
              case 'kenji'
                  [nRows, nclmns] = size(xy);
                  nMarkers = nclmns / 2;
                  xy = [xy(:, 1), xy(:,3), xy(:,2), xy(:,4)]; % the .whl format is inconsistent across dasetTypes
                  genericTrial = genericTrial.ProcessKenji;
                  pos = xy(:,1);
                  pos(~inValidIdx) = 1;
                  t = flipud(pos);
                  t = [-1* t(1);t];
                  pos = [ -1 * pos(1); pos];
                  posBeginIndx = find(diff(pos) == 2);
                  posEndIndx = find(flipud(diff(t) == 2));
                  goodPeriods = [posBeginIndx+1, posEndIndx-1];
                  genericTrial.goodPosPeriods = goodPeriods;
                  xy(inValidIdx,:) = nan;
                  genericTrial.position = reshape(xy, nRows, nMarkers, 2);
                  markerNo = 1;
                  genericTrial.maze.name = char(Beh(strcmp(Beh(:,4), genericTrial.trialName), 5));
                  switch genericTrial.maze.name
                      case 'bigSquare'
                          genericTrial.maze.dimsInCm = [180, 180];
                      case'midSquare'
                          genericTrial.maze.dimsInCm = [120, 120];
                      case 'linear'
                          genericTrial.maze.dimsInCm = [250, 0];
                      otherwise
                          genericTrial.maze.dimsInCm = nan;
                  end
                  xPxRange = [min(genericTrial.position(~inValidIdx, markerNo, 1)), max(genericTrial.position(~inValidIdx, markerNo, 1))];
                  yPxRange = [min(genericTrial.position(~inValidIdx, markerNo, 2)), max(genericTrial.position(~inValidIdx, markerNo, 2))];
                  genericTrial.maze.boundaries = [xPxRange; yPxRange];
                  miscPar.maze.px2CmFactor = genericTrial.maze.dimsInCm ./ [range(xPxRange), range(yPxRange)]; % recenter xy values to zero
                  genericTrial.maze.px2CmFactor = miscPar.maze.px2CmFactor;
                  for mState = 1 : length(state)
                      genericTrial.states{mState} =  GenericState(genericTrial, state{mState});
                  end
                case 'MTA'
                markerNo = 7;
            end
            miscPar.maze = genericTrial.maze;
            save([genericTrial.paths.data, genericTrial.filebase, '.miscPar.mat'], 'miscPar');
            if FileExists(['~/data/', genericTrial.datasetType, '/ElePosition.txt'])
                elPos = importdata(['~/data/', genericTrial.datasetType, '/ElePosition.txt']);
                roi = {'EC2','EC3', 'EC4', 'EC5', 'DG', 'CA1', 'CA3'};
                [~, nElClu] = NClusters(genericTrial);
                cluStartId = cumsum([1, nElClu]); % linear cluster id of the first non-noise cluster on each electrode cluStrtId -by- nElectrodes
                rowId = find(~cellfun(@isempty, regexp(elPos, genericTrial.filebase)));
                %                 cluStartId(nElClu == 0) = [];
                for iRegion = 1 : length(roi)  % look for shanks in regions of interest
                    if ~isempty(rowId)
                        rowCell = regexp(elPos{rowId}, '\s', 'split');
                        genericTrial.elPos(iRegion).shank = find(strcmp(rowCell(:,5:end),roi{iRegion}));
                        genericTrial.elPos(iRegion).region = roi{iRegion};                        
                        if length(genericTrial.elPos(iRegion).shank) > 1 && nElClu(iRegion) > 0 % roi on several shanks, cluster Ids in the shank
                            genericTrial.elPos(iRegion).clu = cluStartId(genericTrial.elPos(iRegion).shank(1)) : cluStartId(genericTrial.elPos(iRegion).shank) + sum(nElClu(genericTrial.elPos(iRegion).shank)) - 1;
                        elseif ~isempty(genericTrial.elPos(iRegion).shank)  && nElClu(iRegion) > 0
                            genericTrial.elPos(iRegion).clu = cluStartId(genericTrial.elPos(iRegion).shank) : cluStartId(genericTrial.elPos(iRegion).shank) + sum(nElClu(genericTrial.elPos(iRegion).shank)) - 1;
                        end
                    end
                end
            end
            if FileExists([genericTrial.paths.analysis, genericTrial.filebase, '.SelectCells.mat'])
                load([genericTrial.paths.analysis, genericTrial.filebase, '.SelectCells.mat']);
                genericTrial.pyrCluIdx = linearPyrCluIdx;
            else
                genericTrial.pyrCluIdx = 1 : NClusters(genericTrial);
            end
        end % END of class constructor 
    %%
    function genericTrial = Convert2Generic(genericTrial,anyTrialObj)
    % subroutine to convert data classes into GenericTrial class
    trClass = class(anyTrialObj);
    switch trClass
        case 'MTATrial'
            genericTrial.datasetType = 'MTA';
            genericTrial.filebase = anyTrialObj.name;
            genericTrial.paths.data = [anyTrialObj.path.nlx, genericTrial.filebase, '/'];
            genericTrial.paths.analysis = [anyTrialObj.path.analysis, genericTrial.filebase, '/'];
            if ~DirExists('~/data/analysis', genericTrial.filebase);
                mkdir([genericTrial.paths.analysis, fileBase]);
            end
            genericTrial.trialName = anyTrialObj.trialName;
            genericTrial.trialSubType = [];
            genericTrial.sampleRate = anyTrialObj.sampleRate;
            genericTrial.lfpSampleRate = anyTrialObj.lfpSampleRate;
            genericTrial.lfp = anyTrialObj.lfp;
            genericTrial.trackingSampleRate = anyTrialObj.xyzSampleRate;
            genericTrial.position = anyTrialObj.xyz;
            genericTrial.clu = anyTrialObj.clu;
            genericTrial.res = anyTrialObj.res;
            genericTrial.pfObject = anyTrialObj.Pfs;
            genericTrial.ccg = anyTrialObj.ccg;
            genericTrial.maze.name = anyTrialObj.Maze.name;
            genericTrial.maze.boundaries = anyTrialObj.Maze.boundaries;
            genericTrial.goodPosPeriods = anyTrialObj.xyzPeriods;
            if ~isempty(anyTrialObj.Bhv)
                for i = 1:length(anyTrialObj.Bhv.States),
                    genericTrial.states{i}.label = anyTrialObj.Bhv.States{i}.label;
                    genericTrial.states{i}.statePeriods = anyTrialObj.Bhv.States{i}.state;
                end
            end
            load([genericTrial.paths.analysis, genericTrial.filebase, '.SelectCells.mat']);
            genericTrial.pyrCluIdx = linearPyrCluIdx;
        case 'GenericTrial'
            genericTrial = anyTrialObj;
        otherwise
            error('not a valid trial class ');
    end
    end
    
    
    
end %END of methods
end % END of classdef







