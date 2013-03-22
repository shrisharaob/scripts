classdef GenericTrial
    % class to structure general recording data
    % use GenericTrial(<filebase>) to add a new filebase / a list of
    % filebases
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
        
        % state
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
        
        % cell array of placefields objects
        pfObject ;
        
        %
        ccg;
        
        
        
    end
    
    
    methods
        function genericTrial = GenericTrial(arg, varargin)
            % arg is the filebase name string
            % genericTrial = GenericTrial(arg.name);
            [trialName, datasetType, mazePar, trackingSampleRate, trialSubType, sampleRate, lfpSampleRate] = DefaultArgs(varargin, { [], 'default', [], [], [], [], []});
            
            fileBase = arg;
            if isempty(arg)
                genericTrial.filebase = '';
                genericTrial.datasetType = '';
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
                                    genericTrial = Generic(MTATrial(fileBase, [], 'all'));
                                else
                                    genericTrial = GenericTrial.Convert2Generic(MTATrial(fileBase, [],trialName));
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
%                         files = dir(genericTrial.paths.data);
%                         nWhlFiles = 0;
%                         for kFile = 1 : length(files)
%                             if ~files(kFile).isdir
%                                 [~,~,ext] = fileparts(files(kFile).name);
%                                 if strcmp(ext, '.whl')
%                                     nWhlFiles = nWhlFiles + 1;
%                                     whlFileNames{nWhlFiles} = files(kFile).name;
%                                 end
%                             end
%                         end
                        if isempty(genericTrial.trialName)
                            fprintf('\n trialName not specified, choose from\n');
                            Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),genericTrial.filebase) , [4, 5, 7])
                            return;
                        else
                            try 
                                xy = importdata([genericTrial.paths.data, genericTrial.trialName '.whl']);
                            catch err
                                fprintf('\n *************************************************** \n ')
                                fprintf([genericTrial.paths.data, genericTrial.trialName '.whl not found']);
                                fprintf('\n *****************************************************\n');
                                return;
                            end
                        end
                end
%                 
%                 for kWhl = 1 : nWhlFiles
%                     
%                     fprintf('\n loading %s \n', [genericTrial.paths.data, whlFileNames{kWhl}]);
%                     sysOut(end) = []; % the last char is a RET
%                     xy = importdata(sysOut);
%                     [~, genericTrial.trialName{kWhl},~] = fileparts(sysOut);
%                 end

                [nRows, nClmns] = size(xy);
                nMarkers = nClmns / 2;
                switch genericTrial.datasetType
                    case 'kenji'
                        xy = [xy(:, 1), xy(:,3), xy(:,2), xy(:,4)]; % the .whl format is inconsistent across dasetTypes
                        genericTrial = genericTrial.ProcessKenji;
                        pos = xy(:,1);
                        pos(pos ~= -1) = 1;
                        t = flipud(pos);
                        t = [-1* t(1);t];
                        pos = [ -1 * pos(1); pos];
                        posBeginIndx = find(diff(pos) == 2);
                        posEndIndx = find(flipud(diff(t) == 2));
                        goodPeriods = [posBeginIndx, posEndIndx];
                        genericTrial.goodPosPeriods = goodPeriods;
                end
                genericTrial.position = reshape(xy, nRows, nMarkers, 2);
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
                case 'MTA'
                    markerNo = 7;
            end
            xPxRange = range(sq(genericTrial.position(:, markerNo, 1))) - 1; %  remove invalid points
            yPxRange = range(sq(genericTrial.position(:, markerNo, 2))) - 1;
            genericTrial.maze.boundaries = [0, xPxRange; 0, yPxRange];
            miscPar.maze.px2CmFactor = genericTrial.maze.dimsInCm ./ [xPxRange, yPxRange];
            genericTrial.maze.px2CmFactor = miscPar.maze.px2CmFactor;
%             if genericTrial
            miscPar.maze = genericTrial.maze;
            save([genericTrial.paths.data, genericTrial.filebase, '.miscPar.mat'], 'miscPar');
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
            for i = 1:length(anyTrialObj.Bhv.States),
                genericTrial.states{i}.label = anyTrialObj.Bhv.States{i}.label;
                genericTrial.states{i}.statePeriods = anyTrialObj.Bhv.States{i}.state;
            end
        case 'GenericTrial'
            genericTrial = anyTrialObj;
        otherwise
            error('not a valid trial class ');
    end
    end
    
    
    
end %END of methods


end % END of classdef







