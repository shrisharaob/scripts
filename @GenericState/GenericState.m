classdef GenericState

    properties
        
        name;

        % statePeriods in original filebase samples 
        statePeriods;

    end

    methods
        function gs = GenericState(arg, varargin)
            [state] = DefaultArgs(varargin, {'RUN'});
            gs.name = state;
                switch arg.datasetType
                    case  'MTA'
                        mtaTrial = MTATrial(arg.filebase, [], arg.trialName);
                        if strcmp(state, 'RUN'),  state = 'walk'; else
                            error('\n no %s state in filebase %s', state, arg.filebase); end
                       gs.statePeriods = mtaTrial.statePeriods(state) + 1; %% state periods starts with 0
                        markerNo = 7;
                    case 'default'
                        gs.statePeriods = load([arg.paths.data, arg.filebase '.sts.', state]);
                    case  'kenji'
                        statePeriods = load([arg.paths.data, arg.filebase '.sts.', state]); % @lfp fs
                        gs.statePeriods = IntersectRanges(statePeriods, arg.trialPeriods);
                end
           end % END OF CLASS CONSTRUCTOR
           
           %%
%            function periods = GetStatePeriods(gs, state);
%                for i = 1:length(trial.states),
%                    if strcmp(trial.states{i}.label,state),
%                        periods = trial.states{i}.statePeriods;
%                        break;
%                    end
%                end
%            end
    end
end
