function timesInPF = GetPFTrajectories(trial, pfPars, pfObject, varargin)
    % find time points when the rat entered the palce fied and when it exited 
    % pfBoundries : cell array - 1 x nPFBoundries

    [stateName, IF_COMPUTE, nSTD] = DefaultArgs(varargin, {'RUN', 1, 3});
    switch trial.datasetType
        case 'MTA'
          stateTimes = trial.Bhv.getState(stateName).state;
          stateTimes = stateTimes + 1;
          trajectories = [];
          markerNo = 7;
          nAcceptedUnits = sum( ~ pfPars.IS_DISCARD);
      case 'kenji'
          stateTimes = 
    end
    if ~IF_COMPUTE
        load(['~/data/analysis/' filebase '/' filebase '.' trial.trialName '.' mfilename '.mat']);
        return;
    end
    if nAcceptedUnits ~= 0
        pfBoundries = cell(1, nAcceptedUnits);
        rateThreshFactor = 0.707
        for mUnit = 1 : nAcceptedUnits
            smoothedRateMap = pfPars.smoothRateMap(:,:,mUnit);
            maxRate = max(smoothedRateMap(:));
            pfBoundries{mUnit} = contourc(pfObject.xbin, pfObject.ybin, smoothedRateMap, [maxRate, maxRate] .* rateThreshFactor);
        end
        for  i = 1 : size(stateTimes,1)
            temp = sq(trial.xyz([stateTimes(i,1) : stateTimes(i,2)], markerNo ,[1,2]));
            trajectories = [ trajectories ; temp];
        end
        allTrajs = sq(trial.xyz(:, markerNo, [1,2]));
        nPlaceCells = length(pfBoundries);
        timesInPF = cell(1, nPlaceCells);
        for kPlaceCell = 1 : nPlaceCells
            fprintf('\n place field %d of %d PFs \n', kPlaceCell, nPlaceCells);
            kPFBoundry = pfBoundries{kPlaceCell};
            % trajectories inside PF boundry
             [TRAJ_IN_PF, TRAJ_ON_PF] = inpolygon(trajectories(:, 1), trajectories(:, 2), kPFBoundry(1, 2:end), kPFBoundry(2, 2:end)) ;
             entryExitPoints = diff(TRAJ_IN_PF); % entry point = 1 , exit = -1
             entryPoint = find(entryExitPoints == 1);
             exitPoint = find(entryExitPoints == -1);
             % if trajectory begins inside a PF
            if TRAJ_IN_PF(1) == 1
                exitPoint(1) = [];
            elseif TRAJ_IN_PF(end) == 1 % if trajectories end inside a PF 
                entryPoint(end) = [];
            end
            entryIdx = find(ismember(allTrajs, trajectories(entryPoint, :), 'rows')); % index of all entry points, index for trial.xyz
            exitIdx = find(ismember(allTrajs, trajectories(exitPoint, :), 'rows'));
            timesInPF{kPlaceCell} = [entryIdx(:), exitIdx(:)];
        end
        filebase = trial.name;
        save(['~/data/analysis/' filebase '/' filebase '.' trial.trialName '.' mfilename '.mat'], 'timesInPF');
end



