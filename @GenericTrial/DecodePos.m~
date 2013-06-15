function [pos, predErr] = DecodePos(gt, varargin)
% [pos, predErr] = DecodePos(gt, varargin)
% [binSize, IF_REPORTFIG, type, state, ratemaps, cluId, binOverlap]
% {200e-3, 0, 'display', 'RUN', defRateMaps, defCluId, 0}
% decode position 
% gt - GenericTrial Object

    if isempty(gt.pfObject), gt = gt.LoadPF; end
    defCluId = gt.pfObject.acceptedUnits;
    defRateMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, defCluId));
    switch gt.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
   
    [binSize, IF_REPORTFIG, type, state, ratemaps, cluId, binOverlap, markerNo] = ...
        DefaultArgs(varargin, {200e-3, 0, 'display', 'RUN', defRateMaps, defCluId, 0, defMarker});

    cluId = cluId(ismember(cluId, gt.pfObject.acceptedUnits));
    ratemaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, cluId));
%     switch gt.datasetType
%       case 'MTA'
%         markerNo = 7;
%         for i = 1:length(gt.states)
%             if strcmp(state, 'RUN')
%                 state = 'walk';
%             end
%             if strcmp(gt.states{i}.label, state),
%                 statePeriods = gt.states{i}.statePeriods;
%                 break;
%             end
%         end
%         % convert 2 samplerate
%         stateXY = SelectPeriods(sq(gt.position(:, markerNo, :)), statePeriods, 'c');
%         statePeriods = ConvertFs(statePeriods, gt.trackingSampleRate, gt.sampleRate);
%       case 'kenji'
%         switch state
%           case 'RUN'
%             for i = 1:length(gt.states)   
%                 if strcmp(gt.states{i}.name, state)
%                     statePeriods = gt.states{i}.statePeriods; % @ lfp fs
%                     break;
%                 end
%             end
%             markerNo = 1;
%             %            stateXY = SelectPeriods(sq(gt.position(:, markerNo, :)), ConvertFs(statePeriods, gt.lfpSampleRate, gt.trackingSampleRate), 'c');
%             statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
%             statePeriods = gt.LoadStatePeriods(state, gt.sampleRate, cluId);
%           case 'SWS'
%             statePeriods = gt.TrajectoryEvents('SWS', cluId);
%         end

    statePeriods = gt.LoadStatePeriods(state, gt.sampleRate, [], 1);
    if isempty(gt.res)
        gt = gt.LoadCR;
    end
    fprintf('\n computing instantaneous firing rate... ')
    keyboard
    [sc, bc, xy] = GetSpikeCounts(gt, binSize, statePeriods, cluId, binOverlap);
    fprintf('\n compution posterior... ');
    options = struct('prior', [], 'bins', size(sc, 2), 'alpha', 5);
    tic, posterior = decode_bayesian_poisson(ratemaps, sc, options);toc
    pos = decodedPosMAP(gt, posterior);
    keyboard;
    % bcidx = round(bc * gt.trackingSampleRate)+1;
    % bcidx(end) = [];
%     bc = round(bc * gt.trackingSampleRate) + 1;
%     for kWin = 1 : size(bc, 1)
%         xy(kWin, :) = nanmean(SelectPeriods(stateXY , bc, 'c'));
%     end
    for kBin = 1 : size(pos,1) - 1
        predErr(kBin) = norm(xy(kBin,:) - pos(kBin,:));
    end
    h = figure;
    plot(predErr);
    title([num2str(binSize * 1e3), 'ms bins']);
    keyboard;
end