function [pos, predErr] = DecodePos(gt, ratemaps, cluId, varargin)
% decode position 
% gt - GenericTrial Object
    [binSize, IF_REPORTFIG, type, state] = DefaultArgs(varargin, {200e-3, 0, 'display', 'RUN'});

    switch gt.datasetType
      case 'MTA'
        markerNo = 7;
        for i = 1:length(gt.states)
            if strcmp(state, 'RUN')
                state = 'walk';
            end
            if strcmp(gt.states{i}.label, state),
                statePeriods = gt.states{i}.statePeriods;
                break;
            end
        end
        % convert 2 samplerate
        statePeriods = ConvertFs(statePeriods, gt.trackingSampleRate, gt.sampleRate);
      case 'kenji'
        for i = 1:length(gt.states)   
            if strcmp(gt.states{i}.name, state)
                statePeriods = gt.states{i}.statePeriods; % @ lfp fs
                break;
            end
        end
        markerNo = 1;
        statePeriods = ConvertFs(statePeriods, gt.lfpSampleRate, gt.sampleRate);
    end

    if isempty(gt.res)
        gt = gt.LoadCR;
    end
    [sc, bc] = GetSpikeCounts(gt, binSize, statePeriods, cluId, 0.6);
    fprintf('compution posterior... ');
    tic, posterior = decode_bayesian_poisson(ratemaps, sc);toc
    pos = decodedPosMAP(gt, posterior);
    % bcidx = round(bc * gt.trackingSampleRate)+1;
    % bcidx(end) = [];
    bc = round(bc * gt.trackingSampleRate) + 1;
    for kWin = 1 : size(bc, 1)
        xy(kWin, :) = nanmean(SelectPeriods(gt.position(:, markerNo , [1, 2]), bc, 'c'));
    end
    for kBin = 1 : size(pos,1) - 1
        predErr(kBin) = norm(xy(kBin,:) - pos(kBin,:));
    end
    h = figure;
    plot(predErr);
    title([num2str(binSize * 1e3), 'ms bins']);
end