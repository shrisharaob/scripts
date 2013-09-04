function pos = LoadPos(gt, state)
% loads position data 

    switch gt.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
    defStatePeriods = gt.LoadStateperiods(state, gt.trackingSampleRate, 0) - ConvertFs(gt.trialPeriods(1), gt.lfpSampleRate, gt.trackingSampleRate); 
    [statePeriods, markerNo, IF_ONLY_VALID_POS] = DefaultArgs(varargin, {defStatePeriods, defMarker});
    
    pos = SelectPeriods(gt.position(:, markerNo, :), statePeriods, 'c');
    
end