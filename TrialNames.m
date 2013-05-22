function trialNames = TrialNames(arg, varargin)
% trialNames = TrialNames(gt)
% [datasetType, roi, arena]
% returns all the trial names in a filebase

    [datasetType, roi, arena] = DefaultArgs(varargin, {[], 'CA3', 'bigSquare'});

    trialNames = [];
    if isa(arg, 'GenericTrial')
        filebase = arg.filebase;
        datasetType = gt.datasetType;
    else
        filebase = arg;
    end

    switch datasetType
      case 'kenji'
        kenjiSearch.roi = roi;
        kenjiSearch.arena = arena;
        matches = SearchKenji(kenjiSearch);
        trialNames = matches(strcmp(matches(:,1), filebase), 2);
      case 'MTA'
 %        try
%             mtaTrial = MTATrial(filebase, [], 'all');
%         catch err
%             return;
%         end
%         trialNames = permute(mtaTrial.list_trialNames, [2, 1]);
        trialNames = MTATrialNames(filebase);
    end
