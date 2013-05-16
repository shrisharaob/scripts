function trialNames = TrialNames(arg, varargin)
  %  trialNames = TrialNames(gt)
   % returns all the trial names in a filebase

  [dataSetType, roi, arena] = DefaultArgs(varargin, {[], 'CA3', 'bigSquare'});

trialName = [];
   if isa(arg, 'GenericTrial')
   filebase = arg.filebase;
 dataSetType = gt.datasetType;
   else
     filebase = arg;
 end

   switch dataSetTyep
   case 'kenji'

    kenjiSearch.roi = roi;
    kenjiSearch.arena = arena;
    matches = SearchKenji(kenjiSearch);
    trialNames = matches(strcmp(matches(:,1), filebase), 2);
case 'MTA'
   try
        mtaTrial = MTATrial(filebase, [], 'all');
    catch err
        return;
    end
    trialnames = permute(mtaTrial.list_trialNames, [2, 1]);
 
end
