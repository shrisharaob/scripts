function trialNames = GetTrialNames(gt, varargin)
    % trialNames = GetTrialNames(gt, varargin)
    % returns trial namespresent in the filebase

    switch gt.datasetType
        case 'zhenya'
         fn = CatStruct(dir(gt.paths.data));
         filenames = fn.names;
         %         idx = find(~cellfun(@isempty, cellfun(@regexp, n, repmat({'.pos$'}, 464, 1), 'uniformoutput', 0)));
         %    idx = find(~cellfun(@isempty, cellfun(@regexp, n, repmat({'.GetPlaceFields2D.mat$'}, length(filenames), 1), 'uniformoutput', 0));
         trialNames = filenames(idx);

         case 'MTA'
           trNames = GenericTrial.MTATrialNames(gt.filebase, [], 'all');
           trialNames = trNames(find(~cellfun(@strcmp, trNames, repmat({'all'}, length(trNames), 1))));
   end
end
