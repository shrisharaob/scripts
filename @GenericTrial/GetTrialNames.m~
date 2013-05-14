function trialNames = GetTrialNames(gt, varargin)

    switch gt.datasetType
        case 'zhenya'
         fn = CatStruct(xdir(gt.paths.data));
         filenames = fn.names;
         %         idx = find(~cellfun(@isempty, cellfun(@regexp, n, repmat({'.pos$'}, 464, 1), 'uniformoutput', 0)));
         find(~cellfun(@isempty, cellfun(@regexp, n, repmat({'.GetPlaceFields2D.mat$'}, length(filenames), 1), 'uniformoutput', 0))
         trialNames = filenames(idx);
   end
end
