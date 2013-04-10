function regionClu = GetRegionClu(trial, roi, varargin)
% regionClu = GetRegionClu(trial, roi, varargin)    
% returns the cluster ids in the regions specified
    % roi - cell array containing the regions or char str 
% IF_CAT - concatinates all the clusters into a singel vector, else
% returns roi clus in different cells 
%[IF_CAT] = DefaultArgs(varargin, {0});


    if ~iscell(roi) == 1
        roi = cellstr(roi);
    end
    regionClu = cell(1, length(roi));
    region = CatStruct(trial.elPos,'region',1); % returns cell array regions
    region = region.region;
    for ii = 1 : length(roi)
        idx = find(~cellfun(@isempty,regexp(region, roi{ii})));
        if ~isempty(idx)
            regionClu{ii} = trial.elPos(idx).clu;
        end
    end   
%     if length(roi) == 1
%         regionClu = regionClu{1};
%     end
%keyboard;
%if IF_CAT
    %for kRoi = 1 : leength(roi)
    %regionClu = cell2mat(regionClu);
    %end
end
 