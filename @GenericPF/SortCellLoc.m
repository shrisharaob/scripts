function sortedCellId = SortCellLoc(gt, varargin)


    rateMap1D = Compute1DRateMap(gt);

    [ clus2Select] = DefaultArgs(varargin, {1 : length(rateMap1D)});
    rateMap1D = rateMap1D(clus2Select)
    VALID_CELLS = ~cellfun(@isempty, rateMap1D);
    validCellId = find(VALID_CELLS);
    validMaps = rateMap1D(VALID_CELLS);
    rateMaps = cell2mat(rateMap1D);
    [~, maxrm] = max(rateMaps, [], 2);
    [~, sortedCellIdx] = sort(maxrm);
    sortedCellId =  validCellId(sortedCellIdx);

    rm =rateMaps;     
    rmm = rm ./ repmat(max(rm , [],2), 1, 100);
    figure;
    imagesc(rmm(sortedCellId, :))
end
