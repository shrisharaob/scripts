function sortedCellId = SortCellLoc(gt, varargin)


    rateMap1D = Compute1DRateMap(gt, 0);

    [ clus2Select] = DefaultArgs(varargin, {1 : length(rateMap1D)});
    rateMap1D = rateMap1D(clus2Select)
    VALID_CELLS = ~cellfun(@isempty, rateMap1D);
    validCellId = find(VALID_CELLS);
    validMaps = rateMap1D(VALID_CELLS);
    rateMaps = cell2mat(rateMap1D);
    [~, maxrm] = max(rateMaps, [], 2);
    [~, sortedCellIdx] = sort(maxrm);
    cluIdx =  validCellId(sortedCellIdx);
    sortedCellId = clus2Select(cluIdx);
    rm =rateMaps;     
    rmm = rm ./ repmat(max(rm , [],2), 1, 100);
    figure;
    imagesc(rmm(cluIdx, :))
    set(gca, 'YTick', [1 : length(clus2Select)] );
    set(gca, 'YTickLabel', sortedCellId);
    ylabel('clu id', 'FontSize', 16);
    set(gca, 'FontSize', 10)
end
