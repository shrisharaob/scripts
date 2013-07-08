function PoolRateMaps(gpf, varargin)
% plot all rate maps

    [IF_SRM, type] =  DefaultArgs(varargin, {1, 'heatmap'});
    filebase = gpf.filebase;
    if FileExists([gpf.paths.analysis, gpf.filebase, '.SelectCells.mat'])
        units  = load(['~/data/analysis/' filebase '/' filebase '.SelectCells.mat']);
        elClu = units.acceptedElClu;
        linearCluIdx = units.linearPyrCluIdx;
    else 
        linearCluIdx = 1 : length(gpf.rateMap);
        elClu = nan(length(linearCluIdx), 2);
    end
    %nCells = length(cellCluIdx);
    
    mapDims = [length(gpf.xBin), length(gpf.yBin)];
    nextSquare = @(x) (floor(sqrt(x)) + 1);
    %    nElements = nextSquare(nCells) .* mapDims;
    %    catMap = nan(nElements);
    figure;
    switch type
      case 'heatmap'
        if IF_SRM     
            nCells = length(gpf.acceptedUnits);
            axHdls = tight_subplot(nextSquare(nCells), nextSquare(nCells), .003);
            for kCell = 1 : nCells
                subplot(axHdls(kCell)), imagescnan(gpf.smoothRateMap(:, :, kCell));
                xAxis = xlim; % recenter to 0 
                yAxis = ylim; 
                set(gca,'YDir','normal');
                text(xAxis(end) - xAxis(end) * 2 / 10 ,yAxis(end) - yAxis(end) * 2 / 10, num2str(max(max(gpf.smoothRateMap(:, :, kCell))), '%.1f'),'color','k','FontWeight', 'bold', 'fontsize', 14);

            end
        else
            validCells = find(~cellfun(@isempty, gpf.rateMap));
            nCells = length(validCells);
            axHdls = tight_subplot(nextSquare(nCells), nextSquare(nCells), .003);
            counter = 0;
            for kCell = 1 : nCells
                subplot(axHdls(kCell)), imagescnan(gpf.rateMap{validCells(kCell)});
                xAxis = xlim; % recenter to 0 
                yAxis = ylim; 
                set(gca,'YDir','normal');
                text(xAxis(end) - xAxis(end) * 2 / 10 ,yAxis(end) - yAxis(end) * 2 / 10, num2str(max(max(gpf.rateMap{validCells(kCell)})), '%.1f'), 'color','k','FontWeight', 'bold', 'fontsize', 14);
            end
        end
    end
    ForAllSubplots('axis off');
end