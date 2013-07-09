function out =  PoolRateMaps(gt, varargin)
% plot all rate maps 
    if isempty(gt.pfObject), gt.LoadPF; end
    gpf = gt.pfObject;
    defClus = gpf.acceptedUnits(gt.pfObject.sparsity < .4);
    [IF_SRM, clus, type, IF_REPORTFIG] =  DefaultArgs(varargin, {1, defClus, 'heatmap', 0});
    out = [];


    filebase = gpf.filebase;
    if FileExists([gpf.paths.analysis, gpf.filebase, '.SelectCells.mat'])
        units  = load(['~/data/analysis/' filebase '/' filebase '.SelectCells.mat']);
        elClu = units.acceptedElClu;
        linearCluIdx = units.linearPyrCluIdx;
    else 
        linearCluIdx = 1 : length(gpf.rateMap);
        elClu = nan(length(linearCluIdx), 2);
    end
    mapDims = [length(gpf.xBin), length(gpf.yBin)];
    nextSquare = @(x) (floor(sqrt(x)) + 1);
    set(gcf, 'position', [1          25        1680         920]);
    switch type
      case 'heatmap'
        if IF_SRM     
            nCells = length(clus);
            axHdls = tight_subplot(nextSquare(nCells), nextSquare(nCells), .003);
            for kCell = 1 : nCells
                subplot(axHdls(kCell)), imagescnan(gpf.smoothRateMap(:, :, ismember(gpf.acceptedUnits, clus(kCell))));
                xAxis = xlim; % recenter to 0 
                yAxis = ylim; 
                set(gca,'YDir','normal');
                text(xAxis(end) - xAxis(end) * 2 / 10 ,yAxis(end) - yAxis(end) * 2 / 10, num2str(max(max(gpf.smoothRateMap(:, :, kCell))), '%.1f'),'color','k','FontWeight', 'bold', 'fontsize', 10);
                axis off;
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
                text(xAxis(end) - xAxis(end) * 2 / 10 ,yAxis(end) - yAxis(end) * 2 / 10, num2str(max(max(gpf.rateMap{validCells(kCell)})), '%.1f'), 'color','k','FontWeight', 'bold', 'fontsize', 10);
            end
        end
    case 'contours'
      % pfff..
    end
    
    for ii = nCells : length(axHdls)
        set(axHdls(ii), 'visible', 'off');
    end

    if IF_REPORTFIG
        [arena, roi] = SearchKenji(gt.trialName);
        arena = arena{2};
        if any(cellfun(@strcmp, roi, repmat({'CA3'}, size(roi))))
            roi = 'CA3';
        else
            roi = 'CA1';
        end
        reportfig(gcf, mfilename, 0, [filebase, '    trialName :   ', gt.pfObject.trialName, '    arena : ' arena, '   roi : ' roi], 300, 0);
        clf;
    end
end