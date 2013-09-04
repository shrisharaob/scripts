function PlotRateMaps(pfObject, varargin)
    % PlotPlaceFields(pfObject, varargin)
    %
    % -------
    % Inputs:
    %     IF_CONTOUR          -  {0, 1, 2}, contour plot ?
    %     IF_WAITFORBTNPRESS  
    %     IF_Srmap            - plot smoothed rate maps ?
    %     nContours           - number of contours
    %     contourColor          
    %     mazeDiameter        - in cm, for circular maze
    %     cellCluIdx          - cluster ids to plot

    if nargin<1, help PlotPlaceFields; return; end
  
    filebase = pfObject.filebase;
    if FileExists([pfObject.paths.analysis, pfObject.filebase, '.SelectCells.mat'])
       units  = load(['~/data/analysis/' filebase '/' filebase '.SelectCells.mat']);
       elClu = units.acceptedElClu;
       linearCluIdx = units.linearPyrCluIdx;
    else 
       linearCluIdx = 1 : length(pfObject.rateMap);
       elClu = nan(length(linearCluIdx), 2);
    end

    [ IF_CONTOUR, IF_WAITFORBTNPRESS, IF_Srmap, nContours, contourColor, mazeDiameter, cellCluIdx] = ...
        DefaultArgs(varargin, { 0, 0, 0, 1, [], 84, pfObject.acceptedUnits});
    mazeDiameter = mazeDiameter * 10;
    nCells = length(cellCluIdx);
    if ~IF_CONTOUR
       
            for mCell = 1 : nCells
                linIdx = ismember(linearCluIdx,cellCluIdx(mCell));
                elCluStr = ['El#' num2str(elClu(linIdx, 1)) ' Clu#' num2str(elClu(linIdx,2))];
                if IF_Srmap
                    try
                        rateMap= pfObject.smoothRateMap(:,:,ismember(pfObject.acceptedUnits, cellCluIdx(mCell)));
                    catch err
                    end
                else
                    rateMap = pfObject.rateMap{linIdx};
                end
                if ~isempty(rateMap)
                    if isfield(pfObject.maze, 'px2CmFactor')
                        xAxis = pfObject.xBin * pfObject.maze.px2CmFactor(1);
                        yAxis = pfObject.yBin * pfObject.maze.px2CmFactor(2);
                        IN_CM = true;
                    else 
                        xAxis = pfObject.xBin;
                        yAxis = pfObject.yBin;
                        IN_CM = false; 
                    end
                    %imagescnan({xAxis, yAxis, rateMap});colormap('hot')
                    imagesc(xAxis, yAxis, rateMap); colormap('hot'); axis square;
                    xAxis = xlim; % recenter to 0 
                    yAxis = ylim; 
                    %                    colorbar;
                    set(gca,'YDir','normal');
                    text(xAxis(end) - xAxis(end) * 3 / 10 ,yAxis(end) - yAxis(end) * 2 / 10, num2str(max(rateMap(:)), '%.1f'),'color','w','FontWeight', 'bold', 'fontsize', 8);
                    hold on;
                    if IN_CM
                        xlabel('cm');
                        ylabel('cm');
                    else
                        xlabel('pixels');
                        ylabel('pixels');
                    end
                    if strcmp(pfObject.maze.name, 'cof')
                        centerX = 0; %range(xlim) / 2;
                        centerY = 0; %range(ylim) / 2;
                        DrawCircle(centerX, centerY, mazeDiameter / 2,'w');
                    end
                    title(['clu #', num2str(cellCluIdx(mCell))]);
                    % title(elCluStr);
                    if IF_WAITFORBTNPRESS, waitforbuttonpress; clf; end
                else
                    title(num2str(cellCluIdx(mCell)));
                end
            end
    else % IF_CONTOUR

        colors = GenColormap(nCells);
        IS_AUTO_COLOR =0;
        if isempty(contourColor), IS_AUTO_COLOR = 1; end
        for mCell = 1 : nCells
            if IS_AUTO_COLOR
                contourColor = colors(mCell,:);
            end
            linIdx = ismember(linearCluIdx,cellCluIdx(mCell));
            idx2 = ismember(pfObject.acceptedUnits, cellCluIdx(mCell)); % idx of selected subset of cells
            if sum(idx2) ~= 0
                smoothedRateMap = pfObject.smoothRateMap(:,:,idx2);
                maxRate = max(smoothedRateMap(:));
                rateThreshPar = 3 * std(smoothedRateMap(:));
                if nContours == 1
                    contour(pfObject.xBin,pfObject.yBin, smoothedRateMap,[1, 1].* rateThreshPar, 'Color', contourColor);
                    hold on;
                else
                    contour(pfObject.xBin,pfObject.yBin, smoothedRateMap,linspace(rateThreshPar  * maxRate, maxRate, nContours), 'Color', contourColor);
                end
                if IF_CONTOUR == 2, hold on; end
                if strcmp(pfObject.maze.name, 'cof')
                    DrawCircle(0,0, mazeDiameter / 2,'k');
                end
                if strcmp(pfObject.datasetType, 'MTA')
                    plot(pfObject.pkLoc(idx2, 1), pfObject.pkLoc(idx2, 2), '*k', 'MarkerSize', 2.5); % x an y are reversed in Justins code
                else
                    plot(pfObject.pkLoc(idx2, 1), pfObject.pkLoc(idx2, 2), '*k', 'MarkerSize', 2.5);
                end
                grid on;
                if IF_WAITFORBTNPRESS, waitforbuttonpress; clf;end
            else
                fprintf('\n *** the specified unit with cluster ID %d was discarded *** \n', cellCluIdx(mCell));
            end
        end
    end
end