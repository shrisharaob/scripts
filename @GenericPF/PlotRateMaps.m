
function PlotRateMaps(pfObject, varargin)
    % PlotPlaceFields(pfObject, varargin)
    % varargin - [cellCluIdx, IS_COUNTOUR, IF_WAITFORBTNPRESS, nContours, contourColor, mazeDiameter]
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

    [ IS_COUNTOUR, IF_WAITFORBTNPRESS, IF_Srmap, nContours, contourColor, mazeDiameter, cellCluIdx] = DefaultArgs(varargin, { 0, 0, 0, 5, [], 84, pfObject.acceptedUnits});
    mazeDiameter = mazeDiameter * 10;
    nCells = length(cellCluIdx);
    if ~IS_COUNTOUR
       
            for mCell = 1 : nCells
                linIdx = ismember(linearCluIdx,cellCluIdx(mCell));
                elCluStr = ['El#' num2str(elClu(linIdx, 1)) ' Clu#' num2str(elClu(linIdx,2))];
%                 pfObject.plot(linearPyrCluIdx(mCell)); %%% FIX INDEXING
%                 smoothedRateMap = pfObject.smoothRateMap(:,:,ismember(pfObject.acceptedUnits, cellCluIdx(mCell)));
                clf;
                if IF_Srmap
                    try
                        rateMap= pfObject.smoothRateMap(:,:,ismember(pfObject.acceptedUnits, cellCluIdx(mCell)));
                    catch err
                    end
                else
                    rateMap = pfObject.rateMap{linIdx};
                end
                
%                 if strcmp(pfObject.maze.name, 'linear')
                    
                if ~isempty(rateMap)
                    if ~isempty(pfObject.maze.px2CmFactor)
                        xAxis = pfObject.xBin * pfObject.maze.px2CmFactor(1);
                        yAxis = pfObject.yBin * pfObject.maze.px2CmFactor(2);
                        IN_CM = true;
                    else 
                        xAxis = pfObject.xBin;
                        yAxis = pfObject.yBin;
                        IN_CM = flase; 
                    end
                    xAxis = xAxis - xAxis(1); % recenter to 0 
                    yAxis = yAxis - yAxis(1); 
                    imagescnan({xAxis, yAxis, rateMap});
                    text(xAxis(1) + 30 ,yAxis(end) - 30, num2str(max(rateMap(:))),'color','w','FontWeight', 'bold');
                    set(gca,'YDir','normal');
                    hold on;
                    if IN_CM
                        xlabel('cm');
                        ylabel('cm');
                    else
                        xlabel('pixels');
                        ylabel('pixels');
                    end
                    if strcmp(pfObject.maze.name, 'cof')
                        DrawCircle(0,0, mazeDiameter / 2,'w');
                    end
                    title(['clu #', num2str(cellCluIdx(mCell))]);
                    %                                 title(elCluStr);
                    if IF_WAITFORBTNPRESS, waitforbuttonpress; end
                else
                    title(num2str(cellCluIdx(mCell)));
                end
            end
    else
        rateThreshPar = 0.707;
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
                if nContours == 1
                    contour(pfObject.xbin,pfObject.ybin, smoothedRateMap,[maxRate, maxRate].*rateThreshPar, 'Color', contourColor, 'LineWidth', 2);
                else
                    contour(pfObject.xbin,pfObject.ybin, smoothedRateMap,linspace(.707 * maxRate, maxRate, nContours), 'Color', contourColor, 'LineWidth', 2);
                end
                hold on;
                if strcmp(pfObject.mazeName, 'cof')
                    DrawCircle(0,0, mazeDiameter / 2,'k');
                end
                plot(pfObject.com(idx2, 2), pfObject.com(idx2, 1), '*k'); % x an y are reversed in Justins code
                grid on
                if IF_WAITFORBTNPRESS, waitforbuttonpress; end
            else
                fprintf('\n *** the specified unit with cluster ID %d was discarded *** \n', cellCluIdx(mCell));
            end
        end
    end
    end