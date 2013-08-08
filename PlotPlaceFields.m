
    function PlotPlaceFields(pfObject, varargin)
    % PlotPlaceFields(pfObject, varargin)
    % varargin - [cellCluIdx, IS_COUNTOUR, IF_WAITFORBTNPRESS, nContours, contourColor, mazeDiameter]
    if nargin<1, help PlotPlaceFields; return; end
    if ~isa(pfObject,'GenericPF')
        if isa(pfObject, 'MTAPlaceField') || isa(pfObject, 'MTATrial')
            pfObject = GenericPF(pfObject);
        end
    end
    nTypes = length(pfObject);
    if nTypes >1
        pfObject = pfObject{end};
    end

    filebase = pfObject.filebase;
    %     posOfDots = regexp(filebase,'\.');
    % %     filebase = filebase(1: posOfDots(1) -1);
    if FileExists([pfObject.paths.analysis, pfObject.filebase, '.SelectCells.mat'])
       units  = load(['~/data/analysis/' filebase '/' filebase '.SelectCells.mat']);
       elClu = units.acceptedElClu;
       linearCluIdx = units.linearPyrCluIdx;
    else 
       linearCluIdx = 1 : length(pfObject.rateMap);
       elClu = nan(length(linearCluIdx), 2);
    end

    [cellCluIdx, IS_COUNTOUR, IF_WAITFORBTNPRESS, nContours, contourColor, mazeDiameter] = ...
        DefaultArgs(varargin, {pfObject.acceptedUnits, 0, 0, 1, [], 84});
    mazeDiameter = mazeDiameter * 10;
    nCells = length(cellCluIdx);
    if ~IS_COUNTOUR
        if(iscell(pfObject)) % when there are multiple types of pfObj
            for mCell = 1 : nCells
                linIdx = ismember(linearCluIdx,cellCluIdx(mCell));
                elCluStr = ['El#' num2str(elClu(linIdx, 1)) ' Clu#' num2str(elClu(linIdx,2))];
                for kType = 1: nTypes
                    pfObject{kType}.plot(cellCluIdx(mCell));
                    subplotfit(kType, nTypes);
                    title(elCluS);
                end
                if IF_WAITFORBTNPRESS, waitforbuttonpress; end
            end
        else
            for mCell = 1 : nCells
                linIdx = ismember(linearCluIdx,cellCluIdx(mCell));
                elCluStr = ['El#' num2str(elClu(linIdx, 1)) ' Clu#' num2str(elClu(linIdx,2))];
%                 pfObject.plot(linearPyrCluIdx(mCell)); %%% FIX INDEXING
%                 smoothedRateMap = pfObject.smoothRateMap(:,:,ismember(pfObject.acceptedUnits, cellCluIdx(mCell)));
                clf;
                subplot(1,2,1)
                try
                    smoothedRateMap = pfObject.smoothRateMap(:,:,ismember(pfObject.acceptedUnits, cellCluIdx(mCell)));
                    imagesc(pfObject.xBin,pfObject.yBin,sq(smoothedRateMap));
                    text(pfObject.xBin(1) + 30 ,pfObject.yBin(end) - 30, num2str(max(smoothedRateMap(:))),'color','w','FontWeight', 'bold');
                    set(gca,'YDir','normal');
                    hold on;
                    if strcmp(pfObject.maze.name, 'cof')
                        DrawCircle(0,0, mazeDiameter / 2,'w');
                    end
                catch err
                end
                subplot(1,2,2)
                rateMap = pfObject.rateMap{linIdx(mCell)};
                imagesc(pfObject.xBin,pfObject.yBin,rateMap); colormap('Hot');
                text(pfObject.xBin(1) + 30 ,pfObject.yBin(end) - 30, num2str(max(rateMap(:))),'color','w','FontWeight', 'bold');
                set(gca,'YDir','normal');
                hold on;
                if strcmp(pfObject.maze.name, 'cof')
                    DrawCircle(0,0, mazeDiameter / 2,'w');
                end
                title(num2str(mCell));
                    
                    
                %                 title(elCluStr);
                if IF_WAITFORBTNPRESS, waitforbuttonpress; end
            end
        end
    else
        %         if(iscell(pfObject))
        rateThreshPar = 0.2;
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
                    contour(pfObject.xBin,pfObject.yBin, smoothedRateMap,[maxRate, maxRate].*rateThreshPar, 'Color', contourColor, 'LineWidth', 2);
                else
                    contour(pfObject.xBin,pfObject.yBin, smoothedRateMap,linspace(.707 * maxRate, maxRate, nContours), 'Color', contourColor, 'LineWidth', 2);
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
        %         end
    end
    end