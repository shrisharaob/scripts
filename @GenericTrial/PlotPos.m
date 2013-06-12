function PlotPos(gt, varargin)
% PlotPos(gt, varargin)

    switch gt.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
    [markerNo, IF_IN_CM] = DefaultArgs(varargin, {defMarker, 0});
    if IF_IN_CM
        if ~isempty(gt.maze), if isfield(gt.maze, 'px2CmFactor'),
                plot(gt.position(:, markerNo, 1) .* gt.maze.px2CmFactor(1), gt.position(:, markerNo, 2) .* gt.maze.px2CmFactor(2));
                xlabel('cm');
                ylabel('cm');
            end; 
        else
            plot(gt.position(:, markerNo, 1), gt.position(:, markerNo, 2));
        end
    else
        plot(gt.position(:, markerNo, 1), gt.position(:, markerNo, 2));
    end
end

