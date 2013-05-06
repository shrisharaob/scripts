function [binnedPos, mask] = Coverage(gt, varargin)
% [selectedPos, mask] = Covarage(gt, varargin)
% select postion samples that lie in the inside the 
% locations covered by recorded place fields

    [roi, arena, markerNo,  IF_PLOT, IF_REPORTFIG, nStd] = DefaultArgs(varargin, {'CA3', 'bigSquare',1, 0, 0, 3});
    gt = gt.LoadPF;
    load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
    srm = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
    nCells = size(srm, 3);
    colors = GenColormap(nCells);
    mask = false(size(srm, 1), size(srm, 2));
    for n = 1 : nCells
        nsrm = srm(:, :,n);
        sd3 = nStd * std(nsrm(:));
        if IF_PLOT
            contour(nsrm, [sd3, sd3], 'Color', colors(n, :));
        end
        mask = mask | nsrm > sd3;
        hold on;
    end
    pos = sq(gt.position(:, markerNo, :));
    binnedPos = BinPos(gt);
    invalidPosIdx = ismember(binnedPos, Ind2Sub([50, 50], find(~mask')), 'rows'); % transpose of mask
    binnedPos(invalidPosIdx, :) = nan;
%     pos(invalidPosIdx, :) = nan;
%     selectedPos = pos;
close(gcf);
    if IF_REPORTFIG  & IF_PLOT
        reportfig(gcf, mfilename, 0, [gt.trialName, '# units :', num2str(sum(ismember(gt.pfObject.acceptedUnits, commonClus)))], [], 0);
        clf;
    end
end
