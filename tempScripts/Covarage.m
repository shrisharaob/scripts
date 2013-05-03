function Covarage(gt, varargin)
    [roi, arena] = DefaultArgs(varargin, {'CA3', 'bigSquare'});
    gt = gt.LoadPF;
    load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
    srm = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, commonClus));
    nCells = size(srm, 3);
    colors = GenColormap(nCells);
    for n = 1 : nCells
        nsrm = srm(:, :,n);
        sd3 = 3 * std(nsrm(:));
        contour(nsrm, [sd3, sd3], 'Color', colors(n, :));
        hold on;
    end
    reportfig(gcf, mfilename, 0, [gt.trialName, '# units :', num2str(sum(ismember(gt.pfObject.acceptedUnits, commonClus)))], [], 0);
    clf;
end
