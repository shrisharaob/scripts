    hFig = figure;
    for kPr = 1 : size(cellPair, 1)
        kPair = ismember(pfObject.selectedPairs, cellPair(kPr, :),'rows');
        ccgIdx  = ismember(out.cellPairs, cellPair(kPr, :), 'rows');
        if sum(kPair) & sum(ccgIdx) 
            subplot(2, 3, 3);
            PlotRateMaps(pfObject, 0, 0, 1, 0, [],[], pfObject.selectedPairs(kPair, 1)); 
            axis square;
            subplot(2, 3, 6);
            PlotRateMaps(pfObject, 0, 0, 1, 0, [],[], pfObject.selectedPairs(kPair,2)); 
            subplot(2, 3, 2);
            PlotRateMaps( pfObject, 1, 0, 1, 0, [],[], pfObject.selectedPairs(kPair, 1)); axis square;
            hold on;
            PlotRateMaps(pfObject, 1, 0, 1, 0, [],[], pfObject.selectedPairs(kPair, 2));
            subplot(2, 3, 1);
            bar(ccg.ccgTimeAx, ccg.ccg{ccgIdx}(:, 1, 1)); axis tight; axis square;
            subplot(2, 3, 4);
            bar(ccg.ccgTimeAx, ccg.ccg{ccgIdx}(:, 1, 2)); axis tight; axis square;
            subplot(2, 3, 5)
            bar(ccg.ccgTimeAx, ccg.ccg{ccgIdx}(:, 2, 2)); axis tight; axis square;
            axis tight;
            subplot(2, 3, 4);
            hold on;
            plot(ccg.smthTAx,ccg.smthCCG(:, ccgIdx),'g-');
            title(['offset: ' num2str(ccg.offset(ccgIdx))]);
            line([0,0], ylim,'Color', 'c','LineWidth', 2);
            line([ccg.firstPeak(ccgIdx), ccg.firstPeak(ccgIdx)], ylim, 'Color', 'm','LineWidth',1.5);
            line([ccg.offset(ccgIdx), ccg.offset(ccgIdx)], ylim, 'Color', 'r','LineWidth', 1.5);
            % subplot(2, 3, 5)
%             bar(ccg.BehTimeAx, ccg.OutBeh(:, 1, 2, ccgIdx), 'b');
%             hold on;
%             plot(ccg.BehSmoothTimeAx,ccg.BehSmooth(:, ccgIdx),'g-');
%             line([0,0], ylim,'Color', 'c','LineWidth', 2);
%             line([ccg.lowPassOffset(ccgIdx), ccg.lowPassOffset(ccgIdx)], ylim, 'Color', 'm','LineWidth',1.5);
%             title(['offset: ' num2str(ccg.lowPassOffset(ccgIdx))]);
%             axis tight;
        else
            fprintf('this pair is not accepted');
        end
        keyboard;
        clf;
    end
