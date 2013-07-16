function PoolEvnts(prePost, roi, varargin)

[minCellInSeq, nResample] = DefaultArgs(varargin, {5, 1e3});


switch prePost
  case 'pre'
    signfEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'preSignfEvntCorr'));
    allEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'preEvntCorrs'));
    surgtCorr = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'preSurrogate'));
    preNCells = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'preNCells'));    be = linspace(-1, 1, 1e2);
    signfEvntCnt = histc(signfEvnts.poolArray, be);
    evntCnt = histc(allEvnts.poolArray, be);
    surCnt = histc(surgtCorr.poolArray, be);
    preNCells = preNCells.poolArray(preNCells.poolArray > minCellInSeq);
    preBinEdg = min(preNCells) - 1 : max(preNCells) + 1; % binEdges for # cells
    preCellsInEvntCnt = histc(preNCells, preBinEdg);
            
    hFig = figure;
    set(hFig, 'position', [1          28        1918         917]);
    hBar1 = bar(be, evntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
    grid on;
    hold on;
    hBar2 = bar(be-.01, surCnt ./ nResample, 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
    hBar3 = bar(be, signfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
    legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events', 'location', 'NorthWest');
    legend('boxoff');
    xlabel('Correlation value');
    ylabel('Number of events');
    ylim([0, max(ylim) + 50]);
    axHdl =  axes; %('position', [.7, .7, .2, .1]);
    bar(axHdl, preBinEdg, log10(preCellsInEvntCnt));
    xlim(axHdl, [min(xlim(axHdl)) + 0.75, max(xlim(axHdl))]);
    set(axHdl, 'XTicks', [minCellInSeq + 1 : round(max(xlim(axHdl)))]);
    grid(axHdl, 'on')
    set(axHdl, 'position', [.65, .8, .25, .1]);
    set(axHdl, 'Box', 'off');
    xlabel('# cells in events', 'FontSize', 5);
    ylabel('# events (log)', 'FontSize', 6);
    set(axHdl, 'FontSize', 6);
    keyboard;    
  case 'post'
    signfEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'postSignfEvntCorr'));
    allEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'postEvntCorrs'));
    surgtCorr = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'postSurrogate'));
    postNCells = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, minCellInSeq, prePost}, 'pool', 0, struct('poolVar', 'postNCells'));
    be = linspace(-1, 1, 1e2);
    signfEvntCnt = histc(signfEvnts.poolArray, be);
    evntCnt = histc(allEvnts.poolArray, be);
    surCnt = histc(surgtCorr.poolArray, be);
    postNCells = postNCells.poolArray(postNCells.poolArray > minCellInSeq); 
    postBinEdg = min(postNCells) : max(postNCells) + 1; % binEdges for # cells
    postCellsInEvntCnt = histc(postNCells, postBinEdg);  
    if ~isempty(allEvnts.poolArray)
        hFig = figure;
        set(hFig, 'position', [1          28        1918         917]);
        hBar1 = bar(be, evntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
        grid on;
        hold on;
        hBar2 = bar(be-.01, surCnt./ nResample , 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
        hBar3 = bar(be, signfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
        legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events', 'Location', 'NorthWest');
        legend('boxoff');
        xlabel('Correlation value');
        ylabel('Number of events');
        ylim([0, max(ylim) + max(ylim) * .25]);
        
        axHdl =  axes; %('position', [.7, .7, .2, .1]);
        bar(axHdl, postBinEdg, postCellsInEvntCnt);
        xlim(axHdl, [min(xlim(axHdl)) + 0.75, max(xlim(axHdl))]);
        set(axHdl, 'XTicks', [minCellInSeq + 1 : round(max(xlim(axHdl)))]);
        grid(axHdl, 'on')
        set(axHdl, 'position', [.65, .8, .25, .1]);
        set(axHdl, 'Box', 'off');
        xlabel(axHdl, '# cells in events', 'FontSize', 5);
        ylabel(axHdl, '# events (log)', 'FontSize', 6);
        set(axHdl, 'FontSize', 6);
    end
    keyboard;
end
