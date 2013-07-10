function PoolEvnts(prePost, roi)

nResample = 1e3;

switch prePost
  case 'pre'
    signfEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'preSignfEvntCorr'));
    allEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'preEvntCorrs'));
    surgtCorr = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'preSurrogate'));
    preNcells = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'preNCells'));    be = linspace(-1, 1, 1e2);
    signfEvntCnt = histc(signfEvnts.poolArray, be);
    evntCnt = histc(allEvnts.poolArray, be);
    surCnt = histc(surgtCorr.poolArray, be);
    preBinEdg = min(preNCells) : max(preNCells) + 1; % binEdges for # cells
    preCellsInEvntCnt = histc(preNCells, preBinEdg);  
            
    hFig = figure;
    set(hFig, 'position', [1          28        1918         917]);
    hBar1 = bar(be, evntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
    grid on;
    hold on;
    hBar2 = bar(be-.01, surCnt , 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
    hBar3 = bar(be, signfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
    legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events');
    legend('boxoff');
    xlabel('Correlation value');
    ylabel('Number of events');
    
    axHdl =  axes; %('position', [.7, .7, .2, .1]);
    bar(axHdl, preBinEdg, preCellsInEvntCnt);
    grid on;
    set(axHdl, 'position', [.7, .8, .2, .1])
    set(axHdl, 'Box', 'off');
    xlabel('# cells in events', 'FontSize', 6);
    ylabel('# events', 'FontSize', 6);
    set(axHdl, 'FontSize', 6);
    
case 'post'
    signfEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'postSignfEvntCorr'));
    allEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'postEvntCorrs'));
    surgtCorr = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'postSurrogate'));
    preNcells = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {0, 0, prePost}, 'pool', 0, struct('poolVar', 'postNCells'));
    be = linspace(-1, 1, 1e2);
    signfEvntCnt = histc(signfEvnts.poolArray, be);
    evntCnt = histc(allEvnts.poolArray, be);
    surCnt = histc(surgtCorr.poolArray, be);

    postBinEdg = min(out.postNCells) : max(out.postNCells) + 1; % binEdges for # cells
    postCellsInEvntCnt = histc(out.postNCells, postBinEdg);  
          
    hFig = figure;
    set(hFig, 'position', [1          28        1918         917]);
    hBar1 = bar(be, evntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
    grid on;
    hold on;
    hBar2 = bar(be-.01, surCnt./ nResample , 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
    hBar3 = bar(be, signfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
    legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events');
    legend('boxoff');
    xlabel('Correlation value');
    ylabel('Number of events');


    axHdl =  axes; %('position', [.7, .7, .2, .1]);
    bar(axHdl, postBinEdg, postCellsInEvntCnt);
    grid on;
    set(axHdl, 'position', [.7, .8, .2, .1])
    set(axHdl, 'Box', 'off');
    xlabel('# cells in events', 'FontSize', 6);
    ylabel('# events', 'FontSize', 6);
    set(axHdl, 'FontSize', 6);
    
end
