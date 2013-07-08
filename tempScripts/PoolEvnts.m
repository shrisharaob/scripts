



preOrPost = 'post';
nResample = 1e3;
roi = 'CA3';
switch preOrPost
  case 'pre'
    signfEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {}, 'pool', 0, struct('poolVar', 'preSignfEvntCorr'));
    allEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {}, 'pool', 0, struct('poolVar', 'preEvntCorrs'));
    surgtCorr = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {}, 'pool', 0, struct('poolVar', 'preSurrogate'));
    be = linspace(-1, 1, 1e2);
    signfEvntCnt = histc(signfEvnts.poolArray, be);
    evntCnt = histc(allEvnts.poolArray, be);
    surCnt = histc(surgtCorr.poolArray, be);

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

case 'post'
    signfEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {}, 'pool', 0, struct('poolVar', 'postSignfEvntCorr'));
    allEvnts = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {}, 'pool', 0, struct('poolVar', 'postEvntCorrs'));
    surgtCorr = BatchProcess(@TemplateMatch, 'kenji', roi, 'linear', 1, {}, 'pool', 0, struct('poolVar', 'postSurrogate'));
    be = linspace(-1, 1, 1e2);
    signfEvntCnt = histc(signfEvnts.poolArray, be);
    evntCnt = histc(allEvnts.poolArray, be);
    surCnt = histc(surgtCorr.poolArray, be);

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
end

% bar(be, allEvntCount, 'facecolor', 'w', 'edgecolor', 'k');
% hold on;
% bar(be, signfEvntCount, 'facecolor', 'r', 'edgecolor', 'none');
% legend( 'All Events', 'Signf Events');
% xlabel('Correlation value');
% ylabel('Number of events');
% title('Preplay');