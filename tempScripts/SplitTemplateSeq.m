function out = SplitTemplateSeq(gt, varargin)
% out = SplitTemplateSeq(gt, varargin)
% script to separate template match according to number of cells in the sequence

    out = [];
    %    [prePost] = DefaultArgs(varargin, {'pre'});
    templateMatch = gt.TemplateMatch;
    for prePost = {'pre', 'post'}
        NCells = eval(['unique(templateMatch.' char(prePost), 'NCells);']);
        NCells(NCells == 0) = [];
        counter = 0;
        for kNCells = 1 : length(NCells)
            kCellSeqCorr = eval(['templateMatch.', char(prePost), 'EvntCorrs(templateMatch.', char(prePost), 'NCells == NCells(kNCells));']);
            binEdg = linspace(-1, 1, 1e2);
            kCount = histc(kCellSeqCorr, binEdg);
            counter = counter + 1;
            subplot(length(NCells), 1, kNCells);
            bar(binEdg, kCount);
            title(['# cells : ' num2str(NCells(kNCells))]);
            grid on;
            eval([char(prePost), 'SplitCorrs(counter, :) = kCount;']);
        end
        reportfig(gcf, [mfilename, '.', char(prePost)] , 0, [gt.filebase, '    ', gt.trialName]);
    end
   
end