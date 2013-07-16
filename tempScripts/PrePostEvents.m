function out = PrePostEvents(gt)
% visualize the time course of  template correlation strength during sleep pre, post trial

   prePost = {'pre', 'post'};
   IF_REPORT = false;
   for k = 1 : 2
       templateMatch = gt.TemplateMatch(0, 5, prePost{k});
       nEvents = eval(['length(templateMatch.IS_SIGNF_', upper(prePost{k}), ');']);
       nSignfEvnts = eval(['sum(templateMatch.IS_SIGNF_', upper(prePost{k}), ');']);
       binEdges = linspace(1, nEvents, ceil(nSignfEvnts / 2) + 1);
       signfEvnts = eval(['find(templateMatch.IS_SIGNF_', upper(prePost{k}), ');']);
       counts = histc(signfEvnts, binEdges);
       if nSignfEvnts > 0
           subplot(1, 2, k);
           bar(binEdges, counts);
           axis square;
           title(prePost{k});
           set(gca, 'XTick', []);
           ylim(gca, [0, max(counts) + 1]);
           set(gca, 'YTick', [0 : 1 :  max(counts) + 1]);
           IF_REPORT = true;
       end
   end
   if IF_REPORT, reportfig(gcf, mfilename , 0, [gt.filebase, '    ', gt.trialName]); end
end   