function out = TemplateTimeCourse(gt, varargin)

    [nChunks, prePost, minCellsInSeq] = DefaultArgs(varargin, {5, 'pre', 5});
    
    templateMatch = gt.TemplateMatch(0, minCellsInSeq, prePost);
    evntperiods = gt.TrajectoryEvents(0, prePost);
    evntperiods = RecenterPeriods(evntperiods);
    bins = 1 : nChunks : max(evntperiods);
    
    [~, chunkId] = histc(find(), bins);
    
  