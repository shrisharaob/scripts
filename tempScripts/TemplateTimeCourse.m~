function out = TemplateTimeCourse(gt, varargin)

    [nChunks, prePost, minCellsInSeq, alpha] = DefaultArgs(varargin, {3, 'pre', 5, .05});
    
    tm = gt.TemplateMatch(0, minCellsInSeq, prePost);
    evntperiods = gt.TrajectoryEvents(0, prePost);
    evntperiods = RecenterPeriods(evntperiods);
    bins = 1 : round(size(evntperiods, 1 )./ nChunks) : size(evntperiods, 1);
    sbins = [1, bins(2:end) * 1e3];
    bins = [bins', circshift(bins', -1)]; bins(end, :) = [];
    sbins = [sbins', circshift(sbins', -1)]; sbins(end, :) = [];
    if ~isempty(tm.preEvntCorrs) & ~isempty(tm.preSurrogate)
        for kChunk = 1 : nChunks
            kIdx = linspace(bins(kChunk, 1), bins(kChunk, 2), diff(bins(kChunk, :)) + 1);
            kSurIdx = linspace(sbins(kChunk, 1), sbins(kChunk, 2), diff(sbins(kChunk, :)) + 1);
            [h(kChunk), p(kChunk), ks2stat(kChunk)] = kstest2(tm.preSurrogate(kSurIdx), tm.preEvntCorrs(kIdx), alpha);
        end
    end
    %[~, chunkId] = histc(find(tm.IS_SIGNF_PRE), bins);
keyboard;    
    
end
    
    
  