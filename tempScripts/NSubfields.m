function out = NSubfields(gt, roi, arena)
% out = NSubfields(gt, roi, arena)
% hist cnt of  number of sub fields 


mpd = gt.MultiPeakPFDistance(roi, arena);
out.counts = histc(cellfun(@length, mpd.cntrVertices), 1:10);

end