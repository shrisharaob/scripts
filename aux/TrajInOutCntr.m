function inoutTime = TrajInOutCntr(pos, cntr)
% inoutTime = TrajInOutCntr(gt, cntr)
% pos nSamples-by-2
% cntr nSamples-by-w, contour vertices
% returns the times when the trajectories lie inside the specified contour

    y = InPolygon(pos, cntr);
    inout = InOut(y);

    inoutTime = inout{1} ;
end
    