function periods = RecenterPeriods(inPeriods)
% periods = RecenterPeriods(inPeriods)

    if size(inPeriods, 1) > 1
        periods = inPeriods - inPeriods(1,1) + 1;
    else
        periods = [1, diff(inPeriods) + 1];
    end
 end
