	function [ccgPeriod, ccgOffset, firstPkOffset] = FindCCGPars(ccgSmooth, ccgTimeAxis)
% finds the peaks of ccg 

	
	ccgPeriod = [];
    ccgOffset = [];

    [pks, locs] = findpeaks(ccgSmooth);
    [maxPk, maxLoc] = max(pks);
    ccgPkTimes = ccgTimeAxis(locs);
    [~, firstPkIndx]=min(abs(ccgPkTimes));
    firstPkOffset = ccgPkTimes(firstPkIndx);
    ccgPeriod = median(diff(ccgPkTimes));
    nNegPks = [];
    nPosPks = [];
    if ~isempty(ccgPkTimes) 
        nNegPks = sum(ccgPkTimes < ccgPkTimes(maxLoc));
        nPosPks = sum(ccgPkTimes > ccgPkTimes(maxLoc));
    end
    
    kcycleNo = [-nNegPks : nPosPks];

   [ccgPars, ~] = polyfit(kcycleNo, ccgPkTimes, 1);  
    ccgPeriod = ccgPars(1);
    ccgOffset = ccgPars(2);
    maxPkTime = ccgTimeAxis(locs(pks == max(pks)));
    if ~isempty(maxPkTime) && ~isempty(ccgOffset)
        if maxPkTime < 0 && ccgOffset >0
            ccgOffset = -ccgOffset;
        end
    end
        ccgOffset = mod(ccgOffset, ccgPeriod);

     end
