function inOutIdx = InOut(y)
     % inoutIdx = InOut(y)
     % returns the in and out indices of a binary time series
     inOutIdx = [];
     if length(y) < 2, return; end
     a = unique(y);
     if length(a) > 2, fprintf('\n input is not a binary time series \n'); return; end
     if size(y, 2) == length(y), y = y';end % clmn vector
     y(y == a(1)) = 0;
     y(y == a(2)) = 1;
   
     inOutIdx{1} = FindIdxInOut(y);
     inOutIdx{2}= FindIdxInOut(~y);
end
 
function idx = FindIdxInOut(y)
    yf = flipud(y);
    y = [~y(1); y];
    yf = [~yf(1); yf];
    inIdx = find(diff(y) == 1);
    outIdx = find(flipud(diff(yf) == 1));
    if length(inIdx) == length(outIdx)
        eqIdx = bsxfun(@eq, inIdx, outIdx);
        inIdx(eqIdx) = [];
        outIdx(eqIdx) = [];
    elseif length(inIdx) > length(outIdx)
        lenDiff = length(inIdx) - length(outIdx);
        gtIdx = bsxFun(@gt, inIdx, [outIdx; zeros(lenDiff, 1)]);
        inIdx(gtIdx) = [];
    elseif length(inIdx) < length(outIdx)
        lenDiff = -1 * length(inIdx) + length(outIdx);
        ltIdx = bsxFun(@lt, outIdx, [inIdx; zeros(lenDiff, 1)]);
        outIdx(ltIdx) = [];
    end
    idx= [inIdx, outIdx];
end