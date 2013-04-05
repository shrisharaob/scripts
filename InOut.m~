function inOutIdx = InOut(y, varargin)
     % inoutIdx = InOut(y)
     % returns the in and out indices of a binary time series
    [alphabets] = DefaultArgs(varargin, {[-1, 1]});
     inOutIdx = [];
     if length(y) < 2, return; end
     a = unique(y);
     if length(a) > 2 , fprintf('\n input is not a binary time series \n'); return; end
     if length(a) == 1  % if the whole series contains only one state
         if a == alphabets(1)
             inOutIdx{1} = [1, length(y)];
             inOutIdx{2}= [];
         else
             inOutIdx{2} = [1, length(y)];
             inOutIdx{1}= [];
         end
         return;
     end
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
    idx= ResolveLenDiff(inIdx, outIdx);
end

function out = ResolveLenDiff(a, b)
    if length(a) == length(b)
        eqIdx = bsxfun(@eq, a, b);
        a(eqIdx) = [];
        b(eqIdx) = [];
    else
        if length(a) > length(b)
            lenDiff = length(a) - length(b);
            b = [b; zeros(lenDiff, 1)];
            gtIdx = bsxfun(@gt, a, b);
            ltIdx = bsxfun(@lt, b, a);
            ltIdx(end: -1 : end-lenDiff) = 0;
        elseif length(a) < length(b)
            lenDiff = -1 * length(a) + length(b);
            a = [a; zeros(lenDiff, 1)];
            ltIdx = bsxfun(@lt, b, a);
            gtIdx = bsxfun(@gt, a, b);
            ltIdx(end: -1 : end-lenDiff) = 0;
        end
        inValidIdx = (gtIdx | ltIdx);
        a(inValidIdx) = [];
        b(inValidIdx) = [];
    end
    out = [a, b];
end