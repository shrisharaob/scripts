function y = RowWiseShift(p, shiftVector)
% y = RowWiseShift(p, shiftVector)
% shift each row of the array individually 

    [nRows, nClmns] = size(p);
    for kRow = 1 : nRows
        y(kRow, :) = circshift(p(kRow, :), [0, shiftVector(kRow)]);
    end
end
