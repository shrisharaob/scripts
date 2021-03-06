function [circlePixels, idx] = GetPxInRadius(a, x, y, radius)
    % returns the elements of the array at distance of radius 
    [rows, columns] = size(a);
    [columnsInImage rowsInImage] = meshgrid(1:columns, 1:rows);
    idx = (rowsInImage - x).^2 + (columnsInImage - y).^2 <= radius.^2; % elements enclosed by circle 
    circlePixels = a(idx);
end