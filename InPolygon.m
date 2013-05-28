function INPOLYGON = InPolygon(a, b)

INPOLYGON = inpolygon(a(:, 1), a(:, 2), b(:,1), b(:, 2));
end