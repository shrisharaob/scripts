function INPOLYGON = InPolygon(a, vertices)
% INPOLYGON = InPolygon(a, vertices)

INPOLYGON = inpolygon(a(:, 1), a(:, 2), vertices(:,1), vertices(:, 2));
end