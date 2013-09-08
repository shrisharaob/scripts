function h = VPlot(x, varargin)

[marker, color] = DefaultArgs(varargin, {'.', 'k'});

h = figure;
plot(x(:,1), x(:, 2), marker, 'Color', color);

end