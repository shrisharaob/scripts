function colors = MakeColorMap(n)

    colors = zeros(n, 3);
    for k = 1:n
        colors(k, :) = hsv2rgb((k-1)/(n), ...
                                1, ...
                                0.6 + 0.2*mod(k, 3));
    end
end