function y = RangeInInterval(rng, interval)
% y = RangeInInterval(rng, interval)

    nRanges = size(rng, 1);
    nIntervals = size(interval, 1);
    for ii = 1 : nRanges
        for jj = 1 : nIntervals
            IN_INTERVAL = Aux(rng(ii, :), interval(jj, :));
            if IN_INTERVAL, y(ii) = jj; end
        end
    end
end

function y = Aux(rng, interval)

    InInterval = @(x, interval) x >= interval(1) & x <= interval(2);
    y = InInterval(rng(1), interval) & InInterval(rng(2), interval);
end
    