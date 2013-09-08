function y = Nan2Zero(a)
    % replace all nan's by yero

    y = a;
    y(isnan(y)) = 0;
end