function IS_EQUAL = IsEqual(a, b, varargin)
% IS_EQUAL = IsEqual(a, b, varargin)
% returns logical for  doe b lie between close interval [a-tol, a+tol] 
% extend for any general matrix

    [tolerence, type] = DefaultArgs(varargin, {.01, 'absolute'});

    switch type
        case 'relative'
            sd = a * tolerence;
            upperBound = a + sd;
            lowerBound = a - sd;
        case 'absolute'
            upperBound = a + tolerence;
            lowerBound = a - tolerence;
        otherwise
            IS_EQUAL = [];
            return;
    end
    b = repmat(b, length(a), 1);
    IS_EQUAL = (b >= lowerBound)  & (b <= upperBound);
end

