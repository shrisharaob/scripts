function IS_EQUAL = IsEqual(a, b, varargin)
% IS_EQUAL = IsEqual(a, b, varargin)
% [tolerence, IF_ALL, type] = [.01, 1, 'absolute'] 
% returns logical for  doe b lie between close interval [a-tol, a+tol] 
% extend for any general matrix

    [tolerence, IF_ALL, type] = DefaultArgs(varargin, {.01, 1, 'absolute'});

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
 
    if ~IF_ALL
        ktt =  double(IS_EQUAL);
        for ii = 1 : length(ktt) - 1
            if IS_EQUAL(ii)
                ktt(ii) = IS_EQUAL(ii+1) == 0;
            end
        end
        IS_EQUAL = ktt;
    end
end

