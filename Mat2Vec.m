function y = Mat2Vec(x, varargin)
    IF_CLM = DefaultArgs(varargin, {1});
    y = x(:);
    if ~IF_CLM
        y = y';
    end
end