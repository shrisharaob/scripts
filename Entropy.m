function bits = Entropy(a, varargin)
% information 

    occ = DefaultArgs(varargin, {ones(size(a))});
    bits = 0;
    a = a ./ sum(a(:));
    a(a==0) = 1;
    occ(isnan(occ)) = 0;
    bits = -1 * sum(occ(:) .* a(:) .* log2(a(:)));
end
