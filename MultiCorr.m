function r = MultiCorr(a, b, dim)

% Remove means 
az = bsxfun(@minus, a, mean(a,dim));
bz = bsxfun(@minus, b, mean(b,dim));
% Standard Pearson correlation coefficient formula
a2 = az .^ 2;
b2 = bz .^ 2;
ab = az .* bz;
r = sq(sum(ab, dim) ./ sqrt(sum(a2, dim) .* sum(b2, dim)));

end