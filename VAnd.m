function out = VAnd(a, varargin)
% perform and operation along dim

[dim] = DefaultArgs(varargin, {[]});
a = logical(a);
if isempty(dim)
    [~, dim] = min(size(a));
end
sizeA = size(a);
out = true(sizeA(dim), 1);
for kk = 1 : sizeA(dim)
    if dim == 1
        out = out & a(:, kk);
    else 
        out = out & a(kk, :)';
    end
end
end