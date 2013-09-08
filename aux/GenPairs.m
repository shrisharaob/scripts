function pairs = GenPairs(a, b)
% pairs = GenPairs(a, b)
% returns all unique pairings of vectors a and b 

    a = a(:)';
    b = b(:)';
    pairs = nchoosek([a, b], 2);
    maxA = max(a);
    maxB = max(b);
    pairs(pairs(:, 1) > maxA, :) = [];
    pairs(pairs(:, 2) > maxB, :) = [];
    pairs = unique(sortrows(pairs), 'rows');
end
    