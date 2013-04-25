function entropy = Entropy(a)
% Shanon entropy

    entropy = 0;
    a = a ./ sum(a(:));
    a(a==0) = 1;
    entropy = -1 * sum(a(:) .* log2(a(:)));
end
