function entropy = Entropy(a)

    entropy = 0;
    a = a ./ sum(a(:));

    entropy = -1 * sum(a(:) .* log2(a(:)));
%     [nRows, nClmns] = size(a);
%     for kRow = 1 : nRows
%         for kClmn = 1 : nClmns
%             entropy = entropy + a(kRow, kClmn) * log2(a(kRow, kClmn));
%         end
%     end
end
