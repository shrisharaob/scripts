function idx = CompareVectors(a, b)
% returns logical of length a for every 
% 
    idx = zeros(length(a), 1);
    if size(a,2) == length(a), a = a'; end
    if size(b,2) == length(b), b = b'; end
    for i = 1 : length(b)
        idx = idx | a == b(i);
    end
end
