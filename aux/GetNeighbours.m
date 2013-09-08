function neighbours = GetNeighbours(A, x, y, varargin)
% neighbours = GetNeighbours(A, x, y)
    [radius] = DefaultArgs(varargin, {1});
    [nRows, nClmns] = size(A);
    displacement = [ 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1];
    neighboursIdx = displacement  + repmat([x, y], [8, 1]);
    sizMat = repmat([nRows, nClmns], [8, 1]);
    %now remove indices out of bounds
    IF_REMOVE = logical(sum(neighboursIdx > sizMat | neighboursIdx < 1, 2));
    neighboursIdx(IF_REMOVE, :) = [];
    nNeighbours = size(neighboursIdx, 1);
    indX = neighboursIdx(:, 1);
    indY = neighboursIdx(:, 2);
    linIndx = sub2ind(size(A), indX, indY);
    neighbours = A(linIndx);
end
     
    