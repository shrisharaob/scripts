function out = SparsifyPV(gt, varargin)
out = [];
[roi, arena, nDims] = DefaultArgs(varargin, {'CA3', 'bigSquare', 2500});

%    load([gt.paths.analysis, gt.filebase, GenFiletag(roi, arena), 'commonClus.mat']);
[oldpv, ~, dotProd] = PopVecTimeCourse(gt, [], [], roi, arena);
[nClus,~,~,nCycles] = size(oldpv);
popVec = [];
save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena),  'PopVecTimeCourse.mat'], 'popVec', '-append','-v7.3');
popVec = sparse(reshape(oldpv, nClus * nDims, nCycles));
save([gt.paths.analysis, gt.filebase, '.', gt.trialName, GenFiletag(roi, arena), 'PopVecTimeCourse.mat'], 'popVec', '-append','-v7.3');  
end