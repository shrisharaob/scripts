function out = PoolPVs(varargin)
% out = ComparePVs(filebase, varargin)
% [datasetType, IF_CHUNKS, nChunks, IF_PLOT, IF_REPORTFIG, roi, arena]
% [], {'CA3'}, {'bigSquare'}, 1, 1
% compute dot product 

[datasetType, IF_CHUNKS, nChunks, IF_PLOT, marker, IF_REPORTFIG, roi, arena, nSpatialBins] = ...
    DefaultArgs(varargin, {[], 0, 3, 1, '*', 1,  {'CA3'}, {'bigSquare'}, 2500});
out = [];
switch datasetType
  case 'kenji'
    analysisFldrPath = '~/data/analysis/kenji/';
    
  case 'MTA'
    roi = 'CA1';
    arena = 'cof';
    analysisFldrPath = ['~/data/analysis/'];

end
filebase = SearchKenji;
filebase = unique(filebase(:, 1));
poolCount = 0;
dpFunc = @(a, b) a' * b ./ (norm(a) * vnorm(b)); % normalized dot product
for kBase = 1 : length(filebase)
    trialNames = TrialNames(filebase{kBase}, datasetType, roi, arena);
    filetag = GenFiletag(roi, arena);
    nTrials = length(trialNames);
    if nTrials > 1
        [trA, trB] = meshgrid(1:nTrials, 1:nTrials);
        trialPairs = [trA(:), trB(:)];
        trialPairs(trA(:) == trB(:), :) = [];
        nTrPairs = size(trialPairs, 1);
    end
    pvNames = genvarname(cellstr(repmat('pv', length(trialNames), 1)));
    rvNames = genvarname(cellstr(repmat('rv', length(trialNames), 1)));
    dpNames = genvarname(cellstr(repmat('dp', length(trialNames), 1)));
    
    for mTr = 1 : nTrials
        if IF_CHUNKS
            if FileExists([analysisFldrPath, filebase{kBase}, '/', filebase{kBase}, '.', trialNames{mTr}, filetag, 'CHUNKS.', num2str(nChunks), '.PopVecTimeCourse.mat'])
                try
                    fprintf(['\n' trialNames{mTr}])
                    load([analysisFldrPath, filebase{kBase}, '/', filebase{kBase}, '.', trialNames{mTr}, filetag, 'CHUNKS.', num2str(nChunks), '.PopVecTimeCourse.mat']);    
                    popVec = out.popVec;
                    dotProd = out.dotProd;
                catch err
                    continue;
                end
            else
                continue; 
            end
        else
            load([analysisFldrPath, filebase{kBase}, '/', filebase{kBase}, '.', trialNames{mTr}, filetag, 'PopVecTimeCourse.mat']);
        end
        popVec = full(popVec);
        eval([pvNames{mTr} '= popVec;']);
        eval([rvNames{mTr} '=mean('  pvNames{mTr} ', 2);']);
        eval([dpNames{mTr} '= dotProd;']);
        eval(['nDims = size(', pvNames{1}, ', 1);']);
        nCells = nDims / nSpatialBins;
        eval(['out.pv{mTr}=' pvNames{mTr} ';']);
        eval(['out.rv{mTr}=' rvNames{mTr} ';']);
        eval(['cdp = dpFunc(Mat2Vec(' rvNames{mTr} '),' pvNames{mTr} ');']);
        poolCount = poolCount + 1;
        pooledDp(poolCount, :) = cdp;
    end
    
    % compare pv across trials with common units
    if ~isempty(out)
        if nTrials > 1
            nCycles = inf;
            for kTrPair = 1 : nTrPairs
                eval(['nCycles = min(nCycles, size(' pvNames{trialPairs(kTrPair, 1)} ', 2));']);
            end
            for mTrPair = 1 : nTrPairs
                eval(['cdp = dpFunc(Mat2Vec(' rvNames{trialPairs(mTrPair, 1)} '),' pvNames{trialPairs(mTrPair, 2)} ');']);
                cdp(isnan(cdp)) = 0;
                poolCount = poolCount + 1;
                pooledDp(poolCount, :) = cdp;
            end
        end
    end
end
keyboard;
scatter(pooledDp(:, 1), pooledDp(:, 2), 'b', marker);
hold on;
scatter(pooledDp(:, 1), pooledDp(:, 3), 'r', marker);
scatter(pooledDp(:, 2), pooledDp(:, 3), 'g', marker);
line([0, 1], [0, 1], 'color', 'k');
legend('1 - 2', '1 - 3', '2 - 3','Location', 'NorthWest');
set(gca, 'FontSize', 14);
keyboard;
end