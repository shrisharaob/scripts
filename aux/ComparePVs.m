function out = ComparePVs(filebase, varargin)
% out = ComparePVs(filebase, varargin)
% [datasetType, IF_CHUNKS, nChunks, IF_PLOT, IF_REPORTFIG, roi, arena]
% [], {'CA3'}, {'bigSquare'}, 1, 1
% compute dot product 

    [datasetType, IF_CHUNKS, nChunks, IF_PLOT, IF_REPORTFIG, roi, arena, nSpatialBins] = ...
        DefaultArgs(varargin, {[], 0, 3, 1, 1,  {'CA3'}, {'bigSquare'}, 2500});
    out = [];
     switch datasetType
       case 'kenji'
         analysisFldrPath = '~/data/analysis/kenji/';
       
       case 'MTA'
         roi = 'CA1';
         arena = 'cof';
         analysisFldrPath = ['~/data/analysis/'];

     end
    trialNames = TrialNames(filebase, datasetType, roi, arena);
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
    if IF_PLOT, figHdl = figure; end
    for mTr = 1 : nTrials
        if IF_CHUNKS
            load([analysisFldrPath, filebase, '/', filebase, '.', trialNames{mTr}, filetag, 'CHUNKS.', num2str(nChunks), '.PopVecTimeCourse.mat']);    
            popVec = out.popVec;
            dotProd = out.dotProd;
        else
            load([analysisFldrPath, filebase, '/', filebase, '.', trialNames{mTr}, filetag, 'PopVecTimeCourse.mat']);
        end
        popVec = full(popVec);
        eval([pvNames{mTr} '= popVec;']);
        eval([rvNames{mTr} '=mean('  pvNames{mTr} ', 2);']);
        eval([dpNames{mTr} '= dotProd;']);
        eval(['nDims = size(', pvNames{1}, ', 1)']);
        nCells = nDims / nSpatialBins;
        eval(['out.pv{mTr}=' pvNames{mTr} ';']);
        eval(['out.rv{mTr}=' rvNames{mTr} ';']);
        if IF_PLOT
            figure(figHdl);
            axHdl = subplot(nTrials, 1, mTr);
            eval(['plot(' dpNames{mTr} ', ''*-'', ''MarkerSize'', 14);']);
            title(trialNames{mTr});
            if mTr == nTrials, 
                if IF_CHUNKS, xlabel('chunk #', 'FontSize', 14);
                else, xlabel('# Theta cycles'); end
                set(axHdl, 'XColor', 'k');
            else
                set(axHdl, 'XTickLabel', '');
                set(axHdl, 'XColor', get(axHdl,'Color'));
            end
            ylabel('dot product');
            ylim([-1, 1]);
        end
    end 
    if IF_REPORTFIG
        filename = ['PopVec',  GenFiletag(roi, arena), datasetType];
        commentString = sprintf(['filebase :::: ' filebase, '<br>'  '# units: ' num2str(nCells)]);
        reportfig(figHdl, filename , 0, commentString, [],0);
    end
    close(figHdl);
    
    % compare pv across trials with common units
    dpFunc = @(a, b) a' * b ./ (norm(a) * vnorm(b)); % normalized dot product
    if IF_PLOT, figHdl = figure; end
    if nTrials > 1
        nCycles = inf;
        for kTrPair = 1 : nTrPairs
            eval(['nCycles = min(nCycles, size(' pvNames{trialPairs(kTrPair, 1)} ', 2));']);
        end
        for mTrPair = 1 : nTrPairs
            if IF_PLOT
                figure(figHdl);
                axHdl = subplot(nTrPairs, 1,  mTrPair);
                %                eval(['cdp = dpFunc(Mat2Vec(' rvNames{trialPairs(mTrPair, 1)} '), reshape(' pvNames{trialPairs(mTrPair, 2)} '(:, :, :, 1 : nCycles), [], nCycles));']);
                eval(['cdp = dpFunc(Mat2Vec(' rvNames{trialPairs(mTrPair, 1)} '),' pvNames{trialPairs(mTrPair, 2)} ');']);
                cdp(isnan(cdp)) = 0;
                plot(cdp, '*-', 'MarkerSize', 14);
                ylim([-1, 1]);
                if mTrPair == nTrPairs, 
                    if IF_CHUNKS, xlabel('chunk #', 'FontSize', 14);
                    else, xlabel('# Theta cycles'); end
                    set(axHdl, 'XColor', 'k');
                else
                    set(axHdl, 'XTickLabel', '');
                    set(axHdl, 'XColor', get(axHdl,'Color'));
                end
                title([trialNames{trialPairs(mTrPair,1)}, ' - ', trialNames{trialPairs(mTrPair,2)} ' ' ' (' char(SearchKenji(trialNames{trialPairs(mTrPair, 1)})), ' , ', char(SearchKenji(trialNames{trialPairs(mTrPair, 2)})) ')']);
            end
        end
        if IF_REPORTFIG
            commentString = sprintf(['filebase :::: ' filebase, '<br>'  '# units: ' num2str(nCells)]);
            reportfig(figHdl, filename , 0, commentString, [],0);
        end

    end
    close all;
end
