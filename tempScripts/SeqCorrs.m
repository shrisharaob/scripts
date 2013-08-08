function out = SeqCorrs(gt, prePost, varargin)
% out = SeqCorrs(gt)
% correlation(rank) mat of sequences 
%
% BatchProcess(@SeqCorrs, 'kenji', 'CA1', 'linear', 1, {'post',[], 'compute'},'allTrials')



    [minCellsInSeq, type, nResample, alpha, nBins] = DefaultArgs(varargin, {5, 'display', 0, 0.05, 5});
    out = [];
    switch type
      case 'compute'
        tm = gt.TemplateMatch(0, minCellsInSeq, prePost); 
        try
            eval(['[out.seqCorr, out.nCells] = SeqCorrAux(tm.', prePost, 'SeqOrder, minCellsInSeq);']);
            eval(['nSeqs = sum(tm.', prePost, 'NCells >  minCellsInSeq);']);
            if nResample
                n = tm.preNCells;
                n(n==0) = [];
                seqLengths = SampleDistri(n, nResample);
                seqNCells = unique(seqLengths);
                seqLnCnts = histc(seqLengths, min(seqNCells) : max(seqNCells));
                sequences = cell(sum(seqLnCnts), 1);
                idxCnt = 0;
                fprintf('\n resampling (M = %d) ... \n', nResample);
                for mSeqLn = 1 : length(seqNCells)
                    for kSeqs = 1 : seqLnCnts(mSeqLn) % generate seqs
                        idxCnt = idxCnt + 1;
                        if ~(mod(idxCnt, floor(nResample / 10))), fprintf('#'); end
                        sequences{idxCnt} = randperm(seqNCells(mSeqLn))';
                    end
                end
                fprintf('  done ... \n cov mat ... ')
                out.shuffledCorMat = SeqCorrAux(sequences, minCellsInSeq);
            end
        catch
            
        end
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat'], 'out'); 
      case 'load'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat']);
        data = out.seqCorr;
        IS_SIGNF = false(size(data));
        inValidIdx = ~(eye(length(out.nCells)) == 0 & out.nCells > 5); 
        data = data(~inValidIdx);
        nCells = out.nCells(~inValidIdx);
        shuffledCorrs = out.shuffledCorMat(setdiff(1:numel(out.shuffledCorMat), linspace(1, numel(out.shuffledCorMat), length(out.shuffledCorMat))));
        fprintf('\n computing pvals...')
        matlabpool open 8
        parfor mSeq = 1 : length(data)
            pval(mSeq) = sum(shuffledCorrs' .* repmat(sign(data(mSeq))', length(shuffledCorrs), 1) >= repmat(abs(data(mSeq))', length(shuffledCorrs), 1)) ./ length(shuffledCorrs); 
        end
        matlabpool close
        fprintf(' done ! \n');
        IS_SIGNF(find(~inValidIdx)) = pval <= alpha;
        out.IS_SIGNF = sparse(IS_SIGNF);
        out.pval(find(~inValidIdx)) = sparse(pval);
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat'], 'out'); 
      case 'display'
        myout = load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat']);
        hFig = figure;
        try
            out.signfTimeCourse = sum(myout.out.IS_SIGNF, 2);
            [count, binCenter] = hist(find(out.signfTimeCourse), nBins);            
            count = count ./ sum(count);
            [maxCnt, cntIdx] = max(count);
            bar(binCenter, count, 'FaceColor', 'k');
            if sum(count) > 0
                [b, ~, ~, ~, stats] = regress(count', [ones(nBins, 1), binCenter(:)]);
                out.regress = [maxCnt, cntIdx(1), b', stats]; % [intercept, slope, R^2, F, pVal, errVariance]
                y = b(1)+ b(2) * (binCenter(1):binCenter(end)) ;
            end
            hold on;
            plot(binCenter(1):binCenter(end), y, 'r');
            %set(gca, 'xtick', []);
            ylabel('Significant Events');
            reportfig(hFig, [mfilename, prePost], 0, [gt.filebase, '    ', gt.trialName, '    ']);
            close(hFig);
        catch err
            fprintf(['\n', err.message, '\n']); 
            close(hFig);
        end
    end
end

%% 
function [seqCorr, nCells] = SeqCorrAux(sequences, minCellsInSeq)
        nSeq = length(sequences);
        seqCorr = nan(nSeq);
        for kSeq = 1 : nSeq
            kSeqOrder = sequences{kSeq};
            for mSeq = kSeq : nSeq
                mSeqOrder = sequences{mSeq};
                mSeqOrder = mSeqOrder(ismember(mSeqOrder, kSeqOrder));
                kmSeqOrder = kSeqOrder(ismember(kSeqOrder, mSeqOrder));
                if ~isempty(mSeqOrder);
                    nCells(kSeq, mSeq) = length(kmSeqOrder);
                    %                    srho = corr([kmSeqOrder, mSeqOrder], 'type', 'spearman');
                    %seqCorr(kSeq, mSeq) = srho(1, 2);
                end
            end
        end
    end
