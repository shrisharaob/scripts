function out = SeqCorrs(gt, prePost, varargin)
% out = SeqCorrs(gt)
% correlation(rank) mat of sequences 

    [minCellsInSeq, type, nResample] = DefaultArgs(varargin, {5, 'load', 0});
    out = [];
    switch type
      case 'compute'
        tm = gt.TemplateMatch(0, minCellsInSeq, prePost); 
        eval(['out.seqCorr = SeqCorrAux(tm.', prePost, 'SeqOrder, minCellsInSeq);']);
        if nResample
            n = tm.preNCells;
            n(n==0) = [];
            seqLengths = SampleDistri(n, nResample);
            seqNCells = unique(seqLengths);
            seqLnCnts = histc(seqLengths, min(seqNCells) : max(seqNCells));
            sequences = cell(sum(seqLnCnts), 1);
            idxCnt = 0;
             for mSeqLn = 1 : length(seqNCells)
                for kSeqs = 1 : seqLnCnts(mSeqLn) % generate seqs
                    idxCnt = idxCnt + 1
                    sequences{idxCnt} = randperm(seqNCells(mSeqLn))';
                end
            end
            shuffledCorMat = SeqCorrAux(sequences, minCellsInSeq);
        end
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat']); 
      case 'load'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat']);
      case 'display'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', prePost, '.', mfilename, '.mat']);
        
    end
end

    function seqCorr = SeqCorrAux(sequences, minCellsInSeq)
    nSeq = length(sequences);
    seqCorr = nan(nSeq);
    for kSeq = 1 : nSeq
        kSeqOrder = sequences{kSeq};
        for mSeq = kSeq : nSeq
            mSeqOrder = sequences{mSeq};
            mSeqOrder = mSeqOrder(ismember(mSeqOrder, kSeqOrder));
            kmSeqOrder = kSeqOrder(ismember(kSeqOrder, mSeqOrder));
            if ~isempty(mSeqOrder);
                %if length(mSeqOrder) > minCellsInSeq
                nCells(kSeq, mSeq) = length(kmSeqOrder);
                srho = corr([kmSeqOrder, mSeqOrder], 'type', 'spearman');
                seqCorr(kSeq, mSeq) = srho(1, 2);
                    %end
            end
        end
    end
    end
