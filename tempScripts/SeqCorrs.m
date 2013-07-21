% corr mat of sequences


    nSeq = length(sequences);
    seqCorr = zeros(nSeq);
    
    for kSeq = 1 : nSeq
        kSeqOrder = sequences{kSeq};
        for mSeq = kSeq : nSeq
            mSeqOrder = sequences{mSeq};
            mSeqOrder = mSeqOrder(ismember(mSeqOrder, kSeqOrder));
            kmSeqOrder = kSeqOrder(ismember(kSeqOrder, mSeqOrder));
            if ~isempty(mSeqOrder);
                %                if length(mSeqOrder) > 2
                    mCells(kSeq, mSeq) = length(kmSeqOrder);
                    srho = corr([kmSeqOrder, mSeqOrder], 'type', 'spearman');
                    seqCorr(kSeq, mSeq) = srho(1, 2);
                    %end
            end
        end
    end