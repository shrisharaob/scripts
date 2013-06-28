function out = TemplateMatch(gt, varargin);

    [nResample, preOrPost] = DefaultArgs(varargin, {1e2, 'pre'});
    
    sqTemplate = SeqTemplate(gt); % returns template struct
    [evntPeriods, params] = gt.TrajectoryEvents(1, preOrPost, [], [], [], [], 0.6);
    [res, clu] = gt.LoadStateRes('SWS');
    fwdPair = [];
    rvrsPair = [];
    fwdCorr = nan(size(evntPeriods, 1), 1);
    rvrsCorr = nan(size(evntPeriods, 1), 1);
    for kEvntPeriod = 1 : size(evntPeriods, 1)
        [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
        c1 = clu(r1i);
        evntSeq = MyUnique(c1);
        seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
        seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
        %  fwdPair =  [fwdPair; evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder]; 
        if length(seqFwdOrder) < 2 & length(seqRvrsOrder) < 2 , continue; end
        fwdPair =  [evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder]; 
        temp = corr(fwdPair, 'type', 'spearman');
        fwdCorr(kEvntPeriod) = temp(1, 2);
        % rvrsPair =  [rvrsPair; evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder]; 
        rvrsPair =  [evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder]; 
        temp = corr(rvrsPair, 'type', 'spearman');
        rvrsCorr(kEvntPeriod) = temp(1, 2);
    end
    keyboard;
    %% surrogate
    for lResample = 1 : nResample
        lResample
        lFwdPair = [];
        lRvrsPair = [];
        for kEvntPeriod = 1 : size(evntPeriods, 1)
            [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
            c1 = clu(r1i);
            evntSeq = MyUnique(c1);
            seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
            seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
            lFwdPair =  [lFwdPair; evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder(randperm(length(seqFwdOrder)))]; 
            lRvrsPair =  [lRvrsPair; evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder(randperm(length(seqRvrsOrder)))]; 
        end
        trs = corr(lRvrsPair, 'type','Spearman');
        lRvrsCorr(lResample) = trs(1, 2);
        trs = corr(lFwdPair, 'type','Spearman');
        lRwdCorr(lResample) = trs(1, 2);
    end
    keyboard;
end