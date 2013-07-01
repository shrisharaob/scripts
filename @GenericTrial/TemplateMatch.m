function out = TemplateMatch(gt, varargin);
% out = TemplateMatch(gt, varargin);
% [nResample, preOrPost]

    [nResample, preOrPost] = DefaultArgs(varargin, {1e2, 'pre'});
    
    sqTemplate = SeqTemplate(gt); % returns template struct
    clus2Select = union(sqTemplate.fwdSortedClu, sqTemplate.rvrsSortedClu);
    [evntPeriods, params] = gt.TrajectoryEvents(1, preOrPost, [], clus2Select, [], [], 0.1);
    [res, clu] = gt.LoadStateRes('SWS');
    fwdPair = [];
    rvrsPair = [];
    fwdCorr = nan(size(evntPeriods, 1), 1);
    rvrsCorr = nan(size(evntPeriods, 1), 1);
    for kEvntPeriod = 1 : size(evntPeriods, 1)
        [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
        c1 = clu(r1i);
        evntSeq = MyUnique(c1); % spike order, only the 1st spike in the window is used to identify the order
        seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
        seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
        %  fwdPair =  [fwdPair; evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder]; 
        if length(seqFwdOrder) < 2 & length(seqRvrsOrder) < 2 , continue; end
        fwdPair =  [evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder]; 
        if ~isempty(fwdPair)
            temp = corr(fwdPair, 'type', 'spearman', 'rows', 'complete');
            fwdCorr(kEvntPeriod) = temp(1, 2);
        end
        % rvrsPair =  [rvrsPair; evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder]; 
        rvrsPair =  [evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder]; 
        if ~isempty(rvrsPair)
            temp = corr(rvrsPair, 'type', 'spearman',  'rows', 'complete');
            rvrsCorr(kEvntPeriod) = temp(1, 2);
        end
    end
    
    %% surrogate
    if nResample
        lfwdCorr = nan(nResample, 1);
        lRvrsCorr = nan(nResample, 1);
        for lResample = 1 : nResample
            lResample
            lFwdPair = nan(1, 2);
            lRvrsPair = nan(1, 2);
            for kEvntPeriod = 1 : size(evntPeriods, 1)
                [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
                c1 = clu(r1i);
                evntSeq = MyUnique(c1);
                seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
                seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
                lFwdPair =  [lFwdPair; evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder(randperm(length(seqFwdOrder)))]; 
                lRvrsPair =  [lRvrsPair; evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder(randperm(length(seqRvrsOrder)))];  
            end
            if ~isempty(lRvrsPair)
                trs = corr(lRvrsPair, 'type','Spearman', 'rows', 'complete');
                lRvrsCorr(lResample) = trs(1, 2);
            end
            if ~isempty(lFwdPair)
                trs = corr(lFwdPair, 'type','Spearman', 'rows', 'complete');
                lFwdCorr(lResample) = trs(1, 2);
            end
        end
        keyboard;
        [fwdNullCount, ~] = histc(lFwdCorr, linspace(-1, 1, 1e2));
        [rvrsNullCount, ~] = histc(lRvrsCorr, linspace(-1, 1, 1e2))
        out.fwdPval = sum(repmat(lFwdCorr(:), 1, size(evntPeriods, 1))   >= repmat(fwdCorr', nResample, 1)) ./  nResample;
        out.rvrsPval = sum(repmat(rvrsCorr(:), 1, size(evntPeriods, 1))   >= repmat(rvrsCorr', nResample, 1)) ./  nResample;
    end
    out.fwdCorr = fwdCorr;
    out.rvrsCorr = rvrsCorr;
    save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' preOrPost, '.', mfilename, '.mat'], 'out');
end