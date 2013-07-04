function out = TemplateMatch(gt, varargin);
% out = Template-Match(gt, ragging);
% [resemble, proposed]

    [nResample, preOrPost, type, overlap, alpha] = DefaultArgs(varargin, {1e3, 'pre', 'display', 0.5, 0.025});
    
    switch type
      case 'compute'
        sqTemplate = SeqTemplate(gt); % returns template stricter
        clus2Select = union(sqTemplate.fwdSortedClu, sqTemplate.rvrsSortedClu);
        [evntPeriods, params] = gt.TrajectoryEvents(1, preOrPost, [], clus2Select, [], [], overlap);
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
            if length(seqFwdOrder) < 2 & length(seqRvrsOrder) < 2 , continue; end
            fwdPair =  [evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder]; 
            if ~isempty(fwdPair)
                temp = corr(fwdPair, 'type', 'spearman', 'rows', 'complete');
                fwdCorr(kEvntPeriod) = temp(1, 2);
            end
            rvrsPair =  [evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder]; 
            if ~isempty(rvrsPair)
                temp = corr(rvrsPair, 'type', 'spearman',  'rows', 'complete');
                rvrsCorr(kEvntPeriod) = temp(1, 2);
            end
        end
        
        %% surrogate
        if nResample
            out.fwdSurrogate = nan(nResample, size(evntPeriods, 1));
            out.rvrsSurrogate = nan(nResample, size(evntPeriods, 1));
            for kEvntPeriod = 1 : size(evntPeriods, 1)
                [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
                c1 = clu(r1i);
                evntSeq = MyUnique(c1);
                seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
                seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
                for lResample = 1 : nResample
                    lFwdPair =  [evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder(randperm(length(seqFwdOrder)))]; 
                    lRvrsPair =  [evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder(randperm(length(seqRvrsOrder)))];  
                    if ~isempty(lRvrsPair)
                        trs = corr(lRvrsPair, 'type','Spearman', 'rows', 'complete');
                        out.rvrsSurrogate(lResample, kEvntPeriod) = trs(1, 2);
                    end
                    if ~isempty(lFwdPair)
                        trs = corr(lFwdPair, 'type','Spearman', 'rows', 'complete');
                        out.fwdSurrogate(lResample, kEvntPeriod) = trs(1, 2);
                    end
                end
            end
        end
        out.fwdCorr = fwdCorr;
        out.rvrsCorr = rvrsCorr;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' preOrPost, '.', mfilename, '.mat'], 'out');
      case 'display'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.pre.', mfilename, '.mat'], 'out');
        preOut =  out;
        nResample = size(preOut.fwdSurrogate, 1);
        %load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.post.', mfilename, '.mat']);
        % postOut = out;
        preFwdPval = sum(preOut.fwdSurrogate .* repmat(sign(preOut.fwdCorr'), nResample, 1) >= repmat(abs(preOut.fwdCorr)', nResample, 1)) ./ length(preOut.fwdSurrogate);
        % postFwdPval = sum(postOut.fwdSurrogate >= repmat(postOut.fwdCorr', size(postOut.fwdSurrogate, 1), 1)) ./ length(postOut.fwdSurrogate);
        preRvrsPval = sum(preOut.rvrsSurrogate  .* repmat(sign(preOut.rvrsCorr'), nResample, 1) >= repmat(preOut.rvrsCorr', nResample, 1)) ./ length(preOut.rvrsSurrogate);
        % postRvrsPval = sum(postOut.rvrsSurrogate >= repmat(postOut.rvrsCorr', size(postOut.rvrsSurrogate, 1), 1)) ./ length(postOut.rvrsSurrogate);
        prePVal = [preFwdPval, preRvrsPval];
        % postPVal = [postFwdPval, postRvrsPval];
        prePVal(prePVal == 0) = nan;
        % postPVal(postPVal == 0) = nan;
        IS_SIGNF_PRE = prePVal < alpha;
        % IS_SIGNF_POST = postPVal < alpha;
        percentPreplay = 100 * nansum(IS_SIGNF_PRE) / length(IS_SIGNF_PRE);
        % percentPost = 100 * nansum(IS_SIGNF_POST) / length(IS_SIGNF_POST);
        %        out.evnts = [nansum(IS_SIGNF_PRE(:)), length(IS_SIGNF_PRE), nansum(IS_SIGNF_POST(:)), length(IS_SIGNF_POST)];
        out.evnts = [nansum(IS_SIGNF_PRE(:)), length(IS_SIGNF_PRE)];
        out.evntCorrs = [preOut.fwdCorr(:); preOut.rvrsCorr(:)];
        out.signfEvntCorr = out.evntCorrs(IS_SIGNF_PRE);
        out.surrogate = [preOut.fwdSurrogate(:); preOut.rvrsSurrogate(:)];
        %        keyboard;
    end
end
% LocalWords:  filebase
