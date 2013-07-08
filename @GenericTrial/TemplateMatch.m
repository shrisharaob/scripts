function out = TemplateMatch(gt, varargin);
% out = Template-Match(gt, ragging);
% [resemble, proposed]

    [IF_PLOT, nResample, preOrPost, type, overlap, alpha, minCellsInSeq] = ...
        DefaultArgs(varargin, {false, 1e3, 'pre', 'load', 0.5, 0.025, 5});
    
    switch type
      case 'compute'
        sqTemplate = SeqTemplate(gt); % returns template stricter
        clus2Select = union(sqTemplate.fwdSortedClu, sqTemplate.rvrsSortedClu);
        [evntPeriods, params] = gt.TrajectoryEvents(0, preOrPost, [], clus2Select, [], [], overlap);
        [res, clu] = gt.LoadStateRes('SWS');
        fwdPair = [];
        rvrsPair = [];
        fwdCorr = nan(size(evntPeriods, 1), 1);
        rvrsCorr = nan(size(evntPeriods, 1), 1);
        matlabpool local 8
        parfor kEvntPeriod = 1 : size(evntPeriods, 1)
            [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
            c1 = clu(r1i);
            evntSeq = MyUnique(c1); % spike order, only the 1st spike in the window is used to identify the order
            seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
            seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
            if length(seqFwdOrder) < minCellsInSeq & length(seqRvrsOrder) < minCellsInSeq , continue; end
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
            if rvrsCorr(kEvntPeriod) == -1 | fwdCorr(kEvntPeriod) == -1 , keyboard; end
        end
        matlabpool close
        keyboard;
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
        out.rvrsCorr = rvrsCorr;keyboard;
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' preOrPost, '.', mfilename, '.mat'], 'out');
      case 'load'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.pre.', mfilename, '.mat'], 'out');
        preOut =  out;
        nResample = size(preOut.fwdSurrogate, 1);
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.post.', mfilename, '.mat']);
        postOut = out;
        clear out;
        preFwdPval = sum(preOut.fwdSurrogate .* repmat(sign(preOut.fwdCorr'), nResample, 1) >= repmat(abs(preOut.fwdCorr)', nResample, 1)) ./ length(preOut.fwdSurrogate);
        postFwdPval = sum(postOut.fwdSurrogate >= repmat(postOut.fwdCorr', size(postOut.fwdSurrogate, 1), 1)) ./ length(postOut.fwdSurrogate);
        preRvrsPval = sum(preOut.rvrsSurrogate  .* repmat(sign(preOut.rvrsCorr'), nResample, 1) >= repmat(preOut.rvrsCorr', nResample, 1)) ./ length(preOut.rvrsSurrogate);
        postRvrsPval = sum(postOut.rvrsSurrogate >= repmat(postOut.rvrsCorr', size(postOut.rvrsSurrogate, 1), 1)) ./ length(postOut.rvrsSurrogate);
        prePVal = [preFwdPval, preRvrsPval];
        postPVal = [postFwdPval, postRvrsPval];
        prePVal(prePVal == 0) = nan;
        postPVal(postPVal == 0) = nan;
        IS_SIGNF_PRE = prePVal < alpha;
        IS_SIGNF_POST = postPVal < alpha;
        percentPreplay = 100 * nansum(IS_SIGNF_PRE) / length(IS_SIGNF_PRE);
        percentPost = 100 * nansum(IS_SIGNF_POST) / length(IS_SIGNF_POST);
        out.evnts = [nansum(IS_SIGNF_PRE(:)), length(IS_SIGNF_PRE), nansum(IS_SIGNF_POST(:)), length(IS_SIGNF_POST)];
        %out.evnts = [nansum(IS_SIGNF_PRE(:)), length(IS_SIGNF_PRE)];
        out.preEvntCorrs = [preOut.fwdCorr(:); preOut.rvrsCorr(:)];
        out.postEvntCorrs = [postOut.fwdCorr(:); postOut.rvrsCorr(:)];
        out.preSignfEvntCorr = out.preEvntCorrs(IS_SIGNF_PRE);
        out.postSignfEvntCorr = out.postEvntCorrs(IS_SIGNF_POST);
        out.preSurrogate = [preOut.fwdSurrogate(:); preOut.rvrsSurrogate(:)];
        out.postSurrogate = [postOut.fwdSurrogate(:); postOut.rvrsSurrogate(:)];
        %        keyboard;
        
        if IF_PLOT
            be = linspace(-1, 1, 1e2);
            %fwdCnt = histc(out.fwdCorr, be);
            %fwdSurr = histc(out.fwdSurrogate, be);
            %rvrsCnt = histc(out.rvrsCorr, be);
            %rvrsSurr = histc(out.rvrsSurrogate, be);
            
            preEvntCnt = histc(out.preEvntCorrs, be);
            preSignfEvntCnt = histc(out.preSignfEvntCorr, be);
            preSurCnt = histc(out.preSurrogate, be);
            
            postEvntCnt = histc(out.postEvntCorrs, be);
            postSignfEvntCnt = histc(out.postSignfEvntCorr, be);
            postSurCnt = histc(out.postSurrogate, be);

            if ~all(isnan(preEvntCnt)) | sum(preEvntCnt) > 10, 
                hFig = figure;
                set(hFig, 'position', [1          28        1918         917]);
                hBar1 = bar(be, preEvntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
                grid on;
                hold on;
                hBar2 = bar(be-.01, preSurCnt / nResample, 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
                hBar3 = bar(be, preSignfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
                legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events');
                legend('boxoff');
                xlabel('Correlation value');
                ylabel('Number of events');
                reportfig(hFig, [mfilename, '.pre'], 0, [gt.filebase, 'trial id     : ' gt.trialName], 200);
                close(hFig);
            end

            
            if ~all(isnan(postEvntCnt)) | sum(postEvntCnt) > 10, 
                hFig = figure;
                set(hFig, 'position', [1          28        1918         917]);
                hBar1 = bar(be, postEvntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
                grid on;
                hold on;
                hBar2 = bar(be-.01, postSurCnt / nResample, 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
                hBar3 = bar(be, postSignfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
                legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events');
                legend('boxoff');
                xlabel('Correlation value');
                ylabel('Number of events');
                reportfig(hFig, [mfilename, '.post'], 0, [gt.filebase, 'trial id     : ' gt.trialName], 200);
                close(hFig);
            end

        end
    end
end
% LocalWords:  filebase
%BatchProcess(@TemplateMatch, 'kenji', 'CA3', 'linear', 1, {1}, 'allTrials', 0)