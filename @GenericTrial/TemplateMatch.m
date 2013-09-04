function out = TemplateMatch(gt, varargin);
    % out = Template-Match(gt, ragging);
    % Template matching analysis for cell sequences on linear track
    % the significance of a sequence is computed by scrambling the cell
    % ids in the template
    % -------
    % Inputs:
    %     IF_PLOT
    %     minCellsInSeq
    %     preOrPost  
    %     nResample
    %     type
    %     overlap,
    %     alpha

    [IF_PLOT, minCellsInSeq, preOrPost, nResample, type, overlap, alpha] = ...
        DefaultArgs(varargin, {false, 2, 'pre', 1e3, 'load', 0.5, 0.025});
    out = [];
    switch type
      case 'compute'
        sqTemplate = SeqTemplate(gt); % returns template struct
        clus2Select = union(sqTemplate.fwdSortedClu, sqTemplate.rvrsSortedClu);
        [evntPeriods, params] = gt.TrajectoryEvents(1, preOrPost, [], clus2Select, [], [], overlap);
        if size(evntPeriods, 1) > 0
            [res, clu] = gt.LoadStateRes('SWS');
            fwdPair = [];
            rvrsPair = [];
            fwdCorr = nan(size(evntPeriods, 1), 1);
            rvrsCorr = nan(size(evntPeriods, 1), 1);
            nCellsFwd = cell(1, size(evntPeriods, 1));
            nCellsRvrs = cell(1, size(evntPeriods, 1));
            for kEvntPeriod = 1 : size(evntPeriods, 1)
                [r1, r1i] = SelectPeriods(res, evntPeriods(kEvntPeriod,:), 'd', 1,1);
                c1 = clu(r1i);
                evntSeq = MyUnique(c1); % spike order, only the 1st spike in the window is used to identify the order
                seqFwdOrder = sqTemplate.fwdSortedClu(ismember(sqTemplate.fwdSortedClu, evntSeq));
                seqRvrsOrder = sqTemplate.rvrsSortedClu(ismember(sqTemplate.rvrsSortedClu, evntSeq));
                fwdPair =  [evntSeq(ismember(evntSeq, seqFwdOrder)), seqFwdOrder]; 
                rvrsPair =  [evntSeq(ismember(evntSeq, seqRvrsOrder)), seqRvrsOrder]; 
                if size(fwdPair, 1) > minCellsInSeq 
                    if ~isempty(fwdPair)
                        temp = corr(fwdPair, 'type', 'spearman', 'rows', 'complete');
                        fwdCorr(kEvntPeriod) = temp(1, 2);
                        nCellsFwd{kEvntPeriod} = fwdPair(:, 1);
                    end
                end
                if size(rvrsPair, 1) > minCellsInSeq
                    if ~isempty(rvrsPair)
                        temp = corr(rvrsPair, 'type', 'spearman',  'rows', 'complete');
                        rvrsCorr(kEvntPeriod) = temp(1, 2);
                        nCellsRvrs{kEvntPeriod} = rvrsPair(:, 1);
                    end
                end
            end
            
            %% surrogate
            if nResample
                fwdSurrogate = nan(nResample, size(evntPeriods, 1));
                rvrsSurrogate = nan(nResample, size(evntPeriods, 1));
                fprintf('resampling \n ');
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
                            rvrsSurrogate(lResample, kEvntPeriod) = trs(1, 2);
                        end
                        if ~isempty(lFwdPair)
                            trs = corr(lFwdPair, 'type','Spearman', 'rows', 'complete');
                            fwdSurrogate(lResample, kEvntPeriod) = trs(1, 2);
                        end
                    end
                end
                out.fwdSurrogate = fwdSurrogate;
                out.rvrsSurrogate = rvrsSurrogate;
            end
            params.nResample = nResample;
            params.overlap = overlap;
            params.minCellsInSeq = minCellsInSeq;
            out.fwdCorr = fwdCorr;
            out.nCellsFwd = nCellsFwd;
            out.nCellsRvrs = nCellsRvrs;
            out.rvrsCorr = rvrsCorr;
            out.params = params;
            save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' preOrPost, '.', mfilename, '.mat'], 'out');
        else
            out = struct('fwdCorr',[], 'nCellsRvrs', [], 'nCellsFwd', [], 'rvrsCorr', [], 'params', []);
        end
      case 'load'
        IF_PLOT_PRE = false;
        IF_PLOT_POST = false;
        if FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.pre.', mfilename, '.mat'])
            temp = load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.pre.', mfilename, '.mat']);
            preOut =  temp.out;
            clear temp;
            nResample = size(preOut.fwdSurrogate, 1);
            preFwdPval = sum(preOut.fwdSurrogate .* repmat(sign(preOut.fwdCorr'), nResample, 1) >= repmat(abs(preOut.fwdCorr)', nResample, 1)) ./ length(preOut.fwdSurrogate);
            preRvrsPval = sum(preOut.rvrsSurrogate  .* repmat(sign(preOut.rvrsCorr'), nResample, 1) >= repmat(preOut.rvrsCorr', nResample, 1)) ./ length(preOut.rvrsSurrogate);        
            prePVal = [preFwdPval, preRvrsPval];
            prePVal(prePVal == 0) = nan;
            IS_SIGNF_PRE = prePVal < alpha;
            percentPreplay = 100 * nansum(IS_SIGNF_PRE) / length(IS_SIGNF_PRE);

            out.preNCells = [cellfun(@length, preOut.nCellsFwd'); cellfun(@length, preOut.nCellsRvrs')];
            preIdx = out.preNCells > minCellsInSeq;
            preSeqOrder = [preOut.nCellsFwd'; preOut.nCellsRvrs']; % se
            out.preSeqOrder = preSeqOrder(preIdx);
            out.preEvntCorrs = [preOut.fwdCorr(:); preOut.rvrsCorr(:)]; out.preEvntCorrs = out.preEvntCorrs(preIdx);
            out.preSignfEvntCorr = out.preEvntCorrs(IS_SIGNF_PRE(preIdx)); 
            out.preSurrogate = [preOut.fwdSurrogate, preOut.rvrsSurrogate]; out.preSurrogate = Mat2Vec( out.preSurrogate(:, preIdx));
            out.IS_SIGNF_PRE = IS_SIGNF_PRE(preIdx);
            IF_PLOT_PRE = true;
        end
        if FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.post.', mfilename, '.mat']);
            temp = load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.post.', mfilename, '.mat']);
            postOut = temp.out;
            clear temp;
            nResample = size(postOut.fwdSurrogate, 1);
            postFwdPval = sum(postOut.fwdSurrogate >= repmat(postOut.fwdCorr', size(postOut.fwdSurrogate, 1), 1)) ./ length(postOut.fwdSurrogate);
            postRvrsPval = sum(postOut.rvrsSurrogate >= repmat(postOut.rvrsCorr', size(postOut.rvrsSurrogate, 1), 1)) ./ length(postOut.rvrsSurrogate);
            postPVal = [postFwdPval, postRvrsPval];
            postPVal(postPVal == 0) = nan;
            IS_SIGNF_POST = postPVal < alpha;
            percentPost = 100 * nansum(IS_SIGNF_POST) / length(IS_SIGNF_POST);
            out.postNCells = [cellfun(@length, postOut.nCellsFwd'); cellfun(@length, postOut.nCellsRvrs')];
            postIdx = out.postNCells > minCellsInSeq;
            postSeqOrder = [postOut.nCellsFwd'; postOut.nCellsRvrs']; % se
            out.postSeqOrder = postSeqOrder(postIdx);
            out.evnts = [nansum(IS_SIGNF_PRE(preIdx)), length(IS_SIGNF_PRE(preIdx)), nansum(Mat2Vec(IS_SIGNF_POST(postIdx))), length(IS_SIGNF_POST(postIdx))];
            out.postEvntCorrs = [postOut.fwdCorr(:); postOut.rvrsCorr(:)]; out.postEvntCorrs = out.postEvntCorrs(postIdx);
            out.postSignfEvntCorr = out.postEvntCorrs(IS_SIGNF_POST(postIdx));
            out.postSurrogate = [postOut.fwdSurrogate, postOut.rvrsSurrogate]; out.postSurrogate = Mat2Vec(out.postSurrogate(:, postIdx));
            out.IS_SIGNF_POST = IS_SIGNF_POST(postIdx);
            IF_PLOT_POST = true;
        end
        if IF_PLOT
            be = linspace(-1, 1, 1e2); % binEdges
            if IF_PLOT_PRE
                preEvntCnt = histc(out.preEvntCorrs, be);
                if ~all(isnan(preEvntCnt)) & sum(preEvntCnt) > 10, 
                    preSignfEvntCnt = histc(out.preSignfEvntCorr, be);
                    preSurCnt = histc(out.preSurrogate, be);
                    preBinEdg = min(out.preNCells) : max(out.preNCells) + 1; % binEdges for # cells
                    preCellsInEvntCnt = histc(out.preNCells, preBinEdg);  
                    %%
                    hFig = figure;
                    set(hFig, 'position', [1          28        1918         917]);
                    hBar1 = bar(be, preEvntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
                    grid on;
                    hold on;
                    hBar2 = bar(be-.01, preSurCnt / nResample, 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
                    hBar3 = bar(be, preSignfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
                    
                    legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events', 'location', 'NorthWest');
                    legend('boxoff');
                    xlabel('Correlation value');
                    ylabel('Number of events');
                    ylim([0, max(ylim) + max(ylim) * .25]);
                    axHdl =  axes; %('position', [.7, .7, .2, .1]);
                    bar(axHdl, preBinEdg, preCellsInEvntCnt, 'FaceColor', 'k');
                    xlim(axHdl, [min(xlim(axHdl)) - 1, max(xlim(axHdl))]);
                    grid on;
                    set(axHdl, 'position', [.7, .8, .2, .1])
                    set(axHdl, 'Box', 'off');
                    xlabel('# template cells in events', 'FontSize', 5);
                    ylabel('# events', 'FontSize', 6);
                    set(axHdl, 'FontSize', 6);
                    reportfig(hFig, [mfilename, '.pre'], 0, [gt.filebase, 'trial id     : ' gt.trialName], [200, 150]);
                    close(hFig);
                end
            end

            if IF_PLOT_POST
                postEvntCnt = histc(out.postEvntCorrs, be);
                if ~all(isnan(postEvntCnt)) & sum(postEvntCnt) > 10, 
                    postSignfEvntCnt = histc(out.postSignfEvntCorr, be);
                    postSurCnt = histc(out.postSurrogate, be);
                    postBinEdg = min(out.postNCells) : max(out.postNCells) + 1; % binEdges for # cells
                    postCellsInEvntCnt = histc(out.postNCells, postBinEdg);  

                    hFig = figure;
                    set(hFig, 'position', [1          28        1918         917]);
                    hBar1 = bar(be, postEvntCnt, 'FaceColor', 'w', 'BarWidth', 0.5, 'LineWidth', 1);
                    grid on;
                    hold on;
                    hBar2 = bar(be-.01, postSurCnt / nResample, 0.1,'FaceColor', 'c', 'BarWidth', .5, 'EdgeColor', 'none');
                    hBar3 = bar(be, postSignfEvntCnt, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 0.5);
                    legend([hBar3, hBar2, hBar1], 'Signf Events', 'surrogate', 'All Events', 'Location', 'NorthWest');
                    legend('boxoff');
                    xlabel('Correlation value');
                    ylabel('Number of events');
                    ylim([0, max(ylim) + max(ylim) * .25]);
                    axHdl =  axes; %('position', [.7, .7, .2, .1]);
                    bar(axHdl, postBinEdg, postCellsInEvntCnt, 'FaceColor', 'k');
                    xlim(axHdl, [min(xlim(axHdl)) - 1, max(xlim(axHdl))]);
                    grid on;
                    set(axHdl, 'position', [.7, .8, .2, .1]);
                    set(axHdl, 'Box', 'off');
                    xlabel('# template cells in events', 'FontSize', 5);
                    ylabel('# events', 'FontSize', 6);
                    set(axHdl, 'FontSize', 6);
                    reportfig(hFig, [mfilename, '.post'], 0, [gt.filebase, 'trial id     : ' gt.trialName], 200);
                    close(hFig);
                end
            end
        end
    end
end
% LocalWords:  filebase
%BatchProcess(@TemplateMatch, 'kenji', 'CA3', 'linear', 1, {1}, 'allTrials', 0)