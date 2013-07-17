function out = TemplateTimeCourse(gt, varargin)

    [nResample, nChunks, prePost, minCellsInSeq, alpha] = DefaultArgs(varargin, {0, 3, 'pre', 5, .05});
    out = [];
    temp = gt.TemplateMatch(0, minCellsInSeq, prePost);
    switch prePost
      case 'pre'
        tm.EvntCorrs = temp.preEvntCorrs;
        tm.Surrogate = temp.preSurrogate;
      case 'post'
        tm.EvntCorrs = temp.postEvntCorrs;
        tm.Surrogate = temp.postSurrogate;
    end
    if ~isempty(tm.EvntCorrs)
        %        evntperiods = gt.TrajectoryEvents(0, Post);
        %evntperiods = RecenterPeriods(evntperiods);
        nEvents = length(tm.EvntCorrs);
        IF_VALID = false;
        if nEvents > nChunks
            %bins = 1 : round(size(evntperiods, 1 )./ nChunks) : size(evntperiods, 1);
            bins = 1 : floor(nEvents / nChunks) : nEvents;
            sbins = [1, bins(2:end) * 1e3];
            bins = [bins', circshift(bins', -1)]; bins(end, :) = [];
            sbins = [sbins', circshift(sbins', -1)]; sbins(end, :) = [];
            colours = GenColormap(nChunks);
            figHdl = figure;
            if size(bins, 1) == nChunks
                if ~isempty(tm.EvntCorrs) & ~isempty(tm.Surrogate)
                    for kChunk = 1 : nChunks
                        kIdx = linspace(bins(kChunk, 1), bins(kChunk, 2), diff(bins(kChunk, :)) + 1);
                        kSurIdx = linspace(sbins(kChunk, 1), sbins(kChunk, 2), diff(sbins(kChunk, :)) + 1);
                        [h(kChunk), p(kChunk), ks2stat(kChunk)] = kstest2(tm.Surrogate(kSurIdx), tm.EvntCorrs(kIdx), alpha);
                        [empCDF, xd] = ecdf(tm.EvntCorrs(kIdx));
                        [xd, xdi] = unique(xd);
                        empCDF = empCDF(xdi);
                        [empSurCDF, xs] = ecdf(tm.Surrogate(kSurIdx));
                        [xs, xsi] = unique(xs);                    
                        empSurCDF = empSurCDF(xsi);
                        empCDF = empCDF(ismember(xd, xs));
                        empSurCDF = empSurCDF(ismember(xs, xd));
                        cdfDiff{kChunk} = empSurCDF - empCDF;
                        if nResample
                            empCDFmat = [empSurCDF, empCDF];
                            for mResample = 1 : nResample
                                shuffleId = rand(size(empCDF, 1), 1) > 0.5;
                                shuffledDiff(:, mResample) = diff(RowWiseShift(empCDFmat, shuffleId), 1, 2);
                            end
                        end
                        figure(figHdl), hold on;
                        linHdl(kChunk) = plot(xd, cdfDiff{kChunk},'o-');
                        set(linHdl(kChunk), 'Color', colours(kChunk, :));
                        grid on;
                        % figure(201),plot(kChunk, ks2stat(kChunk), char(~h(kChunk) * 'r*-' + h(kChunk) * 'g*-'));
                        % hold on;
                    end
                end
                IF_VALID = true;
            end
            %[~, chunkId] = histc(find(tm.IS_SIGNF_), bins);

            %    figure(201),plot(1: nChunks, ks2stat, 'k')
            % set(gca, 'XTick', 1 : nChunks);
            %grid(gca, 'on');
            % xlabel('ks statistic');
            
            [empCDF, xd] = ecdf(tm.EvntCorrs);
            [xd, xdi] = unique(xd);
            empCDF = empCDF(xdi);
            [empSurCDF, xs] = ecdf(tm.Surrogate);
            [xs, xsi] = unique(xs);                    
            empSurCDF = empSurCDF(xsi);
            [~, pVal, ~] = kstest2(tm.EvntCorrs, tm.Surrogate);
            cdfDiffAll = empSurCDF(ismember(xs, unique(xd))) - empCDF(ismember(unique(xd), xs));
            figure(figHdl), linHdl(nChunks + 1) = plot(unique(xd), cdfDiffAll, '*-k'); grid on;
            text(-.9, max(ylim) - .1 * max(ylim) - 0.025 * ~(max(ylim)), ['p = ', num2str(pVal, '%.4f')], 'FontSize', 12);
            legendStr = str2cell(num2str(1 : nChunks));
            legendStr{nChunks + 1} = 'all';
            if IF_VALID
                legend(linHdl, legendStr');
                line(xlim, [0,0], 'color', 'k', 'LineStyle', '--');
                xlabel('Correlation Value');
                try
                    reportfig(figHdl, [mfilename,'.', prePost, '.'], [], [gt.filebase, '    ', gt.trialName]);
                catch, end
                close(figHdl);
            else
                close(figHdl);
            end
        end
    end
end



