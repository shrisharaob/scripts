function out = ComparePVs(filebase, varargin)
% out = ComparePVs(filebase, varargin)
% [datasetType, roi, arena, IF_PLOT, IF_REPORTFIG]
% [], {'CA3'}, {'bigSquare'}, 1, 1
% compute dot product 

    [datasetType, roi, arena, IF_PLOT, IF_REPORTFIG] = ...
        DefaultArgs(varargin, {[], {'CA3'}, {'bigSquare'}, 1, 1});
    out = [];
    switch datasetType
      case 'kenji'
        searchStruct.roi = roi;
        searchStruct.arena = arena;
        matches = SearchKenji(searchStruct);
        matches = matches(strcmp(matches(:, 1), filebase), :);
        trialNames = matches(:, 2);
      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
        trialNames = MTATrialNames(filebase);
      otherwise
        error;
        return;
    end
    filetag = GenFiletag(roi, arena);
    nTrials = length(trialNames);
    if nTrials > 1
        trialPairs = nchoosek(1:nTrials, 2);
        nTrPairs = size(trialPairs, 1);
    end
    pvNames = genvarname(cellstr(repmat('pv', length(trialNames), 1)));
    rvNames = genvarname(cellstr(repmat('rv', length(trialNames), 1)));
    dpNames = genvarname(cellstr(repmat('dp', length(trialNames), 1)));
    for mTr = 1 : length(trialNames)
        load(['~/data/analysis/kenji/' filebase, '/', filebase, '.', trialNames{mTr}, filetag, 'PopVecTimeCourse.mat']);
        eval([pvNames{mTr} '= popVec;']);
        eval([rvNames{mTr} '=mean('  pvNames{mTr} ', 4);']);
        eval([dpNames{mTr} '= dotProd;']);
        eval(['nCells = size(', pvNames{1}, ', 1)']);
        eval(['out.pv{mTr}=' pvNames{mTr} ';']);
        eval(['out.rv{mTr}=' rvNames{mTr} ';']);
        if IF_PLOT
            % figHdl = figure;
            % subplot(nTrials, mTr, 1);
            eval(['plot(' dpNames{mTr} ');']);
            %            title(trialNames{mTr});
            xlabel('# Theta cycles');
            ylabel('dot product');
            if IF_REPORTFIG
                filename = ['PopVec',  GenFiletag(roi, arena), datasetType];
                commentString = sprintf(['filebase :::: ' filebase, '<br>'  '# units: ' num2str(nCells), '<br>','trial - ' trialNames{mTr}]);
                reportfig(gcf, filename , 0, commentString, [],0);
            end
        end
    end
    close(figHdl);
    dpFunc = @(a, b) a' * b ./ (norm(a) * vnorm(b)); % normalized dot product
    if nTrials > 1
        for mTrPair = 1 : nTrPairs
            if IF_PLOT
                %   figHdl = figure;
                % subplot(nTrPairs, mTrPair, 
                eval(['nCycles = size(' pvNames{trialPairs(mTrPair, 1)} ', 4);']);
                eval(['cdp = dpFunc(Mat2Vec(' rvNames{trialPairs(mTrPair, 1)} '), reshape(' pvNames{trialPairs(mTrPair, 1)} ', [], nCycles));']);
                cdp(isnan(cdp)) = 0;
                plot(cdp);
                if IF_REPORTFIG
                    commentString = sprintf(['filebase :::: ' filebase, '<br>'  '# units: ' num2str(nCells), '<br>', trialNames{trialPairs(mTrPair, 1)}, ' -  ' trialNames{trialPairs(mTrPair, 2)}, ' (' char(SearchKenji(trialNames{trialPairs(mTrPair, 1)})), ' , ', char(SearchKenji(trialNames{trialPairs(mTrPair, 2)})) ')' ]);
                    reportfig(gcf, filename , 0, commentString, [],0);
                end
            end
        end
    end

end
