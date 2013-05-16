function out = ComparePVs(filebase, varargin)
% out = ComparePVs(filebase, varargin)
% compare pop vecs 

    [datasetType, roi, arena, IF_PLOT, IF_REPORTFIG] = DefaultArgs(varargin, {{'CA3'}, {'bigSquare'}, 1, 1});
    out = [];
    switch datasetType
      case 'kenji'
        searchStruct.roi = roi;
        searchStruct.arena = arena;
        list = SearchKenji(searchStruct);
        filebases = unique(list(:, 1));
      case 'MTA'
        roi = 'CA1';
        arena = 'cof';
        if ~iscell(filebase), filebase = cellstr(filebase); end
    end
    filetag = GenFiletag(roi, arena);
    sK.roi = roi;
    sK.arena = arena;
    matches = SearchKenji(sK);
    matches = matches(strcmp(matches(:, 1), filebase), :);
    trialNames = matches(:, 2);
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
            eval(['plot(' dpNames{mTr} ');']);
            xlabel('# Theta cycles');
            ylabel('dot product');
            if IF_REPORTFIG
                filename = ['PopVec',  GenFiletag(roi, arena)];
                commentString = [filebase, '- trial id :'  trialNames{mTr}, ' # units: ' num2str(nCells)];
                reportfig(gcf, filename , 0, commentString, [],0);
            end
        end
    end
    
    if nTrials > 1

        for mTrPair = 1 : nTrPairs
            if IF_PLOT
                eval(['nCycles = size(' pvNames{trialPairs(mTrPair, 1)} ', 4);']);
                eval(['plot(transpose(Mat2Vec(' rvNames{trialPairs(mTrPair, 1)} '))* reshape(' pvNames{trialPairs(mTrPair, 1)} ', [], nCycles));']);
                if IF_REPORTFIG
                    filename = ['PopVec',  GenFiletag(roi, arena)];
                    commentString = [trialNames{trialPairs(mTrPair, 1)}, ' -  ' trialNames{trialPairs(mTrPair, 2)}, '::' char(SearchKenji(trialNames{trialPairs(mTrPair, 1)})),' ', char(SearchKenji(trialNames{trialPairs(mTrPair, 2)}))  ' # units: ' num2str(nCells)];
                    reportfig(gcf, filename , 0, commentString, [],0);
                end
            end
        end
    end

end
