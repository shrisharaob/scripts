function CompareSessionPF(filebase, trialTypes, cluIdx,  varargin)
    % sessions - cell array containing sessions to compare
    [filebase, IF_WAITFORBUTTONPRESS, IF_REPORTFIG] = DefaultArgs(varargin, {'jg05-20120315', 0, 0 });
    units  = load(['~/data/analysis/' filebase '/' filebase '.SelectedPyrCells.mat']);
    
    nTypes = length(trialTypes);
    varNames = genvarname(cellstr(repmat('pfObj', nTypes, 1)));
    varParNames = genvarname(cellstr(repmat('pfPars', nTypes, 1)));
    for kType = 1 : nTypes
        eval([varParNames{kType} '=load(''~/data/analysis/' filebase '/' filebase '.FindPFPars.' trialTypes{kType} '.mat'')']);
        eval([varParNames{kType} '=' varParNames{kType} '.pfPars']);
        eval([varNames{kType}  '= LoadPFObject(''' filebase ''',''' trialTypes{kType} ''')']);
    end
    
   
    nCells = length(cluIdx);
    for mCell = 1 : nCells
         idx = ismember(units.linearPyrCluIdx, cluIdx(mCell));
        fprintf('\n %d', units.linearPyrCluIdx(idx));
        clf;
        for kType = 1 : nTypes
            subplot(1, nTypes, kType);
            eval(['PlotPlaceFields(' varNames{kType} ',' varParNames{kType} ',' units.linearPyrCluIdx(mCell) ')']);
%             eval([varNames{kType} '.plot(units.linearPyrCluIdx(mCell))']);
%             eval(['title(''' trialTypes{kType} ''')']);
%             axis off;
        end
        if IF_WAITFORBUTTONPRESS, try waitforbuttonpress; catch err ; end, end
        if IF_REPORTFIG
            MyReportfig(gcf, [filebase '.' mfilename '.selectedPlaceCells.mat'], 0, [ 'clu#:'  num2str(units.linearPyrCluIdx(idx))]);
    end
end
    
        