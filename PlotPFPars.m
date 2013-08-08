function PlotPFPars(pfObject, cellPair, varargin)
% PlotPFPars(pfObj/trial, pfPars, cellPair)
%%
[roi, arena] = DefaultArgs(varargin, {'CA3', 'cof'});
if ~isa(pfObject,'MTAPlaceField')
    if isa(pfObject, 'MTATrial')
        trial = pfObject;
        filebase = trial.name;
        pf_search = MTAPlaceField([]);
        pf_search.mazeName = 'cof';
        pf_search.trialName = 'crt1';
        pf_search.trackingMarker = 'head_front';
        pf_search.stateLabel = 'head.theta';
        pf_search.spk_shuffle = 'n';
        pf_search.pos_shuffle = 0;
        pf_search.numBSiterations = 1;
        pf_search.nbins = 50;
        pf_search.smooth = 0.03;
        pfObject = LoadMTAPFObject('jg05-20120315');
        trialName = trial.trialName;
      elseif isa(pfObject, 'GenericPF')
        trialName = pfObject.trialName;
    end
end
%%
    nTypes = length(pfObject);
    if nTypes >1
        pfObject = pfObject{end};
    end  % load only  one pfObject is
    filebase = pfObject.filebase;
%     posOfDots = regexp(filebase,'\.');
%     filebase = filebase(1: posOfDots(1) -1);
   load([pfObject.paths.analysis, pfObject.filebase, '.', pfObject.trialName, GenFiletag(roi, arena), 'CCGPars.mat'], 'out');

    %    units = load(['~/data/analysis/' filebase '/' filebase '.SelecteCells.mat']);
%     pfPars = FindPFPars(trial,units.linearPyrCluIdx);
%     nPairs = size(pfPars.selectedPairs, 1);
%     for kPair = 1 : nPairs
%       
%     elClus = units.acceptedElClu;
%     elCluStr =[ '(El, Clu) : (' num2str(elClus(1,1)) ',' num2str(elClus(1,2)) ') (' num2str(elClus(2,1)) ',' num2str(elClus(2,2)) ') \n'];
    ccg = out;
    hFig = figure;
    for kPr = 1 : size(cellPair, 1)
        kPair = ismember(pfObject.selectedPairs, cellPair(kPr, :),'rows');
        ccgIdx  = ismember(out.cellPairs, cellPair(kPr, :), 'rows');
        if sum(kPair) & sum(ccgIdx) 
            subplot(2, 3, 3);
            PlotRateMaps(pfObject, 0, 0, 1, 0, [],[], pfObject.selectedPairs(kPair, 1)); 
            axis square;
            subplot(2, 3, 6);
            PlotRateMaps(pfObject, 0, 0, 1, 0, [],[], pfObject.selectedPairs(kPair,2)); 
            subplot(2, 3, 2);
            PlotRateMaps( pfObject, 1, 0, 1, 1, 'r',[], pfObject.selectedPairs(kPair, 1)); axis square;
            hold on;
            PlotRateMaps(pfObject, 1, 0, 1, 1, 'k',[], pfObject.selectedPairs(kPair, 2));
            subplot(2, 3, 1);
            bar(ccg.ccgTimeAx, ccg.ccg{ccgIdx}(:, 1, 1)); axis tight; axis square;
            subplot(2, 3, 5)
            bar(ccg.ccgTimeAx, ccg.ccg{ccgIdx}(:, 2, 2)); axis tight; axis square;
            subplot(2, 3, 4);
            bar(ccg.ccgTimeAx, ccg.ccg{ccgIdx}(:, 1, 2)); axis tight; axis square;
            axis tight;
            hold on;
            plot(ccg.smthTAx,ccg.smthCCG(:, ccgIdx),'r-');
            %    title(['offset: ' num2str(ccg.offset(ccgIdx))]);
            line([0,0], ylim,'Color', 'c');
            line([ccg.firstPeak(ccgIdx), ccg.firstPeak(ccgIdx)], ylim, 'Color', 'm');
            line([ccg.offset(ccgIdx), ccg.offset(ccgIdx)], ylim, 'Color', 'r');
            % subplot(2, 3, 5)
%             bar(ccg.BehTimeAx, ccg.OutBeh(:, 1, 2, ccgIdx), 'b');
%             hold on;
%             plot(ccg.BehSmoothTimeAx,ccg.BehSmooth(:, ccgIdx),'g-');
%             line([0,0], ylim,'Color', 'c','LineWidth', 2);
%             line([ccg.lowPassOffset(ccgIdx), ccg.lowPassOffset(ccgIdx)], ylim, 'Color', 'm','LineWidth',1.5);
%             title(['offset: ' num2str(ccg.lowPassOffset(ccgIdx))]);
%             axis tight;
        else
            fprintf('this pair is not accepted');
        end
        keyboard;
        clf;
    end
    keyboard;
end
        
        
        