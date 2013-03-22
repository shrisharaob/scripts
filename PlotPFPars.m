function PlotPFPars(pfObject, pfPars, cellPair)
% PlotPFPars(pfObj/trial, pfPars, cellPair)
%%
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
        %     pf_search.numZslices = 1;
        pf_search.nbins = 50;
        pf_search.smooth = 0.03;
        
        tr1 = trial.load_Pfs(pf_search); %% loads more than one MTAPF objs
        
        pfObject = tr1.Pfs{2};  %%%%%%% FIX INDEX
        
        trialName = trial.trialName;
    end
else
    trialName = pfObject.trialName;
end
%%
    nTypes = length(pfObject);
    if nTypes >1
        pfObject = pfObject{end};
    end  % this code ensures that only one pfObject is considered

    filebase = pfObject.filebase;
    posOfDots = regexp(filebase,'\.');
    filebase = filebase(1: posOfDots(1) -1);

    load(['~/data/analysis/' filebase '/' filebase  '.ccgAllcells' '.' trialName '.mat'], 'ccg');
    units = load(['~/data/analysis/' filebase '/' filebase '.SelectedPyrCells.mat']);
%     pfPars = FindPFPars(trial,units.linearPyrCluIdx);
%     nPairs = size(pfPars.selectedPairs, 1);
%     for kPair = 1 : nPairs
%       
%     elClus = units.acceptedElClu;
%     elCluStr =[ '(El, Clu) : (' num2str(elClus(1,1)) ',' num2str(elClus(1,2)) ') (' num2str(elClus(2,1)) ',' num2str(elClus(2,2)) ') \n'];
    kPair = ismember(pfPars.selectedPairs, cellPair,'rows');
        if sum(kPair) ~= 0
            subplot(2,3,1);
            PlotPlaceFields(pfObject, pfPars, pfPars.selectedPairs(kPair, 1));
            subplot(2,3,2);
            PlotPlaceFields(pfObject, pfPars, pfPars.selectedPairs(kPair,2));
            subplot(2,3,3);
            PlotPlaceFields(pfObject, pfPars, pfPars.selectedPairs(kPair, 1), 1, 0, [], 'b');
            hold on;
            PlotPlaceFields(pfObject, pfPars, pfPars.selectedPairs(kPair, 2), 1, 0, [], 'r');
            subplot(2, 3, 4)
            bar(ccg.TimeAx, ccg.Out(:, 1, 2, kPair));
            axis tight;
            hold on;
            plot(ccg.SmoothTimeAx,ccg.Smooth(:, kPair),'g-');
            title(['offset: ' num2str(ccg.offset(kPair))]);
            line([0,0], ylim,'Color', 'c','LineWidth', 2);
            line([ccg.lowPassOffset(kPair), ccg.lowPassOffset(kPair)], ylim, 'Color', 'm','LineWidth',1.5);
            line([ccg.offset(kPair), ccg.offset(kPair)], ylim, 'Color', 'r','LineWidth', 1.5);
            subplot(2, 3, 5)
            bar(ccg.BehTimeAx, ccg.OutBeh(:, 1, 2, kPair), 'b');
            hold on;
            plot(ccg.BehSmoothTimeAx,ccg.BehSmooth(:, kPair),'g-');
            line([0,0], ylim,'Color', 'c','LineWidth', 2);
            line([ccg.lowPassOffset(kPair), ccg.lowPassOffset(kPair)], ylim, 'Color', 'm','LineWidth',1.5);
            title(['offset: ' num2str(ccg.lowPassOffset(kPair))]);
            axis tight;
        else
            fprintf('this pair is not accepted');
        end
end
        
        
        