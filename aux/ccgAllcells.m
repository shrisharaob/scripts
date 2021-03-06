function ccgAllcells(trial, pfPars, varargin)
    % ccgAllcells(trial, pfPars, varargin)
    % varargin [IF_REPORT_FIG, pairs2analyse, IF_SAVE]
    [IF_REPORT_FIG, pairs2analyse, IF_SAVE, trialName] = DefaultArgs(varargin, {0, [1:size(pfPars.selectedPairs, 1)],1, 'crt1'});
    cellPairs = pfPars.selectedPairs(pairs2analyse, :);
    load(['~/data/analysis/' trial.name '/' trial.name '.SelectCells.mat']);
    nPairs = size(cellPairs, 1);
    if nPairs ~= 0 & ~isempty(nPairs)

        for kPair = 1 : nPairs
            overlappingPFCells = cellPairs(kPair,:);
            pRes = [];
            pClu = [];
            for kClu = 1:length(overlappingPFCells)
                pRes = [pRes; trial.res(trial.clu == overlappingPFCells(kClu))];
                pClu = [pClu; overlappingPFCells(kClu) * ones(length(trial.res(trial.clu == overlappingPFCells(kClu))),1)];
            end
            %% theta time scale ccg
            binSize = 10e-3; % .5ms
            halfBins = round(500e-3 / binSize);
            binSize = round(binSize * trial.sampleRate);

            [ccgOut, ccgTimeAx, ~] = myCCG(pRes,pClu,binSize,halfBins,trial.sampleRate,overlappingPFCells,'count');
            ccg.Out(:,:,:,kPair) = ccgOut;
            ccg.TimeAx = ccgTimeAx;
            smoothfactor = .02; %  gaussian with std .02ms
            tt = ccgTimeAx(1):.1:ccgTimeAx(end);
            xBins = [-halfBins : halfBins] / halfBins  / 2;
            gw = exp(-power(xBins, 2) ./ (2*smoothfactor ^2));
            gw = gw ./ sum(gw);
            yy = conv(ccgOut(:,1,2), gw, 'same');
            ccgSmooth(:, kPair) = spline(ccgTimeAx',yy,tt);

            %% beh timescale ccg
            binSize = 200e-3; % 100ms
            halfBins = round(2000e-3 / binSize);
            binSize = round(binSize * trial.sampleRate);
            [ccgOutBeh, ccgBehTimeAx, ~] = myCCG(pRes,pClu,binSize,halfBins,trial.sampleRate,overlappingPFCells,'count',[],0);
            ccg.OutBeh(:,:,:,kPair) = ccgOutBeh;
            ccg.BehTimeAx = ccgBehTimeAx;
            smoothfactor = .05; % std 1s
            ttBeh = ccgBehTimeAx(1):.1:ccgBehTimeAx(end);
            xBins = [-halfBins : halfBins] / halfBins  / 2;
            gw = exp(-power(xBins, 2) ./ (2*smoothfactor ^2));
            gw = gw ./ sum(gw);
            yyLarge = conv(ccgOutBeh(:,1,2), gw, 'same');
            ccgLargeSmooth(:, kPair) = spline(ccgBehTimeAx',yyLarge,ttBeh);
            [~ , lowPassPkLoc(kPair)] = max(ccgLargeSmooth(:, kPair));
            lowPassPkTime(kPair) = ttBeh(lowPassPkLoc(kPair));
            %% plot
            subplot(2,2,1);
            hold on;
            plot(ccgTimeAx,yy,'r-');
            plot(tt,ccgSmooth(:, kPair),'g-');
            line([0,0], ylim,'Color', 'c','LineWidth', 2);
            line([lowPassPkTime(kPair), lowPassPkTime(kPair)], ylim, 'Color', 'm','LineWidth',1.5);
            [T(kPair), offset(kPair), firstPeak(kPair)] = FindCCGPars(ccgSmooth(:, kPair), tt);
            line([offset(kPair), offset(kPair)], ylim, 'Color', 'r','LineWidth', 1.5);
            % keyboard
            subplot(224)
            bar(ccgBehTimeAx, ccgOutBeh(:,1,2), 'b');
            hold on
            plot(ttBeh, ccgLargeSmooth(:, kPair),'g-');
            line([0,0], ylim,'Color', 'c','LineWidth', 2);
            line([lowPassPkTime(kPair), lowPassPkTime(kPair)], ylim, 'Color', 'm','LineWidth',1.5);
            axis tight;
%%
            elClus = acceptedElClu([overlappingPFCells],:);

            if IF_REPORT_FIG
                elCluStr =[ '(El, Clu) : (' num2str(elClus(1,1)) ',' num2str(elClus(1,2)) ') (' num2str(elClus(2,1)) ',' num2str(elClus(2,2)) ') \n'];
                reportfig(gcf, [trial.name '_ccgs_' trial.trialName ], 0, ['cell pair ' elCluStr  'T:' num2str(T(kPair))  'tOffset: ' num2str(offset(kPair)) 'firstPeak: ' num2str(firstPeak) ], [],0)
            else
                 waitforbuttonpress
%                disp(['T =' num2str(T(kPair))])
%               disp(['Offset =' num2str(offset(kPair))])
%               lowPassPkTime(kPair)
            end
            if kPair ~= nPairs
                clf;
            end
        end
        ccg.Smooth = ccgSmooth;
        ccg.SmoothTimeAx = tt;
        ccg.BehSmooth = ccgLargeSmooth;
        ccg.BehSmoothTimeAx = ttBeh;
        ccg.T = T;
        ccg.offset = offset;
        ccg.firstPeak = firstPeak;
        ccg.lowPassOffset = lowPassPkTime;
               
        if IF_SAVE
            save(['~/data/analysis/' trial.name '/' trial.name '.' mfilename '.' trial.trialName '.mat'], 'ccg');
        end
end




