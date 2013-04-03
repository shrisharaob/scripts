
    function ccgSegsA2B2A = CCGofSegments(trial, varargin)
    % ccgSegs = CcgOfSegments(trial, pfPars, epochs, varargin)
    % epochs = cell (1 x nPairs) 
    % epochs : dir epochs returned by ProjectAonB
    if nargin<1, help CCGofSegments; return; end
    [IF_REPORT_FIG, pairs2analyse, IF_SAVE, fileTag] = DefaultArgs(varargin, {0, [1:size(trial.pfObject.selectedPairs, 1)],1, []});
    pfObject = trial.pfObject;
    
if FileExists(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat'])
    load(['~/data/analysis/' trial.filebase '/' trial.filebase '.SelectedPyrCells.mat']);
    cellPairs = pfObject.selectedPairs(pairs2analyse, :);
else 
    cellPairs =pfObject.selectedPairs;
    
end
    % load dirEpochs 
    load([trial.paths.analysis, trial.filebase '.ProjectAonB.' trial.trialName '.mat']);
    nPairs = size(cellPairs, 1);

    if nPairs ~= 0 && ~isempty(nPairs)
        ccgSegs = cell(1, nPairs);
%         abDir = {'a2b', 'b2a'};
        for kDir = 1 : 2
            if kDir == 1, fprintf('\n ****** a2b ****** \n'); else fprintf('\n ****** b2a ****** \n'); end
            strNew = '';
            for kPair = 1 : nPairs 
                strOld = strNew;
                str = fprintf([repmat('\b', 1, length(strOld)), 'Pair %d of %d '], kPair, nPairs);
                idx = ismember(pfObject.selectedPairs ,cellPairs(kPair,:), 'rows');
                if sum(idx) == 0; continue; end                
                if ~isempty(dirEpochs{idx})
                    if kDir == 1 % a2b;
                        if isfield(dirEpochs{kPair}, 'a2bPeriods')
                            a2bTimes = dirEpochs{idx}.a2bPeriods;
                        else
                            continue;
                        end
                    else %b2a
                        if isfield(dirEpochs{idx}, 'b2aPeriods')
                            a2bTimes = dirEpochs{idx}.b2aPeriods;
                        else
                            continue;
                        end
                      
                    end
                    overlappingPFCells = cellPairs(kPair,:);
                    a2bTimes = round(a2bTimes .* trial.sampleRate / trial.trackingSampleRate) + 1;
                    % load Clu Res for times when trajectories were from 
                    % A -> B
                    if isempty(trial.res)
                        trial = trial.Load({{'CluRes'}});
                    end
                    [a2bRes,a2nResInd] = SelectPeriods(trial.res, a2bTimes,'d', 1, 1);
                    %             [resB,indB] = SelectPeriods(trial.res,timesCellA,'d',1,1);
                    if isempty(a2bRes),   fprintf(' \n no trajectories from %d to %d \n' , overlappingPFCells(1), overlappingPFCells(2)); end
                    if ~isempty(a2bRes)
                        a2bClus = trial.clu(a2nResInd);
                        pRes = [];
                        pClu = [];
                        for kClu = 1:length(overlappingPFCells)
                            pRes = [pRes; a2bRes(a2bClus == overlappingPFCells(kClu))];
                            pClu = [pClu; overlappingPFCells(kClu) * ones(length(a2bRes(a2bClus == overlappingPFCells(kClu))),1)];
                        end
                        %% theta time scale ccg
                        nSpikesA = sum(a2bClus == overlappingPFCells(1, 1));
                        %                     nSpikesB = sum(cluB == overlappingPFCells(1, 2));
                        if nSpikesA > 10
                            binSize = 10e-3; % .5ms
                            halfBins = round(500e-3 / binSize);
                            binSize = round(binSize * trial.sampleRate);
                            [ccgOut, ccgTimeAx, Pairs] = myCCG(pRes,pClu,binSize,halfBins,trial.sampleRate,overlappingPFCells,'count');
                            ccgSegs{kPair}.Out = ccgOut;
                            ccgSegs{kPair}.TimeAx = ccgTimeAx;
                            smoothfactor = .02; %  gaussian with std .02ms
                            tt = ccgTimeAx(1):.1:ccgTimeAx(end);
                            xBins = [-halfBins : halfBins] / halfBins  / 2;
                            gw = exp(-power(xBins, 2) ./ (2*smoothfactor ^2));
                            gw = gw ./ sum(gw);
                            yy = conv(ccgOut(:,1,2), gw, 'same');
                            ccgSmooth = spline(ccgTimeAx',yy,tt);
                            %% beh timescale ccg
                            
                            binSize = 200e-3; % 100ms
                            halfBins = round(2000e-3 / binSize);
                            binSize = round(binSize * trial.sampleRate);
                            [ccgOutBeh, ccgBehTimeAx, ~] = myCCG(pRes,pClu,binSize,halfBins,trial.sampleRate,overlappingPFCells,'count',[],0);
                            ccgSegs{kPair}.OutBeh = ccgOutBeh;
                            ccgSegs{kPair}.BehTimeAx = ccgBehTimeAx;
                            smoothfactor = .05; % std 1s
                            ttBeh = ccgBehTimeAx(1):.1:ccgBehTimeAx(end);
                            xBins = [-halfBins : halfBins] / halfBins  / 2;
                            gw = exp(-power(xBins, 2) ./ (2*smoothfactor ^2));
                            gw = gw ./ sum(gw);
                            yyLarge = conv(ccgOutBeh(:,1,2), gw, 'same');
                            ccgLargeSmooth = spline(ccgBehTimeAx',yyLarge,ttBeh);
                            [~ , lowPassPkLoc] = max(ccgLargeSmooth);
                            lowPassPkTime = ttBeh(lowPassPkLoc);
                            
                            %%
                            
                            ccgSegs{kPair}.Smooth = ccgSmooth;
                            ccgSegs{kPair}.SmoothTimeAx = tt;
                            ccgSegs{kPair}.BehSmooth = ccgLargeSmooth;
                            ccgSegs{kPair}.BehSmoothTimeAx = ttBeh;
                            ccgSegs{kPair}.lowPassOffset = lowPassPkTime;
                            
                            %% plot
                            subplot(2,2,1);
                            hold on;
                            plot(ccgTimeAx,yy,'r-');
                            plot(tt,ccgSegs{kPair}.Smooth ,'g-');
                            line([0,0], ylim,'Color', 'c','LineWidth', 2);
                            line([ccgSegs{kPair}.lowPassOffset , ccgSegs{kPair}.lowPassOffset ], ylim, 'Color', 'm','LineWidth',1.5);
                            [T, offset, firstPeak] = FindCCGPars(ccgSegs{kPair}.Smooth, tt);
                            line([offset, offset], ylim, 'Color', 'r','LineWidth', 1.5);
                            ccgSegs{kPair}.T = T;
                            ccgSegs{kPair}.offset = offset;
                            ccgSegs{kPair}.firstPeak = firstPeak;
                            subplot(224)
                            bar(ccgBehTimeAx, ccgOutBeh(:,1,2), 'b');
                            hold on
                            plot(ttBeh, ccgLargeSmooth,'g-');
                            line([0,0], ylim,'Color', 'c','LineWidth', 2);
                            line([lowPassPkTime, lowPassPkTime], ylim, 'Color', 'm','LineWidth',1.5);
                            axis tight;
                            %%
%                             elClus = acceptedElClu([overlappingPFCells],:);
                            elClus = [trial.elClu(cellPairs(kPair, 1),:); trial.elClu(cellPairs(kPair, 2), :)];
                            if IF_REPORT_FIG
                                filebase = trial.filebase;
                                elCluStr =sprintf('(El, Clu) : (%d,%d) (%d,%d) \n', elClus(1,1), elClus(1,2), elClus(2,1), elClus(2,2));
                                if kDir == 1
                                    filename = [filebase '.' mfilename '.ccgSegs.a2b.' trial.trialName '.' fileTag ];
                                else
                                    filename = [filebase '.' mfilename '.ccgSegs.b2a.' trial.trialName '.' fileTag];
                                end
                                commentString = sprintf('cell pair (El#, Clu#): (%d,%d) (%d,%d) <br> <br> T: %.1f  &nbsp; &nbsp; offset: %.1f &nbsp; &nbsp; firstPk: %.1f &nbsp; &nbsp behPk: %.1f', ...
                                    elClus(1,1), elClus(1,2), elClus(2,1), elClus(2,2), T, offset, firstPeak, lowPassPkTime);  
                                reportfig(gcf, filename , 0, commentString, [],0)
                            else
                                %                 waitforbuttonpress
                                disp(['T =' num2str(T)])
                                disp(['Offset =' num2str(offset)])
                                lowPassPkTime
                                
                            end
                            
                            if kPair ~= nPairs
                                clf;
                            end
                        end
                    end
                end
            end
            if kDir == 1
                ccgSegsA2B2A.a2b = ccgSegs;
            else
                ccgSegsA2B2A.b2a = ccgSegs;
            end
        end
    end
        if IF_SAVE
            if isempty(trial.trialSubType)
                fileName = [trial.filebase, '.', trial.trialName '.' mfilename, '.mat'];
            else
                fileName = [trial.filebase, '.', trial.trialName '.' trial.trialSubType, '.', mfilename, '.mat'];
            end
            save([trial.paths.analysis, fileName]);
        end
    end
    