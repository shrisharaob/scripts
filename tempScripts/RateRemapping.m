function rr = RateRemapping(varargin)
    %
    % list = importdata('~/data/kenji/list');
    load('~/data/kenji/Beh_time_ind.mat');
    elPos = importdata(['~/data/', 'kenji', '/ElePosition.txt']);
    [roi, IF_PLOT] = DefaultArgs(varargin, {{'CA1'}, 1});
    rowId = find(~cellfun(@isempty, regexp(elPos, roi)));
    list = cell(length(rowId));
    for kk = 1 : length(rowId)
        rowCell = regexp(elPos{rowId(kk)}, '\s', 'split');
        list{kk} = rowCell{2};
    end
    list = list(:);
    list(cellfun(@isempty,list)) = [];
    rmpk =  {};
    for i = 1 : length(list)
        nonSleepTrials{i} = Beh(~strcmp(Beh(:,5) ,'sleep') & strcmp(Beh(:,2),list{i}) & (strcmp(Beh(:,5) ,'bigSquare') | strcmp(Beh(:,5) ,'midSquare')) , 4);
        if ~isempty(nonSleepTrials{i})
            fprintf('\n ********* filebase: %s ************** \n', list{i});
            subTrialNames = nonSleepTrials{i};
            cnt = 1;
            for kSubTr = 1 : length(subTrialNames)
                fprintf('\n subtrial %d of %d \n', kSubTr, length(subTrialNames));
                if FileExists(['~/data/kenji/' list{i} '/' subTrialNames{kSubTr} '.whl'])
                    try
                        gt = GenericTrial(list{i}, subTrialNames{kSubTr});
                        gt = gt.LoadPF;
                        [rmpk{cnt}, kCluIdx{cnt}] = RmPeak(gt, roi);
                        if cnt == 1, commonIdx = kCluIdx{1}; end
                        commonIdx = intersect(commonIdx, kCluIdx{cnt});
                        cnt = cnt + 1;
                    catch err
                        fp = fopen('~/data/analysis/kenji/logFilePFr','a');
                        str = sprintf(['\n' subTrialNames{kSubTr}: 'not done']);
                        fwrite(fp, str);
                        fclose(fp);
                    end
                end
            end
            
            for ktr = 1 : cnt - 1
                   sIdx{ktr} = ismember(kCluIdx{ktr}, commonIdx);
            end
            if ~isempty(rmpk)
                refTrRMPk = rmpk{1};
                refTrRMPk = refTrRMPk(sIdx{1});
                cluid = gt.GetRegionClu(roi);
                r = [];
                if ~isempty(cluid{1})
                    baseno = 1;
                    for mSubTr = 2 : cnt - 1
                        mTrRMPk = rmpk{mSubTr};
                        mTrRMPk = mTrRMPk(sIdx{mSubTr});
                        rho = corrcoef(refTrRMPk(:), mTrRMPk(:), 'rows', 'pairwise');
                        plot(refTrRMPk(:), mTrRMPk(:), '*');
                        hold all
                        r{mSubTr} = atan(rho(1, 2)); % Fisher z transform
                                                     %                        plot(baseno, r{mSubTr}, '*');
                                                     %                        hold on;
                        baseno = baseno + 1;
                    end
                    %             r(1) = []; % corr with the same trial
                end
                xlabel('rate pk arena A', 'FontSize', 14);
                ylabel('rate pk arena B', 'FontSize', 14);
                line([0, 60], [0, 60], 'color','k');
                rr{i} = r;
                rmpk =  {};
            end
        end
        fp = fopen('~/data/analysis/kenji/logFilePFr','a');
        fprintf(fp, ['\n' sprintf( repmat('-+', 1, 30))]);
        fclose(fp);
        fclose('all');
    end
    save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' char(roi), '.' mfilename '.mat'], 'rr');
    keyboard;
    if IF_PLOT
        figure;
        bn = 1;
        for i = 1 : 23
            if ~isempty(rr{i})
                for kk = 1 : length(rr{i})
                    kkr = rr{i};
                    rTr = cell2mat(kkr(kk));
                    if ~isempty(rTr)
                        plot(bn, rTr,'o');
                        hold on;

                    end
                end
                bn = bn + 1;
            end
        end
        grid on;
        line([xlim], [0,0], 'color', 'k');
        ylim([-1 1]) ;
        ylabel('corr coef, z-trans');
        xlabel(['filebases with' char(roi) 'recording']);
        title([char(roi) 'remapping']);
    end
end


function [rmpk, cluIdx] = RmPeak(gt, roi)
    rmpk = nan;
    % roi = 'CA3';
    cluid = gt.GetRegionClu(roi);
    cluid = cluid{1};
    idx = ismember(gt.pfObject.acceptedUnits, cluid);
    cluIdx = gt.pfObject.acceptedUnits(idx);
    smRatemp = gt.pfObject.smoothRateMap(:,:,idx);
    for i = 1 : sum(idx)
        smrm = smRatemp(:,:,i);
        if ~isempty(smrm)
            rmpk(i) = nanmax(smrm(:));
        end
    end
end