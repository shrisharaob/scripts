function pfPars = FindPFPars(arg, varargin)
% pfPars = FindPFPars(trial/pfObject, varargin)
%varargin - [pyrCluIdx, smoothFactor, IF_SAVE]

if nargin < 1, help FindPFPars; return; end
    [pyrCluIdx, trialName, smoothFactor, IF_SAVE, states, absThresh] = DefaultArgs(varargin, {[], 'crt1', 0.01, 0, {'head', 'theta'}, 1});
    
    pfPars.com = [];
    pfPars.smoothRateMap = [];
    pfPars.selectedPairs = [];
    pfPars.idealPFUnits = [];
    pfPars.idealPFPairs = [];
    pfPars.spatialCoherence = [];
    pfPars.sparsity = [];
    pfPars.pkLoc = [];
    pfPars.pkDist = [];
    pfPars.acceptedUnits = [];
    pfPars.comDist = [];
    pfPars.ratePk = [];
    if isa(arg, 'GenericPF')
        pfObject = arg;
%         pyrCluIdx = pfObject.acceptedUnits;
    elseif isa(arg, 'GenericTrial')
        pfObject = GenericPF(arg);
    elseif isa(arg, 'MTATrial')
        filebase = arg.name;
        pfObject =LoadPFObject(filebase, trialName);
        pfObject = GenericPF(pfObject);
    elseif isa(arg, 'MTAPlaceField')
        pfObject = arg;
        filebase = pfObject.filebase;
        posOfDots = regexp(filebase,'\.');
        filebase = filebase(1: posOfDots(1) -1);
        arg = MTATrial(filebase, [],pfObject.trialName);
    end
    nUnits = length(pyrCluIdx);
    nAcceptedUnits = 0;
    if ~isempty(nUnits) && nUnits ~= 0
        for kUnit = 1 : nUnits
            kRateMap = pfObject.rateMap{pyrCluIdx(kUnit)}; %%% FIX INDEXING
            % discard cells with firing rate less then 8
            if isempty(kRateMap)
                IS_DISCARD(kUnit) = 1;
            else
                IS_DISCARD(kUnit) = max(kRateMap(:)) < absThresh;
            end
        end
        nAcceptedUnits = sum(~IS_DISCARD);
        acceptedUnits = pyrCluIdx(~IS_DISCARD);
        for kUnit = 1 : nAcceptedUnits
            kRateMap = pfObject.rateMap{acceptedUnits(kUnit)};
            smoothedRateMap = SmoothSurface(kRateMap, smoothFactor);
            [maxSmoothedRate, linIdx] = max(smoothedRateMap(:));
            pfPars.ratePk(kUnit) = maxSmoothedRate;
            [i,j]=ind2sub(size(smoothedRateMap), linIdx);
            maxRateXLoc = pfObject.xBin(i);
            maxRateYLoc = pfObject.yBin(j);
            rateThresh = .707 * maxSmoothedRate;
            threshMask(:,:,kUnit) = kRateMap > rateThresh;
            maskXY = threshMask(:,:,kUnit);
            z = find(maskXY);
            [i, j] = ind2sub(size(maskXY), z);
            z = [i'; j']; % dims x sampels
            comX = 0;
            comY = 0;
            rateSum = 0;
            for ii = 1 : length(i)
                comX = comX + smoothedRateMap(i(ii),j(ii)) * pfObject.xBin(i(ii));
                comY = comY + smoothedRateMap(i(ii),j(ii)) * pfObject.xBin(j(ii));
                rateSum = rateSum + smoothedRateMap(i(ii),j(ii));
            end
            comX = round(comX / rateSum);
            comY = round(comY /rateSum);
            pfPars.com(kUnit, :) = [comX, comY]; %center of mass
            pfPars.smoothRateMap(:,:,kUnit) = smoothedRateMap;
            pfPars.pkLoc(kUnit, :) = [maxRateXLoc, maxRateYLoc];
        end
    end
    
    %% find cell pairs with overlapping Place fields
    if nAcceptedUnits > 1
        cellPairs = nchoosek(acceptedUnits, 2);
        nPairs = nchoosek(nAcceptedUnits, 2);
        for kPair = 1 : nPairs
            temp1 = threshMask(:,:,acceptedUnits == cellPairs(kPair,1)) & threshMask(:,:,acceptedUnits == cellPairs(kPair,2));
            IS_PF_OVERLAP(kPair) = sum(temp1(:)) ~= 0;
        end
        pfPars.selectedPairs = cellPairs(IS_PF_OVERLAP,:);
        pfPars.IS_DISCARD = IS_DISCARD;
        pfPars.acceptedUnits = acceptedUnits; % units that passed the absolute threshold criterion 
        %%  compute place field overlap params
        nOverlappingPairs = sum(IS_PF_OVERLAP);
        for kPair = 1 : nOverlappingPairs
            units = pfPars.selectedPairs(kPair, :);
            unitA = units(1);
            unitB = units(2);
            comA = pfPars.com(acceptedUnits == unitA, :,:);
            comB = pfPars.com(acceptedUnits == unitB, :,:);
            pfPars.comDist(kPair) = norm(comA - comB);
            pfPars.pkDist(kPair) = norm(pfPars.pkLoc(acceptedUnits == unitA, :,:) - pfPars.pkLoc(acceptedUnits == unitB, :,:));  
        end
        %% spatial coherence
        nAcceptedUnits = sum(~IS_DISCARD);
        for kUnit = 1 : nAcceptedUnits
            rateMap = pfObject.rateMap{acceptedUnits(kUnit)};
            [nRows, nClmns] = size(rateMap);
            rateMap(isnan(rateMap)) = 0;
            for kRow = 1:nRows
                for kClmn = 1:nClmns
                    avgNeigbRate(kRow, kClmn) = mean(GetNeighbours(rateMap, kRow, kClmn));
                end
            end
            avgNeigbRate = avgNeigbRate(:);
            rateMap = rateMap(:);
            rho = corrcoef(rateMap, avgNeigbRate);
            spatialCorretaltion(kUnit) = rho(1,2);
            clear avgNeigbRate;
        end
        
        %fisher z transform - ( var[rEstimator] ~ 1/rTrue) 
        zr = atan(spatialCorretaltion);
        pfPars.spatialCoherence = zr;
        %% compute sparsity
        occupancy = pfObject.occupancy;
        for kUnit = 1 : nAcceptedUnits
            if isempty(occupancy)
                gt = GenericTrial(pfObject.filebase, pfObject.trialName);
                curOccupancy = Occupancy(gt);
            else
                curOccupancy = occupancy{acceptedUnits(kUnit)};
            end
            smoothedRateMap = pfPars.smoothRateMap(:,:,kUnit);
%             smoothedRateMap(isnan(smoothedRateMap)) = 0;
            sparsity(kUnit) = (curOccupancy(:)' * smoothedRateMap(:)) ^2 / (curOccupancy(:)' * (smoothedRateMap(:) .^2)); 
        end
        pfPars.sparsity = 1 ./ sparsity; % 0 is uniform firing
        maxSparsity = 0.3;
        minCoherence = 0.6;
        pfPars.idealPFUnits = ~IS_DISCARD & ismember(pyrCluIdx, acceptedUnits(sparsity < maxSparsity)) ...
            & ismember(pyrCluIdx,acceptedUnits(zr > minCoherence)); 
        if sum(pfPars.idealPFUnits) > 1
            pfPars.idealPFPairs = ismember(pfPars.selectedPairs, nchoosek(pyrCluIdx(pfPars.idealPFUnits),2),'rows');
        else 
            pfPars.idealPFPairs = [];
        end
            %%
        if IF_SAVE
            save(['~/data/analysis/' filebase '/' filebase '.' mfilename '.' trialName '.mat'], 'pfPars');
        end
    end
end


