function out = MultiPeakPFDistance(gpf, varargin)
% out = MultiPeakPFDistance(gpf, varargin)
% 

    [cells2Analyze, IF_SMRM, nSTD, areaThreshFactor] = DefaultArgs(varargin, {[], 1, 3, 0.5});

    nCells = length(cells2Analyze);

    for kCell = 1 : nCells
        if IF_SMRM
            idx = ismember(gpf.acceptedUnits, cells2Analyze(kCell));
            if ~isempty(idx)
                kSmoothRM = gpf.smoothRateMap(:, :, idx);
                rateThresh = nSTD * std(kSmoothRM(:));
                [kContr ] = contour(kSmoothRM, [1, 1] .* rateThresh);
                [nRows, nClmns] = size(kContr);
                nContrs = 0;
                IS_DONE = 0;
                nVals = 0;
                while ~IS_DONE
                    nContrs = nContrs + 1;
                    mClm(nContrs) = 1 + sum(nVals) + (nContrs - 1);
                    if mClm(nContrs) > nClmns, 
                        mClm(nContrs) = nClmns;
                        nContrs = nContrs - 1;
                        IS_DONE = 1; continue; end
                    nVals(nContrs) = kContr(2, mClm(nContrs));
                end
                for mCntr = 1 : nContrs
                    mCntrVertices{kCell, mCntr} = kContr(:, mClm(mCntr) + 1 : mClm(mCntr + 1) - 1 * ~(mCntr == nContrs));
                    mArea(mCntr) = polyarea(mCntrVertices{kCell, mCntr}(1, :), mCntrVertices{kCell, mCntr}(2, :));
                end
                [maxArea, maxCntrId] = max(mArea);
                areaThresh = areaThreshFactor * maxArea;
                IS_VALID_CNTR{kCell} = mArea >= areaThresh;
                clear mArea;
            end
        end
    end
keyboard;
end
