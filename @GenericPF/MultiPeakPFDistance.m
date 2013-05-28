function out = MultiPeakPFDistance(gpf, varargin)
% out = MultiPeakPFDistance(gpf, varargin)
% 
    defClu = gpf.acceptedUnits;
    defCellPairs = nchoosek(defClu, 2);
    [cellPairs, IF_SMRM, nSTD, areaThreshFactor] = DefaultArgs(varargin, {defCellPairs, 1, 3, 0.5});

    clu = unique(cellPairs(:));
    nCells = length(clu);
    smthRateMaps = gpf.smoothRateMap(:, :, nCells);
    for kCell = 1 : nCells
        if IF_SMRM
            idx = ismember(gpf.acceptedUnits, clu(kCell));
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
                    mCntrVertices{mCntr} = kContr(:, mClm(mCntr) + 1 : mClm(mCntr + 1) - 1 * ~(mCntr == nContrs));
                    mArea(mCntr) = polyarea(mCntrVertices{mCntr}(1, :), mCntrVertices{mCntr}(2, :));
                end
                [maxArea, maxCntrId] = max(mArea);
                areaThresh = areaThreshFactor * maxArea; % discard pf subfield if area less than thresh
                IS_VALID_CNTR{kCell} = mArea >= areaThresh;
                cntrVertices{kCell} = mCntrVertices(IS_VALID_CNTR{kCell});
                clear mArea;
                threshMask(:, :, kCell) = kSmoothRM > rateThresh;
            end
        end
    end
    nPairs = size(cellPairs, 1);
    for mPair = 1 : nPairs
      %   temp1 = threshMask(:, :, clu == cellPairs(kPair,1)) & threshMask(:, :, clu == cellPairs(kPair,2));
%         PF_OVERLAP(mPair) = sum(temp1(:)) ~= 0;
      cntrA = cntrVertices{clu == cellPairs(mPair, 1)};
      cntrB = cntrVertices{clu == cellPairs(mPair, 2)};
      nCntrA = length(cntrA);
      nCntrB = length(cntrB);
      
      if nCntrA + nCntrA >= 2, 
          nCntrPairs = nchoosek([1 : nCntrA, 1 : nCntrB], 2);
          nCntrPairs(nCntrPairs(:, 1) > nCntrA, :) = [];
          nCntrPairs(nCntrPairs(:, 2) > nCntrB, :) = [];
          nCntrPairs = sortrows(unique(nCntrPairs, 'rows')); 
          for kCntrPr = 1 : size(nCntrPairs, 1)
              mOVERLAP{kCntrPr} = any(InPolygon(cntrA{nCntrPairs(kCntrPr, 1)}, cntrB{nCntrPairs(kCntrPr, 2)}));
          end
          OVERLAP{mPair} = mOVERLAP;
          clear mOVERLAP;
      else
          PF_OVERLAP(mPair) = false;
      end
      
      
    end
keyboard
end
