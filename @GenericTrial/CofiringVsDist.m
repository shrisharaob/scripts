function out =  CofiringVsDist(gt, varargin)
% CofiringVsDist(gt, varargin)
% script to relate pk distances and cofiring likelihood

    [prePost, type, roi, arena ] = DefaultArgs(varargin, {'pre', 'load','CA3', 'bigSquare'});
    if ~FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat']), type = 'compute'; end
    switch type
        case 'compute'
          out = [];
          disp('computing');
          pairCofiring = PairReactivation(gt, prePost);
          if ~isempty(pairCofiring)
              cluId = unique(pairCofiring.cellPairs(:));

              pks = gt.MultiPeakPFDistance(roi, arena);close all;
              cmnClus = intersect(pks.cluId, cluId);
              cellIds = pks.cluId;

              if cmnClus > 1
                  myPairs = nchoosek(cmnClus, 2);
                  cntrPeaks = pks.cntrPeaks;%(ismember(pks.cluId, cmnClu));
                  cntrVertices = pks.cntrVertices;%(ismember(pks.cluId, cmnClu));
                  validCntrCnt = 0;
                  for mCellPair = 1 : size(myPairs, 1)
                      cntrA = cntrVertices{cellIds == myPairs(mCellPair, 1)}; 
                      cntrB = cntrVertices{cellIds == myPairs(mCellPair, 2)};
                      nCntrA = length(cntrA);
                      nCntrB = length(cntrB);
                      if nCntrA >= 1 && nCntrB >= 1, 
                          cntrPairs = nchoosek([1 : nCntrA, 1 : nCntrB], 2); % all pairs of selected sub contours
                          cntrPairs(cntrPairs(:, 1) > nCntrA, :) = [];
                          cntrPairs(cntrPairs(:, 2) > nCntrB, :) = [];
                          cntrPairs = sortrows(unique(cntrPairs, 'rows')); 
                          pkA = cntrPeaks{cellIds == myPairs(mCellPair, 1)};
                          pkB = cntrPeaks{cellIds == myPairs(mCellPair, 2)};
                          for kCntrPr = 1 : size(cntrPairs, 1)
                              validCntrCnt = validCntrCnt + 1;
                              probCofiring(validCntrCnt) = pairCofiring.cofiringProb(ismember(pairCofiring.cellPairs, myPairs(mCellPair, :), 'rows'));
                              pkDistAB(validCntrCnt) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                              selectedCellpairs(validCntrCnt, :) = myPairs(mCellPair, :);
                              pkAB(validCntrCnt, :) = [pkA(cntrPairs(kCntrPr, 1), :) , pkB(cntrPairs(kCntrPr, 2), :)];
                          end
                      else
                          selectedCellpairs = []; pkDistAB = []; pkAB = [];
                      end
                  end
              end
              out.cfVsPkDist = [probCofiring(:), pkDistAB(:)];
          end          
          save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'out');
      case 'load'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'out');
        plot(cfVsPkDist(:, 2), cfVsPkDist(:, 1));
    end
end



