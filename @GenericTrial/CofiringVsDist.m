function out =  CofiringVsDist(gt, varargin)
% CofiringVsDist(gt, varargin)
% script to relate pk distances and cofiring likelihood

    [prePost, IF_PLOT, type, roi, arena] = ...
        DefaultArgs(varargin, {'pre', 0, 'load','CA3', 'bigSquare'});
    %    if ~FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat']), type = 'compute'; end
    switch type
        case 'compute'
          out = [];
          disp('computing');
          pairCofiring = PairReactivation(gt, prePost);
          if ~isempty(pairCofiring)
              cluId = unique(pairCofiring.cellPairs(:));

              pks = gt.MultiPeakPFDistance(roi, arena); close gcf;
              cmnClus = intersect(pks.cluId, cluId);
              cellIds = pks.cluId;
              if length(cmnClus) > 1
                  myPairs = nchoosek(cmnClus, 2);
                  cntrPeaks = pks.cntrPeaks;%(ismember(pks.cluId, cmnClu));
                  cntrVertices = pks.cntrVertices;%(ismember(pks.cluId, cmnClu));
                  validCntrCnt = 0;
                  for mCellPair = 1 : size(myPairs, 1)
                      cntrA = cntrVertices{cellIds == myPairs(mCellPair, 1)}; 
                      cntrB = cntrVertices{cellIds == myPairs(mCellPair, 2)};
                      nCntrA = length(cntrA);
                      nCntrB = length(cntrB);
                      if nCntrA == 1 && nCntrB == 1, % select only cells with single pf's 
                          cntrPairs = nchoosek([1 : nCntrA, 1 : nCntrB], 2); % all pairs of selected sub contours
                          cntrPairs(cntrPairs(:, 1) > nCntrA, :) = [];
                          cntrPairs(cntrPairs(:, 2) > nCntrB, :) = [];
                          cntrPairs = sortrows(unique(cntrPairs, 'rows')); 
                          pkA = cntrPeaks{cellIds == myPairs(mCellPair, 1)};
                          pkB = cntrPeaks{cellIds == myPairs(mCellPair, 2)};
                          for kCntrPr = 1 : size(cntrPairs, 1)
                              validCntrCnt = validCntrCnt + 1;
                              probCofiringSleep(validCntrCnt) = pairCofiring.sleepCfProb(ismember(pairCofiring.cellPairs, myPairs(mCellPair, :), 'rows'));
                              pkDistAB(validCntrCnt) = norm( pkA(cntrPairs(kCntrPr, 1), :) - pkB(cntrPairs(kCntrPr, 2), :));
                              selectedCellpairs(validCntrCnt, :) = myPairs(mCellPair, :);
                              pkAB(validCntrCnt, :) = [pkA(cntrPairs(kCntrPr, 1), :) , pkB(cntrPairs(kCntrPr, 2), :)];
                          end
                      end
                  end
                  if validCntrCnt > 0
                      out.cfVsPkDist = [pkDistAB(:), probCofiringSleep(:)];
                      save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'out');
                      fprintf('\n DONE... \n');
                  else
                      fprintf('\n no valid cntr pairs \n');
                  end
              else
                  out = [];
                  fprintf('\n no common pairs active both in run and sleep \n ');
              end              
          end          

      case 'load'
        if FileExists([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat']);
            load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.', mfilename, '.mat'], 'out');
            if IF_PLOT
                % hist2(out.cfVsPkDist);
                if ~isempty(out)
                    plot(out.cfVsPkDist(:, 1), out.cfVsPkDist(:, 2), '*');
                    axis square;
                    grid on;
                    title('SWS');
                    xlabel('peak Dist (px)');
                    ylabel('cofiring probablity');
                    reportfig(gcf, [mfilename, '.', char(prePost)] , 0, [gt.filebase, '    ', gt.trialName, '  roi : ' roi, '  arena : ' arena]);
                end
            end
        else, out = [];
        end
    end
end



