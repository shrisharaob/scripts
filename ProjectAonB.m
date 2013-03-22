function dirEpochs = ProjectAonB(trial, pfPars)
    % Projects velocity vector at each point of the trajecory onto the
    % vector from cell A to B 
    % dirEpochs = ProjectAonB(trial, pfPars)
    % History :
    %   Shrisha - Created
    
    if nargin < 1, help ProjectAonB; return; end
    
    nPairs = size(pfPars.selectedPairs,1);
    dirEpochs = cell(1, nPairs);
    for lPair = 1 : nPairs
        cellPair = pfPars.selectedPairs(lPair, :);
%%%*******************
        %         comA = pfPars.com(pfPars.acceptedUnits == cellPair(1), :);
%         comB = pfPars.com(pfPars.acceptedUnits == cellPair(2), :);
%         *************************************

comA = pfPars.pkLoc(pfPars.acceptedUnits == cellPair(1), :);
comB = pfPars.pkLoc(pfPars.acceptedUnits == cellPair(2), :);
        
        a2bVector = comA - comB;
        fprintf('\n')% *** pair %d , %d \n', cellPair(1), cellPair(2));
        resA = trial.res(trial.clu == cellPair(1));
        resB = trial.res(trial.clu == cellPair(2));

        markerNo = 7;
        allTrajs = sq(trial.posistion(:, markerNo, [1,2]));
         % cell array containging the traj segments
        allTrajSegs= PlotTrajectories(trial, [],[], [], 1, [], 0);
        nTrajSegs = size(allTrajSegs, 2);
        a2bPeriods = [];
        b2aPeriods = [];
        strNew = '';
        for kTrajSeg = 1 : nTrajSegs
            strOld = strNew;
            strNew = sprintf('*** pair %d , %d *** \t trajSeg %d of %d', cellPair(1), cellPair(2), kTrajSeg, nTrajSegs);
            fprintf([repmat('\b', 1, length(strOld)), strNew], kTrajSeg, nTrajSegs);
            if ~isempty(allTrajSegs{kTrajSeg})
                velocity = diff(allTrajSegs{kTrajSeg}, 1, 1);
                [nVelRows, nVelClmns] = size(velocity);
                a2bDir = sign(velocity * a2bVector');
                a2bDir(a2bDir == 0) = 1; 
                t = flipud(a2bDir);
                a2bDir = [ -1 * a2bDir(1); a2bDir];
                t = [-1* t(1);t];
                
                a2bBeginIndx = find(diff(a2bDir) == 2);
                a2bEndIndx = find(flipud(diff(t) == 2));
                b2aBeginIndx = find(diff(a2bDir) == -2);
                b2aEndIndx = find(flipud(diff(t) ==-2));
                INVALID_A2B = a2bBeginIndx == a2bEndIndx;
                INVALID_B2A = b2aBeginIndx == b2aEndIndx;
                a2bBeginIndx(INVALID_A2B) = [];
                a2bEndIndx(INVALID_A2B) = [];
                b2aBeginIndx(INVALID_B2A) = [];
                b2aEndIndx(INVALID_B2A) = [];
                                
                a2bEpoch = [a2bBeginIndx, a2bEndIndx];
                b2aEpoch = [b2aBeginIndx, b2aEndIndx];
                
                curTrajSeg = allTrajSegs{kTrajSeg};
                if ~isempty(a2bEpoch)
                    a2bStart = find(ismember(allTrajs, curTrajSeg(a2bEpoch(:, 1), :), 'rows'));
                    a2bEnd = find(ismember(allTrajs, curTrajSeg(a2bEpoch(:, 2), :), 'rows'));
                    tempA2BPeriods = [a2bStart, a2bEnd];
                    a2bPeriods = [a2bPeriods; tempA2BPeriods];
                end
                if ~isempty(b2aEpoch)
                    b2aStart = find(ismember(allTrajs, curTrajSeg(b2aEpoch(:, 1), :), 'rows'));
                    b2aEnd = find(ismember(allTrajs, curTrajSeg(b2aEpoch(:,2), :), 'rows'));
                    tempB2APeriods = [b2aStart, b2aEnd];
                    b2aPeriods  = [b2aPeriods; tempB2APeriods];
                end
                
            end
        end
        dirEpochs{lPair}.a2bPeriods = a2bPeriods;
        dirEpochs{lPair}.b2aPeriods = b2aPeriods;
%         keyboard;
    end
    filebase = trial.name;
    
    save(['~/data/analysis/' filebase '/' filebase '.' mfilename '.' trial.trialName '.mat'], 'dirEpochs');
end