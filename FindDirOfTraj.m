function trajSegs = FindDirOfTraj(trial, timesInPF, cellPairs, pfPars, varargin)
    % returns indices of sections of trajectories from cell A to Cell B and B - A
    % trial     - MTATrial objekt
    % timesInPF - cell array returned by function GetTimesInPF {1
    % ,nCells}(nEntries-by-2) - use empty cells to keep indexing simpler (i.e nEntries = nAcceptedUnits)
    % cellPairs - nPairs-by-2 
    % valid only for continuous trajectory segments 

    [IF_SAVE] = DefaultArgs(varargin, {0});
    
    nPlaceCells = length(timesInPF);
    nPairs = size(cellPairs, 1);
    acceptedUnits = pfPars.acceptedUnits;
    a2bTimes = [];
    b2aTimes = [];
    trajSegs = cell(1, nPairs);
    for kPair = 1 : nPairs
        fprintf(' \n detecting trajectories crossing PF of cells %d %d \n', cellPairs(kPair, 1), cellPairs(kPair, 2));
        inTimesCellA = timesInPF{acceptedUnits == cellPairs(kPair, 1)};  %  [entryIdx, exitIdx]
        inTimesCellB = timesInPF{acceptedUnits == cellPairs(kPair, 2)};
        if ~(isempty(inTimesCellA) && isempty(inTimesCellB))
            nEntriesCellA = size(inTimesCellA, 1);
            nEntriesCellB = size(inTimesCellB, 1);  
            nEntriesBoth = min(nEntriesCellA, nEntriesCellB);
            % number of times both place fields are crossed by a continuous trajectory 
            [maxEntries, higherEntries] = max([nEntriesCellA, nEntriesCellB]); % cell pf crossed more times
            [minEntries, lowerEntries] = min([nEntriesCellA, nEntriesCellB]);
            
            if higherEntries == 2 % convention Cell A has more entries
                temp = inTimesCellA;
                inTimesCellA = inTimesCellB;
                inTimesCellB = temp;
            end
            entryDiff = (maxEntries - minEntries);
            if entryDiff ~= 0
            %discard the last traj samples
            inTimesCellA(end :-1: end - entryDiff + 1 ,: ) = [];
            end
            % indices to select entry and exit times from inTimes
            a2bCellAIdx = logical(zeros(minEntries, 1));
            b2aCellAIdx = logical(zeros(minEntries, 1));
            a2bCellBIdx = logical(zeros(minEntries, 1));
            b2aCellBIdx = logical(zeros(minEntries, 1));
            
            for ii = 1 : minEntries
                % if entryA < entry/exitB  &&
                % exitA > entryB && exitA < exitB 
                entryA = inTimesCellA(:, 1);
                exitA = inTimesCellA(:, 2);
                entryB = inTimesCellB(ii, 1);
                exitB = inTimesCellB(ii, 2);
                tempIdx = (entryB > entryA) & (entryB < exitA) & (entryA < entryB) & (entryA < exitB) &  (exitA < exitB);
                tempIdx(1:ii-1) = 0;
%                 tempIdx(ii+1 : end) = 0;
                a2bCellAIdx = a2bCellAIdx |  tempIdx;
                a2bCellBIdx(ii) = (a2bCellAIdx(ii) ==1);
                
                % B-to-A entry A > entry B && entry A < exit B 
                %  entryB < entry/exitA  && exitB < exitA
                tempIdx2 = (entryA > entryB) & (entryA < exitB) & (entryB < entryA) & (entryB < exitA) &  (exitB < exitA);
                tempIdx2(1:ii-1) = 0;
%                 tempIdx2(ii+1 : end) = 0;
                b2aCellAIdx = b2aCellAIdx | tempIdx2;
                b2aCellBIdx(ii) = (b2aCellAIdx(ii) ==1);
            end
                a2bTimes = [inTimesCellA(a2bCellAIdx, 1), inTimesCellB(a2bCellBIdx, 2)];
                b2aTimes = [inTimesCellB(b2aCellBIdx, 1), inTimesCellA(b2aCellAIdx, 2)];
                
                a2bTrajectories = [];
                b2aTrajectories = [];
                markerNo = 7;
                if ~isempty(a2bTimes) 
                    for  i = 1 : size(a2bTimes,1)
                        temp = sq(trial.xyz([a2bTimes(i,1) : a2bTimes(i,2)], markerNo ,[1,2]));
                        a2bTrajectories = [ a2bTrajectories ; temp];
                    end
                    trajSegs{kPair}.a2bTimes = a2bTimes;
                    trajSegs{kPair}.a2bTrajectories = a2bTrajectories;
                end
                if ~isempty(b2aTimes)
                    for  i = 1 : size(b2aTimes,1)
                        temp = sq(trial.xyz([b2aTimes(i,1) : b2aTimes(i,2)], markerNo ,[1,2]));
                        b2aTrajectories = [ b2aTrajectories ; temp];
                    end
                    trajSegs{kPair}.b2aTimes = b2aTimes;
                    trajSegs{kPair}.b2aTrajectories = b2aTrajectories;
                end      
        end
    end


    if IF_SAVE
        filebase = trial.name;
        save(['~/data/analysis/' filebase '/' filebase '.' mfilename '.' trial.trialName '.mat']);
    end
    
end
    
    


                     