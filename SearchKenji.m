function out = SearchKenji(varargin)
    % returns the filebase and trial names stisfying the search keys
    % args.arena - search struct 
    %   args.arena.arena
    %   args.arena.args.roi
    defArgs = struct('roi', 'CA3', 'arena', 'bigSquare');     
    [args, searchInFilebase, trialName, IF_GET_MAZE] = DefaultArgs(varargin, {defArgs, [], [], false});
    IN_FILEBASE = ~isempty(searchInFilebase);
    load('~/data/kenji/Beh_time_ind.mat'); % loads var Beh
    out = Beh(:, [2, 4]); %[filebase, trialName];
    idx = true(size(Beh,1), 1);
   if ~isstruct(args), IF_GET_MAZE = true; end
    if IF_GET_MAZE
        out = Beh(strcmp(Beh(:, 2), args) & ~strcmp(Beh(:, 5), 'sleep'), [2, 4, 5]);
        if isempty(out) % check if arg is a trial name
             out = Beh(strcmp(Beh(:, 4), args) & ~strcmp(Beh(:, 5), 'sleep'), [5]);
        end
        if isempty(out)
            out{1} = '';
        end
    else
    %% search for matching arenas 
    if isfield(args, 'arena')
        if ~isempty(args.arena)
            if ~iscell(args.arena), args.arena = cellstr(args.arena); end
            
            for kArena = 1 : length(args.arena)
                filebases = unique(Beh(strcmp(Beh(:, 5), args.arena{kArena}), 2));
%                 trialNames = Beh(strcmp(Beh(:, 5), args.arena{kArena}), 4);
                if kArena == 1
                    validFilebases = filebases;
                end
%                 trialNames = trialNames(ismember(filebases, validFilebases));
                validFilebases = filebases(ismember(filebases, validFilebases));
            end
            idx = false(size(Beh,1), 1);
            for nBase = 1 : length(validFilebases)
                for mArena = 1 : length(args.arena)
                    idx = idx | strcmp(Beh(:, 2),validFilebases{nBase}) & strcmp(Beh(:, 5), args.arena{mArena});
                    %                 trialNames = Beh(strcmp(Beh(:, 5), args.arena{kArena}), 4);
                end
            end
%             idx = logical(idx | repmat(sum(idx,2), 1, size(idx,2))); % get all the colums for matches
%             filebase = Beh(idx, 2); % filebase :: trialname
%             trialName = Beh(idx, 4);
            out = Beh(idx, [2, 4]); %[filebase, trialName];
            INVALID = false(1,size(out,1));
            for kMatch = 1 : size(out,1)
                kTrialNames = out(strcmp(out(:, 1), out{kMatch, 1}), 2);
                for ktr = 1 : size(kTrialNames, 1)
                    if ~FileExists(['~/data/kenji/' out{kMatch, 1} '/' kTrialNames{ktr} '.whl'])  || ~FileExists(['~/data/kenji/whl/',out{kMatch, 1}, '.eegTime'])
                        INVALID(kMatch)  = true;
                    end
                end
            end
            out(INVALID, :) = [];
            clear INVALID;
        end
    end
    %% search for args.rois
    if isfield(args, 'roi')
        if ~isempty(args.roi)
            if ~iscell(args.roi), args.roi = cellstr(args.roi); end
            elPos = importdata(['~/data/kenji/ElePosition.txt']);
            roiIdx = false(size(out,1) , 1);
            filebase = unique(out(:,1));
            INVALID = false(1,size(filebase,1));
            for kBase = 1 : size(filebase,1)
                rowId = find(~cellfun(@isempty, regexp(elPos, filebase(kBase, 1))));
                if ~isempty(rowId)
                    for mRoi = 1 : length(args.roi)
                        rowCell = regexp(elPos{rowId}, '\s', 'split');
                        INVALID(kBase) =  ~any(cellfun(@any,regexp(rowCell, args.roi{mRoi})));
                    end
                else
                    INVALID(kBase) = true;
                end            
            end
            filebase(INVALID) = [];
            out(~ismember(out(:,1), filebase), :) = [];
            %             out{:,3} = Beh(ismember(Beh(:, [2, 4]), out, 'rows'), 5);  
        end
    end
    
    if IN_FILEBASE % check if the search obj pars exists in the given filebase
        lst = out;
        out = 0;
        if ~isempty(lst)
            out = any(~cellfun(@isempty, regexp(lst(:, 1), searchInFilebase)));
            if ~isempty(trialName)
                out = out && any(~cellfun(@isempty, regexp(lst(:, 2), trialName)));
            end
        end
    end
end

end

