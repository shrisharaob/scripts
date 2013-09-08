function filetag = GenFiletag(varargin)


    [roi, arena] = DefaultArgs(varargin, {'', ''});
    filetag = '';
    if ~iscell(roi), roi = cellstr(roi); end
    if ~iscell(arena), arena = cellstr(arena); end

    if ~any(cellfun(@isempty, roi))
        for mRoi = 1 : length(roi)
            if mRoi == 1
                filetag = ['.', char(roi{1})];
            end
            if mRoi > 1
                filetag = [filetag, '.', char(roi{mRoi})];
            end
        end
    end
    if ~any(cellfun(@isempty, arena))
        for lArena = 1 : length(arena)
            filetag = [filetag, '.', char(arena{lArena})];
        end
    end
    if  ~isempty(filetag)
        if ~strcmp(filetag(end), '.') 
            filetag = [filetag, '.'];
        end
    end
end 

