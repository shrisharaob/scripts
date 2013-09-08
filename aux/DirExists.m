function IF_DIR_EXISTS =  DirExists(loc,  name)
% check if dirctory <name> exists at location <loc>
 
    IF_DIR_EXISTS = 0;
    files = dir(loc);
    if ~isempty(files)        
        isdir = CatStruct(files, 'isdir', 1);
        IS_DIR = isdir.isdir;
        names = CatStruct(files, 'name', 1);
        names = names.name;
        dirNames = {};
        dirCount = 0;
        dirNames = names(IS_DIR);
        IF_DIR_EXISTS = any(cellfun(@strcmp, dirNames, repmat({name}, length(dirNames), 1)));
    end
end