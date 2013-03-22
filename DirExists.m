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
        for kDir = 1 : length(names)
            if IS_DIR(kDir) == 1
                dirCount = dirCount + 1;
                dirNames{dirCount} = names{dirCount};
            end
        end
        
        for kDir = 1 : sum(IS_DIR)
            
            if strcmp(name, dirNames{kDir})
                IF_DIR_EXISTS = 1;
                continue;
            end
            
        end
    end

end