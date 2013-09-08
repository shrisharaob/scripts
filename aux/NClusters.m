function [nClusters, nElClus] = NClusters(arg)
    % returns number of clusters in a filebase
    %
    % if isa(arg, 'GenericTrial') || isa(arg, 'GenericPF')
    %     arg = fileBase

    files = dir(arg.paths.data);
    nCluFiles = 0;
    for kFile = 1 : length(files)
        if ~files(kFile).isdir
            %             if ~isempty(regexp(files(kFile).name, [arg.filebase '\.clu\.']))
            if FileExists([arg.paths.data, arg.filebase, '.clu.' num2str(nCluFiles + 1)])
                nCluFiles = nCluFiles + 1;
            end
        end
    end
    nClusters = 0;
    for kCluFile = 1 : nCluFiles
        Fp = fopen([arg.paths.data, arg.filebase, '.clu.', num2str(kCluFile)], 'r');
        
        if Fp==-1
            error(['Could not open file ']);
        end
         nElClus(kCluFile) = fscanf(Fp, '%d', 1) -2;
        nClusters = nElClus(kCluFile) + nClusters;
       
        if nElClus(kCluFile) < 0, nElClus(kCluFile) = 0; end
        fclose(Fp);
    end
   
end