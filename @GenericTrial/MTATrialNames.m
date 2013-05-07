function trialnames = MTATrialNames(arg)
% trialnames = MTATrialNames(arg)
% arg : filebase / gt
% returns MTA trial names

    
    trialnames = {};
    if isa(arg, 'GenericTrial')
        filebase = arg.filebase;
    else 
        filebase = arg;
    end
    try
        mtaTrial = MTATrial(filebase, [], 'all');
    catch err
        return;
    end
    trialnames = permute(mtaTrial.list_trialNames, [2, 1]);
    
end