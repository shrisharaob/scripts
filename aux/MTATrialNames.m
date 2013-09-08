function trialnames = MTATrialNames(varargin)
% trialnames = MTATrialNames(arg)
% arg : filebase / gt
% returns MTA trial names

    arg = DefaultArgs(varargin, {[]});
    trialnames = {};
    if isempty(arg), return; end
    if isa(arg, 'GenericTrial')
        filebase = arg.filebase;
    else 
        filebase = arg;
    end
    try
        mtaTrial = MTATrial(filebase, [], 'all');
        trNames = permute(mtaTrial.list_trialNames, [2, 1]);
        trialnames = trNames(find(~cellfun(@strcmp, trNames, repmat({'all'}, length(trNames), 1))));
    catch err
        return;
    end
end