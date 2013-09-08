function pfObject = LoadMTAPFObject(filebase, varargin)
% pfObject = LoadPFObject(filebase, trialName)

if isa(filebase, 'MTATrial')
    filebase = filebase.name;
end

    [trialType] = DefaultArgs(varargin, {'all'});
    trial = MTATrial(filebase,{}, trialType);
    pf_search.filebase = trial.filebase;
    pf_search = MTAPlaceField([]);
    pf_search.mazeName = 'cof';
    pf_search.trialName = trialType;
    pf_search.trackingMarker = 'head_front';
    pf_search.stateLabel = 'head.theta';
    pf_search.spk_shuffle = 'n';
    pf_search.pos_shuffle = 0;
    pf_search.numBSiterations = 1;
    %     pf_search.numZslices = 1;
    pf_search.nbins = 50;
    pf_search.smooth = 0.03;

    tr1 = trial.load_Pfs(pf_search); %% loads more than one MTAPF objs


    pfObject = tr1.Pfs;
    nTrialTypes = length(pfObject);
    idx = [];

    if nTrialTypes ~= 1
        for kTrialType = 1 : nTrialTypes
            trialType =  pfObject{kTrialType}.trialName;
            if strcmp(trialType, trialName)
                idx = kTrialType;
                pfObject = pfObject{idx};
                break
            end
        end
    elseif nTrialTypes == 0
        
        fprintf('trial %s does not exist', trialName);
    end

end
