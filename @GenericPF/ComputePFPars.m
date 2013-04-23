function ComputePFPars(gpf)

    gt = GenericTrial(gpf.filebase, gpf.trialName);
    if isempty(gpf.trialSubType)
        fileName = [gpf.filebase, '.PF.', gpf.trialName '.mat'];
    else
        fileName = [gpf.filebase, '.PF.', gpf.trialName '.' gpf.trialSubType '.mat'];
    end
    fprintf('\n ..... ');
    pfPars = FindPFPars(gpf, gt.pyrCluIdx);
    save([gpf.paths.analysis, fileName], 'pfPars');
end 
    