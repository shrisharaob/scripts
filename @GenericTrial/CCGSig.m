function y = CCGSig(gt, varargin);

    defRandomiz = struct('Type','jitter','nRand',100,'Tau',4,'Alpha',[5 95]);
    [pairs2test, binSiz, halfBins, normalization, randomiz] = DefaultArgs(varargin, {nChoosek(NClusters(gt), 2), 10e-3, 200, 'count', defRandomiz});

    if isempty(gt.res), gt = gt.LoadCR; end
    [res, clu] = gt.LoadStateRes('RUN', 1);

    y = CCGSignif(res, clu, binSiz, halfBins, gt.sampleRate, normalization, pairs2test);
    
end        