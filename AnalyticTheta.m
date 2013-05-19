function [ThPh, ThAmp] = AnalyticTheta(gt, varargin)

    [channel, fRange, IF_COMPUTE] = DefaultArgs(varargin, {1, [4 12], 0});
    
    if ~IF_COMPUTE & FileExists([gt.paths.data, gt.filebase, gt.trialName, '.thpar.mat'])
        load([gt.paths.data, gt.filebase, gt.trialName, '.thpar.mat']);
        return;
    end    
    
    fprintf('\n computing theta phase \n');
    switch gt.datasetType
        case 'MTA'
          mtr = MTATrial(gt.filebase, {{'lfp', channel}}, gt.trialName);
          gt.lfp = mtr.lfp;
    end

    x = ButterBandPass(gt.lfp, gt.lfpSampleRate,fRange, 4);
    xa = hilbert(x);
    ThAmp = abs(xa);
    ThPh = angle(xa);
    save([gt.paths.data, gt.filebase, gt.trialName, '.thpar.mat'], 'ThAmp', 'ThPh');
end
    
    
