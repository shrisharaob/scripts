function [decodedPos, err] = DecodePos(gt, varargin)
    % [pos, predErr] = DecodePos(gt, varargin)
    % Baysean decoding of position from spike times and average rate
    % maps, Ref. Zhang et. al 2013, J. Neurophys
    % -------
    % Inputs: 
    %     binSize      - in seconds, time window within which all spikes are considered
    %     type         - {'display', 'compute'}
    %     cluId        - units to consider
    %     dim          - dimension of the space to decode
    %     state        - {'SWS', 'RUN'}
    %     IF_REPORTFIG 
    %     ratemaps     - nXBins-by-nYBins-by-nClus
    %     binOverlap   - fraction of time bin overlap for moving window
    %     markerNo     - marker to use for computing error 
    % -------
    % Outputs:
    %    decodedPos 
    %    err        - Euclidean distace between decoded position and
    %                 the marker position

    if isempty(gt.pfObject), gt.LoadPF; end
    defCluId = gt.pfObject.acceptedUnits;
    defRateMaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, defCluId));
    switch gt.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
    
    [binSize, type, cluId, dim, state, IF_REPORTFIG, ratemaps,  binOverlap, markerNo] = ...
        DefaultArgs(varargin, {200e-3, 'display', defCluId, 2, 'RUN', false, defRateMaps, 0, defMarker});

    cluId = cluId(ismember(cluId, gt.pfObject.acceptedUnits));
    switch type
      case 'compute'
        if isempty(gt.res), gt.LoadCR; end
        if isempty(gt.pfObject), gt.LoadPF; end
        fprintf('\n computing instantaneous firing rate... ')
        [sc, winEdges, xyInWin] = GetSpikeCounts(gt, binSize,  state, cluId, binOverlap);
        switch dim
          case 2
            ratemaps = gt.pfObject.smoothRateMap(:, :, ismember(gt.pfObject.acceptedUnits, cluId));
          case 1
            ratemaps = cell2mat(Compute1DRateMap(gt, 0))';
            ratemaps = ratemaps(:, ismember(gt.pyrCluIdx, cluId));
        end
        keyboard;
        fprintf('\n compution posterior... ');
        options = struct('prior', [], 'bins', size(sc, 2), 'alpha', 1);
        tic, posterior = decode_bayesian_poisson(ratemaps, sc, options);toc
        [decodedPos, dbin] = decodedPosMAP(gt, posterior);
        switch dim
          case 2 

          case 1
            [~, linearXY] = cart2pol(xyInWin(:,1), xyInWin(:,2));
            decodedPos = linearXY(dbin);
            xyInWin = linearXY;
        end
        err = vnorm(xyInWin - decodedPos, 2);        
        posterior = sparse(posterior);
        pars.binSize = binSize;
        pars.binOverlap = binOverlap;
        pars.clu = cluId;
        keyboard;   
        save([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' num2str(dim), '.' mfilename, '.mat'], 'posterior', 'err', 'pars', 'decodedPos');
   
      case 'display'
        load([gt.paths.analysis, gt.filebase, '.', gt.trialName, '.' num2str(dim), '.' mfilename, '.mat'], 'posterior', 'err', 'pars', 'decodedPos');
        posterior = full(posterior);
        figure;
        imagesc(sum(posterior, 3));
        set(gca, 'YDir', 'normal');
        % set(gca,  'YTickLabel', get(gca, 'XTick') .* par.winSize);
        title(['#cells ' num2str(length(pars.clu))  '  winSize:'  num2str(pars.winSize)])
        figure;
        hist(err, 1e2)
        axis tight
        hold on
        line([mean(err), mean(err)], ylim,'color', 'r');
        title(['error distr   #cells ' num2str(length(pars.clu)) '  winSize:'  num2str(pars.winSize) ' s'])
        xlabel('cm')
    end
end