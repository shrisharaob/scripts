function [pos, predErr] = DecodePos(gt, ratemaps, cluId, varargin)
% decode position 
% gt - GenericTrial Object
    [binSize, IF_REPORTFIG, type] = DefaultArgs(varargin, {200e-3, 0, 'display'});
    [sc,bc] = GetSpikeCounts(gt, binSize, cluId );
    switch gt.datasetType
      case 'kenji'
        markerNo = 1;
      case 'MTA'
        marekerNo = 7;
    end
    fprintf('compution posterior... ');
    tic, posterior = decode_bayesian_poisson(ratemaps, sc);toc
    pos = decodedPosMAP(gt, posterior);
    % figure;
    % plot(pos(:,1), pos(:,2),'.');
    % figure;
    % for kBin = 1 : size(posterior,3)
    %   imagesc(sq(posterior(:,:,kBin)));
    %   pause(.05);
    %   drawnow;
    % end
    bcidx = round(bc * gt.trackingSampleRate)+1;
    bcidx(end) = [];
    xy = gt.position(:, markerNo , [1, 2]);
    xy = xy(bcidx,:);
    for kBin = 1 : size(pos,1) - 1
        predErr(kBin) = norm(xy(kBin,:) - pos(kBin,:));
    end
    h = figure;
    plot(predErr);
    title([num2str(binSize * 1e3), 'ms bins']);
end