function linPos = LinearizePos(gt, varargin);
    % linPos = LinearizePos(gt, varargin);
    % returns linear position
    % -------
    % Inputs:
    %    markerNo
    %    state     

   switch gt.datasetType
      case 'kenji'
        defMarker = 1;
      case 'MTA'
        defMarker = 7;
    end
    
    [markerNo, state] = DefaultArgs(varargin, {defMarker, []});
   
    [~, linPos] = cart2pol(gt.position(:, markerNo, 1), gt.position(:, markerNo, 2));
    linPos = linPos - min(linPos);
    
end
