function smoothedSurface = SmoothSurface(p, smoothfactor)
% smooth the surface p by gaussian with std smoothfactor
% must add other smoothing windows

    [xDim, yDim] = size(p);
    xBins = [-xDim : xDim] ./ xDim /2;
    yBins = [-yDim : yDim] ./ yDim /2;
    smX = exp(-power(xBins, 2) ./ (2*smoothfactor ^2));
    smY = exp(-power(yBins, 2) ./ (2*smoothfactor ^2));
    smX = smX ./ sum(smX);
    smY = smY ./ sum(smY);
    %    p1 = p;
    p(isnan(p)) = 0;
    smoothedSurface = conv2(smX, smY, p, 'same');
end
    
