function [spikeDensity, binEdges] = SpikeDensity(res, varargin)
% spikeDensity = SpikeDensity(res, varargin)
% [sampleRate, binSize, bandWidth]

    [sampleRate, binSize, bandWidth, minSpikes, spikeDensity] = ...
        DefaultArgs(varargin, {20000, 50e-3, 30e-3, 1e2, []});

    res = res ./ sampleRate;
    if length(res) < minSpikes, return; end
    
    res = sort(res);
    binEdges = res(1) : binSize : res(end);
    nBins = length(binEdges) - 1;
    [counts, ~] = histc(res, binEdges);
    counts = counts ./ (nBins * binSize);

%    costFun = @(binSize, mu, sig) (2 * mu - sig) ./ binSize ^ 2;
%    binSize = 1e-3 : 10e-3 : 500e-3;
%    for ii = 1 : length(binSize)
%    binEdges = res(1) : binSize(ii) : res(end);
%    nBins = length(binEdges) - 1;
%    [counts, ~] = histc(res, binEdges);
%    counts = counts ./ (nBins * binSize(ii));
%    mu = mean(counts);
%    sig = var(counts);
%    c(ii) = costFun(binSize(ii), mu, sig);
%    end
%    plot(binSize, c, '*-');
    
    N = 1e2;
    gw = gausswin(N, 1 / bandWidth);
    gw = gw ./ sum(gw);
    spikeDensity = conv(counts, gw, 'same');

    %% ksd with Gaussian kernel
%    nSamples = length(res);
%    maxSamples = max(1 : nSamples + 3 * bandWidth);
%    n = 2 ^ (nextpow2(maxSamples));
%    fftCounts = fft(counts, n);
%    f = (0 : n - 1) / n;
%    f = [-f(1 : n / 2 + 1) f(n / 2 : - 1 : 2)];
%    f = f';
%    kernelFun = exp(- 0.5 * (bandWidth * 2 * pi * f) .^2);
%    spikeDensity = ifft(fftCounts .* kernelFun, n);
%    spikeDensity = spikeDensity(1 : nSamples);
%    keyboard;
end