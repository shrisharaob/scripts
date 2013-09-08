function Plot2DGaussians(mu, sigma, xLimits, yLimits,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mu    : nDims x nGaussians
    % sigma : nDima x nDims x nGaussians
    % xLimits and yLimits specify the plot limits
    % varargin(1) - responsibility : vector of size nGaussians
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % History : shrisha
    
    [~, nGaussians] = size(mu);
    if nargin < 5
        resp = repmat(1/nGaussians, nGaussians, 1);
    else
        resp = varargin{1};
    end

    x = xLimits(1):0.1:xLimits(2);
    y = yLimits(1):0.1:yLimits(2);
    [tX, Y] = meshgrid(x,y);
    tX = tX(:);
    Y = Y(:);
    X = [tX'; Y']';
    for i = 1 :nGaussians
        pix(:,i) = resp(i) * mvnpdf(X, mu(:,i)', sigma(:,:,i));
    end
    pix = sum(pix, 2);
    h = imagesc(x,y,reshape(pix,length(y),[]));
    set(gca, 'Ydir','Normal');
    % mesh(x,y,reshape(pix,length(y),[]));
end