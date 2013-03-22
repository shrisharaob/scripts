function [mu, eCov, piMix, logLike, resp] = MOG(x, kMax)
% x(ptData, kDimention) - data
    % kMax - max number of clusters splitting/merging iterations
    %
    % This code is mainly based on the book of Bishop - "Machine learning..."
    % and the paper of Blekas and Lagaris (2007) "Split–Merge Incremental
    % LEarning (SMILE) of Mixture Models"
    
    K = 2; % starting number of clusters
    [~, nDims] = size(x);
    
    % Init with k-means
    [idx, mu] = kmeans(x, K);
    mu = mu';
    piMix = ones(1, K) / K;
    eCov = zeros(size(x, 2), size(x, 2), K);
    for k = 1:K
        eCov(:, :, k) = cov(x(idx == k, :)) + 1e-5*eye(nDims, nDims);
    end
    
    fprintf('\n Initial EM fitting:');
    [piMix, mu, eCov, resp, logLike] = EM(x, piMix, mu, eCov);
    
    for k = 1:kMax
        bSplitted = 0;
        
        % Splitting step
        divergence = LocalKullbackDivergences(piMix, mu, eCov, x);
        [~, divIndexes] = sort(divergence);
        j = K;
        while (j > 0)
            splitIndex = divIndexes(j);
            clustersLeft = [1:splitIndex-1, splitIndex+1:K, K, K];
            piMixNew = piMix(clustersLeft);
            muNew = mu(:, clustersLeft);
            eCovNew = eCov(:, :, clustersLeft);
            piMixNew(end-1) = 0.5 * piMix(splitIndex);
            piMixNew(end) = 0.5 * piMix(splitIndex);
            [vec, val] = eigs(eCov(:, :, splitIndex), 1);
            muNew(: , end-1) = mu(: , splitIndex) + 0.5*sqrt(val)*vec;
            muNew(: , end) = mu(: , splitIndex) - 0.5*sqrt(val)*vec;
            eCovNew(:, :, end-1) = 0.5 * eCov(:, :, splitIndex);
            eCovNew(:, :, end) = 0.5 * eCov(:, :, splitIndex);
            fprintf('\n Trying to split cluster %d:', splitIndex);
            [piMixNew, muNew, eCovNew, respNew, logLikeNew] = EM(x, piMixNew, muNew, eCovNew);
            bicOld = BIC(logLike, mu, eCov, x);
            bicNew = BIC(logLikeNew, muNew, eCovNew, x);
            if (bicNew < bicOld)
                piMix = piMixNew;
                mu = muNew;
                eCov = eCovNew;
                resp = respNew;
                logLike = logLikeNew;
                K = K + 1;
                bSplitted = 1;
                divIndexes = divIndexes - (divIndexes > splitIndex);
                fprintf('OK, deltaBIC = %f, K = %d', bicNew - bicOld, K);
            else
                fprintf('aborted, deltaBIC = %f, K = %d', bicNew - bicOld, K);
            end
            j = j - 1;
        end
        
        % Merging step
        bMerged = 0;
        while (1) % break if no merging has been done
            J = DistributionDistance(mu, eCov, x);
            [mergeValL, mergeIndexL] = min(J);
            [~, k2] = min(mergeValL);
            k1 = mergeIndexL(k2);
            clustersLeft = [1:k1-1, k1+1:k2-1, k2+1:K, K];
            piMixNew = piMix(clustersLeft);
            piMixNew(end) = piMix(k1) + piMix(k2);
            muNew = mu(:, clustersLeft);
            muNew(: , end) = (piMix(k1)*mu(:, k1) + piMix(k2)*mu(:, k2)) / (piMix(k1) + piMix(k2));
            eCovNew = eCov(:, :, clustersLeft);
            eCovNew(:, :, end) = (piMix(k1)*eCov(:, :, k1) + piMix(k2)*eCov(:, :, k2)) / (piMix(k1) + piMix(k2));
            fprintf('\n Trying to merge clusters %d and %d:', k1, k2);
            [piMixNew, muNew, eCovNew, respNew, logLikeNew] = EM(x, piMixNew, muNew, eCovNew);
            bicOld = BIC(logLike, mu, eCov, x);
            bicNew = BIC(logLikeNew, muNew, eCovNew, x);
            if (bicNew < bicOld)
                piMix = piMixNew;
                mu = muNew;
                eCov = eCovNew;
                resp = respNew;
                logLike = logLikeNew;
                K = K - 1;
                bMerged = 1;
                fprintf('OK, deltaBIC = %f, K = %d', bicNew - bicOld, K);
            else
                fprintf('aborted, deltaBIC = %f, K = %d', bicNew - bicOld, K);
                break;
            end
        end
        
        if (~bSplitted && ~bMerged)
            break;
        end
    end
    fprintf('\n Finished with %d clusters! \n', K);
end

function [piMix, mu, eCov, resp, logLike] = EM(x, piMix, mu, eCov)
[N, nDims] = size(x);
    K = size(eCov, 3);
    maxIterations = 100;
    logLikeOld = -Inf;
    bIterate = 1;
    kIteration = 1;
    logLike = 0;
    resp = nan(N, K);
    
    while (kIteration < maxIterations) && (bIterate)
        % E step
        logLike = 0;
        for n = 1:N
            pt = x(n, :);
            piMixPDF = piMix' .* mvnpdf(repmat(pt, K, 1), mu', eCov);
            denom = sum(piMixPDF);
            resp(n, :) = piMixPDF / denom;
            logLike = logLike + log(sum(piMixPDF));
        end

        % M step
        Nk = sum(resp, 1);
        for k = 1:K
            mu(:, k) = 1/Nk(k) * (resp(:, k)'*x);
            eCov(:, :, k) = zeros(nDims, nDims);
            vpResp = zeros(nDims, nDims);
            for n = 1:N
                vpResp = vpResp + resp(n, k) .* ((x(n, :)' - mu(:, k)) * (x(n, :)' - mu(:, k))');
            end
            eCov(:, :, k) = vpResp;
            eCov(:, :, k) = eCov(:, :, k) / Nk(k) + 1e-5*eye(nDims, nDims);
            piMix(k) = Nk(k) / N;
        end
        
        bIterate = (abs(logLikeOld - logLike) > 0.005*abs(logLike)) || (isinf(logLike));
        logLikeOld = logLike;
        kIteration = kIteration + 1;
        fprintf('.');
    end
end

function res = BIC(logLike, mu, eCov, x)
% % -- THIS PIECE OF CODE IS ADAPTED FROM ALEXANDER ECKER'S SOLUTION
% % -- IT GIVES BETTER CLUSTERING (MORE CLUSTERS)
% nClusters = size(eCov, 3);
% nDims = size(mu, 1);
% nPoints = size(x, 1);
% nParams = nClusters * (0.5 * nDims * (nDims + 1) + nDims) + nClusters;
% res = -2 * logLike + nParams * log(nPoints);
% % -- END OF THE ADAPTED CODE
    
    % Bayessian information criterion estimation
    k = numel(mu) + numel(eCov);
    res = -2*logLike + k*log(size(x, 1));
end

function J = LocalKullbackDivergences(piMix, mu, eCov, x)
% See Blekas and Lagaris (2007)
    
    [N, ~] = size(x);
    K = size(eCov, 3);
    
    % Probabilities for every point
    f = zeros(1, N);
    Pj = zeros(K, N);
    for n = 1:N
        pt = x(n, :);
        vmGauss = mvnpdf(repmat(pt, K, 1), mu', eCov);
        f(n) = f(n) + piMix * vmGauss;
        Pj(:, n) = Pj(:, n) + vmGauss;
    end
    
    J = zeros(1, K);
    for k = 1:K
        for n = 1:N
            J(k) = J(k) + f(n)*log(f(n))/Pj(k, n);
        end
    end
end

function J = DistributionDistance(mu, eCov, x)
% See Blekas and Lagaris (2007)

    [N, ~] = size(x);
    K = size(eCov, 3);
    
    % Probabilities for every point
    Pj = zeros(K, N);
    for n = 1:N
        pt = x(n, :);
        Pj(:, n) = Pj(:, n) + mvnpdf(repmat(pt, K, 1), mu', eCov);
    end
    
    J = nan(K, K);
    for k1 = 1:K
        for k2 = k1+1:K
            s = 0;
            for n = 1:N
                s = s + Pj(k1, n)*log(Pj(k1, n))/Pj(k2, n) + Pj(k2, n)*log(Pj(k2, n))/Pj(k1, n);
            end
            J(k1, k2) = s;
        end
    end
end