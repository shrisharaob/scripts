function out =MOG(x, K)
% x nSamples x nDims
     x = x'; % makecolumns of X as individual obs i.e. x: nDims x nSamples
    [nDims, nSamples] = size(x);
%% initialize
fprintf(1,'\n Running K-means for initialization \n');
%     [muK, ~] = Kmean(x',K,100);
    [cluId, muK] = kmeans(x', 2);
    muK = muK'; %  nDim x K
% muK = rand(2,K);
    sigma = repmat(eye(nDims), [1,1,K]); % nDims x nDims x K
    fractionOfPoints = repmat(1/K, K, 1); % K x 1
    gammaNK = nan(nSamples, K);
    logLik = [];
    runNo = 1;
    IS_CONVERGED = 0;
    convrgTolerence = 0.001;
    maxIter = 10;
    
    %%
%      x = [randn(1,250) , randn(1,250)];
%     x=[x;randn(1,250) + 10 , randn(1,250)+10];
%     plot(x(1,:),x(2,:),'.');
%     drawnow;
%     hold on;
%     for k = 1 : K
%         plot(muK(1,k),muK(2,k),'ro','MarkerSize',(20 - runNo)*((20 -runNo) > 0) + 2);
%         drawnow;
%     end
%% E-step
fprintf(1,'\n ----------------------------------------------------------------------');
fprintf(1,'\n MOG iteration ''#'' ');
while(~ IS_CONVERGED && runNo <= maxIter)
    logLikOld =  LogLiklihood(x, muK, sigma, fractionOfPoints);
    oldMu = muK;
    logLik(runNo) = logLikOld;
    for n = 1 : nSamples
        for k = 1 : K
            denom = 0;
            for lk = 1 : K
                denom = denom + fractionOfPoints(k) * mvnpdf(x(:, n), muK(:,lk), sigma(:,:,lk));
            end
            gammaNK(n, k) = fractionOfPoints(k) * mvnpdf(x(:, n), muK(:,k), sigma(:,:,k))/denom;
        end
    end
  
%% M-step
    Nk = sum(gammaNK, 1)';
    for k = 1 : K
        tempMuK = 0;
        tempSigma = zeros(nDims, nDims);
        for n = 1 : nSamples
            tempMuK = tempMuK + gammaNK(n,k) * x(:, n);
            tempSigma = tempSigma + gammaNK(n,k) * (x(:, n) - muK(:,k)) * (x(:, n) - muK(:,k))';
        end
        
        muK(:,k) = tempMuK/Nk(k);
        sigma(:,:,k) = tempSigma/Nk(k) + 1e-5*eye(nDims, nDims); % add prior for numerical stability
        fractionOfPoints(k) = Nk(k)/nSamples;

    end
    logLikNew = LogLiklihood(x, muK, sigma, fractionOfPoints);
    %%
    fprintf(1,'\b\b%2d',runNo);
%     pause(1);

    for k = 1 : K
        plot(muK(1,k),muK(2,k),'ro','MarkerSize',(20 - runNo)*((20 -runNo) > 0) + 2);
        drawnow;
    end
    %%
    if ((abs(logLikOld - logLikNew) < convrgTolerence) ); %|| sum((abs(oldMu(:) - muK(:)) < convrgTolerence)))
        IS_CONVERGED = 1;
    end
    runNo = runNo + 1;
    
end

clusterIndx = gammaNK(:,1) > gammaNK(:,2);
clusterIndx = clusterIndx + 1; % start numbering clusteris from 1
% figure;
% plot(logLik);
out.mu = muK;
out.sigma = sigma;
out.pk = fractionOfPoints;
out.logLik = logLik;
out.cluIndx = clusterIndx;
end


function logLik = LogLiklihood(x, mu, sigma, piK)
% x : nDims x nSamples
% mu nDim x K
% sigma nDim x nDim x K

    N = length(x);
    [~, K] = size(mu);
    temp = 0;
    logTemp = 0;
    for n = 1 : N
        for k = 1 : K
            temp = temp + piK(k) * mvnpdf(x(:,n), mu(:,k), sigma(:,:,k));
        end
        logTemp = logTemp + log(temp);
    end
    logLik = logTemp;
end