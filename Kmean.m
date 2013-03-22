% function [mu,id] = Kmean(m, K)
function [muK,rnk] = Kmean(x, K, varargin)
% x : nSamples x nDims
% m = 1;
% K = 3;

[nSamples, nDims] = size(x);
rnk = floor(rand(nSamples, 1) * K); % Cluster assignments , nSamples-by-1
muK = [ rand(K,1) * max(x(:,1)), rand(K,1) * max(x(:,2)) ]; % centers, K-by-2
plot(muK(:,1),muK(:,2),'ro','MarkerSize',20);

% colorMap = GenColormap(K);
for k = 1 : K
    plot(muK(k,1),muK(k,2),'ro','MarkerSize',20);
end
if nargin < 3
    maxIter = 100;
else
    maxIter = varargin{1};
end
IS_CONVERGED = 0;
iterNo = 0;
fprintf(1,'\n iteration:  ');
while(~ IS_CONVERGED && iterNo <= maxIter)
    %% E-step
    for n = 1 : nSamples
        rnk(n) = ArgMin(x(n,:), muK);
    end
    
    %% M-step
    oldMu = muK;
    for k = 1 : K
        kSamples = x(rnk == k,:);
        if (~ isempty(kSamples))
            muK(k,:) = mean(kSamples, 1);
        end
    end
    fprintf(1,'\b\b\b%3d',iterNo);
%     pause(1);
    plot(muK(:,1),muK(:,2),'ro','MarkerSize',(20 - iterNo)*((20 -iterNo) > 0) + 2);
        for k = 1 : K
        plot(muK(k,1),muK(k,2),'ro','MarkerSize',(20 - iterNo)*((20 -iterNo) > 0) + 2);
        drawnow;
        end
    
    iterNo = iterNo + 1;
    if(iterNo > maxIter)
        fprintf('\n K-means failed to converge in %d iterations ',maxIter)
    end
    IS_CONVERGED = sum(sum((abs(oldMu - muK)) < repmat(.001,K,nDims)));
end
if(iterNo -1  < maxIter)
fprintf(1,'\n NUMBER OF ITERATIONS = %d',iterNo);
end
end

function j = ArgMin(xi, muK)
for k = 1 : size(muK,1)
    sqNorm(k) = (xi(1) - muK(k,1)) ^ 2 + (xi(2) - muK(k,2)) ^ 2;
end
[~, j] = min(sqNorm);
end
