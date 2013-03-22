function out = FisherLDA(mu, cluSigma, features, assignments)
% mu cluster means (K-by-D)
    % cluSigma cluster covariances (D-by-D-by-K)
    % features features (N-by-D) N spikes and D features
    % assignments cluster index for each spike (1 to K) ; (N-by-K)
    % ---------------------------------------------------------------------
% History:
    % Shrisha - Created

    
    [nClusters, nDims] = size(mu);
    out = figure;
    colors = MakeColorMap(nClusters);
    
    for kCluster = 1 : nClusters -1
        for lCluster = kCluster + 1 : nClusters
        
            Sw = sum(cluSigma,3); % total within class covariance matrix
            kCluMean = mu(kCluster,:);
            lCluMean = mu(lCluster,:);
            
            w=Sw^-1 * (kCluMean' - lCluMean'); % Fischer linear discrimanant
            
            kClu = features((assignments == kCluster), :);
            lClu = features((assignments == lCluster), :);
            kCluProj = kClu * w; % project features onto the line perpendicular to the linear discrimanant
            lCluProj = lClu * w;
            
            subplot(nClusters-1 ,nClusters-1 ,(kCluster - 1) * (nClusters) + lCluster - kCluster)
            hist(kCluProj,100);
%             axis off
            hold on
            h = findobj(gca,'type','patch');
            set(h(1),'facecolor',colors(kCluster,:), 'EdgeColor', colors(kCluster,:));
    
            hist(lCluProj,100);
%             axis off
            h = findobj(gca,'type','patch');
            set(h(1),'facecolor',colors(lCluster,:), 'EdgeColor', colors(lCluster,:));
            %title( [ num2str(kCluster) ' - ' num2str(lCluster) ]);
        end
    end
    out.w = w;
%     delta.Test.Common.ProcessFigure(gcf, ['ClustersLDA'], [nClusters*3, nClusters*2]);
end