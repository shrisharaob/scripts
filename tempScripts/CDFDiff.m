function CDFDiff(gt, varargin)

    [prePost, minCellsInSeq] = DefaultArgs(varargin, {'pre', 5});

    temp = gt.TemplateMatch(0, minCellsInSeq, prePost);
    switch prePost
      case 'pre'
        tm.EvntCorrs = temp.preEvntCorrs;
        tm.Surrogate = temp.preSurrogate;
      case 'post'
        tm.EvntCorrs = temp.postEvntCorrs;
        tm.Surrogate = temp.postSurrogate;
    end
    [dataCDF, xAx] = ecdf(tm.EvntCorrs);
    surCorrs = reshape(tm.Surrogate, 1e3, []);
    [nResample, nEvents] = size(surCorrs);
    shuffleCDF = nan(length(unique(xAx)), nResample);
    for ii = 1 : nResample
        x = surCorrs(ii, :);
        [Fi,xi] = ecdf(x);
        xj = xi(2:end);
        n = length(xj);
        Fj = (Fi(1:end-1)+Fi(2:end))/2;
        xj = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1))); xj; ...
              xj(n)+(1-Fj(n))*((xj(n)-xj(n-1))/(Fj(n)-Fj(n-1)))];
        Fj = [0; Fj; 1];
        [xj, xjidx] = unique(xj);
        F = @(y) interp1(xj,Fj(xjidx),y,'linear','extrap');
        sy = F(xAx(2:end));
        sy(sy > 1) = 1;
        sy(sy < 0) = 0;
        shuffCDF(:, ii) = sy;
    end

    %% 
    figure;
    for kk = 1 : nResample
        %        sHdl = stairs(xAx(2:end), shuffCDF(:, kk), 'Color', [.85, .85, .85]);
        sHdl = plot(xAx(2:end), shuffCDF(:, kk), '.','Color', [.85, .85, .85]);
        hold on
    end
    dHdl =  stairs(xAx(2:end), dataCDF(2:end), 'r');
    mHdl = stairs(xAx(2:end), median(shuffCDF, 2), 'c');
    legend([ dHdl, sHdl, mHdl], {'data', 'shuffle', 'median'}, 'Location', 'NorthWest');
    
    %% pval
    %    be = linspace(0, 1, 3e3); 
    for kk = 1 : length(xAx) - 1
        %        be(:, kk) = linspace(min(shuffCDF(kk, :)), max(shuffCDF(kk, :)), 30);
        cc(:, kk) = histc(shuffCDF(kk, :), be);
        pVal(kk) = sum(cc(be <= dataCDF(kk), kk)) ./ sum(cc(:, kk));
    end
    pVal(pVal == 0) = nan;
    figure;
    h1 = plot(xAx(2:end), pVal, '*-')
    xlabel('Correlation Value');
    ylabel('p-vallue');
    grid on;
    linHdl1 = line(xlim, [.025, .025], 'Color', 'g');
    linHdl2 = line(xlim, [.05, .05], 'Color', 'r');
    legend([linHdl2, linHdl1], {'0.05', '0.025'}, 'location', 'NorthWest');

    %% 
    figure;
    clear cc be
    for kk = 1 : length(xAx) - 1
        be(:, kk) = linspace(min(shuffCDF(kk, :)), max(shuffCDF(kk, :)), 30);
        cc(:, kk) = histc(shuffCDF(kk, :), be(:, kk));
        subplot(10, 10, kk)
        bar(be(:, kk), cc(:, kk) ./ sum(cc(:, kk)))
        hold on
        line([dataCDF(kk+1), dataCDF(kk + 1)], ylim, 'color', 'r')
        pVal(kk) = sum(cc(be(:, kk) <= dataCDF(kk), kk)) ./ sum(cc(:, kk));
    end

    % count density
    figure;
    imagesc(xAx(2:end), be, cc');
    set(gca, 'ydir', 'normal')
    hold on;
    stairs(xAx(2:end), dataCDF(2:end), 'w')

    %% FLIPPED
    [dataCDF, xAx] = ecdf(-tm.EvntCorrs);
    surCorrs = reshape(-tm.Surrogate, 1e3, []);
    [nResample, nEvents] = size(surCorrs);
    shuffleCDF = nan(length(unique(xAx)), nResample);
    for ii = 1 : nResample
        x = surCorrs(ii, :);
        [Fi,xi] = ecdf(x);
        xj = xi(2:end);
        n = length(xj);
        Fj = (Fi(1:end-1)+Fi(2:end))/2;
        xj = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1))); xj; ...
              xj(n)+(1-Fj(n))*((xj(n)-xj(n-1))/(Fj(n)-Fj(n-1)))];
        Fj = [0; Fj; 1];
        [xj, xjidx] = unique(xj);
        F = @(y) interp1(xj,Fj(xjidx),y,'linear','extrap');
        sy = F(xAx(2:end));
        sy(sy > 1) = 1;
        sy(sy < 0) = 0;
        shuffCDF(:, ii) = sy;
    end

    figure;
    for kk = 1 : nResample
        stairs(xAx(2:end), shuffCDF(:, kk), 'Color', [.85, .85, .85]);
        hold on
    end
    stairs(xAx(2:end), dataCDF(2:end), 'r');
    set(gca, 'XTickLabel', flipud(str2num(get(gca, 'XTickLabel'))))
    title('flipped')

    be = linspace(0, 1, 3e3); 
    for kk = 1 : length(xAx) - 1
        cc(:, kk) = histc(shuffCDF(kk, :), be);
        pVal(kk) = sum(cc(be <= dataCDF(kk), kk)) ./ sum(cc(:, kk));
    end

    figure;
    % pval
    pVal(pVal == 0) = nan;
    h1 = plot(xAx(2:end), pVal, '*-')
    xlabel('Correlation Value');
    ylabel('p-vallue');
    grid on;
    linHdl1 = line(xlim, [.025, .025], 'Color', 'g');
    linHdl2 = line(xlim, [.05, .05], 'Color', 'r');
    legend([linHdl2, linHdl1], {'0.05', '0.025'}, 'location', 'NorthWest');
    title('flipped');
    clear cc be

    %%
    figure;
    for kk = 1 : length(xAx) - 1
        be(:, kk) = linspace(min(shuffCDF(kk, :)), max(shuffCDF(kk, :)), 30);
        cc(:, kk) = histc(shuffCDF(kk, :), be(:, kk));
        subplot(10, 10, kk)
        bar(be(:, kk), cc(:, kk) ./ sum(cc(:, kk)))
        hold on
        line([dataCDF(kk+1), dataCDF(kk + 1)], ylim, 'color', 'r')
    end

    figure;
    imagesc(xAx(2:end), be, cc');
    set(gca, 'ydir', 'normal')
    hold on;
    stairs(xAx(2:end), dataCDF(2:end), 'w')


end