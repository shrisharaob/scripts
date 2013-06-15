winSiz = 1000e-3;
%rm = gt.pfObject.smoothRateMap(:,:,ismember(gt.pfObject.acceptedUnits, clus));
rm = cell2mat(Compute1DRateMap(gt, 0))';
rm = rm(:, ismember(1:81, clus));
keyboard;
[sc,winEdges, posinwin] = GetSpikeCounts(gt,  winSiz, 'RUN', clus);
posterior = decode_bayesian_poisson(rm, sc);

keyboard;
figure
imagesc(sum(posterior, 3));
title(['#cells ' num2str(length(clus))  '  winSize:'  num2str(winSiz)])

keyboard;
[dmap, dbin] = decodedPosMAP(gt, posterior);
err = vnorm(posinwin - dmap, 2);
figure
hist(err, 1e2)
hold on
line([mean(err), mean(err)], ylim,'color', 'r');
title(['error distr   #cells ' num2str(length(clus))  '  winSize:'  num2str(winSiz) ' s'])
xlabel('cm')
