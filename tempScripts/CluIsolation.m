gt = GenericTrial('ec013.961_974', 'ec013.970');
gt.LoadPF;

load /data/data/antsiro/blab/kenji/ec013.961_974.clusterquality.mat
eDist = fq.eDist(ismember(gt.pyrCluIdx, gt.pfObject.acceptedUnits));
clus = gt.pfObject.acceptedUnits;
[sedist, edistIdx] = sort(eDist);
spr = gt.pfObject.sparsity(edistIdx);
clus = clus(edistIdx);
goodClu = clus(sedist >= 20 & spr' <= 0.45);
gt.PoolRateMaps('CA3', 'bigSquare', 1, goodClu);

goodClu = gt.pfObject.acceptedUnits(fq.RefracViol(idx) < 2e-3 & gt.pfObject.sparsity' <= .4);
gt.PoolRateMaps('CA3', 'bigSquare', 1, goodClu);


