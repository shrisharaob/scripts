function out = SelectCells(trial, varargin)
    % this script automatically classifies pyramidal cells and interneurons
    % based on spike waveform features and firing rate
    % fits MOG in 3d feature space, could be optimised by adding additional 
    % relavent features

    filebase = trial.filebase;
    %     session = MTASession(filebase,{{'CluRes', 32552},'Pfs'},'cof');
    %     fullpath = ['/data/homes/shrisha/data/nlx/jg05-20120315/' filebase];
    analysisFolderPath = '~/data/analysis/';
    par = LoadPar([trial.paths.data, trial.filebase]);
    [LoadData, Arena, IF_PLOT] = DefaultArgs(varargin, {{'CluRes','Pfs'}, 'cof', 0});
%     session = MTASession(filebase, LoadData, Arena);
    if FileExists([trial.paths.analysis, '/', filebase '.NeuronQuality.mat'])
        load([trial.paths.analysis, '/', filebase '.NeuronQuality.mat']);
    else
        fprintf('\n executing NeuronQuality.m ..... \n ')
        nq = NeuronQuality([trial.paths.data, trial.filebase]);
    end
    %% classify pyr and iterneurons
    %% features
    selectFeatures = {'AmpSym', 'SpkWidthR', 'FirRate'};
    f1 = [];
    f2 = [];
    f3 = [];
    acceptedElClu = [];
    for i = 1 : par.nElecGps
        f1 = [f1; nq(i).(genvarname(selectFeatures{1}))];
        f2 = [f2; nq(i).(genvarname(selectFeatures{2}))];
        f3 = [f3; nq(i).(genvarname(selectFeatures{3}))];
    end
    cnq = CatStruct(nq);
    acceptedElClu = unique([cnq.ElNum,cnq.Clus],'rows');
    linearCluIdx = 1:size(acceptedElClu,1);
    features = [f1, f2, f3]; % nSamples x nDims
    clear f1 f2 f3;
    %% MOG it
    nFittingIterations = 10;
    nClusters = 2;
    logLikOld = -inf;
    for km = 1 : nFittingIterations
        fprintf(1,'\n ************************************************************');
        fprintf(1,'\n Fitting MOG m = \b\b\%2d \n',km);
        MOGout = MOG(features, nClusters);
        logLik(:,km) = MOGout.logLik;
        if(MOGout.logLik(end, km) > logLikOld)
            mu = MOGout.mu;
            sigma= MOGout.sigma;
            pk = MOGout.pk;
            logLikOld = logLik(end,km);
        end
    end
    fprintf(1,'\n DONE !!!');

    out.mu = mu;
    out.sigma = sigma;
    out.pk = pk;
    out.logLik = logLik;
    out.cluIndx = MOGout.cluIndx;
    clear MOGout;
    MOGout = out;
    clear out;
%     save(['~/thesis/analysis/' filebase '.MOG.mat'],'MOGout');
    %% fisher lda
    fLDAOut = FisherLDA(out.mu', out.sigma, features,  out.cluIndx);
    %% figures
    if IF_PLOT
        hold on;
        figure;
        colors = MakeColorMap(nClusters);
        plot(features(out.cluIndx==1,1),features(out.cluIndx==1,2),'.g',...
            features(out.cluIndx==2,1),features(out.cluIndx==2,2),'r.')
        figure;
        plot3(features(out.cluIndx==1,1),features(out.cluIndx==1,2),features(out.cluIndx==1,3),'.g',...
            features(out.cluIndx==2,1),features(out.cluIndx==2,2),features(out.cluIndx==2,3),'r.')

        xLimits = [min(features(:,1)), max(features(:,1))];
        yLimits = [min(features(:,2)), max(features(:,2))];
        figure;
        Plot2DGaussians(mu,sigma,xLimits,yLimits,pk);
        hold on
        for k = 1 : 3
            plot(mu(1,k),mu(2,k),'*','MarkerSize',10,'MarkerEdgeColor',[1 1 1]);
        end
        figure;
        plot(logLik);
    end
    %% identify pyramidal cells
    % cluindx 2 are the pyr cells
    pyrCellClusterIdx = acceptedElClu(MOGout.cluIndx == 1,:);
    linearPyrCluIdx = linearCluIdx(ismember(acceptedElClu,pyrCellClusterIdx,'rows'));
    save([trial.paths.analysis, '/', trial.filebase '.' mfilename '.mat'],'linearPyrCluIdx','acceptedElClu');
    nPairs = nchoosek(length(linearPyrCluIdx), 2);
    cellPairs = nchoosek(linearPyrCluIdx, 2);

