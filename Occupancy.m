function occupancy = Occupancy(trial,varargin)
    % occupancy = Occupancy(trial,varargin)
    % Justin's code

    [Nbin,Smooth,type, states, pos_shuffle] = DefaultArgs(varargin,{50, 0.03,'xy', {'head', 'theta'}, 0});
    
%     [stsp stateLabel] = trial.trialPeriods(states);
    stsp = trial.trialPeriods;
    stsp = round(stsp./trial.lfpSampleRate.*trial.trackingSampleRate)+1;
    if strcmp(trial.datasetType, 'MTA')
        markerNo = 7;
    else 
        markerNo = 1;
    end
    if strcmp(trial.datasetType, 'kenji')
            stsp = [1, diff(stsp)];
    end
    stspos = SelectPeriods(sq(trial.position(:,markerNo,:)),stsp, 'c', 1);
    pos = stspos;
%     shuffled_Pos = @(pos_shuffle,stspos) circshift(stspos,randi([-pos_shuffle,pos_shuffle]));
    %%
    Xmin = trial.maze.boundaries(1,1);
    Xmax = trial.maze.boundaries(1,2);
    Ymin = trial.maze.boundaries(2,1);
    Ymax = trial.maze.boundaries(2,2);

    %% scaling factor for rounding position
    dx = Xmax - Xmin; 
    dy = Ymax - Ymin; 
    k = [Nbin/dx Nbin/dy];

    %% matrix size
%     msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)]); 
    msize = [Nbin, Nbin];
    Bin1 = ([1:msize(1)]-1)/k(1) + Xmin+round(k(1)^-1/2);
    Bin2 = ([1:msize(2)]-1)/k(2) + Ymin+round(k(2)^-1/2);
    
    %% rounded position
X = round((pos(:,1)-Xmin)*k(1))+1;
Y = round((pos(:,2)-Ymin)*k(2))+1;

%% Push back in any stray bins
X(X>Nbin) = Nbin;
Y(Y>Nbin) = Nbin;
X(isnan(X)) =[];
Y(isnan(Y)) = [];

%% Occupancy
occupancy = Accumulate([X Y],1,msize); %./trial.xyzSampleRate;
occupancy = occupancy ./ sum(occupancy(:));

r1 = (-msize(1):msize(1))/msize(1);
r2 = (-msize(2):msize(2))/msize(2);
Smoother1 = exp(-r1.^2/Smooth^2/2);
Smoother1 = Smoother1/sum(Smoother1);
Smoother2 = exp(-r2.^2/Smooth^2/2);
Smoother2 = Smoother2/sum(Smoother2);
occupancy = conv2(Smoother1,Smoother2,occupancy,'same');


end