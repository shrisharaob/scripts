function occupancy = Occupancy(trial,varargin)
    % occupancy = Occupancy(trial,varargin)
    % Justin's code

    [Nbin,Smooth,type, states, pos_shuffle] = DefaultArgs(varargin,{50,[],'xy', {'head', 'theta'}, 0});
    
    [stsp stateLabel] = trial.statePeriods(states);
    stsp = round(stsp./trial.lfpSampleRate.*trial.xyzSampleRate)+1;
    stspos = SelectPeriods(sq(trial.xyz(:,trial.Model.gmi(trial.trackingMarker),:)),stsp, 'c', 1);
    pos = stspos;
%     shuffled_Pos = @(pos_shuffle,stspos) circshift(stspos,randi([-pos_shuffle,pos_shuffle]));
    %%
    Xmin = trial.Maze.boundaries(1,1);
    Xmax = trial.Maze.boundaries(1,2);
    Ymin = trial.Maze.boundaries(2,1);
    Ymax = trial.Maze.boundaries(2,2);

    %% scaling factor for rounding position
    dx = Xmax - Xmin; 
    dy = Ymax - Ymin; 
    k = [Nbin/dx Nbin/dy];

    %% matrix size
    msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)]);

    Bin1 = ([1:msize(1)]-1)/k(1) + Xmin+round(k(1)^-1/2);
    Bin2 = ([1:msize(2)]-1)/k(2) + Ymin+round(k(2)^-1/2);
    
    %% rounded position
X = round((pos(:,1)-Xmin)*k(1))+1;
Y = round((pos(:,2)-Ymin)*k(2))+1;

%% Push back in any stray bins
X(X>Nbin) = Nbin;
Y(Y>Nbin) = Nbin;


%% Occupancy
occupancy = Accumulate([X Y],1,msize); %./trial.xyzSampleRate;
occupancy = occupancy ./ sum(occupancy(:));
end