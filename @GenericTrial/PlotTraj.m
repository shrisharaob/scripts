function varargout = PlotTraj(trial,  varargin)
    % plots trajectories in the specified periods 
    % PlotTraj(trial)
    %[trajPeriods,markerNo, plotColor, IF_PLOT_ALL, state, IF_PLOT]
    [trajPeriods,markerNo, plotColor, IF_PLOT_ALL, state, IF_PLOT] = DefaultArgs(varargin, {[], 7, [],0, 'walk', 1});
    IF_PLOT_ALL = isempty(trajPeriods);
    if IF_PLOT
        figure(gcf);
        hold on;
    end
    varargout{1}=[];
    if IF_PLOT_ALL
        ss =trial.Bhv.getState('walk').state;
        ss = ss + 1; % bug -  xyz index starts from 0
        trajectory = [];
        hold on
        cm = GenColormap(length(ss));
        IS_AUTO_COLOR =0;
        if isempty(plotColor), IS_AUTO_COLOR = 1; end
        for  i = 1 : size(ss,1)
            if IS_AUTO_COLOR, plotColor = cm(i, :);, end
            trajectory = sq(trial.xyz([ss(i,1) : ss(i,2)], markerNo ,[1,2]));
            if IF_PLOT
                plot(trajectory(:,1) , trajectory(:,2), 'Color', plotColor);
            end
            traj{i} = trajectory;
        end
        varargout{1} = traj;
        
    end

    if ~isempty(trajPeriods)
        nPeriods = size(trajPeriods, 1);
        cm = GenColormap(nPeriods);
        IS_AUTO_COLOR =0;
        if isempty(plotColor), IS_AUTO_COLOR = 1; end
        for mPeriod = 1 : nPeriods
            if IS_AUTO_COLOR, plotColor = cm(mPeriod, :); end
            mTraj = sq(trial.xyz(trajPeriods(mPeriod, 1) : trajPeriods(mPeriod,2), markerNo, [1,2]));
                if IF_PLOT
                    plot(mTraj(:, 1), mTraj(:, 2), 'Color', plotColor);
                end
            %         xlim([min(trial.xbin), max(trial.xbin)]);
            %         ylim([min(trial.ybin), max(trial.ybin)]);
        end
        varargout{1} = mTraj;
    end
   
end