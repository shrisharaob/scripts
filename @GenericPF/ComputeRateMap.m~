    function [RateMap, varargout] = ComputeRateMap(Session,spiket,pos,varargin)
    % Justin's code

    [Nbin,Smooth,type] = DefaultArgs(varargin,{50,[],'xy'});


    %% Constraint to maze is forced
    switch type
        case 'xy'
            Xmin = Session.maze.boundaries(1,1);
            Xmax = Session.maze.boundaries(1,2);
            Ymin = Session.maze.boundaries(2,1);
            Ymax = Session.maze.boundaries(2,2);

            %% scaling factor for rounding position
            dx = Xmax - Xmin;
            dy = Ymax - Ymin;
            k = [Nbin/dx Nbin/dy];

            %% matrix size
            %     msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)]);
            msize = [Nbin, Nbin];
            Bin1 = ([1:msize(1)]-1)/k(1) + Xmin+round(k(1)^-1/2);
            Bin2 = ([1:msize(2)]-1)/k(2) + Ymin+round(k(2)^-1/2);

        case 'pfcrz'
            switch Session.maze.shape
                case 'circle'
                    Xmin = 0;
                    Xmax = Session.maze.boundaries(1,2)-Session.maze.boundaries(1,1);
                    Ymin = Session.maze.boundaries(3,1);
                    Ymax = Session.maze.boundaries(3,2);

                    %% scaling factor for rounding position
                    dx = Xmax - Xmin;
                    dy = Ymax - Ymin;
                    k = [Nbin/dx Nbin/dy];

                    %% matrix size
                    msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)]);
                    Bin1 = ([1:msize(1)]-1)/k(1);
                    Bin2 = ([1:msize(2)]-1)/k(2);

                case 'square'
            end
    end


    %% rounded position
    X = round((pos(:,1)-Xmin)*k(1))+1;
    Y = round((pos(:,2)-Ymin)*k(2))+1;

    %% Push back in any stray bins
    X(X>Nbin) = Nbin;
    Y(Y>Nbin) = Nbin;


    %% Occupancy
    Occ = Accumulate([X Y],1,msize)./Session.trackingSampleRate;

    %% spike count
    spiket = spiket;
    indx = spiket(find(spiket>0 & spiket<=length(X)));
    spikep(:,1) = X(indx);
    spikep(:,2) = Y(indx);
    Count = Accumulate(spikep,1,msize);

    %% smooth
    if isempty(Smooth)
        Smooth = Nbin/3000;
    end

    r1 = (-msize(1):msize(1))/msize(1);
    r2 = (-msize(2):msize(2))/msize(2);
    Smoother1 = exp(-r1.^2/Smooth^2/2);
    Smoother1 = Smoother1/sum(Smoother1);
    Smoother2 = exp(-r2.^2/Smooth^2/2);
    Smoother2 = Smoother2/sum(Smoother2);
    SCount = conv2(Smoother1,Smoother2,Count,'same');
    SOcc = conv2(Smoother1,Smoother2,Occ,'same');
    RateMap = NaN(Nbin,Nbin);
    OccThresh = 0.06;
    RateMap(SOcc>OccThresh) = SCount(SOcc>OccThresh)./SOcc(SOcc>OccThresh);
    RateMap(SOcc<OccThresh) = NaN;        rateMap = RateMap;
    varargout{1} = conv2(Smoother1, Smoother2, Occ./sum(Occ(:)),'same');
    varargout{2} = Bin1;
    varargout{3} = Bin2;
end


