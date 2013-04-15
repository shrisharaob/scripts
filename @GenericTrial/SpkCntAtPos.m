    function [spkCnt, varargout] = SpkCntAtPos(Session,spiket,pos,varargin)
    % Justin's code

    [Nbin,Smooth,type] = DefaultArgs(varargin,{50,[],'xy'});
    pos(isnan(pos)) = 0;

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

    end

    %% rounded position
    X = round((pos(:,1)-Xmin)*k(1))+1;
    Y = round((pos(:,2)-Ymin)*k(2))+1;

    %% Push back in any stray bins
    X(X>Nbin) = Nbin;
    Y(Y>Nbin) = Nbin;

    %% spike count
    spkCnt = zeros(msize);
    spiket = spiket;
    indx = spiket(find(spiket>0 & spiket<=length(X)));
   if ~isempty(indx)
        spikep(:,1) = X(indx);
        spikep(:,2) = Y(indx);
        spkCnt = Accumulate(spikep,1,msize);
   end
   spkCnt = spkCnt(:);
end