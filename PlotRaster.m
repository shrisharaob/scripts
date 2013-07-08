function PlotRaster(Res,Clu,varargin)
%[SampleRate, Color, linestyle, cluOrder, ONLY_FIRST_SPK]
[SampleRate, Color, linestyle, cluOrder, ONLY_FIRST_SPK] = DefaultArgs(varargin,{20000, 'k', '-',unique(Clu), 0});

Colors = colormap;
Colors = repmat(Colors,10,1);
%MyClu = unique(Clu, 'stable');
MyClu = cluOrder;

for c=1:length(MyClu)
    if ONLY_FIRST_SPK
        MyRes = Res(find(Clu == MyClu(c), 1, 'first'));
    else
        MyRes = Res(find(Clu==MyClu(c)));
    end
    if isempty(Color)
        plot(repmat(MyRes(:)',2,1)/SampleRate, repmat(c+[0;.6],1,length(MyRes)), 'color', Colors(c,:), 'linestyle', linestyle);
    else
        plot(repmat(MyRes(:)',2,1)/SampleRate, repmat(c+[0;.6],1,length(MyRes)), 'color', Color);
    end
    hold on
end    

set(gca, 'ytick', 1.45:length(MyClu)+.45);
set(gca, 'yticklabel', MyClu);
