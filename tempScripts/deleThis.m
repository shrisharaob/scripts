for i = 1 : 15
icntr = cntrVertices{i};
ipk = cntrPeaks{i};
i
for jj = 1 : length(icntr)
plot(icntr{jj}(:,1), icntr{jj}(:,2),'r.-');
hold on
plot(ipk(jj, 1), ipk(jj, 2), 'k*')
xlim([1 50]);
ylim([1 50]);
end
%waitforbuttonpress
%clf
end
