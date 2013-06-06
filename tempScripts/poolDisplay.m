h = figure;
set(h, 'Position', get(0,'Screensize'));
kkBs = 0;
mmTr = 0;
for iiFb = 1 : length(out)  
    iiFb
    for jjTr = 1 : sum(~cellfun(@isempty, out(iiFb,:))) 
        subplot(1, sum(~cellfun(@isempty, out(iiFb,:))) , jjTr)
        grid on; hold on;
        mmTr = mmTr + 1;
        title(trNames{mmTr});
        cntrVertices = out{iiFb, jjTr}.cntrVertices;
        cntrPeaks = out{iiFb, jjTr}.cntrPeals;
        cc = GenColormap(length(cntrVertices)); 
        for cel = 1 : length(cntrVertices)
            icntr = cntrVertices{cel};
            ipk = cntrPeaks{cel};
            for kkCntr = 1 : length(icntr)
                plot(icntr{kkCntr}(:,1), icntr{kkCntr}(:,2),'.-', 'Color', cc(cel, :));
                hold on;
                plot(ipk(kkCntr, 1), ipk(kkCntr, 2), 'k*')
                xlim([1 50]);
                ylim([1 50]);
                axis square;
            end
        end
    end
kkBs = kkBs + 1;
filebases{kkBs}
 reportfig(h, ['pooled_cntrs_CA3_square'], 0, [filebases{kkBs}])
keyboard;
clf(h);
end