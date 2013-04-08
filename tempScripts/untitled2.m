        bn = 1;
for i = 1 : 23
    if ~isempty(rr{i})

        for kk = 1 : length(rr{i})
            kkr = rr{i};
            rTr = cell2mat(kkr(kk));
            if ~isempty(rTr)
                plot(bn, rTr,'o');
                hold on;
               
            end
        end
         bn = bn + 1;
    end
end
grid on;
line([xlim], [0,0], 'color', 'k');
ylim([-1 1]) ;
ylabel('corr coef, z-trans');
xlabel('filebases with CA1 recording');
title('CA1 remapping');