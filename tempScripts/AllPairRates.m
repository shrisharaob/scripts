function AllPairRates(varargin)


    [cluId, arena, roi, IF_REPORTFIG] = DefaultArgs(varargin, {[], {'bigSquare'}, {'CA3'}, 1});
    load(['~/data/analysis/kenji/GetGoodPFUnits', GenFiletag(arena, roi), 'mat']); % loads goodUnits.mat
    baseCount = 0;
    pkDist = cell(1, length(goodUnits));
    for kBase = 1 : length(goodUnits)
        if ~isempty(goodUnits(kBase).clu)
            fprintf(['\n', repmat('+-', 1, 10), goodUnits(kBase).filebase, repmat('+-', 1, 10), '\n']); 
            baseCount = baseCount + 1;
            if ~isempty(goodUnits(kBase).clu.goodUnits)
                outStruct = PairRemapping(goodUnits(kBase).filebase, goodUnits(kBase).clu.goodUnits, arena, roi,  0, 0);
                if ~isempty(outStruct)
                    pkDist{kBase} = outStruct.pkDist;
                else
                    fprintf('\n none');
                end
            end
        end
    end
keyboard;
end
