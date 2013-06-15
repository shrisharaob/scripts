function CopyToSubplots(figHdlList, subplotDims)

    h = figure;
    for ii = 1 : length(list)
        ax = subplot(subplotDims(1), subplotDims(2), ii);
        copyobj(allchild(get(ll(ii), 'children')), ax);
    end
end