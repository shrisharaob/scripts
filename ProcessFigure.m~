function ProcessFigure(hFigure, fileName, varargin)
% History:
% Dmytro Velychko - Created
% Shrisha - Modified
    
    [fontSiz, paperSize] = DefaultArgs(varargin, {8, [5.5, 5.5]});
    allAxesInFigure = findall(hFigure, 'type', 'axes');
    for hAxis = allAxesInFigure'
        set(hAxis, 'box', 'off');
        set(hAxis, 'TickDir', 'out');
        set(hAxis,'FontSize', 8);
        set(get(hAxis,'XLabel'), 'FontSize', fontSize);
        set(get(hAxis,'YLabel'), 'FontSize', fontSize);
    end
    set(hFigure, 'PaperPosition', [0, 0, paperSize]);
    set(hFigure, 'PaperSize', paperSize);
    set(hFigure, 'Renderer', 'Painters');
    saveas(hFigure, fileName, 'pdf');
end