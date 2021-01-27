function subplot(nrows, ncols, c, FontSize, FontName)
    if nargin < 1
        nrows = nan;
    end
    if nargin < 2
        ncols = nan;
    end
    if nargin < 3
        c = 1;
    end
    if nargin < 4
        FontSize = 18;
    end
    if nargin < 5
        FontName = 'Helvetica';
    end
    if ~isnan(nrows) && ~isnan(ncols)
        subplot(nrows, ncols, c);
    end
    hold on;
    set(gca, 'FontSize', FontSize); 
    set(gca, 'FontName', FontName); % default anyway, but just for fun
    set(gca, 'TickDir', 'out');
    set(gca, 'XColor', 'k', 'YColor', 'k');
end
