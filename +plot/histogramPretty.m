function histogramPretty(counts, bins, clr, lw)
    if nargin < 3
        clr = 'k';
    end
    if nargin < 4
        lw = 2;
    end

    % plot histogram
    ylast = 0;
    assert(numel(bins) == numel(counts)+1, 'expecting histcounts output');
    for kk = 1:numel(counts)
        xcs = bins(kk:kk+1);
        yc = counts(kk);
        plot(xcs, yc*[1 1], '-', ...
            'Color', clr, 'LineWidth', lw, 'HandleVisibility', 'off');
        plot(xcs(1)*[1 1], [ylast yc], '-', ...
            'Color', clr, 'LineWidth', lw, 'HandleVisibility', 'off');
        ylast = yc;
    end
end
