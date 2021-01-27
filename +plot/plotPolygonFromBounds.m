function h = plotPolygonFromBounds(xs, lbs, ubs, clr)
    if nargin < 4
        clr = 0.8*ones(1,3);
    end
    if size(xs,2) == 1
        xs = xs';
        lbs = lbs';
        ubs = ubs';
    end
    assert(size(xs,1) == 1);
    assert(size(lbs,1) == 1);
    assert(size(ubs,1) == 1);
    if isequal(lbs,ubs)
        h = [];
        return;
    end
    if size(lbs,1) ~= size(xs,1)
        lbs = lbs';
        ubs = ubs';
    end

    h = patch([xs xs(end:-1:1)], [lbs ubs(end:-1:1)], clr, ...
        'LineStyle', 'none', 'FaceAlpha', 0.1, ...
        'HandleVisibility', 'off');
    uistack(h, 'bottom');
end
