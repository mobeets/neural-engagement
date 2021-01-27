function [vss,xss] = smoothInPieces(vs, ns, ks, xs, nsteps, fcn)
% smooths vs in pieces defined by the breakpoints in ns with smoothing
%   defined by ks
% i.e., smooths vs(ns(j):ns(j+1),:) using tools.smooth using ks(j)
% 
    if nargin < 2
        % default: no smoothing
        ns = [1];
        ks = [1];
    end
    if nargin < 4 || isempty(xs)
        xs = (1:size(vs,1))';
    end
    if nargin < 5
        nsteps = nan;
    end
    if nargin < 6
        fcn = [];
    end
    if (ns(end) < size(vs,1)) && (numel(ks) == numel(ns))
        ns = [ns inf];
    end
    if size(vs,2) > size(vs,1)
        warning(['default is that size(vs,1) is time,' ...
            ' but vs as provided is wider than it is long']);
    end
    assert(numel(ks) == numel(ns)-1);
    assert(size(vs,1) == size(xs,1));
    
    vss = [];
    xss = [];
    for ll = 1:(numel(ns)-1)
        xMin = ns(ll);
        xMax = min([size(xs,1) ns(ll+1)]);
        ixc = (xs >= xMin) & (xs <= xMax);
        xsc = xs(ixc);
        vsc = vs(ixc,:);
        if size(vsc,1) > 1
            [vsc, xsc] = tools.smooth(xsc, vsc, ks(ll), nsteps, fcn);
        elseif size(vsc,2) > 1
            % tools.smooth will auto-transpose, so we avoid calling it here
            if ks(ll) > 1
                vsc = []; xsc = [];
            end
        end
        if ~isempty(vsc)
            vss = [vss; vsc];
            xss = [xss; xsc];
        end
    end
end
