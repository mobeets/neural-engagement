function [yss, xss, ns] = smooth(xs, ys, K, nsteps, fcn)
% simple boxcar smooth that only includes full bins
%   (i.e., no edge handling)
% assumes ys is [nsamples x ndims]
% 
    if nargin < 3
        K = ys;
        assert(numel(ys) == 1, 'second arg should be integer');
        assert(nargout <= 1, 'only one output when xs is not provided.');
        ys = xs;                
        noX = true;
        if ~exist('nsteps', 'var')
            nsteps = 1;
        end
    else
        noX = false;
    end
    if nargin < 4 || isnan(nsteps)
        nsteps = 1; % number of steps each bin slides over
    end
    if nargin < 5 || isempty(fcn)
        fcn = @(ys) nanmean(ys,1);
    end
    ks = ones(K,1)/K;
    if size(ys,1) < K || (size(ys,1) == 1 && size(ys,2) > 1) || size(xs,1) == 1
        transposeY = true;
        ks = ks';
    else
        transposeY = false;
    end
    if noX
        yss = conv2(ys, ks, 'valid');
        return;
    end
    xsa = unique(xs);
    inds = 1:nsteps:(numel(xsa)-K+1);
    if transposeY
        ys = ys';
    end
    xss = xsa(inds');
    yss = nan(numel(inds), size(ys,2));
    ns = zeros(numel(inds), 1);
    for ii = 1:numel(inds)
        ix = (xsa(inds(ii)) <= xs) & (xs <= xsa(inds(ii)+K-1));
        yss(ii,:) = fcn(ys(ix,:));
        ns(ii) = sum(ix);
    end
    if transposeY
        yss = yss';
        xss = xss';
        ns = ns';
    end
end
