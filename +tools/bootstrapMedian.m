function [meds, ci_lb, ci_ub] = bootstrapMedian(vals, nboots, alpha)
    if nargin < 3
        alpha = 0.95; % for returning CI with level alpha
    end
    if size(vals,1) == 1
        vals = vals';
    end
    assert(size(vals,2) == 1, 'vals must be vector');
    N = numel(vals);    
    meds = nan(nboots,1);
    for ii = 1:nboots
        inds = round(0.5 + N*rand(N,1));
        meds(ii) = nanmedian(vals(inds));
    end
    
    if nargout > 1
        meds_sorted = sort(meds);
        n_lb = round((nboots*(1-alpha))/2);
        n_ub = nboots-n_lb;
        ci_lb = meds_sorted(n_lb);
        ci_ub = meds_sorted(n_ub);
    end
end
