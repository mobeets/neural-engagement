function mdls = fitCosineTuning(X, Y)
    mdls = [];
    for ii = 1:size(Y,2)
        mdl = fitSingleUnit(X, Y(:,ii));
        mdl.cell_id = ii;
        mdls = [mdls mdl];
    end
end

function mdl = fitSingleUnit(X, y, alpha)
% X in radians, y in spike counts
% returns r_0, r_max (min and max mean spike count), and s_max (radians)
    if nargin < 3
        alpha = 0.05; % assess significance at this threshold
    end
    if range(X) > 100
        warning('Is X in radians? The range is large.');
    end

    A = [ones(size(X,1),1), cos(X), sin(X)];
    theta = (A'*A)\(A'*y);

    % To obtain the tuning curve parameters from the x variables, note we have:
    % r_0-(r_max-r_0)*cos(s_i-s_max) = x_1+x_2*cos(s_i)+x_3*sin(s_i).
    % Any cosine with a phase offset can be uniquely written as the weighted linear
    % combination of a cosine and sine with no offset. So we first get that
    % r_0 = x_1. Expanding the LHS using cos(x-y) = cos(x)cos(y)+sin(x)sin(y) and
    % finding like terms, we get x_2 = (r_0-r_max)*cos(s_max) and
    % x_3 = (r_0-r_max)*sin(s_max). Dividing these gets us the next line of code.
    s_max = atan(theta(3)/theta(2));
    r_0 = theta(1); % the offset on both sides of the equation must be same
    r_max = r_0 + theta(2)/cos(s_max);
    
    % if r_0 > r_max, then r_max and s_max are actually r_min and s_min
    if r_0 > r_max
        r_max = r_0 + (r_0 - r_max);
        s_max = s_max + pi;
    end
    
    % assess significance
    eps = y - A*theta; % model residuals
    var_hat = sum(eps.^2)/(size(A,1)-size(A,2));
    se = sqrt(var_hat*diag(inv(A'*A)));
    dof = size(A,1) - size(A,2);
%     mult = tinv(1-alpha/2, dof); % used for CI
    pvals = 2*(1 - tcdf(abs(theta)./se, dof)); % two-sided
    % n.b. pvals matches what's returned by fitlm's mdl.Coefficients.pValue
    isSig = any(pvals(2:end) < alpha); % to check if modulation depth > 0

    mdl.r_0 = r_0;
    mdl.r_max = r_max;
    mdl.s_max = s_max;
    mdl.pvals = pvals;
    mdl.alpha = alpha;
    mdl.isSig = isSig;
    mdl.predict = @(x) r_0 + (r_max - r_0)*cos(x - s_max);
end
