function [engagement_dims_fine, engagement_dims, vexps] = getEngagementDimensions(X, ...
    Y, FactorAnalysisParams, grps_fine, flipByAdjacency, grps)
% 
% for use with tools.getAimingEllipseAndVigorDimensions
% 

    if nargin < 5
        flipByAdjacency = false;
    end
    if nargin < 6
        grps = tools.thetaCenters;
    end
    if isempty(FactorAnalysisParams) || ~isfield(FactorAnalysisParams, 'beta')
        warning(['No beta in FactorAnalysisParams, so dims not flipped.']);
        checkFlip = false;
    else
        checkFlip = true;
    end
        
    engagement_dims = nan(numel(grps), size(Y,2));
    vexps = nan(numel(grps),1);
    for kk = 1:numel(grps)
        ix = X == grps(kk);
        if sum(ix) == 0
            continue;
        end
        
        % handle columns that have no variability by skipping them
        if any(var(Y(ix,:)) == 0)
            indsToSkip = var(Y(ix,:)) == 0;
            Yt = Y(ix,~indsToSkip);
        else
            indsToSkip = false(size(Y,1),1);
            Yt = Y(ix,:);
        end
        [coeff, ~, ~, ~, vexp] = pca(Yt);
        vexps(kk) = vexp(1);
        engagement_dim = coeff(:,1);
        
        % set columns that had no variability to have 0 weight
        if sum(indsToSkip) > 0
            origInds = nan(size(Y,2),1);
            origInds(~indsToSkip) = 1:sum(~indsToSkip);
            new_vdim = zeros(size(origInds));
            for ll = 1:numel(engagement_dim)
                new_vdim(origInds == ll) = engagement_dim(ll);
            end
            engagement_dim = new_vdim;
        end
        
        if checkFlip && doFlipDim(engagement_dim, FactorAnalysisParams)
            engagement_dim = -engagement_dim;
        end
        engagement_dims(kk,:) = engagement_dim;
    end
    
    % flip based on adjacent dims
    if flipByAdjacency
        doFlips = tools.flipEngagementDimsUsingAdjacentDims(engagement_dims);
        for kk = 1:size(engagement_dims,1)
            if doFlips(kk)
                engagement_dims(kk,:) = -engagement_dims(kk,:);
            end
        end
    end

    % interpolate engagement dimensions
    engagement_dims_fine = nan(numel(grps_fine), size(Y,2));
    for kk = 1:numel(grps_fine)
        angs = tools.angdiff(grps_fine(kk), grps);
        [~,ix] = sort(angs);
        ws = angs(ix(1:2));
        ws = 1 - ws/sum(ws);
        vdim1 = engagement_dims(ix(1),:);
        vdim2 = engagement_dims(ix(2),:);
        vdim = ws(1)*vdim1 + ws(2)*vdim2;
        vdim = vdim/norm(vdim);
        engagement_dims_fine(kk,:) = vdim;
    end
%     return;
    % interpolate engagement dimensions with spline
    engagement_dims_fine = tools.interpCircular(engagement_dims, grps, grps_fine);
    nrms = tools.rowwiseNorm(engagement_dims_fine);
    engagement_dims_fine = bsxfun(@times, engagement_dims_fine, 1./nrms);
end

function doFlip = doFlipDim(cdim, FactorAnalysisParams)
    isPosFa1 = mean(sign(FactorAnalysisParams.beta(1,:)) == 1) > 0.5;
    isPosDim = cdim(1) > 0;
    if isPosFa1 ~= isPosDim
        doFlip = true;
    else
        doFlip = false;
    end
end
