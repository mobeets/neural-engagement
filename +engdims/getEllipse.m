function [Ysmu, Ysmu_raw] = getEllipse(X, Y, grps_fine, tol_deg, grps)
    if nargin < 4
        tol_deg = 45/2;
    end
    if nargin < 5
        grps = tools.thetaCenters;
    end
    
    % find repertoire    
    Ysmu = nan(numel(grps), size(Y,2));
    for kk = 1:numel(grps)
        ds = tools.angdiff(X, grps(kk));
        Ysmu(kk,:) = nanmean(Y(ds < tol_deg,:));
    end
    Ysmu_raw = Ysmu;
    if ~isequal(grps, grps_fine)
        Ysmu = tools.interpCircular(Ysmu, grps, grps_fine);
    end
end
