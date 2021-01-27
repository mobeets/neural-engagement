function [behav_vig, behav_base] = targetsHelpedByEngagement(ys_ellipse, ...
    engagement_dims, vfcnwmp, behavNm)
% 
% behavNm is one of: {'progress', 'angular error'}
% 
% returns the behavior of ys_ellipse vs. ys_ellipse + delta*engagement_dim
% 

    if strcmpi(behavNm, 'progress')
        behavFcn = @(vs, vsGoal) tools.getProgress([], [], vsGoal, [], vs);
    elseif strcmpi(behavNm, 'angular error')
        behavFcn = @(vs, vsGoal) tools.getAngularError(vs, vsGoal);
    end
    trgs = tools.thetaCenters(size(ys_ellipse,1));
    trgpos = [cosd(trgs) sind(trgs)];    
    trgsPert = vfcnwmp(ys_ellipse);
    behav_base = behavFcn(trgsPert, trgpos);
    
    if size(engagement_dims,1) == 1
        engagement_dims = repmat(engagement_dims, size(ys_ellipse,1), 1);
    end
    ysv = ys_ellipse + 0.01*engagement_dims;
    trgsPertPlus = vfcnwmp(ysv);
    behav_vig = behavFcn(trgsPertPlus, trgpos);

end
